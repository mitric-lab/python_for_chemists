import torch
from torch import nn
from torchvision import transforms
import os
import torch
import torch.utils.data
from os import mkdir
from torch import optim
from torch.nn import functional as F
import torchvision
from torchvision.utils import save_image
from torchvision.datasets import CelebA
from torch.utils.data import DataLoader

import numpy as np
import matplotlib.pyplot as plt


if torch.cuda.is_available():
    device = 'cuda'
elif torch.backends.mps.is_available():
    device = 'mps'
else:
    device = 'cpu'

print(f'Using device: {device}')


class VAE(nn.Module):
    def __init__(self, IMAGE_SIZE=150, LATENT_DIM=128):
        super(VAE, self).__init__()

        hidden_dims = [32, 64, 128, 256, 512]
        self.final_dim = hidden_dims[-1]
        in_channels = 3
        modules = []

        # Build Encoder
        for h_dim in hidden_dims:
            modules.append(
                nn.Sequential(
                    nn.Conv2d(in_channels, out_channels=h_dim,
                              kernel_size=3, stride=2, padding=1),
                    nn.BatchNorm2d(h_dim),
                    nn.LeakyReLU())
            )
            in_channels = h_dim

        self.encoder = nn.Sequential(*modules)
        out = self.encoder(torch.rand(1, 3, IMAGE_SIZE, IMAGE_SIZE))
        self.size = out.shape[2]
        self.fc_mu = nn.Linear(hidden_dims[-1] * self.size * self.size, LATENT_DIM)
        self.fc_var = nn.Linear(hidden_dims[-1] * self.size * self.size, LATENT_DIM)

        # Build Decoder
        modules = []
        self.decoder_input = nn.Linear(LATENT_DIM, hidden_dims[-1] * self.size * self.size)
        hidden_dims.reverse()

        for i in range(len(hidden_dims) - 1):
            modules.append(
                nn.Sequential(
                    nn.ConvTranspose2d(hidden_dims[i],
                                       hidden_dims[i + 1],
                                       kernel_size=3,
                                       stride=2,
                                       padding=1,
                                       output_padding=1),
                    nn.BatchNorm2d(hidden_dims[i + 1]),
                    nn.LeakyReLU())
            )

        self.decoder = nn.Sequential(*modules)

        self.final_layer = nn.Sequential(
            nn.ConvTranspose2d(hidden_dims[-1],
                               hidden_dims[-1],
                               kernel_size=3,
                               stride=2,
                               padding=1,
                               output_padding=1),
            nn.BatchNorm2d(hidden_dims[-1]),
            nn.LeakyReLU(),
            nn.Conv2d(hidden_dims[-1], out_channels=3,
                      kernel_size=3, padding=1),
            nn.Sigmoid())

    def encode(self, x):
        result = self.encoder(x)
        result = torch.flatten(result, start_dim=1)
        mu = self.fc_mu(result)
        log_var = self.fc_var(result)
        return mu, log_var

    def reparameterize(self, mu, log_var):
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return eps * std + mu

    def decode(self, z):
        result = self.decoder_input(z)
        result = result.view(-1, self.final_dim, self.size, self.size)
        result = self.decoder(result)
        result = self.final_layer(result)
        result = celeb_transform1(result)
        result = torch.flatten(result, start_dim=1)
        result = torch.nan_to_num(result)
        return result

    def forward(self, x):
        mu, log_var = self.encode(x)
        z = self.reparameterize(mu, log_var)
        return self.decode(z), mu, log_var


# Reconstruction + KL divergence losses summed over all elements and batch
def loss_function(recon_x, x, mu, log_var):
    MSE =F.mse_loss(recon_x, x.view(-1, image_dim))
    KLD = -0.5 * torch.mean(1 + log_var - mu.pow(2) - log_var.exp())
    kld_weight = 0.00025
    loss = MSE + kld_weight * KLD  
    return loss


def train(epoch):
    model.train()
    train_loss = 0
    for batch_idx, (data, _) in enumerate(train_loader):
        torch.cuda.empty_cache()
        data = data.to(device)
        optimizer.zero_grad()
        recon_batch, mu, log_var = model(data)
        log_var = torch.clamp_(log_var, -10, 10)
        loss = loss_function(recon_batch, data, mu, log_var)
        loss.backward()
        train_loss += loss.item()
        optimizer.step()
        if batch_idx % 100 == 0:
            print('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                epoch, batch_idx * len(data), len(train_loader.dataset),
                       100. * batch_idx / len(train_loader),
                       loss.item() / len(data)))

    print('====> Epoch: {} Average loss: {:.4f}'.format(
        epoch, train_loss / len(train_loader.dataset)))


def test(epoch):
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for i, (data, _) in enumerate(test_loader):
            data = data.to(device)
            recon_batch, mu, log_var = model(data)
            test_loss += loss_function(recon_batch, data, mu, log_var).item()
            if i == 0:
                n = min(data.size(0), 8)
                comparison = torch.cat([data[:n],
                                        recon_batch.view(BATCH_SIZE, 3, IMAGE_SIZE, IMAGE_SIZE)[:n]])
                save_image(comparison.cpu(),
                           f'{directory}/reconstruction_{str(epoch)}.png', nrow=n)

    test_loss /= len(test_loader.dataset)
    print('====> Test set loss: {:.4f}'.format(test_loss))



if __name__ == "__main__":


    IMAGE_SIZE = 150
    image_dim = 3 * IMAGE_SIZE * IMAGE_SIZE
    LATENT_DIM = 128
    EPOCHS = 20  # number of training epochs
    BATCH_SIZE = 64  # for data loaders
    CELEB_PATH = 'data/'

    # Download the dataset once
    data = torch.utils.data.DataLoader(
            torchvision.datasets.CelebA(CELEB_PATH,
                transform=torchvision.transforms.ToTensor(),
                download=True),
            batch_size=128,
            shuffle=True)

    celeb_transform = transforms.Compose([
        transforms.Resize(IMAGE_SIZE, antialias=True),
        transforms.CenterCrop(IMAGE_SIZE),
        transforms.ToTensor()])  # used when transforming image to tensor

    celeb_transform1 = transforms.Compose([
        transforms.Resize(IMAGE_SIZE, antialias=True),
        transforms.CenterCrop(IMAGE_SIZE)])  # used by decode method to transform final output


    directory = 'models/'
    if not os.path.exists(directory):
        mkdir(directory)

    # download dataset
    train_dataset = CelebA(CELEB_PATH, transform=celeb_transform, download=False, split='train')
    test_dataset = CelebA(CELEB_PATH, transform=celeb_transform, download=False, split='valid') # or 'test'

    # create train and test dataloaders
    train_loader = DataLoader(dataset=train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    test_loader = DataLoader(dataset=test_dataset, batch_size=BATCH_SIZE, shuffle=False)


    model = VAE().to(device)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)


    for epoch in range(1, EPOCHS + 1):
        train(epoch)
        torch.save(model, f'{directory}/vae_model_{epoch}.pth')
        test(epoch)
        with torch.no_grad():
            sample = torch.randn(64, LATENT_DIM).to(device)
            sample = model.decode(sample).cpu()
            save_image(sample.view(64, 3, IMAGE_SIZE, IMAGE_SIZE),
                       f'{directory}/sample_{str(epoch)}.png')
            del sample