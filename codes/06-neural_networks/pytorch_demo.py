#!/usr/bin/env python

import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset

# Set random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

### ANCHOR: pytorch_intro
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from torch.utils.data import DataLoader, TensorDataset

print("PyTorch version:", torch.__version__)
print("CUDA available:", torch.cuda.is_available())

# Check if we're using GPU
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")
### ANCHOR_END: pytorch_intro

### ANCHOR: data_preparation
# Load the same data we used for SLP
path = 'aptamer_xor_data.csv'
df = pd.read_csv(path)
print("Data shape:", df.shape)
print(df.head())

# Prepare data for PyTorch
X = df.drop('labels', axis=1).values.astype(np.float32)
y = df['labels'].values.astype(np.float32)

# Convert to PyTorch tensors
X_tensor = torch.tensor(X, dtype=torch.float32)
y_tensor = torch.tensor(y, dtype=torch.float32).unsqueeze(1)  # Add dimension for batch processing

# Create dataset and dataloader
dataset = TensorDataset(X_tensor, y_tensor)
dataloader = DataLoader(dataset, batch_size=32, shuffle=True)

print(f"Input shape: {X_tensor.shape}")
print(f"Output shape: {y_tensor.shape}")
### ANCHOR_END: data_preparation

### ANCHOR: model_definition
# Define a simple neural network using PyTorch's nn.Module
class SimpleMLP(nn.Module):
    def __init__(self, input_size=2, hidden_size=10, output_size=1):
        super(SimpleMLP, self).__init__()
        # Define layers
        self.layer1 = nn.Linear(input_size, hidden_size)
        self.layer2 = nn.Linear(hidden_size, hidden_size)
        self.layer3 = nn.Linear(hidden_size, output_size)
        
        # Activation function
        self.activation = nn.ReLU()
        
    def forward(self, x):
        # Forward pass through the network
        x = self.activation(self.layer1(x))
        x = self.activation(self.layer2(x))
        x = self.layer3(x)  # No activation on output layer for regression
        return x

# Create model instance
model = SimpleMLP(input_size=2, hidden_size=10, output_size=1)
model = model.to(device)  # Move to GPU if available

print("Model architecture:")
print(model)

# Count parameters
total_params = sum(p.numel() for p in model.parameters())
print(f"Total parameters: {total_params}")
### ANCHOR_END: model_definition

### ANCHOR: training_setup
# Define loss function and optimizer
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.01)

# Training parameters
num_epochs = 100
train_losses = []
### ANCHOR_END: training_setup

### ANCHOR: training_loop
# Training loop
model.train()
for epoch in range(num_epochs):
    epoch_loss = 0.0
    num_batches = 0
    
    for batch_X, batch_y in dataloader:
        # Move data to device
        batch_X = batch_X.to(device)
        batch_y = batch_y.to(device)
        
        # Forward pass
        optimizer.zero_grad()  # Clear gradients
        outputs = model(batch_X)
        loss = criterion(outputs, batch_y)
        
        # Backward pass (automatic differentiation)
        loss.backward()
        optimizer.step()
        
        epoch_loss += loss.item()
        num_batches += 1
    
    # Record average loss for this epoch
    avg_loss = epoch_loss / num_batches
    train_losses.append(avg_loss)
    
    if (epoch + 1) % 20 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {avg_loss:.4f}')

print("Training completed!")
### ANCHOR_END: training_loop

### ANCHOR: evaluation
# Evaluate the model
model.eval()
with torch.no_grad():
    X_eval = X_tensor.to(device)
    y_pred = model(X_eval)
    y_pred = y_pred.cpu().numpy().flatten()

# Calculate accuracy for classification
y_pred_class = np.sign(y_pred)
accuracy = np.mean(y_pred_class == y)
print(f"Classification accuracy: {accuracy:.4f}")

# Calculate MSE for regression
mse = np.mean((y_pred - y) ** 2)
print(f"Mean squared error: {mse:.4f}")
### ANCHOR_END: evaluation

### ANCHOR: visualization
# Create visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Plot decision boundary
x1_min, x1_max = X[:, 0].min() - 1, X[:, 0].max() + 1
x2_min, x2_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx1, xx2 = np.meshgrid(np.linspace(x1_min, x1_max, 100),
                       np.linspace(x2_min, x2_max, 100))

# Get predictions for the grid
grid_X = torch.tensor(np.c_[xx1.ravel(), xx2.ravel()], dtype=torch.float32).to(device)
with torch.no_grad():
    grid_pred = model(grid_X).cpu().numpy().reshape(xx1.shape)

# Plot decision boundary
ax1.contour(xx1, xx2, grid_pred, levels=[0], colors='black', linewidths=2)
ax1.contourf(xx1, xx2, grid_pred, levels=20, alpha=0.3)

# Plot data points
scatter = ax1.scatter(X[:, 0], X[:, 1], c=y, cmap='RdBu', edgecolors='black')
ax1.set_xlabel('Feature 1')
ax1.set_ylabel('Feature 2')
ax1.set_title('Decision Boundary')
plt.colorbar(scatter, ax=ax1)

# Plot training loss
ax2.plot(train_losses)
ax2.set_xlabel('Epoch')
ax2.set_ylabel('Loss')
ax2.set_title('Training Loss')
ax2.grid(True)

plt.tight_layout()
plt.show()
### ANCHOR_END: visualization

# plt.savefig('../../assets/figures/06-neural_networks/pytorch_demo.png')

### ANCHOR: advanced_features
# Demonstrate some advanced PyTorch features
print("\n=== Advanced PyTorch Features ===")

# 1. Different activation functions
class ActivationDemo(nn.Module):
    def __init__(self):
        super(ActivationDemo, self).__init__()
        self.linear = nn.Linear(2, 1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()
        self.tanh = nn.Tanh()
        
    def forward(self, x):
        x = self.linear(x)
        return {
            'linear': x,
            'relu': self.relu(x),
            'sigmoid': self.sigmoid(x),
            'tanh': self.tanh(x)
        }

# 2. Dropout for regularization
class MLPWithDropout(nn.Module):
    def __init__(self, input_size=2, hidden_size=20, output_size=1, dropout_rate=0.2):
        super(MLPWithDropout, self).__init__()
        self.layer1 = nn.Linear(input_size, hidden_size)
        self.layer2 = nn.Linear(hidden_size, hidden_size)
        self.layer3 = nn.Linear(hidden_size, output_size)
        self.dropout = nn.Dropout(dropout_rate)
        self.activation = nn.ReLU()
        
    def forward(self, x):
        x = self.dropout(self.activation(self.layer1(x)))
        x = self.dropout(self.activation(self.layer2(x)))
        x = self.layer3(x)
        return x

# 3. Learning rate scheduling
model_advanced = MLPWithDropout().to(device)
optimizer_advanced = optim.Adam(model_advanced.parameters(), lr=0.01)
scheduler = optim.lr_scheduler.StepLR(optimizer_advanced, step_size=30, gamma=0.5)

print("Learning rate scheduling example:")
for epoch in range(5):
    current_lr = scheduler.get_last_lr()[0]
    print(f"Epoch {epoch+1}, Learning rate: {current_lr:.6f}")
    scheduler.step()

print("\nPyTorch provides automatic differentiation, GPU acceleration, and many other features!")
### ANCHOR_END: advanced_features
