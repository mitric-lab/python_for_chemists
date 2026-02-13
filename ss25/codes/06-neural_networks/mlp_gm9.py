import numpy as np
import matplotlib.pyplot as plt

class Sigmoid():
    def __call__(self, x):
        return 1 / (1 + np.exp(-x))
    
    def gradient(self, x):
        return self(x) * (1 - self(x))

class ReLU():
    def __call__(self, x):
        return np.maximum(x, 0)
    
    def gradient(self, x):
        return np.heaviside(x, 0) 

class Identity():
    def __call__(self, x):
        return x
    
    def gradient(self, x):
        return np.ones_like(x)

class MLP():

    def __init__(self, sizes, tau=0.1, batch_size=5):
        # Initialize the network's parameters
        self.sizes = list(sizes.keys())
        self.activations = list(sizes.values())
        self.num_layers = len(self.sizes)
        self.weights = [np.random.randn(x, y) for x, y in zip(self.sizes[:-1], self.sizes[1:])]
        self.biases = [np.random.randn(y) for y in self.sizes[1:]]
        self.tau = tau
        self.batch_size = batch_size
            
    def feedforward(self, xi):
        # Compute the output of the network for input xi
        h_l = xi
        for b_l, w_l, activation_l in zip(self.biases, self.weights, self.activations[1:]):
            h_l = activation_l(np.dot(w_l.T, h_l) + b_l)
        return h_l

    def train_step(self, X, y):
        # Train the network using mini-batch gradient descent

        # Shuffle the training data
        N = X.shape[0]
        indices = np.arange(N)
        np.random.shuffle(indices)

        # Iterate over mini-batches
        for i in range(0, N, self.batch_size):
            batch_indices = indices[i:i+self.batch_size]
            self.update_mini_batch(X[batch_indices], y[batch_indices])
                
    def update_mini_batch(self, X_batch, y_batch):
        # Update the network's parameters using gradient descent

        # Initialize gradients for this minibatch
        gradient_weights = [np.zeros_like(w) for w in self.weights]
        gradient_bias = [np.zeros_like(b) for b in self.biases]

        # Iterate over all training pairs in this minibatch
        for xi, yi in zip(X_batch, y_batch):
            
            # Compute gradients for this training pair
            gradient_weights_i, gradient_bias_i = self.backprop(xi, yi)

            # Update the gradients for this training pair
            gradient_weights = [dw + dw_i for dw, dw_i in zip(gradient_weights, gradient_weights_i)]
            gradient_bias = [db + db_i for db, db_i in zip(gradient_bias, gradient_bias_i)]

        # Update the network's parameters using the gradients
        self.weights = [w - (self.tau / self.batch_size) * dw for w, dw in zip(self.weights, gradient_weights)]
        self.biases = [b - (self.tau / self.batch_size) * db for b, db in zip(self.biases, gradient_bias)]

    def backprop(self, xi, yi):
        # Compute the gradients of the network's parameters for input xi and target yi

        # Initialize gradients for this input
        gradient_weights_i = [np.zeros_like(w) for w in self.weights]
        gradient_bias_i = [np.zeros_like(b) for b in self.biases]
        
        ### Feedforward pass ###
        
        # First layer is just the input
        h_l = xi
        h_vectors = [h_l] # store hidden layer outputs
        a_vectors = [] # store activations
        
        # Perform forward pass through the network and store activations and outputs
        for w_l, b_l, activation_l in zip(self.weights, self.biases, self.activations[1:]):
            a_l = np.dot(w_l.T, h_l) + b_l
            a_vectors.append(a_l)
            h_l = activation_l(a_l)
            h_vectors.append(h_l)
        
        ### Backward pass ###

        # Compute delta for output layer
        delta_l = (h_vectors[-1] - yi) * self.activations[-1].gradient(a_vectors[-1])
        
        # Compute derivatives of parameters in output layer
        gradient_weights_i[-1] = np.outer(h_vectors[-2], delta_l)
        gradient_bias_i[-1] = delta_l
                
        # Iterate over hidden layers
        for l in range(2, self.num_layers):

            # Compute delta for this layer
            delta_l = self.activations[-l].gradient(a_vectors[-l]) * np.dot(delta_l, self.weights[-(l-1)].T)

            # Compute derivatives of parameters in this layer
            gradient_weights_i[-l] = np.outer(h_vectors[-(l+1)], delta_l)
            gradient_bias_i[-l] = delta_l
        
        return  gradient_weights_i, gradient_bias_i

    def train(self, X, y, epochs=10, validate=False):
        # Train or validate the network for a number of epochs

        for e in range(epochs):
            
            # Train the network
            if not validate:
                self.train_step(X, y)
            
            # Compute accuracy
            N = X.shape[0]
            avg_loss = 0
            for xi, yi in zip(X, y):
                avg_loss += np.abs(self.feedforward(xi)[0] - yi) / N
            
            print(f"Epoch: {e+1}, Loss: {avg_loss:.3f}")


# Load the data
qm9 = np.load("qm9.npy", allow_pickle=True)

X = np.array([x[0] for x in qm9])
y = np.array([x[1][7] for x in qm9])

# split the data into training and test sets
N = X.shape[0]
N_train = int(0.8 * N)
indices = np.arange(N)
np.random.shuffle(indices)
X_train = X[indices[:N_train]]
y_train = y[indices[:N_train]]

#print(y_train[:10])

X_test = X[indices[N_train:]]
y_test = y[indices[N_train:]]

sizes = {25: None,
         12: Sigmoid(),
         6: Sigmoid(),
         1: Identity()}

tau = 0.01
batch_size = 50
epochs = 20

f_hat = MLP(sizes, tau=tau, batch_size=batch_size)

print("Training the network")

f_hat.train(X_train, y_train, epochs=epochs, validate=False)

print("Validating the network")

f_hat.train(X_test, y_test, epochs=1, validate=True)
