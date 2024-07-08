#!/usr/bin/env python

raise ValueError("This script is intended to be run as a Jupyter notebook, not as a standalone script.")

### ANCHOR: slp_interactive
import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import interact, FloatSlider
#%matplotlib inline

# Define the sigmoid activation function
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

# Initialize weights and biases
w = np.random.randn(2)
b = np.random.randn(2)
a = np.random.randn(2)

# Define the neural network forward pass
def neural_network(x, w, b, a):
    h = np.dot(w, x) + b
    return np.dot(a, sigmoid(h))

# Plotting function with sliders
def plot_network(w0, w1, b0, b1, a0, a1):
    # Update weights and biases based on slider values
    w = np.array([w0, w1])
    b = np.array([b0, b1])
    a = np.array([a0, a1])
    
    # Generate input values and compute network output
    x_values = np.linspace(-10, 10, 400)
    y_values = [neural_network(x, w, b, a) for x in x_values]

    # Target function
    f = lambda x: np.exp(-x**2)

    # Plotting
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(x_values, y_values, label='SLP output')
    ax.plot(x_values, f(x_values), label='Target function')
    ax.set_xlabel('Input')
    ax.set_ylabel('Output')
    ax.legend()
    plt.show()

# Create sliders for weights and biases
interact(plot_network,
            w0=FloatSlider(value=w[0], min=-5, max=5, step=0.1),
            w1=FloatSlider(value=w[1], min=-5, max=5, step=0.1),
            b0=FloatSlider(value=b[0], min=-5, max=5, step=0.1),
            b1=FloatSlider(value=b[1], min=-5, max=5, step=0.1),
            a0=FloatSlider(value=a[0], min=-5, max=5, step=0.1),
            a1=FloatSlider(value=a[1], min=-5, max=5, step=0.1))
### ANCHOR_END: slp_interactive