## Regression

You might be wondering why we are revisiting regression. The reason is that linear regression is a fundamental supervised machine learning algorithm used to model the relationship between a set of input features and a continuous output. This section will build upon what you've learned, so don't worry, we'll keep it concise.

### Linear Regression with Multidimensional Data

In the first chapter, we used linear regression to model a relationship between a single input and a single output. In the real world, however, data is rarely described by a single value. Instead, data is often of higher *dimensionality* and is characterized by multiple *features*.

```admonish note title="Note"
In this context, we'll use the term *dimensionality* to refer to the number of features describing a single data point. Think of it as the number of independent variables needed to characterize that point. For example, a 1000x1000 pixel image has a dimensionality of 1,000,000, since changing even one pixel alters the entire image.
```

Let's frame the problem mathematically. We have a dataset $\mathcal{D} = \{(\vec{x}_i, y_i)\}_{i=1, \dots, N}$, where each $\vec{x}_i \in \mathbb{R}^d$ is a $d$-dimensional input vector (our feature vector), and $y_i \in \mathbb{R}$ is the corresponding continuous output (target). Our goal is to find a function $\hat{f}_{\vec{w}}: \mathbb{R}^d \to \mathbb{R}$ with parameters $\vec{w}$ that approximates this relationship, such that $\hat{f}_{\vec{w}}(\vec{x}_i) \approx y_i$ for all data points.

For data with one feature ($d=1$), you'll recall we needed two parameters: $\hat{f}_{w_0, w_1}(x) = w_0 + w_1 x$. For $d$-dimensional data, we need $d+1$ parameters—one for each feature and one bias term:

$$
\hat{f}_{\vec{w}}(\vec{x}) = w_0 + w_1 x_1 + w_2 x_2 + \dots + w_d x_d = w_0 + \sum_{j=1}^{d} w_j x_j
$$

To simplify this, we can augment our input vector $\vec{x}$ by adding a constant value of 1 as its first element, making it $\vec{x}' = (1, x_1, \dots, x_d)$. We can then combine our weights into a single vector $\vec{w} = (w_0, w_1, \dots, w_d)$. Now, the model becomes a simple dot product:

$$
\hat{f}_{\vec{w}}(\vec{x}') = \vec{w} \cdot \vec{x}'
$$

We can extend this to the entire dataset by organizing our inputs into a matrix $\bm{X} \in \mathbb{R}^{N \times (d+1)}$, where each row is an augmented input vector $\vec{x}_i'$. The predictions for all data points can then be computed in a single matrix-vector product:

$$
\hat{\bm{y}} = \bm{X} \vec{w}
$$

where $\hat{\bm{y}} = (\hat{f}_{\vec{w}}(\vec{x}'_1), \dots, \hat{f}_{\vec{w}}(\vec{x}'_N))^T \in \mathbb{R}^N$. The data matrix $\bm{X}$ looks like this:

$$
\bm{X} = \begin{pmatrix}
1 & x_{1,1} & \dots & x_{1,d} \\
1 & x_{2,1} & \dots & x_{2,d} \\
\vdots & \vdots & \ddots & \vdots \\
1 & x_{N,1} & \dots & x_{N,d}
\end{pmatrix}
$$

This compact notation is standard in machine learning and will be used throughout the rest of this chapter.

### Main Components of a Machine Learning Problem

With the above formulation, we have defined the first major component of a supervised machine learning problem: **the model**. Before we can apply it, we need two more components:

The **loss function**, which measures the error between our model's predictions and the true outputs. For linear regression, this is typically the Mean Squared Error (MSE):

$$
\mathcal{L}(\hat{\bm{y}}, \bm{y}) = \frac{1}{2} \sum_{i=1}^N (\hat{f}_{\vec{w}}(\vec{x}_i) - y_i)^2
$$

The **optimization method**, which is the algorithm used to find the model parameters $\hat{\vec{w}} = (w_0, w_1, \dots, w_d)$ that minimize the loss function. While we can always use an iterative method like Gradient Descent, linear regression benefits from having a closed-form analytical solution:

$$
\hat{\vec{w}} = (\bm{X}^T \bm{X})^{-1} \bm{X}^T \bm{y} \,.
$$

where $\bm{y}$ is the vector of true outputs.

```admonish note title="Proof" collapsible=true
The loss function for least squares is given by:

$$
\begin{aligned}
\mathcal{L} &= \frac{1}{2} \sum_{i=1}^N (\hat{f}_{\vec{w}}(\vec{x}_i) - y_i)^2 \\
&= \frac{1}{2} \left\| \bm{X} \vec{w} - \vec{y} \right\|^2_2
\end{aligned}
$$

The gradient of the loss function with respect to $\vec{w}$ is:

$$
\begin{aligned}
\nabla_{\vec{w}} \mathcal{L} &= \nabla_{\vec{w}} \left( \frac{1}{2} (\bm{X}\vec{w} - \vec{y})^T (\bm{X}\vec{w} - \vec{y}) \right) \\
&= \bm{X}^T (\bm{X} \vec{w} - \vec{y})
\end{aligned}
$$

Setting the gradient to zero to find the minimum, we get:

$$
\begin{aligned}
\bm{X}^T (\bm{X} \hat{\vec{w}} - \vec{y}) &= 0 \\
\bm{X}^T \bm{X} \hat{\vec{w}} &= \bm{X}^T \vec{y}
\end{aligned}
$$

Multiplying from the left by $(\bm{X}^T \bm{X})^{-1}$ gives the solution:

$$
\hat{\vec{w}} = (\bm{X}^T \bm{X})^{-1} \bm{X}^T \vec{y}
$$

Linear algebra tells us that the matrix $\bm{X}^T \bm{X}$ is invertible if the columns of $\bm{X}$ are linearly independent. For the common case where $N > d$, this means we have *independent features*.

If the matrix is not invertible (i.e., it is singular), we can use the Moore-Penrose pseudoinverse, defined via the Singular Value Decomposition of $\bm{X}^T \bm{X} = \bm{V} \bm{\Sigma} \bm{V}^T$:

$$
(\bm{X}^T \bm{X})^+ = \bm{V} 
\begin{pmatrix}
1/\sigma_1^2 & & & \\
& \ddots & & \\
& & 1/\sigma_k^2 & \\
& & & 0
\end{pmatrix}
\bm{V}^T
$$

where $k$ is the rank of the matrix. This gives the optimal solution with the minimum possible norm.
```