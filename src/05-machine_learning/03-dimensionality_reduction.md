## Dimensionality Reduction

So far, we have always considered data points as pairs $(\vec{x}_i, y_i)$ for $i = 1, \ldots, N$, where $\vec{x}_i \in \mathbb{R}^d$ is a feature vector and $y_i \in \mathbb{R}$ is a quantitative (regression) or qualitative (classification) label or target. In the case of **unsupervised learning**, this second information is missing and we only have data points $\{(\vec{x}_i)\}_{i=1, \dots, N}$ at our disposal. Here, we are particularly interested in *dimensionality reduction* or *clustering* of the data points.

Dimensionality reduction is the process of transforming high-dimensional data into a lower-dimensional space while preserving as much of the original information as possible. Remember that we use the term *dimensionality* to refer to the number of features in the data. Therefore, our goal is to find a new set of less features that capture as much of the original information as possible. This is particularly useful when dealing with high-dimensional data, such as images or text, where the number of features can be very large.

### Principal Component Analysis (PCA)

In one of the previous chapters, you have learned that by performing a singular value decomposition (SVD) of the data matrix $\bm{X} \in \mathbb{R}^{N \times d}$, we find the set of rank-1 matrices that best approximate the data matrix. In this context, you have also learned that the right singular vectors of the SVD are the so-called *principal components* of the data, and the corresponding singular values denote the contribution of the corresponding principal component to the data, also called *explained variance*. By projecting the data onto the first $k$ principal components, we can therefore reduce the dimensionality of the data from $d$ to $k$ while preserving as much of the original information as possible. This is the idea of *principal component analysis* (PCA).

Before implementing the PCA, we summarize the main steps of the PCA:

1. Center the data matrix $\bm{X}$ by subtracting the mean of each feature. If the features are not on the same scale, you might want to standardize the data by dividing by the standard deviation.
2. Compute the right singular vectors $\bm{V}$ of the SVD of the standardized data matrix $\bm{X}$. We can use the fact that the right singular vectors of the SVD are the eigenvectors of the *covariance matrix* $\bm{C} = \bm{X}^T \bm{X}$:
$$
\bm{C} = \bm{X}^T \bm{X} = (\bm{U} \bm{\Sigma} \bm{V}^T)^T \bm{U} \bm{\Sigma} \bm{V}^T = \bm{V} \bm{\Sigma}^T \bm{U}^T \bm{U} \bm{\Sigma} \bm{V}^T = \bm{V} \bm{\Sigma}^2 \bm{V}^T\,,
$$
where you should recognize the last equation from the eigenvalue decomposition.
3. Sort the eigenvectors by the corresponding eigenvalues in descending order and keep the first $k$ eigenvectors with the largest eigenvalues, given as matrix $\bm{V}_k$.
4. Transform the data matrix $\bm{X}$ into the space of the principal components by $\bm{Z}_k = \bm{X} \bm{V}_k$.

### Implementation of PCA

We will now implement the PCA, of course using our latest skills on object-oriented programming. First, we define the `__init__` method to set the class attributes, such as the number of principal components to keep, the components themselves, and the explained variance (squared singular values). They can be useful for further analysis.

```python
{{#include ../codes/05-machine_learning/pca.py:pca_init}}
```

As with the previous ML models, we will implement a `fit` method to compute the low-dimensional representation of the data according to the steps above. To compute the covariance matrix, we use the numpy function `np.cov`. After computing the eigenvectors and eigenvalues, we sort them by the corresponding eigenvalues in descending order using the [`np.argsort`](https://numpy.org/doc/stable/reference/generated/numpy.argsort.html) function. This function returns the indices of the sorted elements as an array, which we can then use directly to sort the eigenvectors (columns of the eigenvector matrix). In the last step, we store the first $k$ eigenvectors as the principal components. and compute the explained variance.

```python
{{#include ../codes/05-machine_learning/pca.py:pca_fit}}
```

Because we are doing unsupervised learning, we do not need a `predict` method. However, we will implement a `transform` method to project the data onto the principal components that we have computed in the `fit` method. To be sure, we again center the data and then use the `np.dot` function to compute the projection. Because we typically want to perform the projection in the same step as the fitting, we also implement a `fit_transform` method that calls both `fit` and `transform`. Separating the fitting and the projection is useful if we want to perform the projection on new data that we have not seen before.

```python
{{#include ../codes/05-machine_learning/pca.py:pca_transform}}
```

### PCA on the aptamer dataset

We want to apply the PCA to the aptamer dataset to find a low-dimensional representation of the molecules. We will use the `fit_transform` method to compute the principal components and the low-dimensional representation of the data. Instead of the more or less randomly hand-crafted features that we used in the chapter on regression, we will use a more *natural* representation 

