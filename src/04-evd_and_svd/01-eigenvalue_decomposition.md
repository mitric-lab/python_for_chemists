## Eigenvalue Decomposition

To understand the definition and power of the eigenvalue decomposition (EVD), we first need to understand the concept of a vector space. Although this concept is very general, we will restrict ourselves here to $\mathbb{C}^n$, the $n$-dimensional vector space over the complex numbers. Vectors in $\mathbb{C}^n$ can be represented as column vectors in the standard basis. The notation $\vec{v} \in \mathbb{C}^n$ denotes such a vector:
$$
  \vec{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix}
$$
with components $v_i \in \mathbb{C}$. The Hermitian conjugate of a column vector, denoted $\vec{v}^\dag$, is its conjugate transpose (a row vector with complex conjugated components):
$$
  \vec{v}^\dag = \begin{pmatrix} v_1^* & v_2^* & \cdots & v_n^* \end{pmatrix}
$$
with the complex conjugated components $v_i^*$. Two vectors can be mapped to a scalar by the standard scalar product, defined as
$$
  \langle \vec{v}, \vec{w} \rangle = \vec{v}^\dag \vec{w} = \begin{pmatrix} v_1^* & v_2^* & \cdots & v_n^* \end{pmatrix} \begin{pmatrix} w_1 \\ w_2 \\ \vdots \\ w_n \end{pmatrix} = \sum_{i=1}^n v_i^* w_i\,.
$$

### Theoretical Foundations

Now, let's explore the relationship between vectors and matrices. We are particularly interested in vectors that have a special relationship with a given square matrix. For a square matrix $\bm{A} \in \C{n}{n}$, a non-zero vector $\vec{v}_i \in \mathbb{C}^n$ is called an eigenvector of $\bm{A}$ if
$$
  \bm{A} \vec{v}_i = \lambda_i \vec{v}_i
$$
for some scalar $\lambda \in \mathbb{C}$. This equation states that when matrix $\bm{A}$ acts on an eigenvector $\vec{v}_i$, the result is simply the same vector scaled by a number. The scalar $\lambda_i$ is the eigenvalue of $\bm{A}$ corresponding to the eigenvector $\vec{v}_i$. In general, for a matrix of dimension $n \times n$, there are $n$ eigenvalues and corresponding eigenvectors.

Eigenvectors and eigenvalues have powerful properties, especially for the class of [*normal matrices*](https://en.wikipedia.org/wiki/Normal_matrix), which are central to many applications in science and engineering. A matrix $\bm{A}$ is normal if it commutes with its Hermitian conjugate, i.e., $\bm{A}^\dag \bm{A} = \bm{A} \bm{A}^\dag$. This class includes important matrix types like Hermitian (if $\bm{A} = \bm{A}^\dag$) and real symmetric matrices (if $\bm{A} = \bm{A}^T$). For any normal matrix, its eigenvectors are orthogonal to each other, i.e., their scalar product is zero for different eigenvectors. If normalized, they form a complete orthonormal basis, meaning $\vec{v}_i^\dag \vec{v}_j = \delta_{ij}$. This allows any vector $\vec{x} \in \mathbb{C}^n$ to be expressed as a unique linear combination of these eigenvectors:
$$
  \vec{x} = \sum_{i=1}^n \langle \vec{x}, \vec{v}_i \rangle \vec{v}_i\,.
$$
While not all matrices are normal, they are the ones we will focus on for the EVD.

```admonish note title="Example" collapsible=true
Let $\bm{A} = \begin{pmatrix} 1 & 2 \\ 2 & 1 \end{pmatrix}$. This is a real symmetric matrix, so it is normal. Its eigenvectors and eigenvalues are
$$
  \vec{v}_1 = \begin{pmatrix} 1 \\ 1 \end{pmatrix} \quad \text{and} \quad \lambda_1 = 3\,,
$$
$$
  \vec{v}_2 = \begin{pmatrix} 1 \\ -1 \end{pmatrix} \quad \text{and} \quad \lambda_2 = -1\,,
$$
meaning that $\bm{A} \vec{v}_1 = 3 \vec{v}_1$ and $\bm{A} \vec{v}_2 = -1 \vec{v}_2$. You can also verify that the eigenvectors are orthogonal: $\vec{v}_1^\dag \vec{v}_2 = 1 \cdot 1 + 1 \cdot (-1) = 0$.
```

We can collect all eigenvalues into a diagonal matrix $\bm{\Lambda}$:
$$
  \bm{\Lambda} = \begin{pmatrix}
    \lambda_1 & 0 & \cdots & 0 \\
    0 & \lambda_2 & \cdots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \cdots & \lambda_n
  \end{pmatrix}
$$
and all eigenvectors as columns in a matrix $\bm{Q}$:
$$
  \bm{Q} = \begin{pmatrix}
    \vert & \vert &  & \vert \\
    \vec{v}_1 & \vec{v}_2 & \cdots & \vec{v}_n \\
    \vert & \vert &  & \vert
  \end{pmatrix}\,.
$$
Then, the eigenvalue decomposition of $\bm{A}$ can be written compactly as:
$$
  \bm{A} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}\,.
$$

### Eigenvalue Decomposition for Normal Matrices

Computing a matrix inverse $\bm{Q}^{-1}$ is computationally expensive. We can simplify the EVD for normal matrices, where the inverse is easy to compute. As stated before, a key property of normal matrices is that they are unitarily diagonalizable. This means their matrix of eigenvectors, $\bm{Q}$, can be chosen to be a unitary matrix (assuming the eigenvectors have been normalized to unit length). For a unitary matrix, the inverse is simply its Hermitian conjugate: $\bm{Q}^{-1} = \bm{Q}^\dag$. This makes finding the inverse trivial.
$$
  \bm{Q}^{-1} = \bm{Q}^\dag = \begin{pmatrix}
    \text{---}\ \vec{v}_1^\dag\ \text{---} \\
    \text{---}\ \vec{v}_2^\dag\ \text{---} \\
    \vdots \\
    \text{---}\ \vec{v}_n^\dag\ \text{---}
  \end{pmatrix}\,.
$$

For a normal matrix $\bm{A}$, the eigenvalue decomposition can thus be written and expanded as
$$
  \begin{align*}
    \bm{A} &= \bm{Q} \bm{\Lambda} \bm{Q}^\dag \\
    &= \begin{pmatrix}
      \vert & \vert &  & \vert \\
      \vec{v}_1 & \vec{v}_2 & \cdots & \vec{v}_n \\
      \vert & \vert &  & \vert
    \end{pmatrix}
    \begin{pmatrix}
      \lambda_1 & 0 & \cdots & 0 \\
      0 & \lambda_2 & \cdots & 0 \\
      \vdots & \vdots & \ddots & \vdots \\
      0 & 0 & \cdots & \lambda_n
    \end{pmatrix}
    \begin{pmatrix}
      \text{---}\ \vec{v}_1^\dag\ \text{---} \\
      \text{---}\ \vec{v}_2^\dag\ \text{---} \\
      \vdots \\
      \text{---}\ \vec{v}_n^\dag\ \text{---}
    \end{pmatrix} \\
    &= \lambda_1 \vec{v}_1 \vec{v}_1^\dag + \lambda_2 \vec{v}_2 \vec{v}_2^\dag + \cdots + \lambda_n \vec{v}_n \vec{v}_n^\dag \\
    &= \sum_{i=1}^n \lambda_i \vec{v}_i \vec{v}_i^\dag \,.
    {{numeq}}{eq:eigenvalue_decomposition_normal}
  \end{align*}
$$

Note that in the last line, the product of a column vector with a row vector, $\vec{v}_i \vec{v}_i^\dag$, does not produce a scalar like the scalar product. Instead, this is the [*dyadic product*](https://en.wikipedia.org/wiki/Dyadics) (or outer product) and produces a (rank-1) matrix:
$$
  \vec{v} \vec{w}^\dag = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix} \begin{pmatrix} w_1^* & w_2^* & \cdots & w_n^* \end{pmatrix} = \begin{pmatrix} v_1 w_1^* & v_1 w_2^* & \cdots & v_1 w_n^* \\ v_2 w_1^* & v_2 w_2^* & \cdots & v_2 w_n^* \\ \vdots & \vdots & \ddots & \vdots \\ v_n w_1^* & v_n w_2^* & \cdots & v_n w_n^* \end{pmatrix}\,.
$$
This form clearly shows why the EVD is called a decomposition: it breaks the matrix $\bm{A}$ down into a sum of simple rank-1 matrices ($\vec{v}_i \vec{v}_i^\dag$), each scaled by its corresponding eigenvalue $\lambda_i$.

### Eigenvalue Algorithms

While eigenvalues can be found analytically for small matrices (generally for $n < 5$), for larger matrices we must rely on numerical algorithms. We will discuss two such methods below.

#### Power Iteration

Perhaps the simplest eigenvalue algorithm is [*power iteration*](https://en.wikipedia.org/wiki/Power_iteration). It is an iterative method that finds the eigenvector corresponding to the eigenvalue with the largest magnitude (absolute value). For this method to work, the matrix $\bm{A}\in\C{n}{n}$ must be diagonalizable, and its largest magnitude eigenvalue must be unique. Let's assume the eigenvalues of $\bm{A}$ are ordered such that $|\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|$. The iteration proceeds as follows, starting with an arbitrary, non-zero initial vector $\vec{x}^0$:
$$
  \vec{x}^{k+1} = \frac{\bm{A} \vec{x}^k}{\|\bm{A} \vec{x}^k\|}\quad k=0,1,2,\ldots
  {{numeq}}{eq:power_iteration}
$$
The sequence of vectors $(\vec{x}^k)_{k\in\mathbb{N}}$ converges to the eigenvector corresponding to the dominant eigenvalue, $\lambda_1$. The division by the norm in each step prevents the vector's magnitude from growing or shrinking uncontrollably.

Once the eigenvector $\vec{x}^k$ has converged, the corresponding eigenvalue can be found using the Rayleigh quotient:
$$
  \lambda_1 = \frac{x^{k\dag} \bm{A} \vec{x}^k}{x^{k\dag} \vec{x}^k}\,.
  {{numeq}}{eq:power_iteration_eigenvalue}
$$

For those interested in the mathematics, a proof of the convergence of the power iteration can be found here.

```admonish proof title="Proof of Convergence of Power Iteration" collapsible=true
A matrix is diagonalisable if and only if its eigenvectors form a basis for the vector space. Therefore, we can write the arbitrary starting vector $\vec{x}^0$ as a linear
combination of the eigenvectors:
$$
  \vec{x}^0 = \sum_{i=1}^n a_i \vec{v}_i
$$
with coefficients $a_i\in\mathbb{C}$ and eigenvectors $\vec{v}_i$
corresponding to the eigenvalues $\lambda_i$. Applying $\bm{A}$ repeatedly gives:
$$
  \begin{align*}
    \bm{A}^k \vec{x}^0 &= \bm{A}^k \sum_{i=1}^n a_i \vec{v}_i = \sum_{i=1}^n a_i \bm{A}^k \vec{v}_i \\
    &= \sum_{i=1}^n a_i \lambda_i^k \vec{v}_i \\
    &= \lambda_1^k a_1 \vec{v}_1 + \sum_{i=2}^n \lambda_i^k a_i \vec{v}_i \\
    &= \lambda_1^k \left[ 
      a_1 \vec{v}_1 + \sum_{i=2}^n \left(\frac{\lambda_i}{\lambda_1}\right)^k a_i \vec{v}_i
    \right]\,.
  \end{align*}
$$

Because we assumed $|\lambda_i| < |\lambda_1|$ for all $i>1$, the ratio $|\lambda_i/\lambda_1|$ is less than 1. Therefore, as $k\to\infty$, the term $(\lambda_i/\lambda_1)^k$ converges to zero. The sum in the parenthesis vanishes, leaving:
$$
  \lim_{k\to\infty} \bm{A}^k \vec{x}^0 = \lim_{k\to\infty} \lambda_1^k a_1 \vec{v}_1 \propto \vec{v}_1\,.
$$
After a finite number of iterations, the resulting vector is therefore a good approximation of the eigenvector $\vec{v}_1$. We only care about the direction of the vector, as any scalar multiple of an eigenvector is also an eigenvector. The normalization in each step of {{eqref: eq:power_iteration}} prevents the $\lambda_1^k$ term from causing numerical overflow and keeps the vector's length stable.

A generalisation of the power iteration is the so-called *inverse iteration*, which can be used to determine an eigenvalue close to a given but arbitrary number $\mu$ and its eigenvector.
```

#### QR Algorithm

However, we often need all eigenvalues and eigenvectors of a matrix, not just the largest one. A frequently used method that can determine all eigenpairs of a diagonalisable matrix is the [*QR algorithm*](https://en.wikipedia.org/wiki/QR_algorithm).

The key step of the QR algorithm is the decomposition of a matrix $\bm{A}$ into a unitary matrix $\bm{Q}$ and an upper triangular matrix $\bm{R}$:
$$
  \bm{A} = \bm{Q} \bm{R}\,.
$$
Remember that a matrix is *unitary* if its inverse is equal to its Hermitian conjugate ($\bm{Q}^\dag \bm{Q} = \bm{Q} \bm{Q}^\dag = \identity$), and an *upper triangular* matrix has all elements below the diagonal being zero. We will assume that we can compute this QR decomposition for any square matrix.

The QR algorithm is an iterative method. Starting with $\bm{A}_0 = \bm{A}$, each step is defined by:
1. Decompose the current matrix: $\bm{A}_k = \bm{Q}_k \bm{R}_k$.
2. Compute the next matrix by multiplying the factors in reverse order:
$$
  \bm{A}_{k+1} = \bm{R}_k \bm{Q}_k\,.
  {{numeq}}{eq:qr_algorithm}
$$

By applying this transformation iteratively, we can *trace back* the matrix $\bm{A}$ to its original form $\bm{A}_0$:
$$
  \begin{align*}
    \bm{A}_{k+1} &= \bm{R}_k \bm{Q}_k = \textcolor{blue}{\identity} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \textcolor{blue}{\bm{Q}_k^\dag \bm{Q}_k} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{darkcyan}{\bm{A}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\identity} \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\bm{Q}_{k-1}^\dag \bm{Q}_{k-1}} \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\bm{Q}_{k-1}^\dag} \textcolor{brown}{\bm{A}_{k-1}} \textcolor{darkcyan}{\bm{Q}_{k-1}} \bm{Q}_k \\
    &  \phantom{=} \vdots \\
    &= \bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{A}_0} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k \\
    &= \underbrace{\bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{Q}}}_{:=\bm{P}_k} 
      \textcolor{red}{\bm{\Lambda}}
      \underbrace{\textcolor{red}{\bm{Q}^{-1}} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k}_{=\bm{P}_k^{-1}} \\
    &= \bm{P}_k \textcolor{red}{\bm{\Lambda}} \bm{P}_k^{-1} \,,
  \end{align*}
$$
where in the last step, we have inserted the exact but unknown eigenvalue decomposition of the initial matrix $\bm{A}_0 = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$. From this, we can see that all matrices $\bm{A}_k$ have the same eigenvalues as the initial matrix $\bm{A}_0$.

How does this help us? It can be shown that for a general matrix, the sequence $(\bm{A}_k)_{k\in\mathbb{N}_0}$ converges to an upper triangular matrix. The eigenvalues of an upper triangular matrix are simply its diagonal entries.
$$
  \lim_{k\to\infty} \bm{A}_k = 
  \begin{pmatrix}
    \lambda_1 & \cdots & * \\
    \vdots & \ddots & \vdots \\
    0 & \cdots & \lambda_n
  \end{pmatrix}
  \,.
$$
Since all $\bm{A}_k$ have the same eigenvalues as $\bm{A}$, the diagonal of the resulting matrix contains the eigenvalues of our original matrix $\bm{A}$. If the initial matrix $\bm{A}$ is normal, the sequence converges to a fully diagonal matrix.

The matrix of eigenvectors can also be found. It is the product of all the unitary matrices from each step:
$$
  \bm{Q} \approx \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k = \prod_{i=0}^k \bm{Q}_i\,.
  {{numeq}}{eq:qr_algorithm_eigenvectors}
$$
This allows us to approximate the eigenvectors of the matrix $\bm{A}$ as well.

### Implementation

We will now implement the power iteration and QR algorithms. For testing, we will generate random real matrices. To ensure each test matrix is diagonalizable with real eigenvalues, we will make it symmetric by adding it to its transpose. A symmetric matrix is a type of normal matrix.

#### Power Iteration

After importing NumPy
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:import}}
```
we implement the power iteration according to Eq. {{eqref: eq:power_iteration}} as a function:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:power_iteration_function}}
```
This function accepts the matrix `a`, a convergence tolerance `eps`, and a maximum number of iterations `maxiter`. It returns the eigenvalue with the largest magnitude and the corresponding eigenvector.

In the function body, we first determine the matrix dimension `n` from `a.shape[0]`. We then initialize a random starting vector `x` using [`np.random.rand`](https://numpy.org/doc/stable/reference/random/generated/numpy.random.rand.html). This vector is immediately normalized by dividing by its norm, calculated with [`np.linalg.norm`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html). To check for convergence, we need to compare the vector between two consecutive steps, so we initialize a variable `x_prev` to store the vector from the previous iteration.

The iteration then begins with a `for` loop that runs up to `maxiter` times. In each step, the current `x` is saved to `x_prev`. The new `x` is calculated by the matrix-vector product `np.dot(a, x)` and then normalized. We then check if the norm of the difference between `x` and `x_prev` is less than `eps`. If it is, the vector has converged, and the loop is terminated with `break`. If the loop finishes its full course without the `break` statement being reached, it means the convergence criterion was not met within the allowed iterations. A warning is then printed to the console. Finally, the corresponding eigenvalue `lambda_` is calculated using the Rayleigh quotient, and both `lambda_` and the eigenvector `x` are returned.

```admonish tip title="Choosing variable names"
Sometimes a desired variable name is already a reserved keyword in Python, like `lambda`. In this case, you can either choose a different name or, by convention, append an underscore. We used `lambda_` to store the eigenvalue.
```

To test the implementation, we generate a symmetric random matrix `a_mat` and call the `power_iteration` function:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:test_power_iteration}}
```
To verify our results, we use the function [`np.linalg.eigh`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html), which computes the eigenvalues and eigenvectors of a symmetric matrix. The eigenvalues are returned in ascending order in the 1D array `eigvals`, and the normalised eigenvectors as **columns** in the 2D array `eigvecs`. For this reason, we take the last entry of the eigenvalues and the corresponding eigenvector to compare them with the results of the power iteration:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:verify_results}}
```
While we can directly compare the eigenvalues with [`np.isclose`](https://numpy.org/doc/stable/reference/generated/numpy.isclose.html), we must consider that the eigenvectors are only uniquely defined up to a scaling factor. Since our eigenvectors are normalised, as well as those of `np.linalg.eigh`, we can compute the projection of the eigenvectors onto each other using `np.dot`. The absolute value of the projection should be close to 1, which we check with [`np.allclose`](https://numpy.org/doc/stable/reference/generated/numpy.allclose.html).

#### QR Algorithm

Now we will implement the QR algorithm according to Eq. {{eqref: eq:qr_algorithm}} and Eq. {{eqref: eq:qr_algorithm_eigenvectors}}. For this, we define the function `qr_algorithm`, which has the same signature as the function `power_iteration`:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:qr_algorithm_function}}
```
In the body of the function, the input matrix `a` is copied to `a_k`. It's important to make an explicit copy so the original matrix `a` is not modified by the iterative process. To compute the final matrix of eigenvectors by accumulating the products of $\bm{Q}_k$, we initialize a variable `eigvecs` as an identity matrix.

Within the `for` loop, which runs up to `maxiter`, the QR decomposition of `a_k` is first computed using the function [`np.linalg.qr`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.qr.html). Then, `a_k` is overwritten with the product of `r_k` and `q_k`, and `eigvecs` is multiplied by `q_k`.

Since `a_k` approaches an upper triangular matrix (or in this case a diagonal matrix), we can use the largest absolute value of the elements in the strict lower triangle of the matrix as a termination condition. For this, the function [`np.tril`](https://numpy.org/doc/stable/reference/generated/numpy.tril.html) is used with the argument `k=-1` to obtain the lower triangle of the matrix. If this value is smaller than the error tolerance `eps`, the loop is terminated. Afterwards, as with the implementation of the power iteration, it is checked whether the algorithm has converged and a warning is issued if necessary.

The eigenvalues are extracted from the diagonal of the matrix `a_k` using [`np.diag`](https://numpy.org/doc/stable/reference/generated/numpy.diag.html). At the end, the order of the indices of the eigenvalues is computed with [`np.argsort`](https://numpy.org/doc/stable/reference/generated/numpy.argsort.html), and the eigenvalues and eigenvectors are sorted accordingly before being returned.

As with the power iteration, we test the implementation with a symmetric random matrix:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:test_qr_algorithm}}
```
and compare our results with those of `np.linalg.eigh`:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:verify_results}}
```
The pairwise projection of our calculated eigenvectors onto the reference eigenvectors from `np.linalg.eigh` can be performed elegantly with the matrix multiplication `eigvecs.T @ eigvecs_ref`. If the eigenvalues were sorted identically, it should be a diagonal matrix with diagonal elements close to $\pm 1$.
$$
  \begin{pmatrix}
    \text{---}\ v_1^\dag \text{---} \\
    \text{---}\ v_2^\dag \text{---} \\
    \vdots \\
    \text{---}\ v_n^\dag \text{---} \\
  \end{pmatrix}
  \begin{pmatrix}
    \vert & \vert & & \vert \\
    v_{ref,1} & v_{ref,2} & \cdots & v_{ref,n} \\
    \vert & \vert & & \vert
  \end{pmatrix}
  \approx
  \begin{pmatrix}
    \pm 1 & 0 & \cdots & 0 \\
    0 & \pm 1 & \cdots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \cdots & \pm 1
  \end{pmatrix}
  \,.
$$

---

**Self-Study Questions**

1. While the definition of an eigenvector by $\bm{A} \vec{v}_i = \lambda_i \vec{v}_i$ is intuitive, it might not be obvious why we can decompose a normal matrix as $\bm{A} = \bm{Q} \bm{\Lambda} \bm{Q}^\dag$. Try to understand this by inserting this matrix multiplication into $\bm{A} \vec{v}_i = \lambda_i \vec{v}_i$ and use that fact that eigenvectors are orthogonal.

2. The power iteration algorithm uses the Rayleigh quotient to estimate the eigenvalue. Show that if $\vec{x}$ is the exact (and normalized) eigenvector $\vec{v}_i$, the Rayleigh quotient $\lambda = (\vec{x}^\dag \bm{A} \vec{x}) / (\vec{x}^\dag \vec{x})$ simplifies to the corresponding eigenvalue $\lambda_i$.

3. Each step of the QR algorithm computes $\bm{A}_{k+1} = \bm{Q}_k^\dag \bm{A}_k \bm{Q}_k$. Based on this relationship, explain why the final matrix $\bm{A}_{k \to \infty}$ must have the same eigenvalues as the initial matrix $\bm{A}_0$.

**Challenge Questions**

1. Assume a diagonalisable matrix $\bm{A}$ has a single large eigenvalue $\lambda_1$ and all other eigenvalues are much smaller. What does this imply for the number of eigenvectors needed to approximate the matrix according to Eq. {{eqref: eq:eigenvalue_decomposition_normal}}? Think about real-world cases where this might be useful.

2. The power iteration algorithm finds the eigenvector corresponding to the eigenvalue with the largest magnitude. What would happen during the iteration if the initial vector $\vec{x}^0$ happens to be exactly orthogonal to this dominant eigenvector? (Assume no numerical precision errors).
