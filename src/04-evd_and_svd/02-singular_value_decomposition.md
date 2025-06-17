## Singular Value Decomposition

### Definition

In the spirit of the eigenvalue decomposition for square matrices, 
we define the following decomposition for an arbitrary matrix 
$\bm{A} \in \C{m}{n}$:
$$
  \bm{A} = \bm{U} \bm{\Sigma} \bm{V}^\dag
$$
where $\bm{U} \in \C{m}{m}$ and $\bm{V} \in \C{n}{n}$ are unitary matrices,
and $\bm{\Sigma} \in \R{m}{n}$ is a diagonal matrix with non-negative 
diagonal entries $\sigma_1 \geq \sigma_2 \geq \ldots \geq \sigma_p \geq 0$, 
where $p = \min(m,n)$. This decomposition is known as the 
*singular value decomposition* (SVD).

But what does it mean for a rectangular matrix $\bm{\Sigma}$ to be diagonal?
Strictly speaking, this is not possible, but we can bring the matrix into a form
that resembles diagonal form, namely
$$
  \bm{\Sigma} = \left(\begin{array}{c|c}
    \bm{D} & \bm{0} \\ \hline
    \bm{0} & \bm{0}
  \end{array}\right)\,.
$$
Thus, $\bm{\Sigma}$ contains a $p \times p$ diagonal matrix $\bm{D}$ and
possibly zero matrices to fulfill the dimensions.

The diagonal entries $\sigma_i$ are called *singular values*, and the columns of
$\bm{U}$ and $\bm{V}$ are referred to as *left* and *right singular vectors*,
respectively. In analogy to
Equation {{eq: eq:eigenvalue_decomposition_normal}}, the 
singular value decomposition can be written as
$$
  \begin{align*}
    \bm{A} &= \bm{U} \bm{\Sigma} \bm{V}^\dag \\
    &= \begin{pmatrix}
      \vert & \vert &  & \vert \\
      \vec{u}_ 1 & \vec{u}_ 2 & \cdots & \vec{u}_ m \\
      \vert & \vert &  & \vert
    \end{pmatrix}
    \left(\begin{array}{cccc|ccc}
      \sigma_1 & 0 & \cdots & 0 & 0 & \cdots & 0 \\ 
      0 & \sigma_2 & \cdots & 0 & 0 & \cdots & 0 \\
      \vdots & \vdots & \ddots & 0 & \vdots & \ddots & \vdots \\
      0 & 0 & \cdots & \sigma_p & 0 & \cdots & 0 \\ \hline
      0 & 0 & \cdots & 0 & 0 & \cdots & 0 \\
      \vdots & \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\
      0 & 0 & \cdots & 0 & 0 & \cdots & 0
    \end{array}\right)
    \begin{pmatrix}
      \text{---}\ \vec{v}_ 1^\dag\ \text{---} \\
      \text{---}\ \vec{v}_ 2^\dag\ \text{---} \\
      \vdots \\
      \text{---}\ \vec{v}_ n^\dag\ \text{---}
    \end{pmatrix} \\
    &= \sigma_1 \vec{u}_ 1 \vec{v}_ 1^\dag + \sigma_2 \vec{u}_ 2 \vec{v}_ 2^\dag + \cdots + \sigma_p \vec{u}_ p \vec{v}_ p^\dag \\
    &= \sum_{i=1}^p \sigma_i \vec{u}_ i \vec{v}_ i^\dag\,.
    {{numeq}}{eq:singular_value_decomposition}
  \end{align*}
$$

One can see a clear difference to the eigenvalue decomposition:
While in the eigenvalue decomposition all eigenvectors are important, here
$(m-p)$ or $(n-p)$ singular vectors do not contribute to the decomposition.
This means that one of the matrices $\bm{U}$ or $\bm{V}$ can be reduced without
changing the decomposition. Assuming without loss of generality that $m \geq n$,
we can define $\hat{\bm{U}} \in \C{m}{n}$ and $\hat{\bm{\Sigma}} \in \R{n}{n}$
by removing the last $(m-n)$ columns from $\bm{U}$ and the last 
$(m-n)$ rows from $\bm{\Sigma}$. We also define $\hat{\bm{V}} = \bm{V}$.
Then the decomposition
$$
  \bm{A} = \hat{\bm{U}} \hat{\bm{\Sigma}} \hat{\bm{V}}^\dag\,,
$$
which is known as the *economic singular value decomposition*, holds. 
If one of the dimensions $n, m$ is much larger than the other, the 
economic SVD is advantageous because it requires less memory.

### Properties

We now consider the matrix products $\bm{A}^\dag \bm{A}$ and
$\bm{A} \bm{A}^\dag$ using the SVD. We have
$$
  \begin{align*}
    \bm{A}^\dag \bm{A} &= (\bm{U} \bm{\Sigma} \bm{V}^\dag)^\dag \bm{U} \bm{\Sigma} \bm{V}^\dag 
    = (\bm{V}^\dag)^\dag \bm{\Sigma}^\dag \bm{U}^\dag \bm{U} \bm{\Sigma} \bm{V}^\dag
    = \bm{V} \bm{\Sigma}^\dag \bm{\Sigma} \bm{V}^\dag \\
    \bm{A} \bm{A}^\dag &= \bm{U} \bm{\Sigma} \bm{V}^\dag (\bm{U} \bm{\Sigma} \bm{V}^\dag)^\dag
    = \bm{U} \bm{\Sigma} \bm{V}^\dag (\bm{V}^\dag)^\dag \bm{\Sigma}^\dag \bm{U}^\dag
    = \bm{U} \bm{\Sigma} \bm{\Sigma}^\dag \bm{U}^\dag
  \end{align*}\,.
  {{numeq}}{eq:svd_and_evd}
$$

Because the non-zero elements in $\bm{\Sigma}$ are only on the diagonal,
the products $\bm{\Sigma}^\dag \bm{\Sigma}$ and $\bm{\Sigma} \bm{\Sigma}^\dag$ 
are *diagonal matrices* with the diagonal elements $\sigma_i^2$.

The two equations {{eqref: eq:svd_and_evd}} are similar to the 
eigenvalue decomposition for normal matrices 
$\bm{B} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$. 
In fact, the products $\bm{A}^\dag \bm{A}$ and $\bm{A} \bm{A}^\dag$ are even 
hermitian, because
$$
  \begin{align*}
    (\bm{A}^\dag \bm{A})^\dag &= \bm{A}^\dag (\bm{A}^\dag)^\dag = \bm{A}^\dag \bm{A} \\
    (\bm{A} \bm{A}^\dag)^\dag &= (\bm{A}^\dag)^\dag \bm{A}^\dag = \bm{A} \bm{A}^\dag
  \end{align*}\,,
$$
which also makes them normal.

For this reason, the singular values correspond to the square roots
of the eigenvalues of $\bm{A}^\dag \bm{A}$ and $\bm{A} \bm{A}^\dag$. The
left singular vectors are therefore the eigenvectors of $\bm{A} \bm{A}^\dag$,
while the right singular vectors represent the eigenvectors of
$\bm{A}^\dag \bm{A}$. We have thus established a connection between the SVD
and the EVD.


### Use with NumPy

We shall now "implement" the SVD using Eq. {{eqref: eq:svd_and_evd}}.

As for the eigenvalue decomposition, we start by creating a random matrix
after importing NumPy:
```python
{{#include ../codes/04-evd_and_svd/numpy_svd.py:imports}}
```
```python
{{#include ../codes/04-evd_and_svd/numpy_svd.py:random_matrix}}
```
This time, a proper rectangular matrix is created, and no symmetrisation
is performed (which is not even possible for proper rectangular matrices).
We used
[`np.min`](https://numpy.org/doc/stable/reference/generated/numpy.min.html)
to determine the smaller dimension of the matrix, which is stored in the variable
`p`. The SVD is then computed using the double eigenvalue decomposition:
```python
{{#include ../codes/04-evd_and_svd/numpy_svd.py:double_evd}}
```
Since the eigenvalues coming from `np.linalg.eigh` are sorted in ascending order,
and at most `p` eigenvalues are non-zero, the last `p` eigenvalues are
the squared singular values, while the other eigenvalues of the larger matrix
(in this case, $\bm{A} \bm{A}^\dag$) are zero. 
Therefore, we checked for equality of the last `p` eigenvalues coming from
the EVD of $\bm{A} \bm{A}^\dag$ and $\bm{A}^\dag \bm{A}$. The singular
values were then obtained by taking the square root of these eigenvalues.
Since singular values are conventionally sorted in descending order,
we used the slice `[::-1]` to reverse the order of the singular values.
Because of the slicing and reordering of the singular values, the matrices
containing the singular vectors must be adjusted accordingly.

For such a versatile and widely used mathematical tool, there is no surprise 
that the SVD is available in NumPy via the function
[`numpy.linalg.svd`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html).
As an argument, one passes the matrix $\bm{A}$ for which the SVD
should be computed. The return values are the matrices $\bm{U}$,
$\bm{\Sigma}$, and $\bm{V}^\dag$. With the optional argument
`full_matrices=False`, the economic SVD is computed. The default option is
`full_matrices=True`. 

We shall compute the SVD of our random matrix $\bm{A}$ using NumPy:
```python
{{#include ../codes/04-evd_and_svd/numpy_svd.py:svd}}
```

Now, the SVD results obtained *via* the double eigenvalue decomposition
and `numpy.linalg.svd` can be compared:
```python
{{#include ../codes/04-evd_and_svd/numpy_svd.py:comparison_decomposition}}
```
Because the phase (sign) of the eigenvectors is arbitrary, we used the same
trick as in the EVD to compare the singular vectors. 

We now verify that the results of the SVD can be used to reconstruct the
original matrix $\bm{A}$:
```python
{{#include ../codes/04-evd_and_svd/numpy_svd.py:comparison_reconstruction}}
```
While the matrices from `numpy.linalg.svd` managed to reconstruct the original
matrix $\bm{A}$, we get a nasty surprise with the matrices from the
double eigenvalue decomposition. The reconstructed matrix is not equal to
the original matrix $\bm{A}$, not even up to a phase factor.
This happens due to the so-called 

```admonish warning title="Phase ambiguity"
Although Equation {{eqref: eq:svd_and_evd}} is mathematically correct, 
one should not compute the SVD of $\bm{A}$ using the EVD of $\bm{A}^\dag \bm{A}$ 
and $\bm{A} \bm{A}^\dag$.

The eigenvectors of a matrix are only unique up to a (complex) factor.
Even if one normalises the eigenvectors, the phase remains undetermined.
For example, if one applies a minus sign to the first left singular vector
$\vec{u}_1$, the equations {{eqref: eq:svd_and_evd}} remain valid.
However, the SVD in Equation {{eqref: eq:singular_value_decomposition}} is no longer
correct, as the first term $\sigma_1 \vec{u}_1 \vec{v}_1^\dag$ receives a
minus sign. To keep the decomposition valid, the sign of $\vec{v}_1$ would
need to be adjusted. However, since the two EVDs are independent of each other,
such adjustments are not possible.
For this reason, one should not compute the SVD from the double EVD,
unless only the singular values and not the singular vectors are needed.
For the following considerations, it is sufficient to know that
there are special algorithms that compute the left and
right singular vectors simultaneously, which avoids this problem.
```

The SVD is a powerful tool in linear algebra and has numerous applications.
In the following chapters, we will explore some of these applications,
including low-rank approximations, dimensionality reduction, and more.

