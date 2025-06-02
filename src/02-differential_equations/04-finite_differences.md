## Finite Difference Method

The
[Finite Difference Method](https://en.wikipedia.org/wiki/Finite_difference_method)
is another class of numerical methods, besides Runge-Kutta methods,
that can be used to solve ordinary differential equations.
Here too, the domain of the function $y(x)$ is divided into discrete points,
however, the focus of this method is not on calculating
function values, but on approximating the derivative operators.

### Theoretical Foundations

In the implementation of the gradient method, we already
used a formula for the finite difference in Eq. {{eqref: eq:finite_difference_symmetric}}. This formula is a special case of the
general finite difference and, at that time, more or less fell out of the sky.
Here, we strive to find a general recipe for deriving such
formulas.

#### Finite Difference Approximations

Let's assume we want to (approximately) calculate the derivative of a function $y(x)$
at point $x$ from the function values in the vicinity of $x$, such as
$y(x-h)$, $y(x)$, and $y(x+h)$. Accordingly,
$$
  y'(x) \approx c_{-1} y(x-h) + c_0 y(x) + c_1 y(x+h)
$$
must hold, where we want to determine the coefficients $c_{-1}$, $c_0$, and $c_1$.
We expand the right side into a Taylor series around $x$:
$$
  \begin{align}
  y'(x) \approx \
      & c_{-1} y(x) &&+ c_{-1} y'(x) (-h) &&+ \frac{1}{2} c_{-1} y''(x) (-h)^2 + \ldots \\
    +\ & c_0 y(x) \\
    +\ & c_1 y(x) &&+ c_1 y'(x) h &&+ \frac{1}{2} c_1 y''(x) h^2 + \ldots
  \end{align}
$$
The idea now is that if the right side is to be equal to the left side (i.e., $y'(x)$),
then the coefficients in front of the function values, or their derivatives, must also be equal
on both sides. Since we have three unknown coefficients here, we need three equations,
which is why the Taylor series was expanded up to the second order. Comparing the coefficients
on both sides, we obtain the following system of linear equations:
$$
  \begin{align}
    c_{-1} + c_0 + c_1 &= 0 \\
    -c_{-1} h + c_1 h &= 1 \\
    \frac{1}{2} c_{-1} h^2 + \frac{1}{2} c_1 h^2 &= 0
  \end{align}
$$
One could solve this system of equations using Gaussian elimination,
but with only three equations and unknowns, it can also be done by simple rearrangement and substitution.
As a result, we obtain the coefficients $c_{-1} = -1/(2h)$,
$c_0 = 0$, and $c_1 = 1/(2h)$. Substituting these into the formula above, we get
$$
  y'(x) \approx \frac{y(x+h) - y(x-h)}{2h}\,,
  {{numeq}}{eq:finite_difference_symmetric_second_order}
$$
which is identical to Eq. {{eqref: eq:finite_difference_symmetric}}.

This finite difference formula is called the second-order central difference
because it calculates the derivative of the function $y(x)$ at
point $x$ symmetrically from the function values $y(x-h)$ and $y(x+h)$
and uses the Taylor series expansion up to the second order.

Using the grid points $x$ and $x+h$, or $x$ and $x-h$, one can similarly
derive the forward differences
$$
  y'(x) \approx \frac{y(x+h) - y(x)}{h}
  {{numeq}}{eq:finite_difference_forward_first_order}
$$
and the backward differences
$$
  y'(x) \approx \frac{y(x) - y(x-h)}{h}
  {{numeq}}{eq:finite_difference_backward_first_order}
$$
However, these two formulas only consider the first order.

According to this procedure, we can also approximate derivatives of higher orders
using function values at arbitrary points on the grid up to an arbitrary order.

#### Matrix Representation of the Differential Operator

If we discretize a function $y(x)$ on a grid $x_i$ with $i=1,\ldots,N$,
we can (approximately) represent the function as a vector
$\vec{y} = (y_1, \ldots, y_N)^\intercal$ with $y_i = y(x_i)$.
If we want to calculate the function value at an arbitrary grid point $x_i$,
we can represent this as the scalar product of the vector $\vec{y}$ with the
basis vector $\hat{e}^i$, where $\hat{e}^i$ is a vector with
$1$ at position $i$ and $0$ at all other positions. This yields
$$
  \langle \hat{e}^i, \vec{y} \rangle
    = \sum_{j=1}^N \hat{e}^i_ j y_ j
    = \sum_{j=1}^N \delta_{ij} y_ j
    = y_i
$$
with the Kronecker delta $\delta_{ij}$. If we write out the
vector components, the scalar product is given by
$$
  \begin{pmatrix}
    0_1 & \cdots & 0_{i-1} & 1_i & 0_{i+1} & \cdots & 0_N
  \end{pmatrix}
  \cdot
  \begin{pmatrix}
    y_1 \\ \vdots \\ y_{i-1} \\ y_i \\ y_{i+1} \\ \vdots \\ y_N
  \end{pmatrix}
  = y_i
$$
where the indices of the entries in $e^i$ have been written explicitly.

In this way, we can represent $\left(y(x_{i+1}) - y(x_{i-1})\right)$ as
$\langle \hat{e}^{i+1} - \hat{e}^{i-1}, \vec{y} \rangle$,
with the auxiliary vector
$$
  \hat{e}^{i+1} - \hat{e}^{i-1} =
  \begin{pmatrix}
    0_1 & \cdots & 0_{i-2} & -1_{i-1} & 0_i & 1_{i+1} & 0_{i+2} & \cdots & 0_N
  \end{pmatrix}\,.
$$

The derivative at point $x_i$ can thus be written according to Eq. {{eqref: eq:finite_difference_symmetric_second_order}}
as $\frac{1}{2h} \langle \hat{e}^{i+1} - \hat{e}^{i-1}, \vec{y} \rangle$.
Since we want to calculate the derivative at all grid points,
we continue this pattern and obtain a matrix equation
$$
  \vec{y}' = \bm{D} \vec{y}
$$
with
$$
  \vec{y}' = (y'_1, \ldots, y'_N)^\intercal, \quad
  \vec{y} = (y_1, \ldots, y_N)^\intercal
$$
and
$$
  \bm{D} = \frac{1}{2h}
  \begin{pmatrix}
    0 & 1 & 0 & 0 & \cdots & 0 & 0 \\
    -1 & 0 & 1 & 0 & \cdots & 0 & 0 \\
    0 & -1 & 0 & 1 & \cdots & 0 & 0 \\
    \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
    0 & 0 & 0 & 0 & \cdots & -1 & 0 \\
  \end{pmatrix}\,.
  {{numeq}}{eq:finite_difference_symmetric_second_order_matrix}
$$

This *representation of the differential operator* with second-order symmetric
finite differences $\bm{D}$ is therefore a
[*tridiagonal matrix*](https://en.wikipedia.org/wiki/Tridiagonal_matrix)
with $0$ on the main diagonal and $\pm 1$ on the sub-diagonals.

What about the representation of second or higher derivatives?
Of course, one could consider a higher-order differential equation as a
system of first-order differential equations and adapt the
finite difference method for vector-valued functions. But
here we can construct the corresponding matrix directly.

To do this, we first show how not to do it: The form of the operator
for the second derivative $\frac{\du^2}{\du x^2}$ at first glance resembles the
"square" of the derivative operator $\frac{\du}{\du x}$. If we
were to actually square the matrix $\bm{D}$, we would get something that
doesn't look wrong at first glance:
$$
  \bm{D}^2 y = \bm{D} \bm{D} y = \bm{D} y' = y''\,.
$$

This equation would be correct if $\bm{D} y$ were truly the derivative
of $y$. However, this is not the case here, because discretization gives us
finite resolution. Applying the operator $\bm{D}$ to $\bm{D} y$
would therefore amplify the error of the first derivative.

```admonish warning title="Warning"
In general, do **not** use the matrix $\bm{D}^n$ as the
representation of the operator of the $n$-th derivative!
```

Now we show how to do it correctly. Here too, we first consider the second
derivative as the derivative of the first derivative:
$$
  \begin{align}
    y''(x) = \frac{\du}{\du x} y'(x)
      &\approx \frac{1}{h} \left( y'(x) - y'(x-h) \right) \\
      &= \frac{1}{h} \left( \frac{\du}{\du x} y(x) - \frac{\du}{\du x} y(x-h) \right) \\
      &= \frac{1}{h} \left( \frac{y(x+h) - y(x)}{h} - \frac{y(x) - y(x-h)}{h} \right) \\
      &= \frac{1}{h^2} \left( y(x+h) - 2y(x) + y(x-h) \right)\,,
  \end{align}
$$
where we used backward differences {{eqref: eq:finite_difference_backward_first_order}} for $\frac{\du}{\du x} y'$
and forward differences {{eqref: eq:finite_difference_forward_first_order}} for $\frac{\du}{\du x} y$.
This derivation can easily be generalized to the $n$-th derivative.

With this formula, we can construct the matrix representation of the operator for the
second derivative, denoted here as $\bm{D}^{(2)}$:
$$
  \bm{D}^{(2)} = \frac{1}{h^2}
  \begin{pmatrix}
    -2 & 1 & 0 & 0 & \cdots & 0 & 0 \\
    1 & -2 & 1 & 0 & \cdots & 0 & 0 \\
    0 & 1 & -2 & 1 & \cdots & 0 & 0 \\
    \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
    0 & 0 & 0 & 0 & \cdots & 1 & -2 \\
  \end{pmatrix}\,.
  {{numeq}}{eq:second_finite_difference_symmetric_second_order_matrix}
$$
This matrix is again a tridiagonal matrix, but with $-2$ on the
main diagonal and $1$ on the sub-diagonals.

#### Finite Difference Procedure

We now have the matrix representation of the differential operator for the
first and second derivatives. However, in a differential equation,
besides the sought function $y(x)$ and
its derivatives, functions of $x$ also appear. We therefore need a
representation for such functions. A special feature here is that there are
two representation variants for this, which must be chosen according to the context.

If a function of $x$ acts as an inhomogeneity $u(x)$
(cf. Eq. {{eqref: eq:lode_general}}), i.e., as a standalone term without
$y(x)$ or its derivatives, it is discretized exactly as $y(x)$
and represented as a vector $\vec{u} = (u_1, \ldots, u_N)^\intercal$ with $u_i = u(x_i)$.
However, if the function is multiplied by $y(x)$ or its derivatives,
it serves as a coefficient function $a_i(x)$, and must be represented as a
matrix. Since the multiplication of two
functions occurs element-wise, the representation matrix $\bm{A}_ i$
takes a diagonal form, where the diagonal elements are the discretized
function values of $a_i(x)$, i.e.,
$\bm{A}_ i = \text{diag}(a_i(x_1), \ldots, a_i(x_N))$.

Now we know the (approximate) representation of all elements of an
ODE and can convert any **linear** differential equation
(cf. Eq. {{eqref: eq:lode_general}}) into a matrix equation:
$$
  \underbrace{
    \bm{A}_n \bm{D}^{(n)} \vec{y} + \bm{A}_{n-1} \bm{D}^{(n-1)} \vec{y} + \ldots + \bm{A}_1 \bm{D} \vec{y} + \bm{A}_0 \vec{y}
  }_ {\bm{L} \vec{y}} = \bm{B}_ 0 \vec{u}\,
$$
Here $\bm{B}_ 0 = (b_0(x_1), \ldots, b_0(x_N))^\intercal$ serves as the
matrix representation of the coefficient function $b_0(x)$, and in the
above equation, we have combined all linear operators into one operator
$$
  \bm{L} = \bm{A}_n \bm{D}^{(n)} + \bm{A}_{n-1} \bm{D}^{(n-1)} + \ldots + \bm{A}_1 \bm{D} + \bm{A}_0
$$
The solution of the linear system of equations
$\bm{L} \vec{y} = \bm{B}_ 0 \vec{u}$ then yields the discretized function
$\vec{y}$.

```admonish note title="Note on nonlinear ODEs"
Actually, nonlinear ODEs can also be expressed with
finite difference operators. The resulting system of equations,
however, is a nonlinear system of equations, which cannot be solved directly with
methods of linear algebra.
```

There is a crucial difference between solving a linear system of equations and solving a
differential equation:
The initial value. While an initial or
boundary condition is needed to obtain a particular solution for ODEs,
there is no such condition for a linear system of equations.
How then should we incorporate the initial conditions of the ODE within the framework of the discretized
version?

Let's first consider the matrix representation of the differential operator $\bm{D}$
in Eq. {{eqref: eq:finite_difference_symmetric_second_order_matrix}},
particularly the first row $(0, 1, 0, \ldots, 0)$. Multiplied by $\vec{y}$,
this row states
$$
  y'(x_1) = \frac{1}{2h} y(x_2).
$$
However, according to Eq. {{eqref: eq:finite_difference_symmetric_second_order}}, it should be
$$
  y'(x_1) = \frac{1}{2h} \left( y(x_2) - y(x_0) \right)
$$
For both equations to agree, $y(x_0) = 0$ must hold.
The last row of $\bm{D}$ in turn yields the condition $y(x_{N+1}) = 0$.

Thus, by constructing the finite difference operators as matrices, the
boundary conditions are implicitly defined. The finite difference method is
therefore more suitable for boundary value problems than for initial value problems.

```admonish note title="Note on other boundary conditions"
Besides the
[*Dirichlet boundary condition*](https://en.wikipedia.org/wiki/Dirichlet_boundary_condition)
in our case, meaning that the function values outside the grid
must be zero, other boundary conditions, such as
[*periodic boundary conditions*](https://en.wikipedia.org/wiki/Periodic_boundary_conditions)
or
[*Neumann boundary conditions*](https://en.wikipedia.org/wiki/Neumann_boundary_condition),
can be incorporated with the finite difference method. However, we will
not discuss these further here.
```

### Implementation

We now want to implement the finite difference method using the example of the
Schrödinger equation for the harmonic oscillator.

The Schrödinger equation in atomic units is
$$
  -\frac{1}{2} \frac{\du^2}{\du x^2} \psi(x) + \frac{1}{2} k x^2 \psi(x) = E \psi(x)
$$
and is a linear second-order ODE with boundary conditions $\lim_{x \to \pm \infty} \psi(x) = 0$.
The representation matrix $\bm{D}^{(2)}$ in
Eq. {{eqref: eq:second_finite_difference_symmetric_second_order_matrix}} is therefore suitable for the implementation.

After importing the required libraries
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:imports}}
```
we can define the differential operator $\bm{D}^{(2)}$ as follows:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:generate_d2_naive}}
```
While this simple implementation is correct, it is not very
efficient, especially when the number of grid points `n` is large.
A more efficient implementation can be achieved using the function
[`np.diag_indices`](https://numpy.org/doc/stable/reference/generated/numpy.diag_indices.html):
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:generate_d2}}
```
After initializing a zero matrix, we used the `np.diag_indices` function
to get the indices of the main diagonal. This allows the corresponding elements
of the array `d2` to be set to $-2$.
By shifting the indices by $\pm 1$, we can address the entries of the sub-diagonals
and set the corresponding elements to $1$. Avoiding
loops significantly speeds up the calculation.

```admonish note title="Note on Python Loops"
As we mentioned in the introduction, Python is a relatively slow language compared
to compiled languages. Therefore, when efficiency is crucial, one should avoid loops and instead
use functions from the `numpy` library.
```

For applications requiring a large number of grid points,
the use of
[sparse matrices](https://en.wikipedia.org/wiki/Sparse_matrix)
can be advantageous. Since we only need a moderate number of
grid points in this example, we will carry out the implementation with a normal matrix.

```admonish info title="Note on sparse matrices" collapsible=true
For large $N$ due to finer or multi-dimensional grids,
the finite difference matrix can become very large. However, the
majority of the elements in this matrix are zero. Such matrices are called
[*sparse matrices*](https://en.wikipedia.org/wiki/Sparse_matrix)
and can be stored and processed more efficiently using special algorithms. The
finite difference matrix even has entries only on the main and
some sub-diagonals and is therefore also called a
[*band matrix*](https://en.wikipedia.org/wiki/Band_matrix),
which is a generalization of the tridiagonal matrix. Band matrices
can be implemented as sparse matrices using the function
[`scipy.sparse.diags_array`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.diags_array.html).

Then, methods from the submodule
[`scipy.sparse.linalg`](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html)
can be used to solve linear systems of equations or to calculate eigenvalues
and eigenvectors.

If one only wants to calculate the eigenvalues and eigenvectors of band matrices,
the function
[`scipy.linalg.eigh_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_banded.html)
or
[`scipy.linalg.eigh_tridiagonal`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_tridiagonal.html)
can be used, which only expects the occupied diagonals of the band matrix as an argument.
```

We then construct the Hamiltonian operator for the harmonic oscillator
as a matrix:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:build_hamiltonian}}
```
Here we first generated the $D^{(2)}$ matrix with the `generate_d2` function
and then created the potential function $\frac{1}{2} k x^2$ as a
diagonal matrix with `np.diag`. The sum
$-\frac{1}{2} D^{(2)} + \frac{1}{2} k x^2$ then yields the Hamiltonian operator,
or the linear operator $\bm{L}$ (`l_mat`).

Thus, the Schrödinger equation in discrete form is
$$
  \bm{L} \vec{\psi} = E \vec{\psi}\,.
$$
It is easy to see at this point that this is a matrix eigenvalue equation.
We therefore use the function
[`np.linalg.eigh`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html) to calculate the eigenvalues and eigenvectors of the Hamiltonian operator.

For this, we first set $k = 1$ and choose a grid from -5 to 5 with 512 points:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:define_parameters}}
```

Then, the Hamiltonian operator is generated and the eigenvalues and
eigenvectors are calculated:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:solve_matrix_equation}}
```
Since the `np.linalg.eigh` function expects a symmetric matrix, we
used `assert` and
[`np.allclose`](https://numpy.org/doc/stable/reference/generated/numpy.allclose.html)
to check if the `hamiltonian` matrix is identical to its transpose.

We now extract the energies and the corresponding wave functions from the
first 20 eigenvectors and eigenvalues.
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:solve_harmonic_oscillator}}
```
However, the wave functions should be normalized according to
$$
  \int_{-\infty}^{\infty} |\psi(x)|^2 \du x = 1
$$
which in discrete form corresponds to
$$
  \sum_{i=1}^{N} |\psi(x_i)|^2 \Delta x = 1
$$
The eigenvectors from `np.linalg.eigh` are, however, normalized according to
$$
  \sum_{i=1}^{N} |v_i|^2 = 1
$$
Therefore, we divide the eigenvectors by $\sqrt{\Delta x}$, where we
take $\Delta x$ from the first two entries of the grid,
to obtain the normalized wave functions.

Finally, we want to visualize our results. We want to plot the eigenenergies
and the wave functions side by side. Therefore, we pass the following arguments
to the `plt.subplots` function:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:plot_results_subplots}}
```
Here, the arguments `1, 2` mean that we want to create two separate diagrams in one row and
two columns. The axis objects are then stored in the list `axs`.

Then we plot the numerical as well as the analytical eigenenergies
of the harmonic oscillator in the first plot with `axs[0]`:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:plot_eigenenergies}}
```
Subsequently, we plot the harmonic potential
as well as the first 5 numerical wave functions in the second diagram (`axs[1]`):
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:plot_eigenfunctions}}
```
When plotting the potential, we used the argument `lw=2` to increase the
line width.

Finally, we format the plot and display it:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:show_plot}}
```

The diagram should look like this:
![Eigenenergies and wave functions of the harmonic oscillator](../assets/figures/02-differential_equations/fdm_harm_osc.svg)

It can be seen that the numerical eigenenergies of the first approx. 10
states agree very well with the analytical eigenenergies.
The numerically calculated wave functions of the first five states look
sensible, at least at first glance. If greater accuracy is needed for the higher states,
the grid must be chosen finer and also larger, as the
wave functions on the one hand exhibit more oscillations and on the other hand
are spatially more extended.

---

**Self-Study Questions**

1. For the first derivative matrix $\bm{D}$, show how the dot product of row $i$ with vector $\vec{y}$ gives the central difference formula for $y'(x_i)$.

2. Why doesn't $\bm{D}^2$ give the correct second-derivative matrix $\bm{D}^{(2)}$?

3. In the Schrödinger equation, why must the potential $V(x)$ be represented as a diagonal matrix in $\bm{L}\vec{\psi} = E\vec{\psi}$, rather than as a vector?

**Challenge Questions**

1. Starting from Taylor expansions, derive the coefficients for a three-point forward difference approximation of $y'(x)$ using $y(x)$, $y(x+h)$, and $y(x+2h)$. What is the order of accuracy of this approximation?

2. What do 'stability' and 'convergence' mean in the context of Finite Difference Method solutions? How could you empirically test the convergence rate of the chapter's harmonic oscillator solver with respect to the grid spacing $h$?

3. Formulate the finite difference approximation for the 2D Laplacian operator ($\nabla^2 = \frac{\partial^2}{\partial x^2} + \frac{\partial^2}{\partial y^2}$) on a Cartesian grid. Describe the structure of the resulting matrix operator if grid points are ordered lexicographically, and comment on its sparsity.