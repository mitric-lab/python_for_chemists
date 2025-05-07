## Numerical Optimisation

Numerical optimisation is a powerful tool that allows us to find the
minimum or maximum of complex functions, such as the loss functions of
least squares (Eq. {{eqref: eq:least_squares_loss}}), but also functions
in other contexts, such as the energy of a molecule. Since maximising a
function $f$ is equivalent to minimising $-f$, we will only speak of
minimisation in the following. The ability to tackle optimisation
problems for which no closed-form solutions exist, or where the
evaluation of the analytical expression is too expensive, significantly
expands our toolkit in data analysis. In particular, it allows us to
find solutions for models characterised by nonlinearities, high
dimensionalities, or unusual data distributions.

### Theoretical Foundations

A particularly accessible and fundamental approach to numerical
optimisation is the
[gradient descent](https://en.wikipedia.org/wiki/Gradient_descent)
method. This is an *iterative method* that starts from a given initial point
$x^0$ and follows the direction of steepest descent of the function at
each step. Mathematically, this is formulated as
$$
  x^{k+1} = x^k - \alpha \nabla f(x^k)\,,
  {{numeq}}{eq:gradient_descent}
$$
where $x^k$ is the estimate of the minimum at step $k$. The superscript
$k$ has nothing to do with exponentiation, but is merely a notation
to avoid confusion with the $i$-th component of the vector, which is
subscripted here. The gradient of the *objective function* $f$ at $x^k$
can then be noted as
$\nabla f(x^k) = \left(\frac{\partial f}{\partial x^k_1}, \ldots, \frac{\partial f}{\partial x^k_n}\right)^{\intercal}$.

The proportionality constant $\alpha$ is called the *step size* or *learning
rate*. The procedure is repeated until one or more stopping
conditions are met. Typical stopping conditions for iterative
optimisation methods are:
- The change in the function value is smaller than a threshold
- The change in the estimate is smaller than a threshold
- A maximum number of iterations is reached
- The norm of the gradient is smaller than a threshold.

The step size $\alpha$ is the only and at the same time an important
parameter of gradient descent, as it influences the convergence speed
and stability of the method. A too small value for $\alpha$ can lead
to the method converging very slowly, while a too large value can
lead to divergence.

To use gradient descent, we need access to the gradient of the
objective function. If this is not available analytically, a numerical
approximation must be used. A simple approach is the method of
[finite differences](https://en.wikipedia.org/wiki/Finite_difference).
Here, the tangent of the partial derivative is replaced by the secant,
which leads to an approximation of the form
$$
  \frac{\partial f}{\partial x^k_i} \approx 
    \frac{f(x^k + h \hat{e}_i) - f(x^k - h \hat{e}_i)}{2h}\,,
  {{numeq}}{eq:finite_difference_symmetric}
$$
where $\hat{e}_i$ is the $i$-th unit vector and $h$ is a small value
that determines the step size of the approximation. More precisely,
this equation {{eqref: eq:finite_difference_symmetric}} represents the
*central finite difference of order 2*. There are also *one-sided*
approximations and higher-order approximations, which we will not
discuss further here.

### Implementation

#### Finite Differences
The first step is to implement the finite difference. Since we will
need to calculate derivatives multiple times, it is useful to
implement a *function*. Functions in the programming context are
similar to mathematical functions that transform an input into an
output. However, they can do much more.

We first import `numpy`:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:import_numpy}}
```
Then we define the function `finite_difference`, which calculates the
gradient of a function `func` at the point `x0`.
```python
{{#include ../codes/01-regression/numerical_optimisation.py:finite_difference}}
```
The first block of the function is the *signature*, which indicates
which arguments the function expects. In strongly typed languages, 
function signatures also define the types of the arguments and the 
return type. Because this is not the case in Python, the naming
for the arguments becomes important.

The signature of the function `finite_difference` indicates that it
expects the argument `func`, which is probably a function from its name. 
The next argument is `x0`, which is named inline with the mathematical
notation. The next two arguments are `h` and `args`, which are optional. 
This can be seen in the signature, where the default values are defined
as `h=1e-5` and `args=()`. The default value for `h` is a small
positive number, indicating that could be the step size, and the default
value for `args` is an empty tuple, indicating that it could be used
to pass additional arguments to the function `func`. 

In the implementation of the function, we first determine the dimension
of the point `x0` and store it in the integer `n`. The gradient of a
real-valued function with $n$ variables is a vector of length $n$.
Therefore, the variable `grad` is initialised as a numpy array of
length `n` filled with zeros using the function 
[`np.zeros`](https://numpy.org/doc/stable/reference/generated/numpy.zeros.html).

Because we need to apply Eq. {{eqref: eq:finite_difference_symmetric}} to all
$n$ components of the point $x^k$, it is useful to use a *loop*.
Here we use a 
[*for*-loop](https://docs.python.org/3/tutorial/controlflow.html#for-statements)
that iterates over all indices `i` from `0` to `n - 1`. Note that
Python indices start at `0` and the built-in function 
[`range`](https://docs.python.org/3/library/functions.html#func-range)
excludes the end value `n`.

In each iteration of the loop, we first define the unit vector 
$\hat{e}_i$ by first creating a zero array of length `n` with
`np.zeros(n)` and then setting the $i$-th component to `1`. 
Then we can calculate the $i$-th entry of the gradient array `grad`
according to Eq. {{eqref: eq:finite_difference_symmetric}}. The only
difference between our code and this equation is that we pass the
additional argument `args` to the function `func`, with a star `*`
in front of the argument `args`. The use of `*` acts as a unary
operator at this point, which unpacks an object into its components
(see [*unpacking*](https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists)).
This means that the function `func` accepts the further arguments
in the tuple `args` **individually** after the first argument.

After the last iteration, we return the gradient in the variable `grad`.

#### Objective Function

Next, we implement the objective function whose value we want to
minimise. According to Eq. {{eqref: eq:least_squares_loss_linear}}
we can define the function `objective_function` as follows:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:objective_function}}
```

The first argument in this function is `beta`, which are the parameters
to be optimised. The second argument is `args`, which contains the
additional arguments for the objective function. Here we see the use
of `*` in front of the argument `args` again, which makes the function
`objective_function` accept any number of additional arguments after
`beta`.

In the implementation of the function `objective_function`, we first
define the arrays `concenctrations` and `absorbances` from the argument
`args` and then return the value of the loss function of least
squares according to Eq. {{eqref: eq:least_squares_loss_linear}}.
There are no big differences between the mathematical formulation and the
programmatic implementation here.

With the help of the function `finite_difference`, we can implement
the function `objective_function_gradient`, which calculates the
gradient of the objective function:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:objective_function_gradient}}
```
This function takes the same arguments as the function
`objective_function`, but returns the gradient of the objective
function instead of its value. In its implementation, we again
define the arrays `concenctrations` and `absorbances` from the
argument `args`, and then call the function `finite_difference` to calculate
the gradient of the objective function. The result is stored in the
variable `grad` and returned.

#### Gradient Descent

Now we can implement the gradient descent algorithm. According to
Eq. {{eqref: eq:gradient_descent}} we need the gradient of the objective
function `func_grad`, the starting point `x0`, and the step size
`alpha` to start the algorithm. As a stopping condition, we use a
combination of the maximum number of iterations and the norm of
the gradient. This leads to the additional arguments `max_iter` and
`max_norm`. We also need the argument `args`, which contains the
additional arguments for the function `func_grad`.
```python
{{#include ../codes/01-regression/numerical_optimisation.py:gradient_descent}}
```

First, we explicitly copy the array `x0` with the function `np.copy`
and store the copy in the variable `x`. This is necessary because
Python does not copy arrays by default, but only stores references
to them. This means that a change in `x0` will cause a change in `x`,
and *vice versa*. 

Then we use a `for`-loop that iterates the variable `niter` from `0` to
`max_iter - 1`. In each iteration, we calculate the gradient `grad`
by calling the function `func_grad` with the arguments `x` and `args`.
Then we use Eq. {{eqref: eq:gradient_descent}} to update the variable
`x`. After that, we calculate the norm of the gradient with the
function 
[`np.linalg.norm`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html)
and check if it is smaller than `max_norm`. If so, we break the
loop with the `break` command.

After the last iteration, we check if the variable `niter` has reached
`max_iter`. If so, it means that the method may not have converged,
and we issue a warning. At the end, we return the optimum `x` and
the number of iterations `niter` that were actually needed.

### Application

Now we can apply the implemented algorithm to the data from Chapter
[1.2](02-linear_regression.md). First, we define the arrays
`concentrations` and `absorbances` again:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:data_array}}
```
Then we define the initial point `x0`, here `beta_guess`, and call the
function `gradient_descent`:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:gradient_descent_call}}
```

The optimal parameters `beta0` and `beta1` are identical to those of the
analytical solution. On the author's computer, 34683 iterations were
needed to satisfy the stopping condition. The exact number of
iterations may vary slightly depending on the hardware. If you
choose the step size `alpha` slightly larger, fewer iterations are
needed. However, if `alpha` is too large, the method diverges.

```admonish tip title="Tip"
Try changing the step size `alpha` and observe how the number of
iterations changes.
```

Since optimisation is a very general problem, there are many implementations
of various algorithms in libraries such as 
[`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html).
We want to use the function 
[`scipy.optimize.minimize`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)
to find the optimal parameters $\beta$:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:scipy_minimize}}
```
When calling the function `minimize`, we only need to specify the
objective function `objective_function` and the initial point `beta_guess`.
The `minimize` function takes care of the numerical gradient
calculation itself.

~~~admonish note title="Note on the `minimize` function"
The `minimize` function also accepts the argument `jac`
([Jacobian matrix](https://en.wikipedia.org/wiki/Jacobian_matrix)), i.e. a
function that calculates the gradient of the objective function. If you
have the analytical and easily computable gradient available, you can
pass it as the `jac` argument, which can speed up the optimisation
process.
~~~

With the argument `method='CG'`, we have selected the 
[<i>nonlinear **C**onjugate **G**radient method</i>](https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method).
This method is an improvement of the gradient method and reaches the
minimum in only 2 iterations on the author's computer.

In `minimize`, there are a number of other minimisation methods that
can be used. An overview can be found in the
[documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)
of this function. Two important methods among them are:
- `method='Nelder-Mead'`: 
    The [Nelder-Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method)
    is a heuristic method that does not require the calculation of the
    gradient. It is therefore particularly useful when the calculation
    of the gradient is very expensive or the gradient is noisy. It
    is especially suitable for regressions with experimental data.
- `method='BFGS'`:
    The [**B**royden-**F**letcher-**G**oldfarb-**S**hanno method](https://en.wikipedia.org/wiki/BFGS_method)
    is a method that uses the gradient to approximate the Hesse matrix 
    (Hessian) of the objective function. The Hessian contains the
    second derivatives of the function with respect to its parameters.
    The method therefore shows very fast convergence near the minimum.
    In practice, it requires fewer iterations than other optimisation
    methods and is therefore often used.

As you will see in the exercise, the linear regression with the least
squares method has a closed-form solution that directly calculates the optimal
parameters. So why should we bother with numerical optimisation?
In the context of regression, numerical optimisation allows us to
use more complicated models that do not have an analytical solution,
and also to use more sophisticated loss functions, such as the
least absolute deviations (see Eq. {{eqref: eq:least_absolute_deviations_opt}}).
In addition, we can introduce additional control over the parameters
(**regularisation**), which can improve the general performance of the
model.

~~~admonish note title="Function `scipy.optimize.curve_fit`"
There is also the function 
[`scipy.optimize.curve_fit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html),
which performs a (nonlinear) regression of the data directly.
However, it is not as flexible as the general method with the function
`minimize`, as it only has a few optimisation methods and sets the
objective function as the least squares loss function.
It can be applied to simple regressions and usually requires less code.
~~~

