## Numerical Optimisation

In many scientific problems, we are interested in finding the global minimum or maximum of a function $f$, or the parameters $\beta^*$ that yield this minimum or maximum. For example, in the previous chapter, we wanted to find the parameters that minimise the loss function, and in many applications in quantum chemistry, we want to find the coefficients that minimise the energy of a molecule. Since maximising a function $f$ is equivalent to minimising $-f$, we will only discuss minimisation in the following. Because analytical solutions do not exist for most of these problems, we need to find ways to determine the global minimum or maximum of a function numerically.

In the vast majority of cases, finding the global minimum is impossible because the given function $f$ is highly non-convex (i.e., it has many local minima) or is too high-dimensional (i.e., it has many parameters). However, there are several methods that allow us to find a **local minimum** of a function, such as the gradient descent method.

### Theoretical Foundations

[Gradient descent](https://en.wikipedia.org/wiki/Gradient_descent) is an *iterative method* that starts from a given initial set of parameters $\beta^0$ and follows the direction of steepest descent of the function $f$ at each step. Mathematically, one step of the algorithm is formulated as
$$
  \beta^{k+1} = \beta^k - \alpha \nabla f(\beta^k)\,,
  {{numeq}}{eq:gradient_descent}
$$
where $\beta^k$ represents the estimate of the minimum at step $k$. The superscript $k$ is used to denote the iteration step and should not be confused with exponentiation. The gradient of the *objective function* $f$ with respect to the parameters at step $k$ is denoted as $\nabla f(\beta^k) = \left(\frac{\partial f(\beta^k)}{\partial \beta_1}, \ldots, \frac{\partial f(\beta^k)}{\partial \beta_n}\right)^{\intercal}$.

The proportionality constant $\alpha$ is called the *step size* or *learning rate*. The procedure is repeated until one or more stopping conditions are met. Typical stopping conditions for iterative optimisation methods are:
- The change in the function value is smaller than a threshold
- The change in the parameters is smaller than a threshold
- A maximum number of iterations is reached
- The norm of the gradient is smaller than a threshold

The step size $\alpha$ is a crucial parameter of gradient descent, as it influences both the convergence speed and stability of the method. A value that is too small can lead to very slow convergence, while a value that is too large can cause the method to diverge.

To use gradient descent, we need access to the gradient of the objective function with respect to the parameters. If this is not available analytically, we must use a numerical approximation. A simple approach is the method of [finite differences](https://en.wikipedia.org/wiki/Finite_difference). Here, the tangent of the partial derivative is approximated by the secant, leading to an approximation of the form
$$
  \frac{\partial f}{\partial \beta^k_i} \approx 
    \frac{f(\beta^k + h \hat{e}_i) - f(\beta^k - h \hat{e}_i)}{2h}\,,
  {{numeq}}{eq:finite_difference_symmetric}
$$
where $\hat{e}_i$ is the $i$-th unit vector and $h$ is a small value that determines the step size of the approximation. More precisely, this equation {{eqref: eq:finite_difference_symmetric}} represents the *central finite difference of order 2*. There are also *one-sided* approximations and higher-order approximations, which we will not discuss further here.

### Implementation

#### Finite Differences

The first step is to implement the finite difference method. Since we will need to calculate derivatives multiple times, it is useful to implement a *function*. Functions in the programming context are similar to mathematical functions that transform an input into an output, given a set of rules. In contrast to the code blocks that we have seen so far, the code that defines a function is not executed immediately. Instead, it is only executed when the function is called.

We first import `numpy`:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:import_numpy}}
```

Then we define the function `finite_difference`, which calculates the gradient of a function `func` at the point `beta`:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:finite_difference}}
```
The first line of the function is the *signature*, which defines the name of the function and the names and types of the arguments. 

The first line of the function is the *signature*, which defines the name of the function and its arguments. The function takes four arguments:
- `func`: The function to be differentiated
- `beta`: The point at which to calculate the gradient
- `h`: The step size for the finite difference approximation. This is a default argument with value 1e-5, meaning it can be omitted when calling the function
- `args`: Additional arguments to pass to `func`. This is a default argument with value `()`, an empty tuple. Like `h`, it can be omitted when calling the function. Default arguments must always come after required arguments in the function signature

```admonish note title="Dynamic typing in Python"
In strongly typed languages, function signatures also define the types of the arguments and the return type. Since Python is a dynamically typed language, the types of the arguments are only defined when the function is called. We therefore have to be careful in naming the arguments, such that the type of the argument is clear from the name.
```

In the body of the function, we first determine the dimension of the parameters `beta` and store it in the integer `n`. The gradient of a real-valued function with $n$ variables is a vector of length $n$. Therefore, we initialize the variable `grad` as a numpy array of length `n` filled with zeros using the function [`np.zeros`](https://numpy.org/doc/stable/reference/generated/numpy.zeros.html).

We then use a [*for*-loop](https://docs.python.org/3/tutorial/controlflow.html#for-statements) to calculate the gradient component by component. The loop iterates over all dimensions `i` from `0` to `n - 1`. Note that Python indices start at `0` and the built-in function [`range`](https://docs.python.org/3/library/functions.html#func-range) excludes the end value `n`.

In each iteration of the loop, we first define the unit vector $\hat{e}_i$ by first creating a zero array of length `n` with `np.zeros(n)` and then setting the $i$-th component to `1`. We then can calculate the $i$-th entry of the gradient array `grad` according to Eq. {{eqref: eq:finite_difference_symmetric}}. Additionally, we pass the arguments `args` to the function `func`, with a star `*` in front of the variable. The use of `*` acts as a unary operator at this point, which unpacks an object into its components (see [*unpacking*](https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists)). This means that the function `func` accepts the further arguments in the tuple `args` **individually** after the first argument. We will see an example for this in the next section.

After the loop completes, we return the gradient array `grad`.

#### Objective Function

Next, we implement the function whose value we want to minimise. This can be any function, but we choose the least squares loss function {{eqref: eq:least_squares_loss_linear}} of the linear regression model, which we implement as follows:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:least_squares_loss}}
```

Inline with our previous notation, we use the variable `beta` to represent the parameters, which we will later optimise. The second and third arguments are `x` and `y`, which are the data points. 

The gradient of the least squares loss function with respect to the parameters can now easily be calculated using the function `finite_difference` that we implemented earlier. But first, let's implement the gradient descent algorithm.

#### Gradient Descent

According to Eq. {{eqref: eq:gradient_descent}}, the algorithm requires the gradient of the objective function, the initial parameters $\beta^0$, the step size $\alpha$, and the maximum number of iterations. Instead of the gradient, we can also pass the function directly to the algorithm and estimate the gradient numerically using the `finite_difference` function. As a stopping condition, we use a combination of the maximum number of iterations and the norm of the gradient. Therefore, we need the additional arguments `max_iter` and `max_norm`. We again need the argument `args`, which contains the additional arguments for the function that we want to optimise.
```python
{{#include ../codes/01-regression/numerical_optimisation.py:gradient_descent}}
```

In the body of the function, we first explicitly copy the array `beta0` using `np.copy` and store the copy in the variable `beta`. This is necessary because Python does not copy arrays by default, but only stores references to them. Illustrating this with an example, if we simply rename `a` to `b`, a change in `a` will also cause a change in `b`, and *vice versa*. 

Then we use a `for`-loop that iterates the variable `niter` from `0` to `max_iter - 1`. In each iteration, we calculate the gradient `grad` by calling the function `finite_difference` with the arguments `func`, `beta` and `args`. Then we use Eq. {{eqref: eq:gradient_descent}} to update the variable `beta`. Note that Python evaluates the right-hand side of the assignment before assigning the result to `beta`. After that, we calculate the norm of the gradient with the function [`np.linalg.norm`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html) and check if it is smaller than `max_norm`. If so, we break the loop with the `break` command.

After the last iteration, we check if the variable `niter` has reached the maximum number of iterations `max_iter`. If so, it means that the gradient is still large, and the method may not have converged, and we issue a warning. 

At the end, we return the optimum `beta` and the number of iterations `niter` that were actually needed.

### Application

Let's apply our implementation to the data from Chapter [1.2](02-linear_regression.md). First, we define our data arrays:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:data_array}}
```

Then we set an initial guess for the parameters and run gradient descent:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:gradient_descent_call}}
```

The optimal parameters $\beta_0$ and $\beta_1$ are identical to those of the analytical solution. On the author's computer, 34,683 iterations were needed to satisfy the stopping condition. The exact number of iterations may vary slightly depending on the hardware. If you choose the step size $\alpha$ slightly larger, fewer iterations are needed. However, if $\alpha$ is too large, the method diverges.

```admonish tip title="Task"
  Try changing the step size `alpha` and observe how the number of iterations changes.
```

### Built-in Optimisation

Since optimisation is a very general problem, there are many implementations of various algorithms in libraries such as [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html). We want to use the function [`scipy.optimize.minimize`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) to find the optimal parameters $\beta^*$ and compare them with our implementation:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:scipy_minimize}}
```
When calling the function `minimize`, we only need to specify the objective function `objective_function` and the initial point `beta_guess`. The numerical gradient calculation is handled internally.

~~~admonish note title="Gradient Calculation"
If the analytical gradient is available, or you have a good estimate of 
the gradient at hand, you can pass it as the `jac` argument to the 
`minimize` function, which can speed up the optimisation process. 
The `jac` argument expects a function of the form `jac(beta, arg1, arg2, ...)`, 
which takes the current parameters `beta` and the additional arguments 
and returns the gradient. Therefore, we have to use a `lambda` function 
to convert the `finite_difference` function into the required format.

If the `jac` argument is not provided, the `minimize` function will
automatically calculate the gradient using the finite difference method.
Therefore, it is totally fine to not provide the `jac` argument 
in this case.
~~~

With the argument `method='CG'`, we have selected the 
[<i>nonlinear **C**onjugate **G**radient method</i>](https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method). This method is an improvement of the gradient method and reaches the minimum in **only 2 iterations** on the author's computer.

~~~admonish note title="Other Optimisation Methods" collapsible=true
In `minimize`, there are several minimisation methods available. An overview can be found in the [documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) of this function. Two important methods are:
- `method='Nelder-Mead'`: 
    The [Nelder-Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) is a heuristic method that does not require gradient calculation. It is particularly useful when gradient calculation is computationally expensive or when the gradient is noisy. This method is especially suitable for regressions with experimental data.
- `method='BFGS'`:
    The [**B**royden-**F**letcher-**G**oldfarb-**S**hanno method](https://en.wikipedia.org/wiki/BFGS_method) uses the gradient to approximate the Hessian matrix of the objective function. The Hessian contains the second derivatives of the function with respect to its parameters. This method shows very fast convergence near the minimum and typically requires fewer iterations than other optimisation methods, making it a popular choice.
~~~

As you will see in the exercise, linear regression with the least squares method has a closed-form solution that directly calculates the optimal parameters. So why should we use numerical optimisation? In the context of regression, numerical optimisation allows us to:
- Use more complicated models that do not have analytical solutions
- Implement more sophisticated loss functions, such as the least absolute deviations (see Eq. {{eqref: eq:least_absolute_deviations_opt}})
- Introduce additional control over the parameters through **regularisation**, which can improve the model's general performance

~~~admonish note title="Function `scipy.optimize.curve_fit`"
There is also the function [`scipy.optimize.curve_fit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html), which performs nonlinear regression directly. However, it is less flexible than the general `minimize` function, as it only offers a few optimisation methods and uses the least squares loss function. It is suitable for simple regressions and typically requires less code.
~~~

---

**Self-Study Questions**

1. Explain why the step size $\alpha$ in gradient descent is crucial for the algorithm's performance. What happens if $\alpha$ is too small or too large?

2. Illustrate graphically why gradient descent converges to a local minimum, and not the global minimum.

3. Compare and contrast the finite difference method for calculating gradients with analytical gradient calculation. What are the advantages and disadvantages of each approach?

4. To understand why we need the `lambda` function, try to pass the `finite_difference` function to the `jac` argument of the `minimize` function directly. Consult the documentation of `scipy.optimize.minimize` for detailed information.

**Challenge Questions**

1. Show mathematically that one step of the gradient descent method yields a better estimate of the minimum than the initial guess by expanding $f(\beta^{k+1})$ in a Taylor series around $\beta^k$ up to first order, using Eq. {{eqref: eq:gradient_descent}}.

2. Think about how to mitigate the problem of gradient descent converging to a local minimum. How could stochasticity help?

