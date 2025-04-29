## Linear Regression

One of the simplest models in regression analysis is the linear
model, where the estimator $\hat{f}(\beta; x_i)$ is given by a linear
function of the form
$$
  \hat{f}(\beta; x_i) = \beta_0 + \beta_1 x_i
  {{numeq}}{eq:linear_model}
$$
with the scalar parameters $\beta_0$ and $\beta_1$.
The regression analysis with the linear model is called linear
regression.

### Theoretical Background

Let us now insert the linear model (cf. Eq. {{eqref: eq:linear_model}})
into the least squares loss function according to Eq. 
{{eqref: eq:least_squares_loss}}:
$$
  L(\beta; x, y) 
    = \sum_{i=1}^N\, (y_i - \beta_0 - \beta_1 x_i)^2 \,.
  {{numeq}}{eq:least_squares_loss_linear}
$$

Since our goal is to find the parameters $\beta$ that minimize
the loss function (cf. Eq. {{eqref: eq:least_squares_opt}}),
we need to satisfy at least the necessary conditions for this,
i.e. the partial derivatives of $L$ with respect to $\beta_0$ and
$\beta_1$ must vanish:
$$
  \frac{\partial L}{\partial \beta_0} \stackrel{!}{=} 0 
  \quad \text{und} \quad
  \frac{\partial L}{\partial \beta_1} \stackrel{!}{=} 0 \,.
$$

After a somewhat tedious but simple
```admonish derivation title="Derivation" collapsible=true
First, we calculate the partial derivatives of $L$ with respect to $\beta_0$
and $\beta_1$. Since differentiation is linear, we can split the sum in
Eq. {{eqref: eq:least_squares_loss_linear}} over the individual terms:
$$
  \begin{align}
    \frac{\partial L}{\partial \beta_0} &= \frac{\partial}{\partial \beta_0} \sum_{i=1}^N\, (y_i - \beta_0 - \beta_1 x_i)^2 \\
    &= \sum_{i=1}^N\, \frac{\partial}{\partial \beta_0} (y_i - \beta_0 - \beta_1 x_i)^2 \\
    &= \sum_{i=1}^N\, -2 (y_i - \beta_0 - \beta_1 x_i) \,.
  \end{align}
$$
In the same manner, we proceed with the derivative with respect to $\beta_1$ and obtain
$$
  \frac{\partial L}{\partial \beta_1} = \sum_{i=1}^N\, -2 x_i (y_i - \beta_0 - \beta_1 x_i) \,.
$$

Setting the derivatives equal to zero gives us the necessary conditions
for a minimum:
$$
  \begin{align}
    0 &= \sum_{i=1}^N\, -2 (y_i - \beta_0 - \beta_1 x_i) \\
    0 &= \sum_{i=1}^N\, -2 x_i (y_i - \beta_0 - \beta_1 x_i) \,.
  \end{align}
$$
The factor $(-2)$ does not matter, since the expression is set to zero.
Therefore, we can omit it in the following.

By simply rearranging, we can write the above equations as a linear system
of equations in $\beta_0$ and $\beta_1$:
$$
  \begin{align}
    \sum_{i=1}^N\, \beta_0 + \sum_{i=1}^N\, \beta_1 x_i &= \sum_{i=1}^N\, y_i \\
    \sum_{i=1}^N\, \beta_0 x_i + \sum_{i=1}^N\, \beta_1 x_i^2 &= \sum_{i=1}^N\, x_i y_i \,.
  \end{align}
$$
Since the parameters $\beta_0$ and $\beta_1$ are independent of the 
data indices $i$, we can factor them out of the sums and obtain
$$
  \begin{align}
    \beta_0 N + \beta_1 \sum_{i=1}^N\, x_i &= \sum_{i=1}^N\, y_i \\
    \beta_0 \sum_{i=1}^N\, x_i + \beta_1 \sum_{i=1}^N\, x_i^2 &= \sum_{i=1}^N\, x_i y_i \,,
  \end{align}
$$
where we used $\sum_{i=1}^N\, 1 = N$.
This system of equations is equivalent to the matrix equation in
Eq. {{eqref: eq:least_squares_linear_params}}.
```
we obtain the matrix equation for the parameters $\beta_0$ and $\beta_1$:
$$
  \underbrace{
  \begin{pmatrix}
    \displaystyle N & \displaystyle \sum_{i=1}^N\, x_i \\[1.5em]
    \displaystyle \sum_{i=1}^N\, x_i & \displaystyle \sum_{i=1}^N\, x_i^2
  \end{pmatrix}}_{\displaystyle \bm{A}}
  \,
  \underbrace{
  \begin{pmatrix}
    \displaystyle \beta_0 \\[1.5em]
    \displaystyle \beta_1
  \end{pmatrix}
  \vphantom{
    \begin{pmatrix}
      \displaystyle \sum_{i=1}^N\, y_i \\[1.5em]
      \displaystyle \sum_{i=1}^N\, x_i y_i
    \end{pmatrix}
  }
  }_{\displaystyle \vec{x}}
  =
  \underbrace{
  \begin{pmatrix}
    \displaystyle \sum_{i=1}^N\, y_i \\[1.5em]
    \displaystyle \sum_{i=1}^N\, x_i y_i
  \end{pmatrix}}_ {\displaystyle \vec{b}} \,.
  {{numeq}}{eq:least_squares_linear_params}
$$

The solution of this system of equations provides us with the optimal
parameters $\beta_0$ and $\beta_1$ for the linear model.

The (at least formal) solution of the system of equations 
{{eqref: eq:least_squares_linear_params}} is
$\vec{x} = \bm{A}^{-1} \vec{b}$, with the inverse of the matrix
$\bm{A}$. A matrix is invertible if and only if its determinant is non-zero.
The determinant of the system matrix $\bm{A}$ is
$$
  \begin{align}
    \det(\bm{A}) &= N \sum_{i=1}^N\, x_i^2 - \left(\sum_{i=1}^N\, x_i\right)^2 \\
    &= N^2 \left[ 
      \frac{1}{N} \sum_{i=1}^N\, x_i^2 - \left(\frac{1}{N} \sum_{i=1}^N\, x_i\right)^2
    \right] \\
    &= N^2 \left[ 
      \frac{1}{N} \sum_{i=1}^N\, x_i^2 - \bar{x}^2
    \right] \\
    &= N^2 \left[ \frac{1}{N} \sum_{i=1}^N\, (x_i - \bar{x})^2 \right] \,,
  \end{align}
$$
where we have introduced the abbreviation 
$\bar{x} = \frac{1}{N} \sum_{i=1}^N\, x_i$ for the mean of the $x_i$.
If we had only one data point, the mean would be equal to the only data point,
and the determinant of the system matrix would be zero. This means that the
matrix $\bm{A}$ is not invertible, and we do not obtain a unique solution for
the linear model. Only with at least two different data points, a unique
solution of the system of equations {{eqref: eq:least_squares_linear_params}}
can be found.

Now we want to implement the linear regression for a simple example with 
the help of Eq. {{eqref: eq:least_squares_linear_params}}.

### Implementation
Consider the following measurement data[^1] of the Lambert-Beer
relationship for methylene blue in water at different concentrations $c$
and the corresponding absorbances $A$, measured at $610\ \mathrm{nm}$
with a layer thickness $d$ of 1 cm:

| $c$ / &micro;M | $A$      | $c$ / &micro;M | $A$      |
|:--------------:|:--------:|:--------------:|:--------:|
| 2.125          | 0.0572   | 23.38          | 0.8242   |
| 4.250          | 0.1391   | 25.50          | 0.9130   |
| 6.375          | 0.2049   | 27.63          | 1.0043   |
| 8.500          | 0.2754   | 29.75          | 1.0809   |
| 10.63          | 0.3420   | 31.88          | 1.1511   |
| 12.75          | 0.4139   | 34.00          | 1.2483   |
| 14.88          | 0.4956   | 36.13          | 1.3373   |
| 17.00          | 0.5815   | 38.25          | 1.4027   |
| 19.13          | 0.6806   | 40.38          | 1.4927   |
| 21.25          | 0.7481   | 42.50          | 1.5853   |

Before we can proceed, we need to import the data into Python in some way.
The easiest way for such a small dataset is to enter it manually.
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:data_list}}
```
Here we have defined the data in the variables `concentrations` 
and `absorbances` of type *List*. This is indicated by the use of square
brackets `[]` and the comma `,` between the individual values. The
equal sign `=` assigns the values to the variables.

Although Python already provides some mathematical functions for lists
by default, the package `numpy` has even more very useful operations
for such data structures. Therefore, we import `numpy` and convert
the lists into the data type *numpy-array*.
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:data_array}}
```
Since it can become tedious to write `numpy` every time, we define the
alias `np` with `import numpy as np`. From now on, we can refer to
the contents of this package with `np`. In the following lines, we
use the function `np.array` to convert the lists into arrays.
Here it becomes clear that the equal sign `=` in programming has a
different meaning than in mathematics. While in mathematics the
equal sign expresses an equivalence between two sides, in programming
it represents an assignment. Therefore, it is possible to use the
same variable name for the arrays, since the right side of the
equal sign is evaluated first and then the value of the right side
is assigned to the variable on the left side.

Next, we want to calculate the elements of $\bm{A}$ and $\vec{b}$ in 
Eq. {{eqref: eq:least_squares_linear_params}}. First, we ensure that
our data arrays are of equal length. Then we can calculate the
elements using the `np.sum` function.

```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:create_system_array_elements}}
```
Here we used the
[built-in function](https://docs.python.org/3/library/functions.html)
`len` to determine the length of the arrays. The function `np.sum` calculates 
the sum of all elements in an array. The asterisk `*` stands for 
multiplication. The use of `*` between two arrays leads to an element-wise
multiplication. The double asterisk `**` means exponentiation in Python, 
meaning that the array `concentrations` is squared element-wise.

```admonish info title="Info for Advanced Learners" collapsible=true
The operators `*` and `**` can also be used as unary operators, i.e. 
operators with only one argument, contrary to binary operators which
require two arguments, such as multiplication. As a unary operator,
they have a different meaning. Interested readers with some programming
experience can read the details, for example, 
[here](https://book.pythontips.com/en/latest/args_and_kwargs.html).
```

We can now calculate the elements of the system matrix $\bm{A}$ and the
vector $\vec{b}$ in Eq. {{eqref: eq:least_squares_linear_params}}:
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:create_system_arrays}}
```

Notice that a list of lists is used to create a matrix 
(or a so-called *2D-array*). The inner list corresponds to a row of the matrix, 
while the outer list contains these rows. Now we can finally solve the
system of equations {{eqref: eq:least_squares_linear_params}} to obtain
the optimal parameters $\beta_0$ and $\beta_1$:
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:solve_system}}
```
The function 
[`np.linalg.solve`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.solve.html)
solves a linear system of equations. We then output the
solution with the `print` function.

If one already knows which values $\beta_0$ and $\beta_1$ should have and
only wants to test the algorithm, one can perform the verification as follows:
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:beta_verification}}
```
Here we first *indexed* the array `beta` with square brackets `[]`. Be aware 
that indexing in Python starts at 0. This means that $\beta_0$ corresponds 
to the 0-th entry and $\beta_1$ to the first entry of the array `beta`.
Then we compared the values of $\beta_0$ and $\beta_1$ with reference values.
```admonish note title="Hinweis"
It was not the exact values that were compared with `==`, since 
`float` numbers 
([floating point numbers](https://en.wikipedia.org/wiki/Floating-point_arithmetic))
cannot be represented exactly in general. Therefore, the function
[`np.isclose`](https://numpy.org/doc/stable/reference/generated/numpy.isclose.html)
was used to compare the values with a (in this case preset) tolerance.
```

The Lambert-Beer law states that the absorbance $A$ is linearly dependent
on the concentration $c$ with the proportionality constant $\varepsilon d$, 
i.e.
$$
  A = \varepsilon d c \,.
$$
The molar extinction coefficient $\varepsilon$ can be calculated from the
linear parameter $\beta_1$ as
$$
  \varepsilon = \beta_1 / d = 0.038\ \mathrm{{\mu M}^{-1}\cdot cm^{-1}} 
    = 38000\ \mathrm{M^{-1}\cdot cm^{-1}} \,.
$$

Now we have calculated the optimal parameters $\beta_0$ and $\beta_1$ for the
above dataset. How can we know how good the linear regression is? 
There are various quality measures, which we will not consider here for the
time being. Instead, we want to graphically represent the optimal parameters
and assess the quality of the regression visually. 
This is probably less rigorous than the mathematically defined quality 
measures, but certainly more fun.

### Visualisation
For the visualisation of the results of the linear regression, we will use
the Python package `matplotlib`. Here we import the submodule `pyplot`
with the alias `plt`:
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:import_mpl}}
```

For the graphical representation with `matplotlib`, we always need
the objects `Figure` and `Axes`. We can create these with the function
[`plt.subplots`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html):
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:make_axes}}
```
Here, the argument `figsize=(8, 6)` is passed, which specifies the size
of the figure in inches. The size of the figure is not set in stone
and can be adjusted as needed.

Now we want to plot the measurement data as points and the linear regression
as a line in the diagram. For this, we use the function
[`plot`](https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.plot.html):
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:plot_data}}
```
For the measurement data, we used the argument `o` so that they are
displayed as points. Other
[markers](https://matplotlib.org/stable/api/markers_api.html)
like `s` for squares or `x` for crosses can also be used. 
Without specifying this argument, the data will be displayed as a line.
In both calls to `plot`, we used the `label` argument. These are used
to name the lines in the legend.

```admonish note title="Note for Advanced Learners" collapsible=true
The `plot` function is actually a method of the `Axes` class.
This means that it is bound to instances of the `Axes` class.
We will learn more about classes and methods in the following chapters.
For now, it is sufficient to know that the `plot` function
must be called on instances of the `Axes` class.
```

For the completeness of the diagram, we will add the axes labels
and a legend:
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:customize_ax}}
```

The functions `set_xlabel` and `set_ylabel` accept a *str* (string)
as an argument, which serves as the axis label. The function `legend`
adds a legend. When called without arguments, the legend is automatically
generated from the `label` arguments of the `plot` functions.

After all elements of the diagram have been added, we can display
the figure with the function `plt.show`:
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:show_plot}}
```

If the program is executed successfully, a diagram like the following
should appear:
![Lineare Regression der Lambert-Beer-Beziehung](../assets/figures/01-regression/linreg_lambert_beer.svg)

This diagram shows the results of the linear regression, even if one can
argue about the aesthetic design of the diagram. Throughout this course,
we will learn more functionalities of `matplotlib` that help us
to design elements of the diagram according to personal taste.

At first glance, the regression line seems to fit the measurement data
very well. To make the difference between the measurement data and the
regression line more apparent, we can represent the residuals, 
i.e. the difference between the measurement data and the regression line, 
in another diagram.
```python
{{#include ../codes/01-regression/linreg_lambert_beer.py:plot_residuals}}
```
After calculating and storing the residuals in the variable `residuals`,
another `Figure` and `Axes` object were created. The residuals were 
represented in a bar diagram, which can be created with the method
`bar`. The diagram looks like this:
![Residuals of the linear regression of the Lambert-Beer relationship](../assets/figures/01-regression/linreg_lambert_beer_residuals.svg)

Now it can be seen that the deviation is positive at low and high 
concentrations, while it is negative at medium concentrations. This could 
indicate a slight positive curvature of the data, which is not captured by
the linear model. However, assigning this behaviour to a systematic or random
deviation usually requires a deeper analysis. In fact, there is a slight
non-linearity in the data, which demonstrates the
[deviation from the Lambert-Beer law](https://chem.libretexts.org/Bookshelves/Analytical_Chemistry/Instrumental_Analysis_(LibreTexts)/13%3A_Introduction_to_Ultraviolet_Visible_Absorption_Spectrometry/13.02%3A_Beer%27s_Law)
especially at high concentrations. Those interested in the exact explanation
of this deviation of methylene blue in water can read the publication
[A. Fernández-Pérez, T. Valdés-Solís, G. Marbán, *Dyes and Pigments* **2019**, *161*, 448&ndash;456](https://www.sciencedirect.com/science/article/pii/S014372081831845X).
Whether the linear regression is reasonable in this case depends on the
desired accuracy of the modelling.

~~~admonish tip title="Tip for Programming Style"
Here we placed the `import` statements in the respective sections to 
emphasise the dependencies of the different parts of the script.
The resulting code can be interpreted without errors.
However, it is common in Python to place the imports at the beginning of
the script, i.e.
```python
import numpy as np
import matplotlib.pyplot as plt

# Rest of your code
```
~~~

---
[^1]: The authors would like to thank Dr. Hans-Christian Schmitt for providing the
lab equipments for the measurements.

### Exercise

#### Problem 1.1: Linear and Quadratic Regression

{{#include ../psets/01.md:aufgabe_1}}

