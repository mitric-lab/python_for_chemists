## Nonlineare Regression

In the previous sections, we learned about linear regression, which allows us 
to model linear relationships between variables. While many physical 
relationships are linear or can be formulated as such, there are also many 
nonlinear relationships that we want to model. In the context of the 
least squares method, we simply replace the model $\hat{f(\beta; x)}$ in 
Eq. {{eqref: eq:least_squares_opt}} with a nonlinear function. However, 
in this case, an analytical solution like 
Eq. {{eqref: eq:least_squares_linear_params}} is not always possible, 
which is why numerical optimisation methods must be used.

### Application

#### Reaction Kinetics

In the Physical Chemistry lab, you probably performed the experiment
"Determination of the rate constant and activation energy of the
Manganese(III)-trioxalate decomposition reaction", also known as 
"Mn-Zerfall", by measuring the absorbance $A$ as a function of time $t$ and
fitting the measurement data to determine the rate constant $k$. 
The underlying relationship is exponential:
$$
  A(t) = A_0\, \eu^{-k t}
$$
with the parameters $A_0$ and $k$, i.e. $\beta = (A_0, k)^\intercal$.

Here, we have $A(t) = \hat{f}(\beta; t)$, and we can formulate the regression
problem using the least squares loss function as the following
optimisation problem:
$$
  \beta^{* } = \argmin{\beta\in\mathbb{R}^2} \sum_{i=1}^N\, (A_i - A_0\, \eu^{-k t_i})^2
  {{numeq}}{eq:least_squares_exp_opt}
$$

For the implementation, we first import the required modules and libraries
as always:
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:import}}
```

Next, we need to provide the data in the form of arrays. Since we have
quite a few data points here, manually typing them in would be quite tedious.
Therefore, we use the function 
[`np.loadtxt`](https://numpy.org/doc/stable/reference/generated/numpy.loadtxt.html)
to read the data from a text file. The text file `mn_decay.txt` 
(<a href="../codes/01-regression/mn_decay.txt" download>download here</a>)
has two columns containing the values for time $t$ and
absorbance $A$. The first few lines of the file look like this:
```txt
{{#include ../codes/01-regression/mn_decay.txt::10}}
```

To read the data, the function `np.loadtxt` requires the filename:
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:read_data}}
```
In this case, the text file is located in the same directory as the 
executing script. If you have the text file in a different directory, you 
need to adjust the path accordingly. The optional argument `unpack=True` 
makes the columns of the file output separately and stores them as arrays
`time` and `absorbance`. If we set `unpack=False`, which is also the default
value, the function would return a 2D array in which the columns are
combined. The first line of this file starts with a comment character `#`,
which causes `np.loadtxt` to ignore this line.

```admonish note title="Note"
With the optional argument `comments`, we can change the character used for comments.
```


Now we can define the model:
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:exp_model}}
```

Although `t` is an array and `k` is a scalar, the multiplication
`-k * t` works element-wise. In this case, `k` is also interpreted as an array,
which is called
[broadcasting](https://numpy.org/doc/stable/user/basics.broadcasting.html).
The function `np.exp` calculates the element-wise exponential
of the array. The result is again an array, which we multiply with 
the scalar `a0`.

Next, we define the objective function in
Eq. {{eqref: eq:least_squares_exp_opt}}:
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:objective_function}}
```

Now we use the `minimize` function to solve this
optimisation problem:
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:optimise}}
```

Here, we applied the Nelder-Mead method with the initial parameters
$A_0^0 = 1$ and $k^0 = 0.01$. To print the result with the
`print` command, we placed an `f` before the string.
This indicates that the string is a so-called
[f-string](https://realpython.com/python-f-strings/#doing-string-interpolation-with-f-strings-in-python)
in which we can embed variables with curly braces `{}`.
In fact, f-strings can do much more, which we will see in the future.

The optimised parameters should have the following values:
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:verification}}
```

Finally, we can plot the results:
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:plot}}
```

You should be familiar with most of the functions in the above code block from
Chapter [1.2](02-linear_regression.md). One difference is the use of the function
[`np.linspace`](https://numpy.org/doc/stable/reference/generated/numpy.linspace.html),
which generates times between the measurement points so that we can
use the regression model for interpolation.
This function takes three arguments: the start value, the end value, and
the number of points to be generated. It then produces an array
with evenly distributed values between the start and end value.

Another new function is
`fig.tight_layout()`, which automatically adjusts the layout of the plot.
We can see from the following diagram that the exponential function describes the data
very well.
![Exponential Regression of the Mn Decay](../assets/figures/01-regression/nonlinreg_mn.svg)

Some of you might ask why we didn't linearise the function to use
linear regression, which is a valid question. In fact, it is possible 
to linearise the function by
taking the logarithm of both sides of the equation:
$$
  \ln(A(t)) = \ln(A_0) - k t
$$
A linear regression using this equation, however, does not yield the
same results, as you can see in the following diagram:
![Exponential Regression of the Mn Decay with Linearised Function](../assets/figures/01-regression/nonlinreg_mn_wlin.svg)

As a voluntary exercise, you can try to reproduce the diagram above.
It is easy to see that the linearised fit does not match the data as well.
This is because linear regression does not treat the errors in absorbance
equally when taking the logarithm. The resulting parameters
```python
{{#include ../codes/01-regression/nonlinreg_mn.py:verification_lin}}
```
are indeed quite different from the previous ones.
Therefore, it is often necessary to perform nonlinear regressions on the original
data rather than using linear models with linearised data.

#### Titration Curve

In the analytical chemistry lab, you probably performed a titration of 
a strong base against a strong acid with a pH meter. At that time, 
you probably had to plot the values on millimeter paper and determine the
equivalence point based on the position of the pH jump. This is tedious
and inaccurate, as only the few data points near the steep rise are
considered.

Since the pH curve is a function of the amount of base added, we can
model it using nonlinear regression and determine the equivalence point
with much higher accuracy.

The $\mathrm{H^+}$ concentration during the titration of a strong
base against a strong acid is (under certain approximations) given by:
$$
  [\mathrm{H^+}] = \frac{\Delta + \sqrt{\Delta^2 + 4K_w}}{2}\,,
  {{numeq}}{eq:titration_sasb_hplus}
$$
where $\Delta = [\mathrm{A^-}] - [\mathrm{B^+}]$ is the concentration difference
between the counterions of the acid $\mathrm{A^-}$ and the base $\mathrm{B^+}$.
$K_w$ is the ion product of water.

Since strong acids and bases dissociate completely, the concentrations of their
counterions can be expressed as follows:
$$
  \begin{align}
    [\mathrm{A^-}] &= \frac{n_\mathrm{A}}{V} 
      = \frac{c_\mathrm{A}^0 V^0}{V^0 + V_\mathrm{B}} \\
    [\mathrm{B^+}] &= \frac{n_\mathrm{B}}{V}
      = \frac{c_\mathrm{B}^0 V_\mathrm{B}}{V^0 + V_\mathrm{B}}\,,
  \end{align}
$$
where $c_\mathrm{A}^0$ and $c_\mathrm{B}^0$ are the concentrations of the acid
and base to be analysed, $V^0$ is the initial volume of the sample solution,
and $V_\mathrm{B}$ is the volume of the added base.

The pH value can be calculated from the $\mathrm{H^+}$ concentration:
$$
  \mathrm{pH} = -\lg \left( \frac{[\mathrm{H^+}]}{1\ \mathrm{M}} \right)\,.
  {{numeq}}{eq:titration_sasb_ph}
$$

By combining the equations {{eqref: eq:titration_sasb_hplus}} and
{{eqref: eq:titration_sasb_ph}} into a single function, we obtain the model
$f(\beta; V_\mathrm{B})$, where
$\beta = (c_\mathrm{A}^0, V^0)^\intercal$.

As before, we first implement the model and the objective function:
```python
{{#include ../codes/01-regression/nonlinreg_titration.py:titration_model}}
```
```python
{{#include ../codes/01-regression/nonlinreg_titration.py:objective_function}}
```

Although the function of the pH value is relatively complicated, we can
simplify the implementation of the function in Python significantly by
defining intermediate variables, as in
Eq. {{eqref: eq:titration_sasb_hplus}}. The objective function is almost 
identical to that of the "Mn-Zerfall", with the main difference being 
the replacement of our model.

Just like in the previous example, we read the data from a text file
(<a href="../codes/01-regression/titration_sasb.txt" download>download here</a>):
```python
{{#include ../codes/01-regression/nonlinreg_titration.py:read_data}}
```
Additionally, we have defined the concentration of the titrant `C0_B`.
According to the general convention, all constants in Python should be written 
in uppercase letters.
After that, we can perform the nonlinear regression
and plot the results:
```python
{{#include ../codes/01-regression/nonlinreg_titration.py:optimise}}
```
```python
{{#include ../codes/01-regression/nonlinreg_titration.py:plot}}
```

From the optimised parameters, we determined the amount of the analyte
$n_\mathrm{A}^0 = 1.7698\ \mathrm{mmol}$. The corresponding diagram should 
look like this:
![Nonlinear Regression of the Titration Curve](../assets/figures/01-regression/nonlinreg_titration.svg)
Although the fit is not perfect at the beginning and end of the curve,
the agreement near the equivalence point is excellent.

You may have wondered why the optimised parameters $c_\mathrm{A}^0$
and $V^0$ are not explicitly listed here, but only their product.
```admonish warning title="Warning: correlated parameters"
The model parameters $\beta$ can be correlated, i.e., changing
two or more parameters leads to a similar change in the
objective function. This circumstance may indicate so-called 
*overfitting* of the data by the model.
In this case, it is advisable to check whether the model can also
be simplified with fewer parameters.

In our case, the parameters $c_\mathrm{A}^0$ and $V^0$ are correlated,
as it is the amount of the analyte $n_\mathrm{A}^0$ that has a significant
impact on the titration curve. The initial volume $V^0$, on the other hand,
plays only a minor role, which is why the product $c_\mathrm{A}^0 V^0$ serves
as an approximate single parameter. Nevertheless, $V^0$ is an important
parameter here, as it is on the same order of magnitude as $V_\mathrm{B}$ and
the dilution cannot be neglected. However, since we only need the product
$c_\mathrm{A}^0 V^0$ for the titration analysis,
the correlation of the parameters does not bother us in this case.
```

```admonish tip title="Tip"
Change the initial guesses and observe how the
optimised parameters change, but their product remains nearly constant.
```

Finally, it should be emphasised again that it is often important to
perform the regression on the original data and not to use
transformed data. In the case of the acid-base titration, one might think
of fitting $[\mathrm{H^+}]$ instead of the pH value. For the same
reasons as before, this is not sensible: Since the same error on the pH scale 
leads to larger
errors at higher $\mathrm{H^+}$ concentrations and smaller errors
at lower $\mathrm{H^+}$ concentrations, the
data points at higher $\mathrm{H^+}$ concentrations are fitted better.
This leads to a distortion of the results, as can be seen in the following
diagram:
![Nonlinear Regression of the Titration Curve with Regression on H+ Concentration](../assets/figures/01-regression/nonlinreg_titration_ch.svg)
It can be seen that the regression on $[\mathrm{H^+}]$
favours the earlier data points, resulting in a completely
incorrect determination of the equivalence point.

As a voluntary exercise, you can try to reproduce the diagram above.

