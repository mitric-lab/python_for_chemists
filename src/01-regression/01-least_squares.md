## Least Squares

In all regression problems, we aim to find a model that best fits the data, 
which can hopefully be used to predict future data points. 
Assuming that our data only has one independent variable $x$ and 
one dependent variable $y$, we can write
$$
  y_i = f(x_i) + \epsilon_i \,, {{numeq}}{eq:regression_true}
$$
where the function $f$ describes the **true** relationship in the data 
between $x$ and $y$, and $\epsilon$ is the statistical error
term. Note that the variables $x$, $y$, and $\epsilon$ are 
vectors of length $n$, where $n$ is the number of data points. 
In eq. {{eqref: eq:regression_true}}, we refer to the $i$-th entry of 
these vectors using the subscript $i$.

In most cases, we do not know the true relationship $f$, and thus we
must use a model to estimate it. In this case, we can write
$$
  y_i = \hat{f}(\beta; x_i) + \hat{\epsilon}_ i\,, {{numeq}}{eq:regression_estimated}
$$
where $\hat{f}(\beta; x_i)$ is the parameterised model 
estimator of $f$ with parameters $\beta$, and $\hat{\epsilon}_ i$
is a combination of the statistical error and the unmodeled part of the
relationship. The goal of regression is to find the parameters $\beta$
that best fit the data.

But what does "best fit" mean? On common approach is to use the
**least squares** method, which aims to minimise the sum of the squared
residuals, which are the difference between the observed and predicted 
values. Written mathematically, we aim to find the parameters $\beta^* $ 
using
$$
  \beta^* = \underset{\beta}{\mathrm{arg\,min}}
  \sum_{i=1}^n\, (y_i - \hat{f}(\beta; x_i))^2 \,. {{numeq}}{eq:least_squares}
$$


