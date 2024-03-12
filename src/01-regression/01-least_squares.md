## Least Squares

In all regression problems, we aim to find a model that best fits the data, 
which can hopefully be used to predict future data points. 
For simplicity, we shall assume that our data only has one independent
variable ($x$) and one dependent variable ($y$). A generalisation to
multiple independent variables is straightforward, and will be
discussed in a later chapter. In addition, we shall assume that the 
data is real-valued. For a total of $N$ data points, we write

<!-- i18n:skip --> $$
  \begin{align*}
    x &= (x_1, x_2, \ldots, x_N) \in \mathbb{R}^N \\
    y &= (y_1, y_2, \ldots, y_N) \in \mathbb{R}^N \,,
  \end{align*}
$$ 

$$$$ i.e., $x$ and $y$ are $N$-dimensional vectors with real-valued components 
$x_i$ and $y_i$. 

We can think of the relationship between $x_i$ and $y_i$ as a function $f$.
In the ideal case, we have

<!-- i18n:skip --> $$
  y_i = f(x_i) + \epsilon_i \,, {{numeq}}{eq:regression_true}
$$

$$$$ where the function $f$ describes the *true* relationship between the 
data points $x_i$ and $y_i$, while $\epsilon_i$ describtes the random
statistical error. 

In the real world, however, the relationship $f$ is often unknown, and
we must use a model to estimate it. In this case, we can write

<!-- i18n:skip --> $$
  y_i = \hat{f}(\beta; x_i) + \hat{\epsilon}_ i\,, {{numeq}}{eq:regression_estimated}
$$

$$$$ where $\hat{f}(\beta; x_i)$ is the model estimator of $f$ parametrised by
$\beta \in \mathbb{R}^n$, with $n$ being the number of parameters in the
model. The error term $\hat{\epsilon}_i$ now is a combination of the 
statistical error and the unmodeled part of the relationship. The goal of 
the regression analysis is to find the parameters $\beta$ that best fit 
the data. 

But what does "best fit" mean? On common approach is to use the
*least squares* method, which aims to minimise the sum of the squared
errors $\hat{\epsilon}_i$, i.e. the least squares loss function

<!-- i18n:skip --> $$
  L(\beta; x, y) 
    = \sum_{i=1}^N\, \hat{\epsilon}_i^2
    = \sum_{i=1}^N\, (y_i - \hat{f}(\beta; x_i))^2 \,. {{numeq}}{eq:least_squares_loss}
$$

$$$$ is minimised. Written mathematically, we want to find the parameters 
$\beta^* $ using

<!-- i18n:skip --> $$
  \beta^* 
    = \underset{\beta\in\mathbb{R}^n}{\mathrm{arg\,min}}\ L(\beta; x, y)
    = \underset{\beta\in\mathbb{R}^n}{\mathrm{arg\,min}} \sum_{i=1}^N\, (y_i - \hat{f}(\beta; x_i))^2 \,. {{numeq}}{eq:least_squares_opt}
$$

$$$$ If we define the vector of predicted values $\hat{y}$ as

<!-- i18n:skip --> $$
  \hat{y}_i = \hat{f}(\beta; x_i)\,,
$$

$$$$ we can write eq. {{eqref: eq:least_squares_opt}} as

<!-- i18n:skip --> $$
  \beta^* 
    = \underset{\beta\in\mathbb{R}^n}{\mathrm{arg\,min}}\ \sum_{i=1}^N\, (y_i - \hat{y}_i)^2 
    = \underset{\beta\in\mathbb{R}^n}{\mathrm{arg\,min}}\ \|y - \hat{y}\|_2^2 \,, {{numeq}}{eq:least_squares_opt2}
$$

$$$$ where $\| v \|_2$ denotes the Euclidean norm or the $\ell_2$-norm of 
a vector $v$, defined as

<!-- i18n:skip --> $$
  \| v \|_2 = \sqrt{\sum_{i=1}^N\, v_i^2} \,.
$$

$$$$ The least squares method is a popular choice for regression analysis
because it has closed-form solutions for some simple but important
models, and because it can be solved efficiently using numerical
linear algebra methods for more complex models. However, it has some
drawbacks, such as being sensitive to outliers and prone to overfitting. 
An alternative to the least squares method is the least absolute
deviations method, which uses the $\ell_1$-norm instead of the
$\ell_2$-norm, defined for a vector $v$ as

<!-- i18n:skip --> $$
  \| v \|_1 = \sum_{i=1}^N\, |v_i| \,.
$$

$$$$ The optimisation problem is then formulated as:

<!-- i18n:skip --> $$
  \beta^* 
    = \underset{\beta\in\mathbb{R}^n}{\mathrm{arg\,min}}\ \|y - \hat{y}\|_1 \,. {{numeq}}{eq:least_absolute_deviations_opt}
$$

$$$$ There are many other loss functions that can be used for regression, 
e.g. the [Huber loss](https://en.wikipedia.org/wiki/Huber_loss),
the log-cosh loss, 
the [quantile loss](https://en.wikipedia.org/wiki/Quantile_regression#Quantile_of_a_random_variable),
etc. While the least squares method works well in many cases, it is
important to be aware of the limitations of the method and to consider
other loss functions when appropriate.

