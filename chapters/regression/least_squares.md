(sec:least_squares)=

# Least Squares

Regression analysis aims to identify and quantify relationships between variables in a dataset. The goal is to construct a mathematical model that not only explains observed data but also enables predictions for new observations. In this introduction, we focus on the simplest case with one independent variable ($x$) that influences one dependent variable ($y$). While extensions to multiple independent variables are possible and will be covered later, this basic scenario helps establish fundamental concepts. We consider real-valued data consisting of $N$ paired observations, which we can write as vectors:
:::{math}
:enumerated: false
  \begin{aligned}
    x &= (x_1, x_2, \ldots, x_N) \in \mathbb{R}^N \\
    y &= (y_1, y_2, \ldots, y_N) \in \mathbb{R}^N \,,
  \end{aligned}
:::
where $x$ and $y$ are $N$-dimensional vectors containing the real-valued measurements $x_i$ and $y_i$ for each observation $i$.

We can think of the relationship between $x_i$ and $y_i$ as a function $f$. Under idealized conditions, we have
:::{math}
:label: eq:regression_true
  y_i = f(x_i) + \epsilon_i \,,
:::
where the function $f$ describes the *true* relationship between the data points $x_i$ and $y_i$, while $\epsilon_i$ represents the random statistical error.

In the real world, however, the relationship $f$ is often unknown. Therefore, we need to use a model to estimate it. In this case, we write
:::{math}
:label: eq:regression_estimated

  y_i = \hat{f}(\beta; x_i) + \hat{\epsilon}_ i\,.
:::
where $\hat{f}(\beta; x_i)$ is the model's estimator of $f$, parameterized by $\beta \in \mathbb{R}^n$ with $n$ parameters in the model. The error term $\hat{\epsilon}_i$ is now a combination of the statistical error and the part of the true relationship not captured by the model. The goal of regression analysis is to find the parameters $\beta$ such that the model best represents the data.

But what does "best" mean in this context and how can we find the best parameters $\beta$? A common approach is the *least squares method*, which aims to minimize the sum of the squared residuals ($\hat{\epsilon}_i$), representing the differences between observed data points and the model's predictions. We can write the least squares loss function as
:::{math}
:label: eq:least_squares_loss

  L(\beta; x, y) 
    = \sum_{i=1}^N\, \hat{\epsilon}_i^2
    = \sum_{i=1}^N\, (y_i - \hat{f}(\beta; x_i))^2 \,.
:::
Mathematically speaking, we want to find the parameters $\beta^*$ by solving the optimization problem
:::{math}
:label: eq:least_squares_opt

  \beta^* 
    = \argmin_{\beta\in\mathbb{R}^n} L(\beta; x, y)
    = \argmin_{\beta\in\mathbb{R}^n} \sum_{i=1}^N\, (y_i - \hat{f}(\beta; x_i))^2 \,.
:::
If we write the vector of predicted values $\hat{y}$ as
:::{math}
:enumerated: false

  \hat{y}_ i = \hat{f}(\beta; x_i)\,,
:::
we can rewrite Eq. {eq}`eq:least_squares_opt` as
:::{math}
:label: eq:least_squares_opt2

  \beta^{* } 
    = \argmin_{\beta\in\mathbb{R}^n} \sum_{i=1}^N\, (y_ i - \hat{y}_ i)^2 
    = \argmin_{\beta\in\mathbb{R}^n} \|y - \hat{y}\|_ 2^2 \,.
:::
where $\| v \|_2$ denotes the Euclidean norm or the $\ell_2$-norm of a
vector $v$, defined as
:::{math}
:enumerated: false

  \| v \|_2 = \sqrt{\sum_{i=1}^N\, v_i^2} \,.
:::

The least squares method is a popular choice for regression analysis because it provides closed-form solutions for some simple but important models and can be solved efficiently with numerical linear algebra methods. However, it has some disadvantages, such as sensitivity to outliers and vulnerability to overfitting. An alternative to the least squares method is the *least absolute deviations method*, which uses the $\ell_1$-norm instead of the $\ell_2$-norm. This is defined for a vector $v$ as
:::{math}
:enumerated: false

  \| v \|_1 = \sum_{i=1}^N\, |v_i| \,.
:::
Because the $\ell_1$-norm sums absolute errors rather than squared errors, it gives less weight to large deviations, making this method less sensitive to outliers compared to the least squares ($\ell_2$-norm) approach. The optimization problem is then formulated as
:::{math}
:label: eq:least_absolute_deviations_opt

  \beta^* 
    = \argmin_{\beta\in\mathbb{R}^n} \|y - \hat{y}\|_1 \,.
:::

<!-- There are many other loss functions that can be used for regression,
e.g. the
- [Huber loss function](https://en.wikipedia.org/wiki/Huber_loss),
- log-cosh loss function,
- [quantile loss function](https://en.wikipedia.org/wiki/Quantile_regression#Quantile_of_a_random_variable)
etc.
Although the least squares method works well in many cases, it is important
to be aware of its limitations and to consider other loss functions 
if necessary. -->
