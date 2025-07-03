## Classification

From the previous section, we have seen that finding a line (or hyperplane) through your data, consisting of features and targets, already is the simplest form of supervised learning. You should, however, apply regression only if the targets are continuous, e.g. the intensity of a fluorescence signal. For discrete targets that originate from a finite set of labels, we call this *classification* and we will discuss in this chapter how to effectively tackle such a problem.

Formally, the abstract goal of classification is to find a function $\hat{f}_{\vec{w}} : \mathbb{R}^d \to \{1,2,\dots,K\}$ that maps the input $\vec{x} \in \mathbb{R}^d$ to one of $K$ classes, such that most of the $N$ data points are assigned to the correct class $y_i \in \{1,2,\dots,K\}$. We will focus on the simplest case, binary classification, where the data can be assigned to two possible classes $y_i \in \{-1,+1\}$. An example of such a problem is the prediction of whether a given image contains a cat or a dog, or whether a given molecule is toxic or not.

### Rosenblatts Perceptron

Let us tackle this problem with the same concepts in mind that we defined in the previous section, namely identifying a **model**, a **loss function** and an **optimisation algorithm**. 

We define the model as a linear function of the input features:

$$
\hat{f}_{\vec{w}}(\vec{x}) = \text{sign}(\left\langle \vec{w}, \vec{x} \right\rangle + b)\,,
$$

which we can again simplify to 

$$
\hat{f}_{\vec{w}}(\vec{x}) = \text{sign}(\left\langle \vec{w}, \vec{x} \right\rangle)\,,
$$

by including the bias $b$ as the last element of the weight vector $\vec{w}$ and adding 1 to the feature vector $\vec{x}$. The $\text{sign}(z)$ function returns 1 if $z > 0$, -1 otherwise. Therefore, the model assigns the class 1 to all data points $\vec{x}$ for which $\left\langle \vec{w}, \vec{x} \right\rangle > 0$ and the class -1 to all data points for which $\left\langle \vec{w}, \vec{x} \right\rangle \leq 0$.

Secondly, we need to define a loss function that quantifies how well the model performs. A natural choice is

$$
\mathcal{L} = \sum_{i=1}^N \max(0, - y_i \left\langle \vec{w}, \vec{x}_i \right\rangle)\,,
$$

where $\max(0,z)$ returns $0$ if $z \leq 0$ and $z$ otherwise. Let us contemplate for a moment why this loss function is a good choice. Consider the case where we have a input $\vec{x}_i$ that is misclassified. This can happen in two ways:

1. $\left\langle \vec{w}, \vec{x}_i \right\rangle > 0$ and $y_i = -1$.
2. $\left\langle \vec{w}, \vec{x}_i \right\rangle \leq 0$ and $y_i = 1$.

In both cases, the $i$-th part of the loss function is greater than 0. On the other hand, if the data point is correctly classified, the loss function for this data point is 0, which is exactly what we want. Consequently, we can only achieve a minimal loss of 0 if the data points are *linearly separable*, i.e. if there exists a hyperplane that completely separates the two classes.

Finally, we need to define an optimisation algorithm that finds the weights $\vec{w}$ that minimise the loss function. The simplest algorithm is the *gradient descent*, which we have already seen in the first section. However, we have to be careful here, as the loss function is not differentiable when $\left\langle \vec{w}, \vec{x}_i \right\rangle = 0$. We can, however, compute *subgradients* of the non-zero parts of the loss function:

$$
\nabla_{\vec{w}} \left( -y_i \left\langle \vec{w}, \vec{x}_i \right\rangle \right) = -y_i \vec{x}_i\,,
$$

that give rise to the following update rule:

$$
\vec{w} \leftarrow \vec{w} + \eta y_i \vec{x}_i\,,
$$

where $\eta$ is the learning rate. Put together, we need to iterate over all misclassified data points (loss is greater than 0) and update the weights according to the update rule. Remember that we can easily check if a data point $\vec{x}_i$ is misclassified by comparing the prediction $\hat{f}_{\vec{w}}(\vec{x}_i)$ with the true label $y_i$. 

### Inspecting the Data

Before we implement the algorithm, we need to load the data that we want to classify. We will use data from the same study as in the previous section, which you can <a href="../codes/05-machine_learning/aptamer_classification_data.csv" download>download here</a>. We again use the [`pandas`](https://pandas.pydata.org/) library to import and inspect the data.

```python
{{#include ../codes/05-machine_learning/classification.py:load_data_from_csv}}
```

| Index | PC1 | PC2 | lambda_abs_continuous | lambda_abs_class |
|-------|-----|-----|----------------------|------------------|
| 0 | 0.0515 | -1.7994 | 400.0 | -1 |
| 1 | 0.7790 | -1.7741 | 400.0 | -1 |
| 2 | 0.7006 | -1.4832 | 400.0 | -1 |
| 3 | 0.8047 | -1.4448 | 400.0 | -1 |
| 4 | 0.6912 | -1.5818 | 400.0 | -1 |

You will probably not believe it at first, but this data set represents the same set of molecules as in the previous section. However, now each molecule (row) is not characterized by a set of physical properties (columns), like molecular weight, number of rotatable bonds, and number of aromatic rings, but have been projected onto two abstract features, `PC1` and `PC2`. You have seen in the chapter on [Principal Component Analysis](../04-evd_and_svd/04-principal_component_analysis.md) what the background of this projection is, and we will revisit it in more detail in the next chapter. The last column, `lambda_abs_class`, is the binary label that denotes whether the molecule absorbs light in the visible spectrum ($\lambda_{\text{abs}} > 415$ nm), and is the target variable that we want to predict. To use the data from the dataframe for our machine learning model, we again transform the features and labels into NumPy arrays.

```python
{{#include ../codes/05-machine_learning/classification.py:process_data}}
```

From visualizing the data in the two-dimensional feature space, and coloring the data points according to the class label (by setting the `c` parameter of the `scatter` function to the labels `y`), we can see that the data is indeed linearly separable. We can therefore expect to be able to train a Rosenblatt Perceptron that successfully separates molecules that absorb light in the visible spectrum from those that do not.

```python
{{#include ../codes/05-machine_learning/classification.py:plot_data}}
```

![Data points in the feature space](../assets/figures/05-machine_learning/classification_data.svg)

### Implementing Rosenblatts Perceptron

```python
{{#include ../codes/05-machine_learning/classification.py:classification_class_init}}
```

```python
{{#include ../codes/05-machine_learning/classification.py:classification_class_fit}}
```

```python
{{#include ../codes/05-machine_learning/classification.py:classification_class_predict}}
```

```python
{{#include ../codes/05-machine_learning/classification.py:fit_model}}
```

```python
{{#include ../codes/05-machine_learning/classification.py:calculate_accuracy}}
```

```python
{{#include ../codes/05-machine_learning/classification.py:plot_decision_boundary}}
```

![Decision boundary](../assets/figures/05-machine_learning/classification_decision_boundary.svg)

