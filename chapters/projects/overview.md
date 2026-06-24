---
numbering:
  headings: false
---

# Overview 

Below is an overview of the projects, which are designed to help you apply and deepen the knowledge gained in the course. Please note that **these project outlines are preliminary** and meant to assist you as you begin your work. Stay adaptable and be ready to adjust your approach as you gain new understanding during the project.

## Fitting absorbance data for the decay of tris(oxalato)manganate(III)

In the physical chemistry practical course, you used a linearised model of the absorbance
decay of tris(oxalato)manganate(III) to determine the reaction rate constant. In this
project, you will return to the original data and fit it directly with a more realistic
nonlinear model, including a background correction.

**Related knowledge:**
- Numerical optimisation (Sec.&nbsp;{ref}`sec:computational_optimisation`)
- Nonlinear least-squares fitting (Sec.&nbsp;{ref}`sec:regression`)

**Suggested project guide:**

A possible starting point is the dataset below. You are also welcome to use your own
measurements from the practical course if they are available.

| Temperature | Dataset |
| --- | --- |
| 20 °C | {download}`mn_decay_20.txt </assets/data/regression/mn_decay_20.txt>` |
| 30 °C | {download}`mn_decay_30.txt </assets/data/regression/mn_decay_30.txt>` |
| 40 °C | {download}`mn_decay_40.txt </assets/data/regression/mn_decay_40.txt>` |
| 50 °C | {download}`mn_decay_50.txt </assets/data/regression/mn_decay_50.txt>` |

One useful baseline would be the exponential fit from Exercise 1.2A,

$$
A(t) = A_0 \exp(-kt),
$$

applied separately to each temperature. As in Exercise 1.2A, you can formulate the
regression as an optimisation problem by defining a loss function and minimising it
with respect to the model parameters. This gives one rate constant $k$ for each
measurement series.

Afterwards, you could compare this baseline with models that account for an
absorbance offset. For example, if the background measurement was not done perfectly,
the data may be better described by

$$
A(t) = (A_0 - A_\infty) \exp(-kt) + A_\infty,
$$

where $A_0$ is the fitted absorbance at $t = 0$ and $A_\infty$ is fitted separately
for each temperature. Another possible correction model is a double-exponential decay,
which can describe cases where the absorbance change contains two kinetic
contributions:

$$
A(t) = C_1 \exp(-k_1 t) + C_2 \exp(-k_2 t).
$$

If a residual background seems relevant, the same idea can be extended by adding an
optional offset term $A_\infty$.

You are free to propose and test your own correction models if you can motivate them
chemically or experimentally. For example, one could test a common-background model
with one shared $A_\infty$ for all temperatures, or fit all curves globally while
constraining the rate constants by the Arrhenius relation
$k(T) = A \exp(-E_\mathrm{a}/RT)$.

For the final analysis, it would be interesting to estimate the activation energy from
the Arrhenius equation,

$$
k = A \exp\left(-\frac{E_\mathrm{a}}{RT}\right),
$$

using the rate constants obtained from the different fitting approaches. Comparing the
resulting activation energies can help you discuss how sensitive the kinetic
interpretation is to the chosen correction model. If you use the double-exponential
model, it is worth discussing whether one of the two fitted rate constants can be
assigned to the main decay before including it in an Arrhenius analysis.

## Predicting TLC mobile phase compositions for organic molecules

Choosing a suitable mobile phase for thin-layer chromatography (TLC) is essential for
successful separations, but it can be difficult without practical experience. In this
project, you will curate a small dataset of organic molecules, train a simple regression
model to predict suitable TLC mobile phase compositions, and evaluate its performance
on a test set.

**Related knowledge:**
- Regression models (Sec.&nbsp;{ref}`sec:regression`)
- Neural networks (Sec.&nbsp;{ref}`sec:multi_layer_neural_networks`)
- Molecular representations (Sec.&nbsp;{ref}`sec:molecular_representations`)

**Suggested project guide:**

Start by treating mobile-phase prediction as a supervised machine learning task and define the problem statement carefully, since the data collection strategy depends on this choice. In particular, specify which solvent systems are considered and what the input $\vec{x}$ and target $\vec{y}$ should represent. Once the predictive task is fixed, search the literature for relevant datasets and protocols. Extracting and postprocessing these data can be time consuming, so use tools such as `pandas` for data management and `rdkit` for chemical data handling.

Similar studies show that both molecular representation and model architecture can strongly affect predictive performance.[^1] To make a meaningful comparison, evaluate all models on data that was not used for training, for example using nested cross-validation. Metrics such as $R^2$ help assess how much of the target variance is captured by the predictions.

[^1]: https://doi.org/10.3390/molecules26092474

## Classifying chemical samples from spectra

Spectra often contain enough information to distinguish between different types of
chemical samples. In this project, you will train a simple classification model to assign
spectra to predefined classes of molecules and evaluate how reliably the model recognises new samples.

**Related knowledge:**
- Classification models (Sec.&nbsp;{ref}`sec:classification`)
- Neural networks (Sec.&nbsp;{ref}`sec:multi_layer_neural_networks`)
- Molecular representations (Sec.&nbsp;{ref}`sec:molecular_representations`)

**Suggested project guide:**

Start by defining the classification task, including the type of spectra you want to use, such as IR, NMR, simulated spectra, or experimental spectra. A rich data source is the [NIST Chemistry WebBook](https://webbook.nist.gov/chemistry/). Choose the class labels carefully, for example molecular classes, functional groups, compound identities, or sample types, because this choice determines the structure of the dataset.

Next, decide how the spectra should be represented numerically, for example as raw intensities on a common grid, binned spectra, extracted peak positions, or peak intensities. Apply preprocessing steps such as baseline correction, normalisation, smoothing, or standardisation consistently, and avoid leaking information from the test set into the training process. Suitable models include support vector machines, random forests, and neural networks with a softmax output layer. Evaluate the final model on unseen data using metrics such as accuracy, precision, recall, F1 score, and, where useful, a confusion matrix. Discuss possible pitfalls such as overfitting, underfitting, class imbalance, and overly optimistic results caused by data leakage.

## Solubility prediction for organic molecules

The solubility of organic molecules is important for synthesis, purification, and
formulation, but it is often hard to estimate from structure alone. In this project,
you will collect or use a small molecular dataset, build a regression model for
aqueous solubility, and evaluate which molecular features are most useful for
prediction.

**Related knowledge:**
- Regression models (Sec.&nbsp;{ref}`sec:regression`)
- Neural networks (Sec.&nbsp;{ref}`sec:multi_layer_neural_networks`)
- Molecular representations (Sec.&nbsp;{ref}`sec:molecular_representations`)

**Suggested project guide:**

Treat solubility prediction as a supervised regression task in which molecular structure is used to predict an experimentally measured solubility value. Use the literature as a starting point, since numerous studies and benchmark datasets exist for aqueous solubility prediction of organic molecules.[^2] Choose suitable molecular representations, such as hand-crafted descriptors, fingerprints, or learned embeddings, and compare how strongly this choice affects model performance.

Train and compare regression models such as linear regression, random forests, kernel methods, or neural networks, while keeping a separate test set for unbiased evaluation. Assess predictive performance with metrics such as mean squared error, mean absolute error, and $R^2$. In addition to prediction accuracy, focus on explainability: analyse which descriptors, substructures, or functional groups contribute most strongly to the predicted solubility, and discuss whether the model decisions are chemically plausible. For example, polar groups, hydrogen-bond donors and acceptors, charge, molecular size, and hydrophobic fragments should influence the prediction in chemically meaningful ways.

[^2]: https://doi.org/10.1021/acs.jcim.4c02399

## Inferring ionic solution composition

Conductivity curves for dilution series contain information about the ions present in
a solution and their concentrations. In this project, you will generate synthetic
conductivity curves and use them to infer a plausible composition of an unknown ionic
solution.

**Related knowledge:**
- Numerical optimisation (Sec.&nbsp;{ref}`sec:computational_optimisation`)
- Electrolyte conductivity models

**Suggested project guide:**

As a starting point, it is useful to collect simple conductivity models for strong and
weak electrolytes. For a strong electrolyte, the molar conductivity often follows
Kohlrausch's law at low concentrations,

$$
\mathit{\Lambda}_\mathrm{m}(c) = \mathit{\Lambda}_\mathrm{m}^0 - K \sqrt{c},
$$

where $\mathit{\Lambda}_\mathrm{m}^0$ is the limiting molar conductivity and $K$ is an
empirical constant. The conductivity is then

$$
\kappa(c) = c \mathit{\Lambda}_\mathrm{m}(c).
$$

For a weak electrolyte, the degree of dissociation $\alpha$ changes with concentration.
Using Ostwald's dilution law,

$$
K_\mathrm{a} = \frac{c \alpha^2}{1 - \alpha},
$$

and

$$
\mathit{\Lambda}_\mathrm{m}(c) = \alpha \mathit{\Lambda}_\mathrm{m}^0,
\qquad
\kappa(c) = c \mathit{\Lambda}_\mathrm{m}(c).
$$

In an identification task, however, the molar mass of the unknown compound is not
known. Therefore, it is often more practical to generate and analyse conductivity as a
function of mass concentration $\rho$, for example $\kappa(\rho)$ with $\rho$ in g/L.

One possible first goal is to develop a method that identifies a single unknown
compound from its conductivity-mass concentration curve. This could be approached by
generating synthetic reference curves for a list of candidate compounds, interpolating
them onto a common concentration grid, and comparing a synthetic unknown curve to each
reference curve using a distance such as a mean squared error. You could also add a
small amount of random noise to make the identification task more realistic.
Good starting candidates might include common, water-soluble salts such as NaCl, KCl,
LiCl, NaNO$_3$, KNO$_3$, NH$_4$Cl, MgCl$_2$, CaCl$_2$, Na$_2$SO$_4$, and MgSO$_4$.
If you also want weak or medium-strength electrolytes in the comparison, acetic acid,
ammonia, or ammonium acetate could provide useful contrast.

For sufficiently dilute solutions, conductivities can often be treated as
approximately additive. This suggests an extension to unknown mixtures containing only
a few compounds. If the total mass concentration is $\rho$ and $w_i$ is the mass
fraction of compound $i$, a simple mixture model is

$$
\kappa_\mathrm{mix}(\rho) \approx \sum_i \kappa_i(w_i \rho),
\qquad
w_i \ge 0,
\qquad
\sum_i w_i = 1.
$$

The inverse problem is then to find a small set of compounds and mass fractions that
reproduces the synthetic target curve. This could be formulated as an optimisation
problem. Because the mixture is expected to contain only a few compounds, some form of
sparsity constraint or penalty is useful. A one-norm penalty is often used for
sparsity, but under the constraint $\sum_i w_i = 1$ it becomes constant and therefore
does not help much. Alternatives include explicitly limiting the number of non-zero
components, testing all one-, two-, or three-compound combinations, or using a
non-convex $p$-pseudonorm penalty with $0 < p < 1$.

As an additional data-driven extension, you could train a classification model that
uses conductivity curves to distinguish strong, medium, and weak electrolytes. This
could be done by sampling each curve on a fixed mass-concentration grid and using the
sampled conductivities, slopes, or fitted model parameters as features.
