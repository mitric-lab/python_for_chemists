---
downloads: []
---

# Data-driven Modelling

In this chapter, we shift our focus from strictly rule-based modelling to a data-driven modelling approach. Within *machine learning*, flexible models featuring many adjustable parameters are used to uncover patterns and relationships in large datasets. The complexity of these models ranges continuously from basic linear regression to deep neural networks; the latter have made a significant impact in many applications in recent years.

We cover supervised learning techniques (regression and classification) and unsupervised learning techniques (dimensionality reduction and clustering), applied to real-world chemical data sets. We further introduce powerful libraries for handling and structuring large datasets and explore new concepts for object-oriented programming.

The accompanying problem set for this chapter is available at {ref}`sec:pset_3`.

| Section | Covered Examples | New Concepts and Tools |
| ------- | ---------------- | ---------------------- |
| Regression | Lambert-Beer's law | Loss functions, `pandas`, `np.polyfit`, cross-validation, overfitting |
| Classification | Classification of wines | Logistic regression, Support Vector Machines, object-oriented programming, Stochastic gradient descent, `scikit-learn` |
| Dimensionality Reduction | PAH-solvent conformers | Covariance matrix, explained variance |
| Clustering | PAH-solvent conformers | `np.random.choice`, `np.argmin`, `np.linalg.norm`, adding an extra axis to an array |