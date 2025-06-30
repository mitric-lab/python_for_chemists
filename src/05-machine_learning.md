# Machine Learning

We have finally gained enough coding skills to turn to the fascinating topic of machine learning. Although we will not cover all aspects of machine learning, we will learn enough to be able to apply it to real-world problems in chemical research.

## What is Machine Learning?

Let us first clarify some basic terminology. In recent years, **Artificial Intelligence** (AI) has become a buzzword to describe a range of new technologies that have revolutionized many fields of science and industry or simply made our lives easier. In fact, AI has been around for a long time and is a collective term for techniques that solve problems that typically require human intelligence. For example, the ability to recognize objects in images, drive a car autonomously, or understand natural language are all forms of AI. In this course, we will focus on a subset of AI called **machine learning** (ML), which tackles these problems in a very specific way. Later in the course, we will also touch on **deep learning**, a subfield of ML that uses neural networks and is mainly responsible for the recent breakthroughs in AI.

<figure>
    <center>
    <img src="./assets/figures/05-machine_learning/AI_ML_DL.svg"
         alt="AI, ML, and DL"
         width="400"\>
    <figcaption>Illustration of the relationship between AI, ML, and Deep Learning.</figcaption>
    </center>
</figure>

To understand why ML is so powerful, let's first consider how we would solve a problem without it. Take an example from chemistry: we want to estimate the solubility of a molecule in a given solvent without performing any lab experiments. Since solubility depends on factors like the presence of polar or nonpolar groups, we could write a program that checks for these groups (e.g., using `for` loops and `if/else` statements) and then increases or decreases the predicted solubility accordingly. This approach does not require examples of molecules with known solubilities, but it does require us to know the underlying rules of solubility. However, in chemistry and many other fields, we often face problems where these rules are unknown, too complex for simple logic, or too difficult to implement.

<figure>
    <center>
    <img src="./assets/figures/05-machine_learning/ML_classical_programming.svg"
         alt="Classical Algorithms vs. ML"
         width="500"\>
    <figcaption>Illustration of classical programming vs. ML.</figcaption>
    </center>
</figure>

ML, on the other hand, does not rely on a predefined set of rules. Instead, it *learns* to make predictions from a set of examples where the answers are already known. This means the ML algorithm independently extracts the underlying patterns from the data, rather than being explicitly programmed with a set of rules. The learned relationships, which map the input (e.g., the molecule) to the output (e.g., the solubility), are captured in *models*. These models have many parameters that can be flexibly adapted to the data. In the next sections, we will explore a programming style for creating and applying such models.

From the above example, we see that ML does not require us to know the underlying rules of solubility. However, it does require a (preferably large) set of examples where the input and output are known. This is why ML has made such breakthroughs in fields like natural language processing and image recognition, where vast amounts of data are readily available on the internet.

## Types of Machine Learning

One typically distinguishes between two types of ML:

- **Supervised Learning**: This is the scenario described above with the solubility example. Here, the algorithm is given data consisting of inputs and their corresponding outputs. The algorithm then learns to map the inputs to the outputs.

- **Unsupervised Learning**: Here, the algorithm is given data consisting of inputs only, without any corresponding outputs. The algorithm then learns to group the inputs into meaningful clusters. An example of this is identifying groups of molecules with similar properties.

Sometimes, other types of ML are distinguished, like **reinforcement learning** or **semi-supervised learning**. We will not cover these in this course.

```admonish tip title="Tip"
Try to think about the algorithms we have already covered in this course and whether they fit the classical programming approach or the machine learning approach (supervised or unsupervised). You will find that some techniques can be classified as machine learning, so you can already claim to have some knowledge of ML. Congratulations!
```

