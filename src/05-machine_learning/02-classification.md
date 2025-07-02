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

### Implementating Rosenblatts Perceptron


<!-- 

Was Sie bei genauerer Betrachtung des Plots sicherlich bemerkt haben, ist, dass die Labels der Datenpunkte diskrete
Werte zwischen 3 und 8 annehmen. Im obigen Beispiel haben wir diese Labels als kontinuierliche Werte
interpretiert, um die Regression durchzuführen. In der Praxis ist es jedoch oft sinnvoller, solche 
Probleme als *Klassifikation* zu betrachten, bei denen die Labels diskrete Klassen repräsentieren.
Labels, die zwar diskret, aber geordnet sind, wie im obigen Beispiel, bezeichnet man übrigens als *ordinal*.

Das abstrakte Ziel der Klassifikation besteht darin, eine Funktion 
$\hat{f}_{\theta} : \mathbb{R}^n \to \{1,2,\dots,K\}$ zu finden, die die Eingabe auf eine von $K$ Klassen 
abbildet. Wir werden uns im Folgenden auf den einfachsten Fall der Klassifikation beschränken, die *binäre 
Klassifikation*, bei denen die Daten zwei möglichen Klassen zugeordnet werden können. Ein Beispiel für ein 
solches Problem ist die Vorhersage, ob ein Molekül synthetisierbar ist oder nicht.

```admonish note title="Binäre Klassifikation mit K > 2 Klassen"
Für den Fall, dass die Daten mehr als zwei Klassen aufweisen, kann binäre Klassifikation ebenfalls 
angewendet werden. Dazu gibt es zwei gängige Verfahren:

1. **One-vs-All**: Hierbei wird für jede Klasse ein binäres Modell trainiert, das die Daten dieser Klasse 
   von den anderen Klassen unterscheidet. Die Vorhersage erfolgt dann durch das Modell, das die höchste 
   Wahrscheinlichkeit für die gegebene Eingabe liefert.

2. **One-vs-One**: Hierbei wird für jede mögliche Kombination von zwei Klassen ein binäres Modell trainiert,
    das die Daten dieser beiden Klassen voneinander unterscheidet. Die Vorhersage erfolgt dann durch das Modell,
    das die meisten Stimmen für die gegebene Eingabe erhält.

```

In [Übung 4](../psets/04.md), als wir die Gesichter von zwei Personen auf zwei Hauptkomponenten
(Eigenfaces) projiziert haben, ist uns bereits eine Darstellung der Daten in $\mathbb{R}^2$ begegnet, 
die durch eine Gerade getrennt werden könnte. Basierend auf einer solchen *Entscheidungsgrenze*
(engl. *decision boundary*) wollen nun wir nun ein Modell trainieren, welches möglichst alle Datenpunkte korrekt
klassifiziert und auch für neue, unbekannte Datenpunkte eine korrekte Vorhersage treffen kann. Sie können dazu Ihre 
Implementation aus der Übung verwenden oder die Daten 
der Gesichter <a href="../codes/05-machine_learning/eigenfaces_pca.csv" download>hier</a>
herunterladen, wobei die dritte Spalte die Labels $y \in \{-1, 1\}$ der Personen enthält.

Wir zeigen zunächst anhand eines Negativbeispiels, wie eine solche lieare Entschiedungsgrenze zustande kommen kann.
Dazu interpretieren wir die Labels wie im obigen Beispiel als kontinuierliche Werte und führen lineare
Regression durch. Die (kontinulierlichen) Vorhersagen des Modells könnten dann als Klassen interpretiert werden, indem
wir die Werte auf die nächstgelegene Klasse abbilden. Da die Klassen hier durch -1 und 1 repräsentiert werden,
bilden wir die Vorhersagen auf die Klasse ab, die dem Vorzeichen der Vorhersage entspricht. Das Modell hat also die Form

$$
\hat{f}_ {\theta}(\vec{x}_ i) = 
\begin{cases}
    1 & \text{falls } \left\langle \vec{w}, \vec{x}_ i \right\rangle + b > 0 \\
    -1 & \text{sonst}
\end{cases}\,.
{{numeq}}{eq:linear_classification_model}
$$

Wir nutzen dazu im Grunde den gleichen Code wie im obigen Beispiel:

~~~admonish note title="Code" collapsible=true
```python
{{#include ../codes/05-machine_learning/eigenfaces_regression.py:eigenfaces_regression}}
```
~~~

Zusätzlich haben wir die Entscheidungsgrenze des Modells 

$$
    \left\langle \vec{w}, \vec{x} \right\rangle + b = 0
$$

in den Plot eingefügt, die durch die gestrichelte Linie dargestellt wird: 

![Eigenfaces Regression](../assets/figures/05-machine_learning/eigenfaces_regression.svg)

Wenn Sie allerdings den Plot von *oben* betrachten, werden Sie feststellen, dass die Entscheidungsgrenze nicht 
optimal ist, da sie nicht alle Datenpunkte korrekt klassifiziert. Die Methode der linearen Regression durch 
Minimierung der quadratischen Fehler ist also nicht geeignet, um Klassifikationsprobleme zu lösen.

#### Rosenblatt-Perzeptron

Anstatt die Labels als kontinuierliche Werte zu interpretieren und eine lineare Regression durchzuführen,
wäre es sinnvoller, die Gewichte des Models zu lernen, indem die Anzahl der falsch klassifizierten Datenpunkte 
minimiert wird. Ein Modell, welches die Daten nach einer linearen Projektion auf das Vorzeichen der
Vorhersage abbildet, wird auch als *Perzeptron* bezeichnet.
Dem Perzeptron liegt ein einfacher Algorithmus zugrunde, der die Datenpunkte iterativ durchgeht und die Gewichte 
$\vec{w}$ und $b$ anpasst, wenn ein Datenpunkt falsch klassifiziert wurde. Betrachten wir dazu die folgende 
Verlustfunktion:

$$
\mathcal{L} = - \sum_{i \in \mathcal{M}} y_i (\left\langle \vec{w}, \vec{x}_ i \right\rangle + b) \,,
{{numeq}}{eq:rosenblatt_loss}
$$

wobei die Summe über die Menge $\mathcal{M}$ der **falsch** klassifizierten Datenpunkte läuft. Dabei erinnern wir uns 
daran, dass wir überprüfen können, ob ein Datenpunkt falsch klassifiziert wurde, indem wir die Vorhergesage gemäß
Gleichung {{eqref: eq:linear_classification_model}} mit dem tatsächlichen Label vergleichen. Werden alle Datenpunkte 
von unserem Modell $\hat{f}_{\theta}$ korrekt klassifiziert, so ist $\mathcal{L} = 0$ und damit minimal. Für den
Fall, dass ein (oder mehrere) Datenpunkt(e) falsch klassifiziert wurde(n), gibt es zwei Mögkichkeiten:

1. $\left\langle \vec{w}, \vec{x}_i \right\rangle + b > 0$ und $y_i = -1$.

2. $\left\langle \vec{w}, \vec{x}_i \right\rangle + b < 0$ und $y_i = 1$.

In beiden Fällen ist die Verlustfunktion größer als Null, was bedeutet, dass die Gewichte $\vec{w}$ und $b$ so
angepasst werden müssen, dass der Fehler minimiert wird. Mit ein wenig linearer Algebra kann zudem gezeigt werden, 
dass $\mathcal{L}$ dann proportional zur Distanz des falsch klassifizierten Datenpunkts zur Entscheidungsgrenze ist.
Das einfachste, aber auch effektivste Verfahren, um die Verlustfunktion zu minimieren, ist das *Gradientenabstiegsverfahren*, welches Sie bereits in Abschnitt [(1.3)](../01-regression/03-numerical_optimisation.md) kennengelernt haben. Dazu benötigen wir den Gradienten der Verlustfunktion nach den Gewichten $\vec{w}$ und dem Bias $b$:

$$
\begin{aligned}
\nabla_{\vec{w}} \mathcal{L} &= - \sum_{i \in \mathcal{M}} y_i \vec{x}_i\,, \\
\nabla_b \mathcal{L} &= - \sum_{i \in \mathcal{M}} y_i\,.
\end{aligned}
$$

Im Gegensatz zu bisherigen Verfahren, verwenden wir zur Aktualisierung der Gewichte und des Bias jedoch nicht 
den gesamten Gradienten, sondern lediglich den Gradienten für einen einzelnen Datenpunkt.
Dies wird als *Stochastisches Gradientenabstiegsverfahren* (engl. *Stochastic Gradient Descent*, SGD) bezeichnet,
und kann insbesondere bei hochdimensionalen Daten die Rechenzeit erheblich reduzieren, sowie die Konvergenz
beschleunigen.
Das heißt, dass die Gewichte und der Bias in jedem Schritt des Gradientenabstiegs 
für jeden falsch klassifizierten Datenpunkt angepasst werden:
$$
\begin{aligned}
\vec{w} &\leftarrow \vec{w} + \frac{\tau}{N} y_i \vec{x}_ i\,, \\
b &\leftarrow b + \frac{\tau}{N} y_i\,,
\end{aligned}
$$

wobei $\tau$ die Lernrate ist, die die Schrittweite des Gradientenabstiegs bestimmt. Hat das Modell alle Datenpunkte 
einmal durchlaufen, so nennen wir dies eine *Epoche*. Dieser Algorithmus, 
der auch als *Rosenblatt Perzeptron* bekannt
ist, wird dann für eine festgelegte Anzahl von Epochen durchgeführt, oder bis alle Datenpunkte korrekt klassifiziert
wurden. Für Daten, die durch eine Gerade bzw. eine (Hyper-)Ebene linear separierbar sind, kann bewiesen werden, dass 
der Algorithmus konvergiert und eine Entscheidungsgrenze findet, die die Daten korrekt klassifiziert. -->

