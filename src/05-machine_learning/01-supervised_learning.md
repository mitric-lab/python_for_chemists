## Überwachtes Lernen

Liegen uns Datenpaare $\{(\vec{x}_i, y_i)\}_{i=1,2,\dots,N}$ vor, wobei $\vec{x}_i \in \mathbb{R}^n$ 
die Eingabe und $y_i \in \mathbb{R}$ die Ausgabe ist, wollen wir eine approximative
Abbildung $\hat{f}_{\theta} : \mathbb{R}^n \to \mathbb{R}$ finden, die die Eingabe auf die Ausgabe abbildet, 
d.h. $\hat{f}_{\theta}(\vec{x}_i) \approx y_i$ für alle $i$. Man bezeichnet $y_i$ auch als *Labels* oder 
*Targets*. Unser Modell $\hat{f}$ ist in 
der Regel parametrisiert, d.h. es gibt eine Menge von Parametern $\theta$, die die Abbildung definieren. 
Demnach ist unser Ziel, die Parameter $\theta$ so zu wählen, dass der Fehler zwischen der
tatsächlichen Ausgabe $y_i$ und der approximierten Ausgabe $\hat{f}_{\theta}(\vec{x}_i)$ minimiert wird.
Dieses Lernen wird als **überwachtes Lernen** bezeichnet, da wir $y_i$ für jede Eingabe 
$\vec{x}_i$ kennen und somit den Fehler direkt berechnen können.

Gemäß den obigen Überlegungen können wir festhalten, dass jedes ML-Modell, das auf überwachtem Lernen basiert,
die folgenden Komponenten enthält:

1. **Modell**: Die Funktion $\hat{f}_{\theta}$, die die Eingabe auf die Ausgabe abbildet.
2. **Verlustfunktion**: Die Funktion, die den Fehler zwischen der tatsächlichen Ausgabe und der approximierten 
Ausgabe misst.
3. **Optimierungsverfahren**: Der Algorithmus, der die Parameter $\theta$ des Modells so anpasst, dass der Fehler 
minimiert wird.

Dieses Setting lässt sich tatsächlich auch auf das unüberwachte Lernen übertragen, wobei die Verlustfunktion
durch eine allgemeine Objektivfunktion ersetzt wird, die es zu maximieren oder minimieren gilt.

### Regression

Für den Fall, dass die Labels $y_i$ kontinuierlich sind, also reelle Zahlen $y_i \in \mathbb{R}$ annnehmen, 
spricht man von **Regression**. Den einfachsten Fall, die **lineare Regression**, und mit $x_i \in \mathbb{R}$,
haben Sie bereits ausführlich in Kapitel [(1.2)](../01-regression/02-linear_regression.md) kennengelernt und 
diskutiert. Wir werden daher nur kurz darauf eingehen und die lineare Regression in den Kontext des überwachten 
Lernens einordnen.

Zunächst möchten wir anmerken, dass wir bisher lediglich Inputs $x_i$ betrachtet haben, die nur eine Dimension 
besitzen, also $x_i \in \mathbb{R}$ für $i = 1,2,\dots,N$. Dies stellten beispielsweise die Konzentrationen von 
Methylenblau, welche auf die Absorption abgebildet wurden, dar. In der Praxis haben wir jedoch oft mit 
mehrdimensionalen Inputs $\vec{x}_i \in \mathbb{R}^n$ zu tun, welche $n$ *Features* besitzen. In den vorherigen
Abschnitten haben wir zur Darstellung dieser Daten die Datenmatrix $\bm{X}$ kennengelernt, welche die Inputs 
$\vec{x}_i$ in den Zeilen
und die Features in den Spalten speichert. Die lineare Regression für mehrdimensionale Inputs lautet dann:

$$
\hat{f}(\vec{x}_i) = w_1 x_{i1} + w_2 x_{i2} + \dots + w_n x_{in} + b = \sum_{j=1}^n w_j x_{ij} + b\,,
$$

was wir als Skalarprodukt zwischen dem Vektor $\vec{w} = (w_1, w_2, \dots, w_n)^T$ und dem Inputvektor $\vec{x}_i$
schreiben können:

$$
\hat{f}_{\theta}(\vec{x}_i) = \left\langle \vec{w}, \vec{x}_i \right\rangle + b =\vec{w}^T \vec{x}_i + b\,.
$$

Unter Berücksichtigung aller Datenpunkte durch die Matrix $\mathbf{X}$ ergeben sich die Vorhersagen 
$\hat{\vec{y}}$ aller $N$ Inputs durch eine Matrix-Vektor-Multiplikation:

$$
\hat{\vec{y}} = \mathbf{X} \vec{w} + b\,.
$$

Da wir die Vorhersage des Modells als gewichtete Summe der Features berechnen, bezeichnet man die 
Parameter $\theta = (\vec{w}, b)$ des Modells auch als Gewichte $\vec{w}$ und Bias $b$. 
Das $j$-te Gewicht $w_j$ gibt dabei an, wie stark das $j$-te Feature $x_{ij}$ in die Vorhersage eingeht.
Als Verlustfunktion wählt man in der Regel die Summe der quadratischen Fehler, welche Sie bereits kennen.

```admonish note title="Hinweis"
Um die Modellparameter $\theta = (\vec{w}, b)$ zusammenzufassen, kann $b$ auch als ein zusätzliches Gewicht 
$w_0$ eingeführt werden, sodass $\vec{w} = (w_0, w_1, w_2, \dots, w_n)^T$. Dazu muss der Inputvektor um
einen konstanten Wert $1$ erweitert werden, sodass $\vec{x}_i = (1, x_{i1}, x_{i2}, \dots, x_{in})^T$. 
Überprüfen Sie, dass die lineare Regression in diesem Fall äquivalent zur obigen Formulierung ist.
```

Gemäß den oben diskutierten Komponenten eines ML-Modells, fehlt nun noch die Angabe eines Optimierungsverfahrens, 
um die Modellparameter $\theta$ so zu wählen, dass der Fehler zwischen den Vorhersagen $\hat{\vec{y}}$ und den
tatsächlichen Labels $\vec{y}$ minimiert wird. In der ersten Übung haben Sie gesehen, dass dieses Problem 
eine analytische Lösung besitzt. Auch für den multi-dimensionalen Fall existiert eine analytische Lösung,

$$
\hat{\theta} = (\mathbf{X}^T \mathbf{X})^{-1} \mathbf{X}^T \vec{y}\,,
$$

was Sie leicht durch Berechnen und Nullsetzen des Gradienten der Verlustfunktion nach $\theta$ überprüfen können. 

```admonish note title="Eindeutigkeit der Lösung"
Die lineare Algebra besagt, dass die Matrix $\mathbf{X}^T \mathbf{X}$ invertierbar ist, wenn die Spalten von
$\mathbf{X}$ linear unabhängig sind. Für den Fall $N > n$ ist dies sehr wahrscheinlich, und wir sprechen von
*unabhängigen Features*.

Ist die Matrix nicht invertierbar, so liefert uns die Moore-Penrose-Pseudoinverse 

$$
(\mathbf{X}^T \mathbf{X})^+ = \mathbf{V} 
\left(\begin{array}{cccc}
\frac{1}{\sigma_1^2} & & & \\
& \ddots & & \\
& & \frac{1}{\sigma_k^2} & \\
& & & 0_{n-k, n-k}
\end{array}\right)
\mathbf{V}^T\,,
$$

die Lösung mit minimaler Norm. 
```

### Regression am Wine Quality Dataset
