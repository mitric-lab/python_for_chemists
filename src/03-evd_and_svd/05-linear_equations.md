## Lineare Gleichungssysteme

Viele Probleme in der Naturwissenschaft, Technik und darüber hinaus lassen
sich durch lineare Gleichungssysteme beschreiben. Solche Systeme lassen sich
kompakt in Matrixschreibweise darstellen als
$$
  \bm{A} \vec{x} = \vec{b}
$$
mit der Koeffizientenmatrix $\bm{A} \in \R{m}{n}$, dem Vektor der Unbekannten 
$\vec{x} \in \mathbb{R}^n$ und dem Vektor der rechten Seite 
$\vec{b} \in \mathbb{R}^m$, welcher auch als Inhomogenität bezeichnet wird. 
Die Anzahl der Gleichungen $m$ entspricht der Anzahl der Zeilen von $\bm{A}$, 
während die Anzahl der Unbekannten $n$ der Anzahl der Spalten von $\bm{A}$
entspricht.

Das Standardverfahren zur Lösung linearer Gleichungssysteme ist der
Gauß-Algorithmus. Die übliche Implementierung funktioniert allerdings 
nur für Gleichungssysteme mit eindeutiger Lösung. Für unterbestimmte
Systeme mit unendlich vielen Lösungen ist der Gauß-Algorithmus numerisch
instabil und für überbestimmte Systeme ohne Lösungen liefert er uns
keine nützlichen Informationen.

Die formale Lösung eines eindeutigen Systems mit $m = n$ ist
$\vec{x} = \bm{A}^{-1} \vec{b}$, wobei $\bm{A}^{-1}$ die Inverse von $\bm{A}$
bezeichnet. Mit exakter Arithmetik ist diese Lösung identisch zur der
aus dem Gauß-Algorithmus. Für singuläre oder gar nicht-quadratische Matrizen
ist die Inverse aber nicht definiert und diese Lösung daher nicht anwendbar.
Es wäre doch schön, wenn wir eine Operation hätten, die einige Eigenschaften
der Inversen erfüllt, aber für alle Matrizen definiert ist. Solche Operationen
nennen sich Pseudoinversen und ein bekannter Vertreter ist die
[*Moore-Penrose-Pseudoinverse*](https://de.wikipedia.org/wiki/Pseudoinverse#Die_Moore-Penrose-Inverse),
die wir im Folgenden kennenlernen wollen.

### Theoretische Grundlagen

Ist eine invertierbare Matrix $\bm{A}$ auch diagonalisierbar mit der
Eigenwertzerlegung $\bm{A} = \bm{U} \bm{\Lambda} \bm{U}^{-1}$, so kann die
ihre Inverse durch $\bm{A}^{-1} = \bm{U} \bm{\Lambda}^{-1} \bm{U}^{-1}$
berechnet werden, wobei $\bm{\Lambda}^{-1}$ die Diagonalmatrix der 
invertierten Eigenwerte ist. Die Richtigkeit dieser Formel lässt sich leicht
durch
$$
  \bm{A} \bm{A}^{-1} = \bm{U} \bm{\Lambda} \bm{U}^{-1} \bm{U} \bm{\Lambda}^{-1} \bm{U}^{-1} = \identity
$$
and
$$
  \bm{A}^{-1} \bm{A} = \bm{U} \bm{\Lambda}^{-1} \bm{U}^{-1} \bm{U} \bm{\Lambda} \bm{U}^{-1} = \identity
$$
zeigen. Da die Singulärwertzerlegung einer Verallgemeinerung der 
Eigenwertzerlegung für beliebige Matrizen ist, können wir überlegen, eine 
Verallgemeinerung der Inversen über die Singulärwertzerlegung zu definieren.

Wir definieren die *Moore-Penrose-Pseudoinverse* (MP-Pseudoinverse)
$\bm{A}^+$ einer Matrix $\bm{A} \in \R{m}{n}$ als
$$
  \bm{A}^+ = \bm{V} \bm{\Sigma}^+ \bm{U}^T
  {{numeq}}{eq:mp_pseudoinverse}
$$
mit der pseudoinversen Diagonalmatrix $\bm{\Sigma}^+$, die durch Transponieren
und Invertieren der nicht-verschwindenden Singulärwerte von $\bm{\Sigma}$
entsteht. Es kann gezeigt werden, dass die MP-Pseudoinverse die 
sog.
```admonish note title="Moore-Penrose-Bedingungen"
$$
  \begin{alignat}{2}
    \text{1.}\ \ &\mathbf{A} \mathbf{A}^+ \mathbf{A} &&= \mathbf{A} {{numeq}}{eq:mp_conditions_general_inverse} \\
    \text{2.}\ \ &\mathbf{A}^+ \mathbf{A} \mathbf{A}^+ &&= \mathbf{A}^+ {{numeq}}{eq:mp_conditions_weak_inverse} \\
    \text{3.}\ \ &\left(\mathbf{A} \mathbf{A}^+\right)^\intercal &&= \mathbf{A} \mathbf{A}^+ {{numeq}}{eq:mp_conditions_symmetry_aapt} \\
    \text{4.}\ \ &\left(\mathbf{A}^+ \mathbf{A}\right)^\intercal &&= \mathbf{A}^+ \mathbf{A} {{numeq}}{eq:mp_conditions_symmetry_apat}
  \end{alignat}
$$
```
erfüllt. Umgekehrt definieren diese Bedingungen die MP-Pseudoinverse eindeutig.
Ist die Matrix $\bm{A}$ invertierbar, so ist die MP-Pseudoinverse
identisch zur Inversen, also $\bm{A}^+ = \bm{A}^{-1}$.

Betrachten wir nun ein überbestimmtes Gleichungssystem 
$\bm{A} \vec{x} = \vec{b}$ mit $\bm{A}\in \R{m}{n}$, $m > n$ 
und $\vec{b} \in \mathbb{R}^m$, welches i.A. keine Lösung besitzt.
Dann gilt das Theorem
$$
  x_0 := \bm{A}^+ \vec{b} = \argmin{x\in \mathbb{R}^n} \| \bm{A} \vec{x} - \vec{b} \|_2\,,
$$
also liefert $x_0$ die Lösung, die den quadratischen Fehler beider Seite
des Gleichungssystems minimiert.
~~~admonish proof title="Beweis" collapsible=true
Wir starten mit der Differenzvektor $\bm{A} \vec{x} - \vec{b}$ 
für ein beliebiges $\vec{x} \in \mathbb{R}^n$ und addieren 0:

$$
  \bm{A} \vec{x} - \vec{b} = \bm{A} \vec{x} - \bm{A} \vec{x}_0 + \bm{A} \vec{x}_0 - \vec{b} = \bm{A} (\vec{x} - \vec{x}_0) + (\bm{A} \vec{x}_0 - \vec{b})
$$
und berechnen dann ihre euklidische Norm:
$$
  \begin{align}
    \|\bm{A}x - b\|_2^2 &= \left(\bm{A}(x - x_0) + \left(\bm{A}x_0 - b\right)\right)^\intercal \left(\bm{A}(x - x_0) + \left(\bm{A}x_0 - b\right)\right) \\
    &= \left(\bm{A}(x - x_0)\right)^\intercal \left(\bm{A}(x - x_0)\right) + \left(\bm{A}x_0 - b\right)^\intercal \left(\bm{A}x_0 - b\right) \\
    &\hphantom{=} + \left(\bm{A}(x - x_0)\right)^\intercal \left(\bm{A}x_0 - b\right) + \left(\bm{A}x_0 - b\right)^\intercal \left(\bm{A}(x - x_0)\right) \\
    &=: \|\bm{A}(x - x_0)\|_2^2 + \|\bm{A}x_0 - b\|_2^2 + p_1 + p_2\,,
  \end{align}
$$
wo $p_1$ und $p_2$ die Skalarprodukte zwischen den beiden Vektoren
$\bm{A}(x - x_0)$ und $(\bm{A}x_0 - b)$ sind. Diese können durch Einsetzen von
$x_0 = \bm{A}^+ b$ berechnet werden:
$$
  \begin{align}
    p_1 &= \left(\bm{A}(x - x_0)\right)^\intercal \left(\bm{A}x_0 - b\right) = \left(\left(x - \bm{A}^+b\right)^\intercal \bm{A}^\intercal\right) \left(\bm{A}\bm{A}^+b - b\right) \\
    &= x^\intercal\bm{A}^\intercal\bm{A}\bm{A}^+b 
    - x^\intercal\bm{A}^\intercal b 
    - b^\intercal\left(\bm{A}^+\right)^\intercal\bm{A}^\intercal \bm{A}\bm{A}^+b 
    + b^\intercal\left(\bm{A}^+\right)^\intercal\bm{A}^\intercal b \\
    &= x^\intercal\bm{A}^\intercal\left(\bm{A}\bm{A}^+\right)^\intercal b 
    - x^\intercal\bm{A}^\intercal b 
    - b^\intercal\left(\bm{A}^+\right)^\intercal \bm{A}^\intercal \left(\bm{A}\bm{A}^+\right)^\intercal b 
    + b^\intercal\left(\bm{A}^+\right)^\intercal\bm{A}^\intercal b \\
    &= x^\intercal\left(\bm{A}\bm{A}^+\bm{A}\right)^\intercal b 
    - x^\intercal\bm{A}^\intercal b 
    - b^\intercal \left(\bm{A}^+\right)^\intercal \left(\bm{A}\bm{A}^+\bm{A}\right)^\intercal b 
    + b^\intercal\left(\bm{A}^+\right)^\intercal\bm{A}^\intercal b \\
    &= x^\intercal\bm{A}^\intercal b 
    - x^\intercal\bm{A}^\intercal b 
    - b^\intercal \left(\bm{A}^+\right)^\intercal \bm{A}^\intercal b 
    + b^\intercal\left(\bm{A}^+\right)^\intercal \bm{A}^\intercal b \\
    &= 0\,,
  \end{align}
$$
wobei wir für die 3. Zeile Gl. {{eqref: eq:mp_conditions_symmetry_aapt}} und 
für die 5. Zeile Gl. {{eqref: eq:mp_conditions_general_inverse}} eingesetzt
wurden. 

Aus der Symmetrie des Skalarprodukts folgert man $p_2 = 0$. Der ursprüngliche
Ausdruck vereinfacht sich nun zu
$$
  \begin{align}
    \|\mathbf{A}x - b\|_2^2 &= \|\mathbf{A}(x - x_0)\|_2^2 + \|\mathbf{A}x_0 - b\|_2^2 + 0 \\
    &= \|\mathbf{A}(x - x_0)\|_2^2 + \|\mathbf{A}x_0 - b\|_2^2 \geq \|\mathbf{A}x_0 - b\|_2^2.
  \end{align}
$$
Also ist die euklidische Norm des Differenzvektors für beliebiges 
$\vec{x} \in \mathbb{R}^n$ immer größer als $\|\bm{A} \vec{x_0} - \vec{b}\|$.
Damit minimiert $x_0 = \bm{A}^+ \vec{b}$ den quadratischen Fehler beider 
Seite des Gleichungssystems.
~~~

Für ein unterbestimmtes Gleichungssystem $\bm{A} \vec{x} = \vec{b}$ 
mit $\bm{A}\in \R{m}{n}$, $m < n$ und $\vec{b} \in \mathbb{R}^m$,
die unendlich viele Lösungen besitzt, gilt das Theorem
$$
  x_0 := \bm{A}^+ \vec{b} = \argmin{\bm{A}\vec{x}=\vec{b}} \| x \|_2\,,
$$
also liefert $x_0$ die Lösung mit der kleinsten euklidischen Norm.

### Implementierung

Wir implementieren die MP-Inverse am Beispiel eines überbestimmten 
Gleichungssystems:
$$
\bm{A} = \begin{pmatrix}
  1 & 1 \\
  2 & 1 \\
  3 & 1 \\
  4 & 1 \\
\end{pmatrix}
\quad \text{und} \quad
\vec{b} = \begin{pmatrix}
  4.0 \\
  3.5 \\
  5.0 \\
  6.5 \\
\end{pmatrix}\,.
$$

Nach dem Importieren von NumPy
```python
{{#include ../codes/04-evd_and_svd/pinv_lineq.py:imports}}
```
definieren wir die Funktion `pinv`:
```python
{{#include ../codes/04-evd_and_svd/pinv_lineq.py:pinv_function}}
```
Nach der SVD der Ausgangsmatrix initiieren wir die invertierte Matrix
der Singulärwerte `s_inv` mit der Dimension der Transposition von `s`.
Danach iterieren wir über $p = \min(m,n)$ und invertieren die 
Singulärwerte, die ungleich null sind. Da die Fließkommazahlen nicht exakt
sind, setzen wir einen Schwellenwert `rcond`. Singulärwerte, die kleiner als
diese Schwelle liegen, werden als null behandelt und nur größere 
Singulärwerte werden invertiert. Am Ende gibt die Funktion die MP-Inverse
gemäß Gl. {{eqref: eq:mp_pseudoinverse}} zurück.

Nun können wir das angegebene Gleichungssystem definieren:
```python
{{#include ../codes/04-evd_and_svd/pinv_lineq.py:define_lineq}}
```
und diese mit der MP-Inversen lösen:
```python
{{#include ../codes/04-evd_and_svd/pinv_lineq.py:solve_lineq}}
```
Wir erhalten die folgenden Ergebnisse:
```python
{{#include ../codes/04-evd_and_svd/pinv_lineq.py:verify_results}}
```

Einige Einträge der Matrix $\bm{A}$ dürfte etwas komisch erscheinen, wie 
z.B. die zweite Spalte, die nur Einsen enthält. Wir schreiben diese
Matrixgleichung nun aus:
$$
  \begin{align}
    x_1 + x_2 &= 4.0 \\
    2x_1 + x_2 &= 3.5 \\
    3x_1 + x_2 &= 5.0 \\
    4x_1 + x_2 &= 6.5
  \end{align}
$$
oder mit familiären Variablennamen
$$
  \begin{align}
    m + t &= 4.0 \\
    2m + t &= 3.5 \\
    3m + t &= 5.0 \\
    4m + t &= 6.5 \,.
  \end{align}
$$

Es soll jetzt unschwer erkennbar sein, dass dieses Gleichungssystem
das Problem darstellt, eine Gerade zu finden, die durch die 4 Punkte
$(1.0, 4.0)$, $(2.0, 3.5)$, $(3.0, 5.0)$ und $(4.0, 6.5)$ verläuft.
Mit einwenig Vorstellung wird es einem klar, dass so eine Gerade nicht
existiert. Die Lösung durch die MP-Inverse ist daher diejenige Gerade,
die den quadratischen Fehler zwischen den Punkten und der Geraden
minimiert. Falls einem diese Idee bekannt vorkommt, dann ist das kein
Zufall. Das ist gerade die lineare Regression mit Methode der kleinsten
Quadrate, die wir bereits in 
Kapitel [1.2](../01-regression/02-linear_regression.md)
kenngelernt haben. 

Die MP-Inverse liefert uns eine weitere Möglichkeit,
die Parameter der (multi-)lineare Regression zu berechnen: Für die
Datenpunkte $\{(\text{---}\ \vec{x}_ i \text{---}, y_i)\}_ {i=1,\cdots,N}$
mit $\vec{x}_ i \in \mathbb{R}^n$, $y_i \in \mathbb{R}$
sowie das Modell $y = \sum_{j=1}^n w_j x_j + b = \vec{w}^\intercal \vec{x} + b$,
sind die optimalen Parameter $\vec{w} \in \mathbb{R}^n$ und $b \in \mathbb{R}$
durch die Least-Squares-Lösung des Gleichungssystems 
$$
  \underbrace{\begin{pmatrix}
    \text{---}\ \vec{x}_ 1 \text{---} & 1 \\
    \vdots & \vdots \\
    \text{---}\ \vec{x}_ N \text{---} & 1 \\
  \end{pmatrix}}_ {\bm{A}\in \R{N}{(n+1)}}
  \underbrace{\begin{pmatrix}
    \vert \\
    \vec{w} \\
    \vert \\
    b
  \end{pmatrix}}_ {\vec{\beta} \in \mathbb{R}^{n+1}}
  = \underbrace{\begin{pmatrix}
    y_1 \\
    \vdots \\
    y_N
  \end{pmatrix}}_ {\vec{y} \in \mathbb{R}^N}
$$
gegeben, also 
$$
\vec{\beta} = \bm{A}^+ \vec{y}\,.
{{numeq}}{eq:multilinear_regression_pinv}
$$
