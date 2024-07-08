## Singulärwertzerlegung



### Definition

Im Sinne der Eigenwertzerlegung für quadratische Matrizen definieren wir die 
folgende Zerlegung für eine beliebige Matrix $\bm{A} \in \C{m}{n}$:
$$
  \bm{A} = \bm{U} \bm{\Sigma} \bm{V}^\dag
$$
wobei $\bm{U} \in \C{m}{m}$ und $\bm{V} \in \C{n}{n}$ unitäre Matrizen sind 
und $\bm{\Sigma} \in \R{m}{n}$ eine Diagonalmatrix mit nicht-negativen
Diagonaleinträgen $\sigma_1 \geq \sigma_2 \geq \ldots \geq \sigma_p \geq 0$ ist,
wobei $p = \min(m,n)$. Diese Zerlegung wird als *Singulärwertzerlegung*
(engl. *singular value decomposition*, SVD) bezeichnet.

Aber was bedeutet es, dass eine rechteckige Matrix $\bm{\Sigma}$
diagonal ist? Streng genommen ist das nicht möglich, aber wir können die
Matrix in eine Form bringen, die der Diagonalform ähnelt, nämlich
$$
  \bm{\Sigma} = \left(\begin{array}{c|c}
    \bm{D} & \bm{0} \\ \hline
    \bm{0} & \bm{0}
  \end{array}\right)\,.
$$
Demnach enthält $\bm{\Sigma}$ eine $p \times p$ Diagonalmatrix $\bm{D}$ und
möglicherweise Nullmatrizen, um die Dimensionen zu erfüllen.

Die Diagonaleinträge $\sigma_i$ werden als *Singulärwerte* bezeichnet und 
die Spalten von $\bm{U}$ und $\bm{V}$ heißen *Links-* bzw. 
*Rechts-Singulärvektoren*. In Analogie zu 
Gl. {{eqref: eq:eigenvalue_decomposition_normal}} lässt sich die
Singulärwertzerlegung als
$$
  \begin{align}
    \bm{A} &= \bm{U} \bm{\Sigma} \bm{V}^\dag \\
    &= \begin{pmatrix}
      \vert & \vert &  & \vert \\
      \vec{u}_ 1 & \vec{u}_ 2 & \cdots & \vec{u}_ m \\
      \vert & \vert &  & \vert
    \end{pmatrix}
    \left(\begin{array}{cccc|ccc}
      \sigma_1 & 0 & \cdots & 0 & 0 & \cdots & 0 \\ 
      0 & \sigma_2 & \cdots & 0 & 0 & \cdots & 0 \\
      \vdots & \vdots & \ddots & 0 & \vdots & \ddots & \vdots \\
      0 & 0 & \cdots & \sigma_p & 0 & \cdots & 0 \\ \hline
      0 & 0 & \cdots & 0 & 0 & \cdots & 0 \\
      \vdots & \vdots & \ddots & \vdots & \vdots & \ddots & \vdots \\
      0 & 0 & \cdots & 0 & 0 & \cdots & 0
    \end{array}\right)
    \begin{pmatrix}
      \text{---}\ \vec{v}_ 1^\dag\ \text{---} \\
      \text{---}\ \vec{v}_ 2^\dag\ \text{---} \\
      \vdots \\
      \text{---}\ \vec{v}_ n^\dag\ \text{---}
    \end{pmatrix} \\
    &= \sigma_1 \vec{u}_ 1 \vec{v}_ 1^\dag + \sigma_2 \vec{u}_ 2 \vec{v}_ 2^\dag + \cdots + \sigma_p \vec{u}_ p \vec{v}_ p^\dag \\
    &= \sum_{i=1}^p \sigma_i \vec{u}_ i \vec{v}_ i^\dag 
    {{numeq}}{eq:singular_value_decomposition}
  \end{align}
$$
schreiben. 

Man erkennt hier einen deutlichen Unterschied zur Eigenwertzerlegung:
Während dort alle Eigenvektoren wichtig sind, tragen hier $(m-p)$ bzw. $(n-p)$
Singulärvektoren nicht zur Zerlegung bei. Das bedeutet, dass eine der 
Matrizen $\bm{U}$ oder $\bm{V}$ reduziert werden kann, ohne die Zerlegung zu
verändern. Nehmen wir ohne Einschränkung der Allgemeinheit $m \geq n$ an, so können 
wir $\hat{\bm{U}} \in \C{m}{n}$ und $\hat{\bm{\Sigma}} \in \R{n}{n}$
durch Entfernen der letzten $(m-n)$ Spalten von $\bm{U}$ und der letzten 
$(m-n)$ Zeilen von $\bm{\Sigma}$ definieren. Zudem definieren wir 
$\hat{\bm{V}} = \bm{V}$. Dann gilt die Zerlegung
$$
  \bm{A} = \hat{\bm{U}} \hat{\bm{\Sigma}} \hat{\bm{V}}^\dag\,,
$$
welche als *econiomic singular value decomposition* bekannt ist. Sollte eine der
Dimensionen $n, m$ viel größer sein als die andere, so ist die economic SVD
von Vorteil, weil sie weniger Speicherplatz benötigt.

### Eigenschaften
Wir Betrachten nun die Matrixprodukte $\bm{A}^\dag \bm{A}$ und
$\bm{A} \bm{A}^\dag$ unter Verwendung der SVD. Es gilt
$$
  \begin{align}
    \bm{A}^\dag \bm{A} &= (\bm{U} \bm{\Sigma} \bm{V}^\dag)^\dag \bm{U} \bm{\Sigma} \bm{V}^\dag 
    = (\bm{V}^\dag)^\dag \bm{\Sigma}^\dag \bm{U}^\dag \bm{U} \bm{\Sigma} \bm{V}^\dag
    = \bm{V} \bm{\Sigma}^\dag \bm{\Sigma} \bm{V}^\dag \\
    \bm{A} \bm{A}^\dag &= \bm{U} \bm{\Sigma} \bm{V}^\dag (\bm{U} \bm{\Sigma} \bm{V}^\dag)^\dag
    = \bm{U} \bm{\Sigma} \bm{V}^\dag (\bm{V}^\dag)^\dag \bm{\Sigma}^\dag \bm{U}^\dag
    = \bm{U} \bm{\Sigma} \bm{\Sigma}^\dag \bm{U}^\dag
  \end{align}\,.
  {{numeq}}{eq:svd_and_evd}
$$
Weil die von Null verschiedenen Elemente in $\bm{\Sigma}$ nur auf der 
Diagonalen stehen, sind die Produkte $\bm{\Sigma}^\dag \bm{\Sigma}$ und
$\bm{\Sigma} \bm{\Sigma}^\dag$ *Diagonalmatrizen* mit den Diagonalelementen
$\sigma_i^2$. 

Die zwei Gleichungen {{eqref: eq:svd_and_evd}} ähneln der Eigenwertzerlegung für normale Matrizen
$\bm{B} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$ sehr. Tatsächlich sind die Produkte $\bm{A}^\dag \bm{A}$ und 
$\bm{A} \bm{A}^\dag$ sogar hermitesch, weil
$$
  \begin{align}
    (\bm{A}^\dag \bm{A})^\dag &= \bm{A}^\dag (\bm{A}^\dag)^\dag = \bm{A}^\dag \bm{A} \\
    (\bm{A} \bm{A}^\dag)^\dag &= (\bm{A}^\dag)^\dag \bm{A}^\dag = \bm{A} \bm{A}^\dag
  \end{align}\,,
$$
wodurch sie außerdem [normal](https://de.wikipedia.org/wiki/Normale_Matrix) sind. 
Aus diesem Grund entsprechen die Singulärwerte den Wurzeln
der Eigenwerte von $\bm{A}^\dag \bm{A}$ und $\bm{A} \bm{A}^\dag$. Die 
Links-Singulärvektoren sind demnach die Eigenvektoren von $\bm{A} \bm{A}^\dag$, 
während die Rechts-Singulärvektoren die Eigenvektoren von
$\bm{A}^\dag \bm{A}$ darstellen. Somit haben wir eine Verbindung zwischen der SVD
und der EVD hergestellt.

Obwohl Gl. {{eqref: eq:svd_and_evd}} mathematisch korrekt ist, sollte man die
SVD von $\bm{A}$ nicht über die EVD von $\bm{A}^\dag \bm{A}$ und
$\bm{A} \bm{A}^\dag$ berechnen:

```admonish warning title="Phasenmehrdeutigkeit"
Die Eigenvektoren einer Matrix sind nur bis auf einen (komplexen) Faktor
eindeutig. Selbst wenn man die Eigenvektoren normiert, bleibt die Phase
unbestimmt. Versieht man z.B. den ersten Links-Singulärvektor $\vec{u}_1$ mit einem
Minuszeichen, bleiben die Gleichungen {{eqref: eq:svd_and_evd}} gültig.
Die SVD in Gl. {{eqref: eq:singular_value_decomposition}} ist dann jedoch nicht mehr
korrekt, da der ersten Term $\sigma_1 \vec{u}_1 \vec{v}_1^\dag$ ein 
Minuszeichen erhält. Damit die Zerlegung weiterhin stimmt, müsste das 
Vorzeichen von $\vec{v}_1$ angepasst werden. Da die zwei EVDs aber 
unabhängig von einander sind, sind solche Anpassungen nicht möglich.

Aus diesem Grund sollte man die SVD nicht aus der doppleten EVD berechnen, 
es sei denn, es werden lediglich die Singulärwerte und nicht die Singulärvektoren
benötigt. Für die folgenden Überlegungen ist es ausreichend zu wissen, dass
es gibt spezielle Algorithmen gibt, die die Links- und
Rechts-Singulärvektoren simultan berechnen.
```

### Verwendung mit NumPy
Die Singulärwertzerlegung ist in NumPy über die Funktion 
[`numpy.linalg.svd`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.svd.html)
verfügbar. Als Argument übergibt man die Matrix $\bm{A}$, von der die SVD
berechnet werden soll. Als Rückgabewerte erhält man die Matrizen $\bm{U}$,
$\bm{\Sigma}$ und $\bm{V}^\dag$. Mit dem optionalen Argument 
`full_matrices=False` wird die economic SVD berechnet. Die Standardoption ist
`full_matrices=True`. Ein Beispiel-Aufruf sieht wie folgt aus:
```python
import numpy as np

a = np.random.rand(9, 6)
u, s, vh = np.linalg.svd(a, full_matrices=False)
print(u.shape) # (9, 6)
print(s.shape) # (6,)
print(vh.shape) # (6, 6)
```

---

### Übung

#### Aufgabe 2: SVD und Eigengesichter

{{#include ../psets/04.md:aufgabe_2}}

