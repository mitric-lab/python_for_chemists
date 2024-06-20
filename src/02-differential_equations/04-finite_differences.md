## Finite-Differenzen-Verfahren

Die 
[finite-Differenzen-Methode](https://de.wikipedia.org/wiki/Finite-Differenzen-Methode)
ist eine weitere Klasse an numerischen Verfahren neben Runge-Kutta-Verfahren,
die zur Lösung von gewöhnlichen Differentialgleichungen eingesetzt werden 
kann. Auch hier wird der Definitionsbereich der Funktion $y(x)$ in diskrete Punkte unterteilt,
allerdings liegt der Fokus bei diesem Verfahren nicht auf der Berechnung
von Funktionswerten, sondern auf der Approximation der Ableitungsoperatoren.

### Theoretische Grundlagen

Bei der Implementierung des Gradientenverfahrens haben wir bereits 
in Gl. {{eqref: eq:finite_difference_symmetric}} eine Formel für die
finite Differenz verwendet. Diese Formel stellt ein Spezialfall der
allgemeinen finite Differenz dar und fiel damals mehr oder weniger vom Himmel.
Hier bemühen wir uns, ein allgemeines Rezept für die Herleitung solcher 
Formeln zu finden.

#### Finite-Differenzen-Approximationen

Gehen wir davon aus, dass wir die Ableitung einer Funktion $y(x)$
am Punkt $x$ aus den Funktionswerten in der Umgebung von $x$, wie z.B.
$y(x-h)$, $y(x)$ und $y(x+h)$ (nährungsweise) berechnen wollen. Demnach muss
$$
  y'(x) \approx c_{-1} y(x-h) + c_0 y(x) + c_1 y(x+h)
$$
gelten, wobei wir die Koeffizienten $c_{-1}$, $c_0$ und $c_1$ bestimmen wollen.
Wir entwickeln die rechte Seite in eine Taylor-Reihe um $x$:
$$
  \begin{align}
  y'(x) \approx \ 
      & c_{-1} y(x) &&+ c_{-1} y'(x) (-h) &&+ \frac{1}{2} c_{-1} y''(x) (-h)^2 + \ldots \\
    +\ & c_0 y(x) \\
    +\ & c_1 y(x) &&+ c_1 y'(x) h &&+ \frac{1}{2} c_1 y''(x) h^2 + \ldots
  \end{align}
$$
Die Idee ist nun, dass falls die rechte Seite gleich der linken Seite (also $y'(x)$) sein soll,
auch die Koeffizienten vor den Funktionswerten, bzw. den Ableitungen, auf beiden Seiten gleich
sein müssen. Da wir hier drei unbekannte Koeffizienten haben, benötigen wir drei Gleichungen,
weshalb die Taylor-Reihe bis zur zweiten Ordnung entwickelt wurde. Vergleichen wir die Koeffizienten
auf beiden Seiten, erhalten wir das folgende lineare Gleichungssystem:
$$
  \begin{align}
    c_{-1} + c_0 + c_1 &= 0 \\
    -c_{-1} h + c_1 h &= 1 \\
    \frac{1}{2} c_{-1} h^2 + \frac{1}{2} c_1 h^2 &= 0
  \end{align}
$$
Man könnte dieses Gleichungssystem mit dem Gauß-Algorithmus lösen,
bei nur drei Gleichungen und Unbekannten geht es auch durch einfaches Umformen und Einsetzen. 
Als Ergebnis erhalten wir die Koeffizienten $c_{-1} = -1/(2h)$,
$c_0 = 0$ und $c_1 = 1/(2h)$. Nach Einsetzen in die obige Formel erhalten wir
$$
  y'(x) \approx \frac{y(x+h) - y(x-h)}{2h}\,,
  {{numeq}}{eq:finite_difference_symmetric_second_order}
$$
was identisch zu Gl. {{eqref: eq:finite_difference_symmetric}} ist.

Diese Finite-Differenzen-Formel wird als symmetrische Differenz zweiter
Ordnung bezeichnet, da sie die Ableitung der Funktion $y(x)$ an der
Stelle $x$ symmetrisch aus den Funktionswerten $y(x-h)$ und $y(x+h)$
berechnet und dabei die Entwicklung der Taylor-Reihe bis zur zweiten
Ordnung verwendet.

Mit den Stützpunkten $x$ und $x+h$ bzw. $x$ und $x-h$ kann man auf 
gleiche Weise die Forwärtsdifferenzen
$$
  y'(x) \approx \frac{y(x+h) - y(x)}{h}
  {{numeq}}{eq:finite_difference_forward_first_order}
$$
und die Rückwärtsdifferenzen
$$
  y'(x) \approx \frac{y(x) - y(x-h)}{h}
  {{numeq}}{eq:finite_difference_backward_first_order}
$$
herleiten. Diese beiden Formeln berücksichtigen allerdings nur die erste Ordung.

Gemäß dieses Verfahrens können wir auch die Ableitungen höherer Ordnungen
mithilfe der Funktionswerte an beliebigen Stellen des Grid bis zur beliebigen Ordnung
approximieren. 

#### Matrixdarstellung des Differentialoperators

Wenn wir eine Funktion $y(x)$ auf einem Gitter $x_i$ mit $i=1,\ldots,N$
diskretisieren, können wir die Funktion als einen Vektor
$\vec{y} = (y_1, \ldots, y_N)^\intercal$ mit $y_i = y(x_i)$ (näherungsweise) darstellen. 
Wollen wir den Funktionswert an einem beliebigen Gridpunkt $x_i$ berechnen, 
können wir das als Skalarprodukt des Vektors $\vec{y}$ mit dem
Basisvektor $\hat{e}^i$ darstellen, wobei $\hat{e}^i$ ein Vektor mit
$1$ an der Stelle $i$ und $0$ an allen anderen Stellen ist. Das ergibt
$$
  \langle \hat{e}^i, \vec{y} \rangle 
    = \sum_{j=1}^N \hat{e}^i_ j y_ j 
    = \sum_{j=1}^N \delta_{ij} y_ j
    = y_i
$$
mit dem Kroenecker-Delta $\delta_{ij}$. Schreiben wir die 
Vektorkomponenten aus, ist das Skalarprodukt durch
$$
  \begin{pmatrix}
    0_1 & \cdots & 0_{i-1} & 1_i & 0_{i+1} & \cdots & 0_N
  \end{pmatrix}
  \cdot
  \begin{pmatrix}
    y_1 \\ \vdots \\ y_{i-1} \\ y_i \\ y_{i+1} \\ \vdots \\ y_N
  \end{pmatrix}
  = y_i
$$
gegeben, wobei die Indizes der Einträge in $e^i$ explizit geschrieben wurden.

Auf diese Weise können wir $\left(y(x_{i+1}) - y(x_{i-1})\right)$ als
$\langle \hat{e}^{i+1} - \hat{e}^{i-1}, \vec{y} \rangle$ darstellen,
mit dem Hilfsvektor 
$$
  \hat{e}^{i+1} - \hat{e}^{i-1} = 
  \begin{pmatrix}
    0_1 & \cdots & 0_{i-2} & -1_{i-1} & 0_i & 1_{i+1} & 0_{i+2} & \cdots & 0_N
  \end{pmatrix}\,.
$$

Die Ableitung an der Stelle $x_i$ kann gemäß Gl. {{eqref: eq:finite_difference_symmetric_second_order}} 
also als  $\frac{1}{2h} \langle \hat{e}^{i+1} - \hat{e}^{i-1}, \vec{y} \rangle$
geschrieben werden. Da wir die Ableitung an allen Gridpunkten berechnen wollen,
setzen wir das Muster fort und erhalten eine Matrixgleichung
$$
  \vec{y}' = \bm{D} \vec{y}
$$
mit
$$
  \vec{y}' = (y'_1, \ldots, y'_N)^\intercal, \quad
  \vec{y} = (y_1, \ldots, y_N)^\intercal
$$
und
$$
  \bm{D} = \frac{1}{2h}
  \begin{pmatrix}
    0 & 1 & 0 & 0 & \cdots & 0 & 0 \\
    -1 & 0 & 1 & 0 & \cdots & 0 & 0 \\
    0 & -1 & 0 & 1 & \cdots & 0 & 0 \\
    \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
    0 & 0 & 0 & 0 & \cdots & -1 & 0 \\
  \end{pmatrix}\,.
  {{numeq}}{eq:finite_difference_symmetric_second_order_matrix}
$$

Diese *Darstellung des Differentialoperators* mit symmetrischen 
Finite-Differenzen zweiter Ordnung $\bm{D}$ ist also eine 
[*Tridiagonalmatrix*](https://de.wikipedia.org/wiki/Tridiagonalmatrix)
mit $0$ auf der Hauptdiagonale und $\pm 1$ auf den Nebendiagonalen.

Wie sieht es mit der Darstellung der zweiten oder höheren Ableitungen aus?
Natürlich könnte man eine Differentialgleichung höherer Ordnung als ein
System von Differentialgleichungen erster Ordnung auffassen und die 
Finite-Differenzen-Methode für vektorwertige Funktionen anpassen. Aber 
hier können wir die entsprechende Matrix auf einem direkten Weg konstruieren.

Dazu zeigen wir zunächst, wie man es nicht machen sollte: Die Form des Operators
für die zweite Ableitung $\frac{\du^2}{\du x^2}$ erinnert auf den ersten Blick an das
"Quadrat" des Ableitungsoperators $\frac{\du}{\du x}$. Sollten wir
tatsächlich die Matrix $\bm{D}$ quadrieren, erhalten wir etwas, was
auf den ersten Blick nicht falsch aussieht:
$$
  \bm{D}^2 y = \bm{D} \bm{D} y = \bm{D} y' = y''\,.
$$

Diese Gleichung wäre dann korrekt, wenn $\bm{D} y$ wirklich die Ableitung
von $y$ wäre. Das ist hier aber nicht der Fall, weil wir durch die Diskretisierung 
eine endliche Auflösung haben. Das Anwenden des Operators $\bm{D}$ auf $\bm{D} y$
würde daher den Fehler der ersten Ableitung verstärken.

```admonish warning title="Warnung"
Benutzen Sie im Allgemeinen **nicht** die Matrix $\bm{D}^n$ als die 
Darstellung des Operators der $n$-ten Ableitung!
```

Jetzt zeigen wir, wie es richtig gemacht wird. Auch hier fassen wir zunächst die zweite
Ableitung als die Ableitung der ersten Ableitung auf:
$$
  \begin{align}
    y''(x) = \frac{\du}{\du x} y'(x) 
      &\approx \frac{1}{h} \left( y'(x) - y'(x-h) \right) \\
      &= \frac{1}{h} \left( \frac{\du}{\du x} y(x) - \frac{\du}{\du x} y(x-h) \right) \\
      &= \frac{1}{h} \left( \frac{y(x+h) - y(x)}{h} - \frac{y(x) - y(x-h)}{h} \right) \\
      &= \frac{1}{h^2} \left( y(x+h) - 2y(x) + y(x-h) \right)\,,
  \end{align}
$$
wobei wir für $\frac{\du}{\du x} y'$ die Rückwärtsdifferenzen
{{eqref: eq:finite_difference_backward_first_order}} und
für $\frac{\du}{\du x} y$ die Vorwärtsdifferenzen
{{eqref: eq:finite_difference_forward_first_order}} benutzt haben.
Diese Herleitung lässt sich leicht auf die $n$-te Ableitung
verallgemeinern.

Mit dieser Formel können wir die Matrixdarstellung des Operators für die
zweite Ableitung, hier als $\bm{D}^{(2)}$ bezeichnet, konstruieren:
$$
  \bm{D}^{(2)} = \frac{1}{h^2}
  \begin{pmatrix}
    -2 & 1 & 0 & 0 & \cdots & 0 & 0 \\
    1 & -2 & 1 & 0 & \cdots & 0 & 0 \\
    0 & 1 & -2 & 1 & \cdots & 0 & 0 \\
    \vdots & \vdots & \vdots & \vdots & \ddots & \vdots & \vdots \\
    0 & 0 & 0 & 0 & \cdots & 1 & -2 \\
  \end{pmatrix}\,.
  {{numeq}}{eq:second_finite_difference_symmetric_second_order_matrix}
$$
Diese Matrix ist wieder eine Tridiagonalmatrix, aber mit $-2$ auf der
Hauptdiagonale und $1$ auf den Nebendiagonalen.

#### Finite-Differenzen-Verfahren

Wir haben nun die Matrixdarstellung des Differentialoperators für die
ersten und zweiten Ableitungen. In einer Differentialgleichung tauchen 
neben der gesuchten Funktion $y(x)$ und
ihren Ableitungen allerdings auch Funktionen von $x$ auf. Wir benötigen also eine
Darstellung für solche Funktionen. Eine Besonderheit hierbei ist, dass es
es dafür zwei Darstellungsvarianten gibt, welche man gemäß dem Kontext wählen
muss. 

Fungiert eine Funktion von $x$ als eine Inhomogenität $u(x)$
(vgl. Gl. {{eqref: eq:lode_general}}), also als ein alleinstehender Term ohne 
$y(x)$ oder ihre Ableitungen, wird sie genau so wie bei $y(x)$ diskretisiert 
und als ein Vektor $\vec{u} = (u_1, \ldots, u_N)^\intercal$ mit $u_i = u(x_i)$ dargestellt. 
Wird die Funktion aber mit $y(x)$ oder ihren Ableitungen multipliziert, 
dient sie als eine Koeffizientenfunktion $a_i(x)$, und muss als eine 
Matrix darstellt werden. Da die Multiplikation zweier 
Funktionen elementweise erfolgt, nimmt die Darstellungsmatrix $\bm{A}_ i$ 
Diagonalform an, wobei die Diagonalelemente die diskretisierten 
Funktionswerte von $a_i(x)$ sind, also 
$\bm{A}_ i = \text{diag}(a_i(x_1), \ldots, a_i(x_N))$.

Nun kennen wir die (approximative) Darstellung aller Elemente einer 
DGL und können eine beliebige **lineare** Differentialgleichung
(vgl. Gl. {{eqref: eq:lode_general}}) in eine Matrixgleichung umwandeln:
$$
  \underbrace{
    \bm{A}_n \bm{D}^{(n)} \vec{y} + \bm{A}_{n-1} \bm{D}^{(n-1)} \vec{y} + \ldots + \bm{A}_1 \bm{D} \vec{y} + \bm{A}_0 \vec{y}
  }_ {\bm{L} \vec{y}} = \bm{B}_ 0 \vec{u}\,
$$
Hier dient $\bm{B}_ 0 = (b_0(x_1), \ldots, b_0(x_N))^\intercal$ als die 
Matrixdarstellung der Koeffizientenfunktion $b_0(x)$ und wir können
in der obigen Gleichung haben wir alle lineare Operatoren zu einem Operator
$$
  \bm{L} = \bm{A}_n \bm{D}^{(n)} + \bm{A}_{n-1} \bm{D}^{(n-1)} + \ldots + \bm{A}_1 \bm{D} + \bm{A}_0
$$
zusammengefasst. Die Lösung des linearen Gleichungssystems 
$\bm{L} \vec{y} = \bm{B}_ 0 \vec{u}$ liefert demnach die diskretisierte Funktion
$\vec{y}$.

```admonish note title="Anmerkung zu nichtlinearen DGLs"
Tatsächlich lassen sich auch nichtlineare DGLs mit 
Finiten-Differenzen-Operatoren ausdrücken. Das resultierte Gleichungssystem
ist allerdings ein nichtlineares Gleichungssystem, welches nicht direkt mit
Methoden der linearen Algebra gelöst werden kann. 
```

Zwischen dem Lösen des linearen Gleichungssystems und dem Lösen einer 
Differentialgleichung gibt es einen entscheidenden Unterschied: 
Der Anfangswert. Während bei der Lösung der DGLs eine Anfangs- oder
Randbedingung benötigt wird, um eine spezielle Lösung zu erhalten, 
gibt es eine solche Bedingung bei einem linearen Gleichungssystem nicht. 
Wie sollen wir dann die Anfangsbedingungen der DGL im Rahmen der diskretisierten 
Version berücksichtigen?

Betrachten wir dazu zunächst die Matrixdarstellung des Differentialoperators $\bm{D}$
in Gl. {{eqref: eq:finite_difference_symmetric_second_order_matrix}}, 
insbesondere die erste Zeile $(0, 1, 0, \ldots, 0)$. Multipliziert mit $\vec{y}$ 
besagt diese Zeile
$$
  y'(x_1) = \frac{1}{2h} y(x_2).
$$
Gemäß Gl. {{eqref: eq:finite_difference_symmetric_second_order}} sollte aber
doch
$$
  y'(x_1) = \frac{1}{2h} \left( y(x_2) - y(x_0) \right)
$$
gelten. Damit beide Gleichungen übereinstimmen, muss $y(x_0) = 0$ gelten. 
Die letzte Zeile von $\bm{D}$ liefert wiederum die Bedingung $y(x_{N+1}) = 0$.

Durch die Konstruktion der Finite-Differenzen-Operatoren als Matrizen werden die 
Randbedingungen also implizit festgelegt. Das Finiten-Differenzen-Verfahren ist 
deshalb eher geeignet für Randwertprobleme als für Anfangswertprobleme.

```admonish note title="Anmerkung zu anderen Randbedingungen"
Neben der 
[*Dirichlet-Randbedingung*](https://de.wikipedia.org/wiki/Dirichlet-Randbedingung)
in unserem Fall, also dass die Funktionswerte außerhalb des Grids 
Null sein müssen, können noch andere Randbedingungen, wie z.B. 
[*periodische Randbedingungen*](https://de.wikipedia.org/wiki/Periodische_Randbedingung)
oder 
[*Neumann-Randbedingungen*](https://de.wikipedia.org/wiki/Neumann-Randbedingung)
mit dem Finite-Differenzen-Verfahren berücksichtigt werden. Diese werden wir
aber hier nicht weiter behandeln.
```


### Implementierung

Wir wollen nun das Finite-Differenzen-Verfahren am Beispiel der 
Schrödingergleichung für den harmonischen Oszillator implementieren.

Die Schrödingergleichung in atomaren Einheiten lautet
$$
  -\frac{1}{2} \frac{\du^2}{\du x^2} \psi(x) + \frac{1}{2} k x^2 \psi(x) = E \psi(x)
$$
und ist eine lineare DGL zweiter Ordung mit Randbedingungen $\lim_{x \to \pm \infty} \psi(x) = 0$.
Für die Implementierung eignet sich also die Darstellungsmatrix $\bm{D}^{(2)}$ in 
Gl. {{eqref: eq:second_finite_difference_symmetric_second_order_matrix}}.

Nach dem Importieren der benötigten Libraries
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:imports}}
```
können wir den Differentialoperator $\bm{D}^{(2)}$ wie folgt definieren:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:generate_d2_naive}}
```
Während diese einfache Implementierung zwar korrekt ist, ist sie allerdings nicht sehr 
effizient, insbesondere wenn die Anzahl der Gridpunkte `n` groß ist.
Eine effizientere Implementierung kann mithilfe der Funktion
[`np.diag_indices`](https://numpy.org/doc/stable/reference/generated/numpy.diag_indices.html)
erzielt werden:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:generate_d2}}
```
Nach der Initialisierung einer Nullmatrix haben wir uns mit Hilfe der Funktion `np.diag_indices`
die Indizes der Hauptdiagonalen ausgeben lassen. Damit können die entsprechenden Elemente 
des Arrays `d2` auf $-2$ gesetzt werden.
Durch Verschieben der Indizes um $\pm 1$ können wir die Einträge der Nebendiagonalen ansprechen
und die entsprechenden Elemente auf $1$ setzen. Der Verzicht auf
Schleifen beschleunigt die Berechnung erheblich.

```admonish note title="Anmerkung zu Python-Loops"
Wie wir schon in der Einführung erwähnt haben, ist Python im Vergleich 
zu kompilierten Sprachen eine eher langsame Sprache. Deshalb sollte man, 
wenn Effizienz von entscheidender Bedeutung ist, Schleifen vermeiden und stattdessen
Funktionen aus der Bibliothek `numpy` benutzen.
```

Bei Anwendungen, bei denen eine große Anzahl von Gitterpunkten benötigt wird,
kann die Verwendung von
[dünnbesetzten Matrizen](https://de.wikipedia.org/wiki/Dünnbesetzte_Matrix)
(engl. *sparse matrices*) von Vorteil sein. Da wir in diesem Beispiel nur eine moderate Anzahl von
Gitterpunkten benötigen, werden wir die Implementierung mit einer normalen Matrix durchführen.

```admonish info title="Anmerkung zu dünnbesetzten Matrizen" collapsible=true
Bei großen $N$ durch feinere oder auch mehrdimensionale Grids kann
die Finite-Differenzen-Matrix sehr groß werden. Allerdings ist der
Großteil der Elemente dieser Matrix Null. Solche Matrizen werden als
[*dünnbesetzte Matrizen*](https://de.wikipedia.org/wiki/Dünnbesetzte_Matrix)
(engl. *sparse matrices*) bezeichnet und können durch spezielle Algorithmen
effizienter gespeichert und verarbeitet werden. Die 
Fininte-Differenzen-Matrix hat sogar nur Einträge auf der Haupt- und
einigen Nebendiagonalen und wird deshalb auch als 
[*Bandmatrix*](https://de.wikipedia.org/wiki/Bandmatrix) bezeichnet,
was eine verallgemeinerung der Tridiagonalmatrix ist. Bandmatrizen
können mit Hilfe der Funktion
[`scipy.sparse.diags_array`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.diags_array.html)
als dünnbesetzte Matrizen implementiert werden.

Dann können Methoden aus dem Untermodul
[`scipy.sparse.linalg`](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html)
eingesetzt werden, um lineare Gleichungssysteme zu lösen oder Eigenwerte
und Eigenvektoren zu berechnen.

Möchte man nur die Eigenwerte und Eigenvektoren von Bandmatrizen berechnen,
kann die Funktion
[`scipy.linalg.eigh_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_banded.html)
oder
[`scipy.linalg.eigh_tridiagonal`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_tridiagonal.html)
eingesetz werden, die als Argument lediglich die besetzten Diagonalen der Bandmatrix erwartet.
```

Wir konstruieren anschließend den Hamiltonoperator für den harmonischen Oszillator
als Matrix:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:build_hamiltonian}}
```
Hier haben wir zuerst die $D^{(2)}$-Matrix mit der Funktion `generate_d2`
erzeugt und dann die Potentialfunktion $\frac{1}{2} k x^2$ als eine
Diagonalmatrix mit `np.diag` erstellt. Die Summe 
$-\frac{1}{2} D^{(2)} + \frac{1}{2} k x^2$ ergibt dann den Hamiltonoperator,
bzw. den linearen Operator $\bm{L}$ (`l_mat`).

Damit lautet die Schrödingergleichung in diskreter Form
$$
  \bm{L} \vec{\psi} = E \vec{\psi}\,.
$$
Man erkennt an dieser Stelle leicht, dass es sich um eine Matrix-Eigenwertgleichung handelt.
Wir verwenden deshalb die Funktion 
[`np.linalg.eigh`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html), um die Eigenwerte und Eigenvektoren des Hamiltonoperators zu
berechnen. 

Wir setzen dafür zunächst $k = 1$ und wählen ein Grid von -5 bis 5 mit 512 Punkten:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:define_parameters}}
```

Danach wird der Hamiltonoperator erzeugt und die Eigenwerte und
Eigenvektoren berechnet:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:solve_matrix_equation}}
```
Da die Funktion `np.linalg.eigh` eine symmetrische Matrix erwartet, haben
wir mit `assert` und
[`np.allclose`](https://numpy.org/doc/stable/reference/generated/numpy.allclose.html)
überprüft, ob die Matrix `hamiltonian` identisch zu ihrer Transponierten ist.

Wir entnehmen nun die Energien und die zugehörigen Wellenfunktionen aus den
ersten 20 Eigenvektoren und Eigenwerten. 
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:solve_harmonic_oscillator}}
```
Allerdings sollten die Wellenfunktionen gemäß 
$$
  \int_{-\infty}^{\infty} |\psi(x)|^2 \du x = 1
$$
normiert sein, was in der diskreten Form
$$
  \sum_{i=1}^{N} |\psi(x_i)|^2 \Delta x = 1
$$
entspricht. Die Eigenvektoren aus `np.linalg.eigh` sind aber gemäß
$$
  \sum_{i=1}^{N} |v_i|^2 = 1
$$
normiert. Deshalb teilen wir die Eigenvektoren durch $\sqrt{\Delta x}$, wobei wir
$\Delta x$ aus den ersten beiden Einträgen des Grids entnehmen,
um die normierten Wellenfunktionen zu erhalten.

Zuletzt wollen wir unsere Ergebnisse visualisieren. Wir wollen dabei die Eigenenergien 
und die Wellenfunktionen nebeneinander plotten. Deshalb übergeben wir die folgenden Argumente
an die Funktion `plt.subplots`:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:plot_results_subplots}}
```
Hier bedeuten die Argumente `1, 2`, dass wir zwei separate Diagramme in einer Zeile und
zwei Spalten erzeugen möchten. Die Achsen-Objekte werden dann in der Liste `axs` gespeichert.

Danach plotten wir die numerischen sowie die analytischen Eigenenergien
des harmonischen Oszillators in dem ersten Plot mit `axs[0]`:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:plot_eigenenergies}}
```
Anschließend plotten wir in dem zweiten Diagramm (`axs[1]`) das harmonische
Potential sowie die ersten 5 numerischen Wellenfunktionen:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:plot_eigenfunctions}}
```
Beim Plotten des Potentials haben wir das Argument `lw=2` benutzt, um die
Linienbreite zu erhöhen. 

Zuletzt formatieren wir den Plot und zeigen ihn an:
```python
{{#include ../codes/02-differential_equations/fdm_harm_osc.py:show_plot}}
```

Das Diagramm sollte wie folgt aussehen:
![Eigenenergien und Wellenfunktionen des harmonischen Oszillators](../assets/figures/02-differential_equations/fdm_harm_osc.svg)

Man erkennt, dass die numerischen Eigenenergien der ersten ca. 10
Zuständen sehr gut mit den analytischen Eigenenergien übereinstimmen.
Die numerisch berechneten Wellenfunktionen der ersten fünf Zustände sehen zumindest 
auf den ersten Blick sinnvoll aus. Sollte man größere Genauigkeit für die höheren Zustände 
benötingen, muss das Grid feiner und auch größer gewählt werden, da die 
Wellenfunktionen einerseits mehr Oszillationen aufweisen und andererseits
räumlich ausgedehnter sind.


---

### Übung

#### Aufgabe 2.3: Lösen der Schrödingergleichung mit dem Schießverfahren

{{#include ../psets/02.md:aufgabe_3}}

