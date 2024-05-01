## Runge-Kutta-Verfahren

Die Runge-Kutta-Verfahren sind eine
nach Carl Runge und Martin Wilhelm Kutta benannte
Familie von Methoden zur numerischen 
Lösung von Anfangswertproblemen für gewöhnliche Differentialgleichungen.
Diese Methoden berechnen iterativ die Lösung im nächsten Zeitschritt aus 
einer Linearkombination von dem Funktionswert und den Steigungen an
verschiedenen Stellen. Sprich man jedoch von **dem** Runge-Kutta-Verfahren,
ist das 
[klassische Runge-Kutta-Verfahren](https://de.wikipedia.org/wiki/Klassisches_Runge-Kutta-Verfahren)
gemeint, was aber nur ein Spezialfall der Runge-Kutta-Verfahren darstellt.

### Theoretische Grundlagen
Genau so wie beim Euler-Verfahren, wählen wir zunächst ein gleichmäßiges
Grid $\{x_n\}_ {n=0,\cdots,N}$ mit $x_n = x_0 + n\cdot h$ und $h$ als
Schrittweite, sowie eine Anfangsbedingung $y_0 = y(x_0)$.
Dann entwickeln wir $y(x_{n+1})$ wieder in eine Taylor-Reihe um $y(x_n)$:
$$
  y(x_{n+1}) = y(x_n) 
  + h y'(x_n) 
  + \frac{h^2}{2} y''(x_n) 
  + \frac{h^3}{6} y'''(x_n)
  + \frac{h^4}{24} y^{(4)}(x_n)
  + \cdots\,.
  {{numeq}}{eq:ode_taylor_series}
$$
Man spricht von einem Runge-Kutta-Verfahren der *(Konsistenz-)Ordnung* $p$,
wenn das Taylorpolynom bis zum $p$-ten Grad im Verfahren berücksichtigt wird.

Ein allgemeines Runge-Kutta-Verfahren für das AWP in 
Gl. {{eqref: eq:ivp_first_order}} ist dann durch
$$
  y_{n+1} = y_n + h \sum_{i=1}^s b_i k_i
  {{numeq}}{eq:runge_kutta_general}
$$
gegeben, wobei
$$
  k_i = f(x_n + h c_i, y_n + h \sum_{j=1}^{s} a_{ij} k_j)\,.
$$
Dabei bezeichnet $s$ die *Stufe* des Verfahrens. Die Koeffizienten $a_{ij}$,
$b_i$ und $c_i$ sind charakteristische Parameter des Verfahrens.
Hier werden wir nur *explizite* Runge-Kutta-Verfahren betrachten, bei denen
die Koeffizienten $c_1 = 0$ und $a_{ij} = 0$ für $j \geq i$ gelten.

Weil die allgemeine Form doch sehr unhandlich ist, betrachten wir zuerst ein
zwei-stufiges Verfahren zweiter Ordnung.

#### Runge-Kutta-Verfahren zweiter Ordnung (RK2)
Nach Gl. {{eqref: eq:runge_kutta_general}} können wir ein explizites
zwei-stufiges Runge-Kutta-Verfahren zweiter Ordnung formulieren als
$$
  \begin{align}
    y_{n+1} &= y_n + h(b_1 k_1 + b_2 k_2) \\
    k_1 &= f(x_n, y_n) \\
    k_2 &= f(x_n + h c_2, y_n + h a_{21} k_1)\,.
  \end{align}
  {{numeq}}{eq:runge_kutta_2_stage}
$$

Für ein Verfahren zweiter Ordnung benötigen wir ein Taylorpolynom zweiten
Grades, in dem die zweite Ableitung von $y$ vorkommt. Da wir aber $y$ nicht kennen,
haben wir keinen direkten Zugriff auf die zweite Ableitung. Wir können jedoch
mit Hilfe der Kettenregel $y''$ wie folgt ausdrücken:
$$
  \frac{\du^2 y}{\du x^2} = \frac{\du}{\du x} y' = \frac{\du f(x, y(x))}{\du x}
    = \frac{\partial f}{\partial x} + \frac{\partial f}{\partial y} y'
    = \frac{\partial f}{\partial x} + f \frac{\partial f}{\partial y} \,.
$$ 
Das liefert und das Taylorpolynom zweiten Grades:
$$
  y(x_{n+1}) = y(x_n) + h f(x_n, y_n) + \frac{h^2}{2} \left[
    \frac{\partial f}{\partial x} + f \frac{\partial f}{\partial y}
  \right](x_n, y_n) + \mathcal{O}(h^3)\,.
  {{numeq}}{eq:ode_taylor_polynomial_2nd_order}
$$

Anschließend entwickeln wir $k_2$ in Gl. {{eqref: eq:runge_kutta_2_stage}}
linear um $(x_n, y_n)^\intercal$:
$$
  k_2 = f(x_n, y_n) + h c_2 \frac{\partial f}{\partial x}(x_n, y_n)
    + h a_{21} k_1 \frac{\partial f}{\partial y}(x_n, y_n) + \mathcal{O}(h^2)\,.
$$
Setzen wir dies in Gl. {{eqref: eq:runge_kutta_2_stage}} ein, erhalten wir
$$
  y_{n+1} = y_n + h(b_1 + b_2) f(x_n, y_n) + h^2 b_2 \left[
    c_2 \frac{\partial f}{\partial x}(x_n, y_n) + a_{21} \left(
      f \frac{\partial f}{\partial y}
    \right)(x_n, y_n)
  \right] + \mathcal{O}(h^3)\,.\\
  {{numeq}}{eq:runge_kutta_2_stage_2nd_order}
$$

Damit das Verfahren tatsächlich eine Konsistenzordnung von $p=2$ hat, müssen die
Koeffizienten vor den Funktionen $f$ und ihren Ableitungen in Gl.
{{eqref: eq:ode_taylor_polynomial_2nd_order}} und
{{eqref: eq:runge_kutta_2_stage_2nd_order}} übereinstimmen, da diese für beliebige
Funktionen $f$ gelten müssen. Das führt zu den Bedingungen
$$
  \begin{align}
    b_1 + b_2 &= 1 \\
    b_2 c_2 &= \frac{1}{2} \\
    b_2 a_{21} &= \frac{1}{2}\,.
  \end{align}
  {{numeq}}{eq:runge_kutta_2_stage_conditions}
$$
Das ist ein unterbestimmtes Gleichungssystem mit drei Gleichungen für vier
Unbekannte. Wir können demnach einen der Koeffizienten frei wählen und
erhalten so eine Familie von konsistenten
Runge-Kutta-Verfahren zweiter Ordnung mit zwei Stufen.

Wählen wir $b_1 = b_2 = \frac{1}{2}$, $c_2 = 1$ und $a_{21} = 1$, erhalten wir
das sog. [Heun-Verfahren](https://de.wikipedia.org/wiki/Heun-Verfahren).
Wählen wir dagegen $b_1 = 0$, $b_2 = 1$, $c_2 = \frac{1}{2}$ und 
$a_{21} = \frac{1}{2}$, erhalten wir das sog. 
[Mittelpunktsverfahren](https://en.wikipedia.org/wiki/Midpoint_method).
Neben diesen beiden geläufigen Verfahren kann man natürlich auch andere
Kombinationen der Koeffizienten wählen, solange die Bedingungen in Gl.
{{eqref: eq:runge_kutta_2_stage_conditions}} erfüllt sind. Eine allgemeine
Parametrisierung der Koeffizienten lautet
$$
  \begin{align}
    a_{21} &= \alpha \\
    b_1 &= 1 - \frac{1}{2 \alpha} \\
    b_2 &= \frac{1}{2 \alpha} \\
    c_2 &= \alpha \,.
  \end{align}
$$

Die Angabe der einzelnen Koeffizienten auf diese Weise ist jedoch nicht sehr übersichtlich,
wenn man bedenkt, dass es auch Verfahren mit mehr als zwei Stufen gibt. Daher gibt es eine
kompaktere Schreibweise, das sog. *Butcher-Tableau*.

#### Butcher-Tableau
Das Butcher-Tableau ist eine übersichtliche Darstellung der Koeffizienten
$a_{ij}$, $b_i$ und $c_i$ eines Runge-Kutta-Verfahrens. Das Tableau für ein
allgemeines $s$-stufiges Verfahren lautet
$$
  \begin{array}{c|c}
    \vec{c} & \bm{A} \\ \hline
            & \vec{b}^\intercal
  \end{array}
  =
  \begin{array}{c|cccc}
    c_1 & a_{11} & a_{12} & \cdots & a_{1s} \\
    c_2 & a_{21} & a_{22} & \cdots & a_{2s} \\
    \vdots & \vdots & \vdots & \ddots & \vdots \\
    c_s & a_{s1} & a_{s2} & \cdots & a_{ss} \\ \hline
        & b_1 & b_2 & \cdots & b_s \\
  \end{array}\,,
  {{numeq}}{eq:butcher_tableau}
$$

wobei die Koeffizienten mit einem Index ($b_i$ und $c_i$) als Vektoren
und die Koeffizienten mit zwei Indices ($a_{ij}$) als eine Matrix
dargestellt werden.

Das Heun-Verfahren kann damit als
$$
  \begin{array}{c|cc}
    0 & 0 & 0 \\
    1 & 1 & 0 \\ \hline
      & 1/2 & 1/2
  \end{array}
$$
und das Mittelpunktsverfahren als
$$
  \begin{array}{c|c}
    0 & 0 & 0 \\
    1/2 & 1/2 & 0 \\ \hline
        & 0 & 1
  \end{array}
$$
dargestellt werden. Die allgemein parametrisierte Form eines 2-stufigen Verfahrens
zweiter Ordnung lautet dann
$$
  \begin{array}{c|cc}
    0 & 0 & 0 \\
    \alpha & \alpha & 0 \\ \hline
      & 1 - \frac{1}{2\alpha} & \frac{1}{2\alpha}
  \end{array}\,.
$$

```admonish note title="Hinweis"
Wenn man möchte, kann man das Euler-Verfahren als ein einstufiges
Runge-Kutta-Verfahren erster Ordnung auffassen, was als 
$$
  \begin{array}{c|c}
    0 & 0 \\ \hline
      & 1
  \end{array}
$$
dargestellt werden kann. Verifizieren Sie dies durch Einsetzen der Koeffizienten
in Gl. {{eqref: eq:runge_kutta_general}}, was zu Gl. {{eqref: eq:euler_method}}
führen sollte.
```

Als eine letzte Bemerkung sei noch gesagt, dass die Bedingung $a_{ij} = 0$ 
für $j \geq i$ bei expliziten Runge-Kutta-Verfahren in dieser Darstellung
bedeutet, dass die Matrix $\bm{A}$ eine strikte untere Dreiecksmatrix ist.

#### Klassisches Runge-Kutta-Verfahren (RK4)
Das am häufigsten verwendete Runge-Kutta-Verfahren ist das klassische
Runge-Kutta-Verfahren. Diese ist eine vierstufige Methode vierter Ordnung,
hat also die Form
$$
  \begin{align}
    y_{n+1} &= y_n + h(b_1 k_1 + b_2 k_2 + b_3 k_3 + b_4 k_4) \\
    k_1 &= f(x_n, y_n) \\
    k_2 &= f(x_n + h c_2, y_n + h a_{21} k_1) \\
    k_3 &= f(x_n + h c_3, y_n + h(a_{31} k_1 + a_{32} k_2)) \\
    k_4 &= f(x_n + h c_4, y_n + h(a_{41} k_1 + a_{42} k_2 + a_{43} k_3)) \,.
  \end{align}
  {{numeq}}{eq:runge_kutta_4_stage}
$$
Die Bedingungen für die Koeffizienten kann auf die gleiche Weise wie beim
RK2 hergeleitet werden; die Rechnungen sind jedoch deutlich aufwendiger, 
weswegen wir sie hier nicht durchführen.

Beim klassischen Runge-Kutta-Verfahren lauten die Koeffizienten
$$
  \begin{array}{c|cccc}
    0 & 0 & 0 & 0 & 0 \\
    1/2 & 1/2 & 0 & 0 & 0 \\
    1/2 & 0 & 1/2 & 0 & 0 \\
    1 & 0 & 0 & 1 & 0 \\ \hline
      & 1/6 & 1/3 & 1/3 & 1/6
  \end{array}\,.
$$
Eine weiteres, in dem gleichen Paper<sup>[</sup>[^kutta1901]<sup>]</sup> 
wie das klassische RK4-Verfahren vorgestellte, aber bei weitem nicht so 
bekanntes Verfahren nutzt die Koeffizienten
$$
  \begin{array}{c|cccc}
    0 & 0 & 0 & 0 & 0 \\
    1/3 & 1/3 & 0 & 0 & 0 \\
    2/3 & -1/3 & 1 & 0 & 0 \\
    1 & 1 & -1 & 1 & 0 \\ \hline
      & 1/8 & 3/8 & 3/8 & 1/8
  \end{array}\,.
$$

Sie wundern sich vielleicht, warum das klassische RK4-Verfahren so beliebt ist,
obwohl es Verfahren mit höherer Ordnung gibt. Dabei sollte man allerdings das
folgende bedenken:

```admonish warning title="Ordnung &ne; Stufe!"
Die Konsistenzordnung $p$ und die Stufe $s$ des Verfahrens sind zwei 
verschiedene Dinge, obwohl wir bis jetzt nur Verfahren mit
$s = p$ betrachtet haben. Tatsächlich gilt für die minimale 
Stufenzahl $s_{\mathrm{min}}$ zum Erreichen einer Konsistenzordnung $p$
bei expliziten Runge-Kutta-Verfahren stets
$s_{\mathrm{min}} \geq p$.
```

Man kann sogar zeigen, dass für $p \geq 5$ die strikte Ungleichung
$s_{\mathrm{min}} > p$ gilt.<sup>[</sup>[^butcher1987]<sup>]</sup>
In anderen Worten: Die Verbesserung der Genauigkeit von $p = 4$
auf $p = 5$ unter Verwendung von expliziten Runge-Kutta-Verfahren
ist mit einer Erhöhung der Stufezahl um mindestens 2
verbunden. Dies erklärt, warum das klassische RK4-Verfahren so beliebt ist.

Der Zusammenhang zwischen $p$ und $s_{\mathrm{min}}$ für einige Ordnungen
von expliziten Runge-Kutta-Verfahren
ist in der folgenden Tabelle zusammengefasst:<sup>[</sup>[^butcher1987]<sup>]</sup>
|                    |    |    |    |    |    |    |    |    |
|--------------------|---:|---:|---:|---:|---:|---:|---:|---:|
| $p$                |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |
| $s_{\mathrm{min}}$ |  1 |  2 |  3 |  4 |  6 |  7 |  9 | 11 |

Die Zahlen $s_{\mathrm{min}}$ sind auch als *Butcher-Schranken* bekannt.

[^kutta1901]: M. W. Kutta, *Z. Math. Phys.* **1901**, *46*, 435&ndash;453.
[^butcher1987]: J. C. Butcher, in *The Numerical Analysis of Ordinary Differential Equations*, John Wiley & Sons, Chichester, **1987**, pp. 185&ndash;194.

### Implementierung
Wir verwenden als Beispiel erneut die Dynamik der Belousov-Zhabotinsky-Reaktion.
Genau wie im Abschnitt [2.2](02-euler_method.md#belousov-zhabotinsky-reaktion)
importieren wir zunächst die notwendigen Libraries 

```python
{{#include ../codes/02-differential_equations/rk4_bz.py:imports}}
```

und kopieren die Implementierung der Funktion `dydx`.

```python
{{#include ../codes/02-differential_equations/rk4_bz.py:dydx}}
```

#### RK4-Verfahren

Dann implementieren wir die Funktion `rk4_step`, die den Funktionswert
$y_{n+1}$ mit Hilfe des RK4-Verfahrens gemäß 
Gl. {{eqref: eq:runge_kutta_4_stage}} berechnet.

```python
{{#include ../codes/02-differential_equations/rk4_bz.py:rk4_step}}
```
Die Funktion sieht zwar auf den ersten Blick kompliziert aus, ein 
Großteil der Zeilen ist allerdings nur für die Definition der Koeffizienten
des RK4-Verfahrens belegt, wobei wir die weniger bekannte Variante des RK4-
Verfahrens verwendet haben. Danach werden die
vier Stufen $k_i$ berechnet und zum Schluss die Lösung $y_{n+1}$ entsprechend
Gl. {{eqref: eq:runge_kutta_4_stage}} ausgegeben.

Als nächstes implementieren wir die Funktion `rk4_method`:
```python
{{#include ../codes/02-differential_equations/rk4_bz.py:rk4_method}}
```
Diese Funktion ist tatsächlich identisch mit der Funktion `euler_method` aus
dem Abschnitt [2.2](02-euler_method.md##mn-zerfall), nur dass wir hier
`rk4_step` statt `euler_step` aufrufen. Man könnte auch eine allgemeingültige
Funktion `rk_method` schreiben, die als Argument `rk_step` sowohl `euler_step` als auch `rk4_step` akzeptiert.

Zum Schluss lösen wir das AWP mit dem RK4-Verfahren:
```python
{{#include ../codes/02-differential_equations/rk4_bz.py:solve_ode}}
```
Leider wird hier immer noch eine Schrittweite von `h = 0.001` benötigt, um
eine stabile Lösung zu erhalten. Dann plotten wir wieder die Lösung:
```python
{{#include ../codes/02-differential_equations/rk4_bz.py:plot}}
```
Optisch sollte das Ergebnis identisch wie das des Euler-Verfahrens sein:
![RK4 für Oregonator-Modell](../assets/figures/02-differential_equations/rk4_bz.svg)

Ein möglicher Grund, warum das RK4-Verfahren nur eine geringfügig größere
Schrittweite als das Euler-Verfahren verwenden kann, ist, dass das AWP
[steif](https://de.wikipedia.org/wiki/Steifes_Anfangswertproblem) ist,
also dass explizite Verfahren erhebliche Schwierigkeiten haben, eine stabile
Lösung zu finden. 

Um das AWP mit weniger Schritten zu lösen, können wir
z.B. adaptive Schrittweitenverfahren verwenden, welche die Schrittweite
an schwierigen Stellen automatisch verkleinern. Alternativ können wir
implizite Verfahren verwenden, die stabiler sind als explizite Verfahren.
Wir wollen uns hier jedoch nicht mit Details dieser Verfahren beschäftigen,
sondern lediglich diskutieren, wie und wann man sie einsetzen sollte. 
Deshalb werden wir im folgenden Abschnitt die Funktion `solve_ivp` aus der 
Bibliothek `scipy` verwenden, die eine Vielzahl von Verfahren zur Lösung 
von AWP bereitstellt.

#### Lösen von AWP mit `scipy.integrate.solve_ivp`
Die Funktion 
[`scipy.integrate.solve_ivp`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
bietet ein universelles Interface für eine Vielzahl von Verfahren zur Lösung
von AWP. Wir importieren diese Funktion sowie andere notwendige Libraries wie
gewohnt:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:imports}}
```
Wir nutzen hier wieder die Funktion `dydx` vom Oregonator-Modell.
Danach definieren wir die Anfangsbedingungen und die Parameter für den
DGL-Solver, genau so wie wir bisher gemacht haben:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:solve_setup}}
```
Ein Unterschied hier ist, dass wir die Konstante `STEP` durch `MAXSTEP` 
ersetzt haben, da der Algorithmus von `solve_ivp` die Schrittweite selbst
anpasst und wir nur eine obere Schranke setzen können.

Danach rufen wir die Funktion `solve_ivp` auf mit der Methode `RK45`:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:solve_rk45}} 
```
Als Ergebnis erhalten wir ein Objekt mit verschiedenen nützlichen Attributen.
Das Grid ${x_n}_{n=0,\cdots,N}$ ist im Attribut `t` gespeichert und die
Lösung $y_n$ im Attribut `y`. Die Anzahl der Schritte berechnen wir aus der
Länge des Grids abzüglich 1 (Anfangsbedingung). 

`RK45` ist ein adaptives Runge-Kutta-Verfahren mit einer Konsistenzordnung
von 5, weshalb die Schrittweite $h$ nicht konstant ist. Trotzdem können wir
die kleinste Schrittweite berechnen, indem wir zuerst die Differenz zwischen 
allen Gridpunkten mit
[np.diff](https://numpy.org/doc/stable/reference/generated/numpy.diff.html)
berechnen und dann das Minimum davon mit
[np.min](https://numpy.org/doc/stable/reference/generated/numpy.min.html)
bestimmen. Die minimale Schrittweite gibt an, mit welcher Präzision das Verfahren die
Lösung an den schwierigsten Stellen berechnet hat und bietet einen guten
Vergleich gegenüber Verfahren mit konstanter Schrittweite.

In diesem Fall werden nur knapp 30.000 Schritte benötigt, um eine stabile
Lösung des AWPs zu erhalten. Das ist deutlich weniger als die 200.000 Schritte, die unsere
Implementierung des klassischen RK4-Verfahrens benötigt hat. Die minimale Schrittweite beträgt dabei 
ca. 0.0012, also nur unwesentlich größer als die Schrittweite beim RK4-Verfahren.
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:verification_rk45}}
```
Die Erhöhung der Ordnung von 4 auf 5 hat dem Lösungsverfahren demnach
kaum geholfen; die adaptive Schrittweitensteuerung dagegen schon.

```admonish tip title="Tipp"
Testen Sie die Funktion `solve_ivp` mit den Argumenten `method='DOP853'` 
und `max_step=0.02`. DOP8(5,3) ist ein adaptives Runge-Kutta-Verfahren
mit einer Konsistenzordnung von 8. Sie werden feststellen, dass die
minimale Schrittweite bei ca. 0.002 liegt, also trotz der hohen Ordnung
immer noch sehr klein ist. Das bestätigt die Aussage, dass das AWP steif ist.
```

Nun probieren wir ein implizites Verfahren aus, z.B. mit `method='Radau'`
(und wieder `MAXSTEP=0.1`):
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:solve_radau}}
```
Dabei sollten Sie ungefähr die folgenden Werte für `nsteps` und `minstep` erhalten:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:verification_radau}}
```
Hier gibt es einen erheblicher Unterschied: Das Radau-Verfahren benötigt
nur knapp über 2000 Schritte und die minimale Schrittweite ist ca. 0.015.
Dass implizite Verfahren die Stabilität der Lösung signifikant verbessern,
ist eine weitere Eigenschaft von steifen AWPs.

Wir könnten die Lösung des AWPs der Belousov-Zhabotinsky-Reaktion an dieser Stelle erneut plotten, 
was uns allerdings keine neuen Erkenntnisse liefern wird. Stattdessen widmen wir uns zwei weiteren 
Visualierungsmethoden von Lösungen von AWPs: Konfigurationsraum- und Phasenraumtrajektorien.

#### Konfigurationsraum und Phasenraum

Der 
[*Konfigurationsraum*](https://de.wikipedia.org/wiki/Konfigurationsraum)
ist der Raum der Freiheitsgrade eines Systems. Für das
Oregonator-Modell sind das die Konzentrationen der drei Spezies, also 
$[\mathrm{X}]$, $[\mathrm{Y}]$ und $[\mathrm{Z}]$. Die Lösung des DGL-Systems 
zum Zeitpunkt $t_n$ ist also durch den Punkt 
$([\mathrm{X}](t_n), [\mathrm{Y}](t_n), [\mathrm{Z}](t_n))^\intercal$
im Konfigurationsraum gegeben. Die Zeitentwicklung des Systems kann demnach durch
eine Reihe von Punkten im Konfigurationsraum beschrieben. Die Ansammlung 
dieser Punkte wird dann als *Konfigurationsraumtrajektorie* bezeichnet.

Da die Information, zu welchem Zeitpunkt welcher Punkt im Konfigurationsraum durchlaufen
wird, verloren geht, wollen wir festlegen, dass der zeitliche Abstand zwischen zwei Punkten
in der Trajektorie konstant bleibt. Somit können wir wenigsten eine grobe Vorstellung
der Zeitentwicklung des Systems erhalten, da Punkte mit größerem Abstand im Konfigurationsraum
*schneller* durchlaufen werden. Das können wir erreichen, indem
wir das Argument `dense_output=True` an die Funktion `solve_ivp` übergeben.
Damit wird das Attribut `sol` des Rückgabewerts `res`ein
[`scipy.integrate.OdeSolution`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.OdeSolution.html)-Objekt
sein, was wie eine Funktion behandelt werden kann.
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:solve_dense}}
```
Hier definieren wir nach dem Lösen des AWPs ein gleichmäßiges Grid mit 
[`np.linspace`](https://numpy.org/doc/stable/reference/generated/numpy.linspace.html),
wobei wir 5000 gleichmäßig verteilte Punkte zwischen `T0` und `TMAX` wählen.
Die Lösung des AWPs an diesen Punkten erhalten wir dann durch Aufrufen der
Funktion `res.sol` mit dem Grid als Argument.

```admonish note title="Hinweis"
Alternativ kann mit dem Argument `t_eval` eine Liste von Zeitpunkten an `solve_ivp` 
übergeben werden, an denen die Lösung berechnet werden soll.
```

Nun können wir die Konfigurationsraumtrajektorie der Lösung mit den zeitlich gleichmäßig
verteilten Punkten plotten:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:configuration_space_plot}}
```

Da der Konfigurationsraum dreidimensional ist, müssen wir beim Aufrufen der 
Funktion `plt.subplots` das Argument `subplot_kw={'projection': '3d'}`
übergeben. Wir verwenden hier die Methode `scatter` anstatt von `plot`, um
die Punkte einzeln darzustellen. Mit dem Argument
`s` kann die Größe der Punkte eingestellt werden und mit `alpha` 
die Transparenz, wobei wir `alpha=0.1` (d.h. 10%) gewählt haben.
Die Methode
[`tight_layout`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tight_layout.html)
hat leider Schwierigkeiten mit 3D-Plots, weshalb wir hier den 
gewünschten Bereich des Plots mit dem Argument `rect=[0, 0, 0.95, 1.00]`
manuell angepasst haben.

Das Ergebnis sollte wie folgt aussehen:
<p align="center">
  <img src="../assets/figures/02-differential_equations/bz_configuration_space.svg" alt="Konfigurationsraumtrajektorie des Oregonator-Modells">
</p>

Durch die Einstellung `alpha=0.1` können wir jetzt erkennen, welche Bereiche
des Konfigurationsraums mit welcher Frequenz besucht werden.
Am Anfang der Reaktion ist nur $\mathrm{Br^-}$ vorhanden (hintere Ecke). Danach nimmt seine 
Konzentration ab, während die Konzentrationen der anderen Spezies erstmal auf einem Niveau nahe Null bleiben. 
Da dieses Ereignis nur einmal stattfindet, sind die Punkte auch nur sehr
schwach sichtbar. Danach beginnt die Oszillation des Systems. Auch hier kann 
man anhand der Farbstärke erkennen, dass die Änderung von $[\mathrm{HBrO_2}]$ (oberer Bogen)
schneller verläuft als die von $[\mathrm{Br^-}]$ (rechter Bogen).

Der Konfigurationsraum allein ist allerdings nicht ausreichend für eine vollständige
Beschreibung des Systems. Wir wissen z.B. bei einer gegebenen
$\mathrm{Ce^{4+}}$-Konzentration nicht, ob diese gerade steigt oder fällt, d.h. in welcher
Richtung das System sich auf der geschlossenen Kurve bewegt.
Um das zu bestimmen, benötigen wir zusätzlich die "Geschwindigkeiten" 
oder "Impulse" der Koordinaten. Ein Raum, der sowohl die Koordinaten als
auch die Geschwindigkeiten enthält, wird als 
[*Phasenraum*](https://de.wikipedia.org/wiki/Phasenraum) bezeichnet.

Der Phasenraum des Oregonator-Modells ist also sechsdimensional, was über
die Grenzen des menschlichen Vorstellungsvermögens hinausgeht. Deshalb 
plotten wir hier einen zweidimensionalen Schnitt durch den Phasenraum, indem wir nur 
die Konzentrationen von $[\mathrm{Ce^{4+}}]$ und ihre Ableitungen zeigen. Auch hier verwenden wir 
die gleichmäßig ausgewertete Lösung und berechnen die Ableitung mit der 
Funktion `dydx`:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:phase_space_plot}}
```
Hier plotten wir die Ableitung `dzdt` gegen die Konzentration `c_z`, ebenfalls
mit der Methode `scatter` und dem Argument `alpha=0.1`. Das Diagramm sollte
wie folgt aussehen:

<p align="center">
  <img src="../assets/figures/02-differential_equations/bz_phase_space.svg" alt="Phasenraum des Oregonator-Modells">
</p>

Die Trajektorie startet hier im Ursprung, verläuft dann im Uhrzeigersinn, 
also sowohl die Konzentration als auch ihre Ableitung steigen zunächst. Dann wird 
der Punkt der maximalen Zunahme erreicht, während die Konzentration weiter steigt.
Zu einem späteren Punkt wird die Ableitung schließlich negativ und die Konzentration nimmt leicht ab.
Besonders auffällig ist das letzte Stück der Trajektorie:
Die Konzentration nimmt stetig ab, während die Ableitung von stark negativ zu
null wird. Dieser Bereich im Phasenraum ähnelt sehr stark einer Graden, was einer Kinetik erster Ordnung entspricht.
In diesem Bereich gilt $\frac{\du [\mathrm{Z}]}{\du t} \propto [\mathrm{Z}]$, was einen exponentiellen Zerfall 
darstellt. Die intensivste Färbung in diesem
Plot ist um den Ursprung, was bedeutet, dass für die meisten Zeitpunkte
die Konzentration von $\mathrm{Ce^{4+}}$ und ihre Ableitung sehr klein sind.

```admonish tip title="Tipp"
Zoomen Sie in den interaktiven Plot des Phasenraums hinein, um die Details
der Trajektorie besser zu erkennen. Nehmen Sie dazu das Diagramm des zeitlichen Verlaufes
der Konzentration als Vergleich und versuchen Sie, die Merkmale in der
Phasenraumtrajektorie wiederzuerkennen.
```

---

### Übung

#### Aufgabe 2.2: Lösen des klassischen harmonischen Oszillators mit Runge-Kutta-Verfahren

{{#include ../psets/02.md:aufgabe_2}}

