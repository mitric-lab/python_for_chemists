## Runge-Kutta-Verfahren

Die Runge-Kutta-Verfahren sind eine
nach Carl Runge und Martin Wilhelm Kutta benannte
Familie von Methoden zur numerischen 
Lösung von Anfangswertproblemen für gewöhnliche Differentialgleichungen.
Diese Methoden berechnen iterativ die Lösung im nächsten Zeitschritt aus 
einer Linearkombination von dem Funktionswert und den Steigungen an
verschiedenen Stellen. Sprich man jedoch vom **dem** Runge-Kutta-Verfahren,
ist das 
[klassische Runge-Kutta-Verfahren](https://de.wikipedia.org/wiki/Klassisches_Runge-Kutta-Verfahren)
gemeint, was aber nur ein Spezialfall der Runge-Kutta-Verfahren darstellt.

### Theoretische Grundlagen
Genau so wie beim Euler-Verfahren, wählen wir zunächst ein gleichmäßiges
Grid $\{x_n\}_ {n=0,\cdots,N}$ mit $x_n = x_0 + n\cdot h$ und $h$ als
Schrittweite, sowie eine Anfangsbedingung $y_0 = y(x_0)$.
Dann entwickeln wir wieder $y(x_{n+1})$ in eine Taylor-Reihe um $y(x_n)$:
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
  k_i = f(x_n + h c_i, y_n + h \sum_{j=1}^{s} a_{ij} k_j)
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
    k_2 &= f(x_n + h c_2, y_n + h a_{21} k_1)
  \end{align}
  {{numeq}}{eq:runge_kutta_2_stage}
$$

Für ein Verfahren zweiter Ordnung benötigen wir ein Taylorpolynom zweiten
Grades, wo die Ableitung von $y$ vorkommt. Da wir aber $y$ nicht kennen,
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
Setzen wir dies in Gl. {{eq: eq:runge_kutta_2_stage}} ein, erhalten wir
$$
  y_{n+1} = y_n + h(b_1 + b_2) f(x_n, y_n) + h^2 b_2 \left[
    c_2 \frac{\partial f}{\partial x}(x_n, y_n) + a_{21} \left(
      f \frac{\partial f}{\partial y}
    \right)(x_n, y_n)
  \right] + \mathcal{O}(h^3)\,.\\
  {{numeq}}{eq:runge_kutta_2_stage_2nd_order}
$$

Damit das Verfahren tatsächlich eine Konsistenzordnung von 2 hat, müssen die
Koeffizienten vor den Funktionen $f$ und ihren Ableitungen in Gl.
{{eqref: eq:ode_taylor_polynomial_2nd_order}} und
{{eqref: eq:runge_kutta_2_stage_2nd_order}} übereinstimmen, da diese für beliebige
Funktionen $f$ gelten müssen. Das führt zu den Bedingungen
$$
  \begin{align}
    b_1 + b_2 &= 1 \\
    b_2 c_2 &= \frac{1}{2} \\
    b_2 a_{21} &= \frac{1}{2}
  \end{align}
  {{numeq}}{eq:runge_kutta_2_stage_conditions}
$$
Das ist ein unterbestimmtes Gleichungssystem mit drei Gleichungen für vier
Unbekannte. Wir erhalten also eine Familie von konsistenten
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
    c_2 &= \alpha
  \end{align}\,.
$$

Die Angabe der einzelnen Koeffizienten ist jedoch nicht sehr übersichtlich,
wenn man bedenkt, dass es Verfahren mit mehr Stufen gibt. Daher gibt es eine
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

wobei die Koeffizienten mit einer Index ($b_i$ und $c_i$) als ein Vektor
und die Koeffizienten mit zwei Indices ($a_{ij}$) als eine Matrix
dargestellt werden.

Damit kann man das Heun-Verfahren mit
$$
  \begin{array}{c|cc}
    0 & 0 & 0 \\
    1 & 1 & 0 \\ \hline
      & 1/2 & 1/2
  \end{array}
$$
und das Mittelpunktsverfahren mit
$$
  \begin{array}{c|c}
    0 & 0 & 0 \\
    1/2 & 1/2 & 0 \\ \hline
        & 0 & 1
  \end{array}
$$
darstellen. Die allgemein parametrisierte Form eines 2-stufigen Verfahrens
zweiter Ordnung lautet dann
$$
  \begin{array}{c|cc}
    0 & 0 & 0 \\
    \alpha & \alpha & 0 \\ \hline
      & 1 - \frac{1}{2\alpha} & \frac{1}{2\alpha}
  \end{array}\,.
$$

Wenn man es möchte, kann man das Euler-Verfahren als ein einstufiges
Runge-Kutta-Verfahren erster Ordnung auffassen mit dem Tableau
$$
  \begin{array}{c|c}
    0 & 0 \\ \hline
      & 1
  \end{array}\,.
$$


Als eine letzte Bemerkung sei noch gesagt, dass die Bedingung $a_{ij} = 0$ 
für $j \geq i$ bei expliziten Runge-Kutta-Verfahren in dieser Darstellung
zu der Bedingung übersetzt wird, dass die Matrix $\bm{A}$ eine strikte
untere Dreiecksmatrix ist.

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
    k_4 &= f(x_n + h c_4, y_n + h(a_{41} k_1 + a_{42} k_2 + a_{43} k_3))
  \end{align}
  {{numeq}}{eq:runge_kutta_4_stage}
$$
Die Bedingungen für die Koeffizienten kann auf die gleiche Weise wie beim
RK2 hergeleitet werden, die Rechnungen sind jedoch deutlich aufwendiger, 
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

Der Erfolg von RK4 mag auf den ersten Blick überraschen, da es prinzipiell
möglich ist, Verfahren mit beliebig hoher Ordnung zu konstruieren. Allerdings
soll man nicht vergessen: 
```admonish warning title="Ordnung &ne; Stufe!"
Die Konsistenzordnung $p$ und die Stufe $s$ des Verfahrens sind zwei 
verschiedene Dinge, obwohl wir bis jetzt nur Verfahren mit
$s = p$ betrachtet haben. Tatsächlich gilt für minimale 
Stufenzahl $s_{\mathrm{min}}$ zum Erreichen einer Konsistenzordnung $p$
bei expliziten Runge-Kutta-Verfahren stets
$s_{\mathrm{min}} \geq p$.
```

Man kann sogar zeigen, dass für $p \geq 5$ die strikte Ungleichung
$s_{\mathrm{min}} > p$ gilt.<sup>[</sup>[^butcher1987]<sup>]</sup>
In anderen Worten: Die Verbesserung der Genauigkeit von $p = 4$
auf $p = 5$ unter Verwendung von expliziten Runge-Kutta-Verfahren
ist mit einer Erhöhung der Stufezahl um mindestens 2
verbunden. Das ist ein möglicher Grund, warum RK4 so beliebt ist.

Die Korrespondenz zwischen $p$ und $s_{\mathrm{min}}$ für einige Ordnungen
bei expliziten Runge-Kutta-Verfahren
ist in der folgenden Tabelle zusammengefasst:<sup>[</sup>[^butcher1987]<sup>]</sup>
|                    |    |    |    |    |    |    |    |    |
|--------------------|---:|---:|---:|---:|---:|---:|---:|---:|
| $p$                |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |
| $s_{\mathrm{min}}$ |  1 |  2 |  3 |  4 |  6 |  7 |  9 | 11 |

Die Zahlen $s_{\mathrm{min}}$ sind als *Butcher-Schranken* bekannt.

[^kutta1901]: M. W. Kutta, *Z. Math. Phys.* **1901**, *46*, 435&ndash;453.
[^butcher1987]: J. C. Butcher, in *The Numerical Analysis of Ordinary Differential Equations*, John Wiley & Sons, Chichester, **1987**, pp. 185&ndash;194.

### Implementierung
Wir nehmen als Beispiel die Dynamik der Belousov-Zhabotinsky-Reaktion her.
Genau wie im Abschnitt [2.2](02-euler_method.md#belousov-zhabotinsky-reaktion)
importieren wir die notwendigen Libraries und kopieren die Implementierung
der Funktion `dydx`

```python
{{#include ../codes/02-differential_equations/rk4_bz.py:imports}}
```
```python
{{#include ../codes/02-differential_equations/rk4_bz.py:dydx}}
```

#### RK4-Verfahren

Dann implementieren wir die Funktion `rk4_step`, die den Funktionswert
$y_{n+1}$ mit Hilfe des RK4-Verfahrens in 
Gl. {{eqref: eq:runge_kutta_4_stage}} berechnet.

```python
{{#include ../codes/02-differential_equations/rk4_bz.py:rk4_step}}
```
Die Funktion sieht zwar auf den ersten Blick kompliziert aus, aber ein 
Großteil der Zeilen sind nur für die Definition der Koeffizienten
des RK4-Verfahrens. Hier haben wir die weniger bekannte Variante des RK4-
Verfahrens verwendet. Nach der Definition der Koeffizienten werden die
vier Stufen $k_i$ berechnet und zum Schluss die Lösung $y_{n+1}$ nach
Gl. {{eqref: eq:runge_kutta_4_stage}} ausgegeben.

Als nächstes implementieren wir die Funktion `rk4_method`:
```python
{{#include ../codes/02-differential_equations/rk4_bz.py:rk4_method}}
```
Diese Funktion ist tatsächlich identisch mit der Funktion `euler_method` aus
dem Abschnitt [2.2](02-euler_method.md##mn-zerfall), nur dass wir hier
`rk4_step` statt `euler_step` aufrufen. Man könnte auch eine allgemeine
Funktion `rk_method` schreiben, die als Argument die Funktion `rk_step`,
wo man sowohl `euler_step` als auch `rk4_step` übergeben kann.

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
z.B. adaptive Schrittweitenverfahren verwenden, die die Schrittweite
bei schwierigen Stellen automatisch verkleinern und bei unproblematischen
Stellen vergrößern, oder implizite Verfahren verwenden.
Wir wollen uns hier jedoch nicht mit Details dieser Verfahren beschäftigen,
sondern nur wie und wann man sie einsetzen sollte. Deshalb werden wir im
folgenden Abschnitt die Funktion `solve_ivp` aus der Bibliothek `scipy`
verwenden, die eine Vielzahl von Verfahren zur Lösung von AWP bereitstellt.

#### `scipy.integrate.solve_ivp`
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
Solver, genau so wie wir bis jetzt gemacht haben:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:solve_setup}}
```
Ein Unterschied hier ist, dass wir die Konstante `STEP` durch `MAXSTEP` 
ersetzt haben, da die Algorithmus von `solve_ivp` die Schrittweite selbst
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
von 5, deshalb ist die Schrittweite $h$ nicht konstant. Trotzdem können wir
die kleinste Schrittweite berechnen, indem wir zuerst die Differenz zwischen 
allen Gridpunkten mit
[np.diff](https://numpy.org/doc/stable/reference/generated/numpy.diff.html)
berechnen und dann das Minimum davon mit
[np.min](https://numpy.org/doc/stable/reference/generated/numpy.min.html)
bestimmen. Die minimale Schrittweite gibt an, wie fein das Verfahren die
Lösung an den schwierigsten Stellen berechnet hat und stellt einen guten
Vergleich mit der Schrittweite bei Verfahren mit konstanter Schrittweite dar.

In diesem Fall werden nur knapp 30.000 Schritte benötigt, um eine stabile
Lösung des AWPs zu erhalten, im Vergleich zu den 200.000 Schritten beim
RK4-Verfahren. Die minimale Schrittweite beträgt dabei aber ca. 0.0012,
also nur unwesentlich größer als die Schrittweite beim RK4-Verfahren.
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:verification_rk45}}
```
Also eine Erhöhung der Ordnung von 4 auf 5 hat dem Lösungsverfahren
kaum geholfen, die adaptive Schrittweitensteuerung dagegen schon.

```admonish tip title="Tipp"
Probieren Sie die Funktion `solve_ivp` mit dem Argument `method='DOP853'` 
und `max_step=0.02` aus. DOP8(5,3) ist ein adaptives Runge-Kutta-Verfahren
mit einer Konsistenzordnung von 8. Sie werden aber feststellen, dass die
minimale Schrittweite bei ca. 0.002 liegt, also trotz der hohen Ordnung
immer noch sehr klein ist. Das bestätigt die Aussage, dass das AWP steif ist.
```

Nun probieren wir ein implizites Verfahren aus, z.B. mit `method='Radau'`
(und wieder `MAXSTEP=0.1`):
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:solve_radau}}
```
Sie sollten ungefähr die folgenden Werte für `nsteps` und `minstep` erhalten:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:verification_radau}}
```
Hier gibt es ein erheblicher Unterschied: Das Radau-Verfahren benötigt
nur knapp über 2000 Schritte und die minimale Schrittweite ist ca. 0.015.
Das ist ein weitere Bestätigung für die Steifheit des AWPs, da durch die
eines impliziten Verfahrens die Stabilität der Lösung signifikant verbessert
wird.

Wir können hier natürlich auch die Lösung plotten, aber das haben wir
schon oft genug gesehen. Stattdessen wollen wir uns zwei anderen Plots
anschauen: Konfigurationsraumtrajektorie und Phasenraumtraektorie.

#### Konfigurationsraum und Phasenraum

Der 
[*Konfigurationsraum*](https://de.wikipedia.org/wiki/Konfigurationsraum)
ist der Raum der Koordinaten eines Systems. Für das
Oregonator-Modell sind das die Konzentrationen der drei Spezies, also 
$[\mathrm{X}]$, $[\mathrm{Y}]$ und $[\mathrm{Z}]$. Die Lösung des DGL-Systems 
zum Zeitpunkt $t_n$ ist also durch den Punkt 
$([\mathrm{X}](t_n), [\mathrm{Y}](t_n), [\mathrm{Z}](t_n))^\intercal$
im Konfigurationsraum gegeben. Die Zeitentwicklung des Systems wird also durch
eine Reihe von Punkten im Konfigurationsraum beschrieben. Die Ansammlung 
dieser Punkte wird dann als *Konfigurationsraumtrajektorie* bezeichnet.

Da die Zeitinformation in der Konfigurationsraumtrajektorie verloren geht,
wollen wir feststellen, dass der zeitliche Abstand zwischen zwei Punkten
in der Trajektorie konstant bleibt, damit wir wenigsten eine grobe Vorstellung
von der Zeitentwicklung des Systems bekommen. Das können wir erreichen, indem
wir das Argument `dense_output=True` an die Funktion `solve_ivp` übergeben.
Dann wird das Attribut `sol` des Rückgabewerts ein
[`scipy.integrate.OdeSolution`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.OdeSolution.html)-Objekt
zurück, was wie eine Funktion aufgerufen werden kann:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:solve_dense}}
```
Nachdem Lösen des AWPs definieren wir ein gleichmäßiges Grid mit 
[`np.linspace`](https://numpy.org/doc/stable/reference/generated/numpy.linspace.html),
wobei wir 5000 gleichmäßig verteilte Punkte zwischen `T0` und `TMAX` wählen.
Die Lösung des AWPs an diesen Punkten erhalten wir durch Aufrufen der
Funktion `res.sol` mit dem Grid als Argument.

```admonish note title="Hinweis"
Mit dem Argument `t_eval` kann eine Liste von Zeitpunkten übergeben werden,
an denen die Lösung berechnet werden soll. Damit könnten wir das gleiche 
erreichen. 
```

Nun können wir die Konfigurationsraumtrajektorie mit der zeitlich gleichmäßig
ausgewerteten Lösung plotten:
```python
{{#include ../codes/02-differential_equations/scipy_bz.py:configuration_space_plot}}
```

Da der Konfigurationsraum dreidimensional ist, müssen wir beim Aufrufen der 
Funktion `plt.subplots` das Argument `subplot_kw={'projection': '3d'}`
übergeben. Danach wurde die Methode `scatter` verwendet anstatt `plot`, um
die Punkte einzeln darzustellen ohne sie zu verbinden. Mit dem Argument
`s` kann die Größe der Punkte eingestellt werden und mit `alpha` 
die Transparenz. Hier haben wir `alpha=0.1` gewählt, damit die Punkte
nur schwach sichtbar sind.
Die Methode
[`tight_layout`](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.tight_layout.html)
hat bisschen Schwierigkeiten mit 3D-Plots, deshalb haben wir hier den 
gewünschten Bereich des Plots mit dem Argument `rect=[0, 0, 0.95, 1.00]`
manuell angepasst.

Das Ergebnis sollte wie folgt aussehen:
<p align="center">
  <img src="../assets/figures/02-differential_equations/bz_configuration_space.svg" alt="Konfigurationsraumtrajektorie des Oregonator-Modells">
</p>

Durch die Einstellung `alpha=0.1` können wir jetzt erkennen, welche Bereiche
des Konfigurationsraums häufiger besucht werden und welche weniger.
Am Anfang der Reaktion ist nur $\mathrm{Br^-}$ vorhanden. Danach nimmt seine 
Konzentration ab, während die Konzentration anderer Spezies erstmal konstant
bleiben. Da dieses Ereignis nur einmal stattfindet, sind die Punkte nur sehr
schwach sichtbar. Danach beginnt die Oszillation des Systems. Auch hier kann 
man an der Farbstärke erkennen, dass die Änderung von $[\mathrm{HBrO_2}]$
schneller verläuft als die von $[\mathrm{Br^-}]$.

Der Konfigurationsraum ist aber nicht ausreichend für eine vollständige
Beschreibung des Systems. Wir wissen z.B. nicht, bei einer gegebenen
$\mathrm{Ce^{4+}}$-Konzentration, ob diese gerade steigt oder fällt.
Um das zu bestimmen, benötigen wir zusätzlich die "Geschwindigkeiten" 
oder "Impulse" der Koordinaten. Ein Raum, der sowohl die Koordinaten als
auch die Geschwindigkeiten enthält, wird als 
[*Phasenraum*](https://de.wikipedia.org/wiki/Phasenraum) bezeichnet.

Der Phasenraum des Oregonator-Modells ist also sechsdimensional, was über
die Grenze des menschlichen Vorstellungsvermögens hinausgeht. Deshalb plotten
wir hier einen 2D-Schnitt durch den Phasenraum, indem wir die Konzentrationen
von $[\mathrm{Ce^{4+}}]$ und ihre Ableitungen plotten. Auch hier verwenden wir 
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
also sowohl die Konzentration als auch ihre Ableitung steigen. Danach wird 
ein Maximum der Ableitung erreicht, während die Konzentration weiter steigt.
Danach wird die Ableitung negativ und die Konzentration nimmt leicht ab.
Besonders auffällig ist das letzte Stück der Trajektorie in der Periode,
wo die Konzentration absteigt, während die Ableitung von stark negativ zu
null wird. Dieses Stück im Phasenraum ähnelt sehr stark einer Grade, was
$\frac{\du [\mathrm{Z}]}{\du t} \propto [\mathrm{Z}]$ bedeutet. Das stellt 
garede einen exponentiellen Zerfall dar. Die intensivsten Farben in diesem
Plot ist um den Ursprung, was bedeutet, dass die für die meisten Zeitpunkte
die Konzentration von $\mathrm{Ce^{4+}}$ und ihre Ableitung sehr klein sind.

```admonish tip title="Tipp"
Zoomen Sie in den interaktiven Plot des Phasenraums hinein, um die Details
der Trajektorie besser zu erkennen. Nehmen Sie dazu das Diagramm der Lösung
des AWPs als Vergleich und versuchen Sie, die Merkmale der Lösung in der
Phasenraumtrajektorie wiederzuerkennen.
```

---

### Übung

#### Aufgabe 2.2
<!--
{{#include ../psets/02.md:aufgabe_2}}
-->
