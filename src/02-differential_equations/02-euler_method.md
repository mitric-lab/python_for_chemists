## Euler-Verfahren

Wir betrachten das Anfangswertproblem
$$
\begin{cases}
  y'(x) = f(x, y(x))\\ 
  y(x_0) = y_0
\end{cases}\,,
{{numeq}}{eq:ivp_first_order}
$$
wobei $f(x, y(x))$ eine gegebene Funktion von der unabhängigen Variablen $x$ 
und der gesuchten Funktion $y(x)$ ist. Gl. {{eqref: eq:ivp_first_order}} ist ein
AWP erster Ordnung, mit einer DGL erster Ordnung wie in 
Gl. {{eqref: eq:ode_first_order}} und einer Anfangsbedingung.

Unser Ziel ist es, die spezielle Lösung $y(x)$ des AWP numerisch zu finden.
Wie bei vielen numerischen Verfahren, beginnen wir mit einer Diskretisierung
der Funktion $y(x)$, d.h. wir wählen eine Menge von Punkten $x_i$ und
betrachten die Funktion nur an diesen Punkten anstatt auf dem gesamten
Definitionsbereich. Die konzeptionell wohl einfachste Wahl der Punkte ist
ein gleichmäßiges Gitter (Grid), was bedeutet, dass man einen Anfangspunkt
$x_0$ und eine Schrittweite $h$ wählt, so dass die weiteren
Punkte durch $x_n = x_0 + n h$ für $n = 0, 1, 2, \ldots$ festgelegt werden.
Unser Ziel in den folgenden Abschnitten wird es sein, die Funktionswerte von
$y(x)$ an diesen Punkten $x_1, x_2, \ldots$ zu berechnen, bzw. zu approximieren.

Der Funktionswert von $y$ an dem Punkt $x_{n+1}$ kann (unter bestimmten
Voraussetzung an $y(x)$) durch eine Taylor-Entwicklung um den Punkt
$x_n$ ersetzt werden, also
$$
  y(x_{n+1}) = y(x_n) 
  + h y'(x_n) 
  + \frac{h^2}{2} y''(x_n) 
  + \frac{h^3}{6} y'''(x_n)
  + \cdots\,.
$$

```admonish info title="Info für Mathematik-Interessierte" collapsible=true
Die Voraussetzung ist, dass $y(x)$ im Punkt $x_n$ analytisch ist, d.h. dass es
eine Potenzreihe
$$
  \sum_{k=0}^\infty a_k (x - x_n)^k
$$
gibt, die für alle $\|x - x_n\| < R$ konvergiert, wobei $R \in \mathbb{R}^+$ 
der Konvergenzradius der Potenzreihe ist. Zudem muss der Konvergenzradius
größer als die Schrittweite $h$ sein, also $h < R$.
```

### Theoretische Grundlagen

Gehen wir nun davon aus, dass die Funktion $y(x)$ gut durch ihr
Taylor-Polynom 1. Ordnung
$$
  y(x_{n+1}) = y(x_n) + h y'(x_n) + \mathcal{O}(h^2)
$$
approximiert werden kann, also dass der Fehler $\mathcal{O}(h^2)$, der proportional zu $h^2$ ist,
klein ist. Mathematisch unsauber, aber praktisch für die Implementierung, schreiben wir 
im Folgenden 
$$
  y(x_{n+1}) = y(x_n) + h y'(x_n)\,,
  {{numeq}}{eq:ode_taylor_first_order}
$$
wobei wir den Fehler vernachlässigen. 
Gl. {{eqref: eq:ode_taylor_first_order}} zeigt, dass wir, wenn wir den
Funktionswert $y$ und die Ableitung $y'$ an dem Punkt $x_n$ kennen,
den Funktionswert an dem nächsten Punkt $x_{n+1}$
berechnen können. Mit Hilfe der DGL in Gl. {{eqref: eq:ivp_first_order}} können
wir dabei die Ableitung $y'(x_n)$ durch $f(x_n, y(x_n))$ ersetzen, was zu
$$
  y(x_{n+1}) = y(x_n) + h f(x_n, y(x_n))
  {{numeq}}{eq:euler_method}
$$
führt. Gl. {{eqref: eq:euler_method}} beschreibt das *explizite Euler-Verfahren* 
für die numerische Lösung eines AWP erster Ordnung. Durch die Anfangsbedingung
kennen wir $y(x_0)$ und können dann, Schritt für Schritt, alle weiteren
$y(x_{n})$ mit Hilfe von Gl. {{eqref: eq:euler_method}} berechnen.

### Implementierung

#### Mn-Zerfall

Wir nehmen wieder den Mn-Zerfall als Beispiel. Dort haben wir eine 
Reaktion 1. Ordnung, welche durch die DGL
$$
  c'(t) = -k c(t)
  {{numeq}}{eq:mn_decay}
$$
beschrieben wird, wobei $c(t)$ die Konzentration des Mn-Komplexes zum Zeitpunk $t$ beschreibt. 
Wir können nach dem Importieren der benötigten Libraries
```python
{{#include ../codes/02-differential_equations/euler_mn.py:imports}}
```
die Funktion `dydx` implementieren, die $f(x_n, y(x_n))$ berechnet:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:dydx}}
```
Diese Funktion akzeptiert die Argumente `x` ($t_n$) und `y` ($c(t_n)$) und
gibt die Ableitung $c'(t_n) = -k c(t_n)$ zurück. Hier haben wir die gefittete
Geschwindigkeitskonstante $k$ aus Abschnitt 
[1.4](../01-regression/04-nonlinear_regression.md#reaktionskinetik) verwendet.
Anschließend können wir gemäß Gl. {{eqref: eq:euler_method}} die Funktion `euler_step`
implementieren, die den Funktionswert $y(x_{n+1})$ berechnet:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:euler_step}}
```

Wir können jetzt das Euler-Verfahren implementieren:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:euler_method}}
```

Diese Funktion akzeptiert neben der Anfangsbedingungen `x0` und `y0` die
Schrittweite `h`, die Ableitungsfunktion `dydx` und die Anzahl der Schritte
`n`. Wir erstellen zunächst das Grid `x` mit der Funktion 
[`np.arange`](https://numpy.org/doc/stable/reference/generated/numpy.arange.html)
und initialisieren das Nullarray `y`, um später die Lösung zu speichern. Dann wird 
der erste Eintrag `y[0]` dieses Arrays mit 
dem Anfangswert `y0` überschrieben. Anschließend verwenden wir eine for-Schleife
über die Anzahl der Schritte und rufen in jedem Schritt die Funktion
`euler_step` auf, die den Funktionswert an dem jeweils nächsten Punkt berechnet, 
und speichern diesen in `y[i + 1]`. Am Ende wird das Grid `x` und die Lösung `y`
zurückgegeben.

Nun wenden wir das Euler-Verfahren auf Gl. {{eqref: eq:mn_decay}} an:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:solve_ode}}
```
Wir setzen dazu zunächst die Anfangsbedingungen `C0 = 1.0` und `T0 = 0.0`, die Schrittweite
`h = 1.0` sowie die maximale Zeit `MAXTIME = 600.0`. Die Anzahl der Schritte
`nsteps` wird durch `int(MAXTIME / h)` berechnet. Die `int`-Funktion rundet
dabei das Ergebnis der Division ab und konvertiert es in eine Ganzzahl. Anschließend rufen
wir die Funktion `euler_method` auf und speichern das Ergebnis in `x` und `y`.

Zum Schluss können wir das numerische Ergebnis mit der analytischen Lösung
vergleichen. Dabei sei angemerkt, dass die analytische Lösung von Gl. 
{{eqref: eq:mn_decay}} leicht zu berechnen ist und durch $c(t) = c_0 \exp(-k t)$
gegeben ist.
```python
{{#include ../codes/02-differential_equations/euler_mn.py:plot}}
```
Dabei sollte dieses Diagramm erscheinen:
![Euler-Verfahren für Mn-Zerfall](../assets/figures/02-differential_equations/euler_mn.svg)

Erst durch Hineinzoomen können wir den Unterschied zwischen der analytischen
und der numerischen Lösung erkennen. Das Euler-Verfahren mit $h = 1$ liefert
demnach eine sehr gute Approximation.

#### Belousov-Zhabotinsky-Reaktion

Nun wollen wir etwas kompliziertere Kinetik betrachten. Die 
[Belousov-Zhabotinsky-Reaktion](https://de.wikipedia.org/wiki/Belousov-Zhabotinsky-Reaktion)
ist eine klassisches Beispiel für eine oszillierende Reaktion. 
Dabei wird ein Redox-System ($\mathrm{Ce^{3+}/Ce^{4+}}$ in der Originalreaktion)
zum Oszillieren zwischen der oxidierten und der reduzierten Form gebracht.
In dem untenstehende Video können Sie die Reaktion beobachten:

<iframe width="750" height="422" src="https://www.youtube.com/embed/kw9wF-GNjqs?start=141&end=154" title="Everything about the Belousov Zhabotinsky reaction" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
<!--
<div onclick="this.nextElementSibling.style.display='block'; this.style.display='none'">
   <img src="../assets/figures/02-differential_equations/bz_thumbnail.png" style="cursor:pointer" />
</div>
<div style="display:none">
<iframe width="750" height="422" src="https://www.youtube.com/embed/kw9wF-GNjqs?start=141&end=154" title="Everything about the Belousov Zhabotinsky reaction" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
</div>
-->

*Hier wird als Redox-System [Ferroin](https://de.wikipedia.org/wiki/Ferroin) 
verwendet, welches eine stärkere Farbänderung zeigt.*

Der Mechanismus dieser Reaktion ist sehr
kompliziert, weshalb wir hier nur eine vereinfachte Version, das sog.
[Oregonator-Modell](https://de.wikipedia.org/wiki/Oregonator), betrachten.
Ein häufig verwendetes Oregonator-Modell besteht aus fünf gekoppelten
Reaktionen mit sechs Spezies:<sup>[</sup>[^field1986]<sup>, </sup>[^schneider1996]<sup>]</sup>

$$
\begin{align}
  \mathrm{A} + \mathrm{Y} &\xrightarrow{k_1} \mathrm{X} + \mathrm{P} \\
  \mathrm{X} + \mathrm{Y} &\xrightarrow{k_2} 2\ \mathrm{P} \\
  \mathrm{A} + \mathrm{X} &\xrightarrow{k_3} 2\ \mathrm{X} + \mathrm{Z} \\
  2\ \mathrm{X} &\xrightarrow{k_4} \mathrm{A} + \mathrm{P} \\
  \mathrm{B} + \mathrm{Z} &\xrightarrow{k_5} \mathrm{Y}
\end{align}
$$
mit der Identifikation
$\mathrm{A} = \mathrm{BrO_ 3^-}$, 
$\mathrm{B} = \mathrm{Brommalonsäure}$, 
$\mathrm{X} = \mathrm{HBrO_ 2}$, 
$\mathrm{Y} = \mathrm{Br^-}$, 
$\mathrm{Z} = \mathrm{2\ Ce^{4+}}$ 
und $\mathrm{P} = \mathrm{HOBr}$. 
Möglicherweise fällt Ihnen auf, dass die
obigen Reaktionsgleichungen nicht ausbalanciert sind. Das liegt daran, dass
dieses annimmt, dass die Konzentrationen dieser Spezies entweder aufgrund
ihrer hohen Konzentration ($\mathrm{H^+}$, Malonsäure, etc.) oder ihrer schnellen
Reaktionen ($\mathrm{BrO_ 2^\cdot}$) als konstant angenommen werden können
und in die Geschwindigkeitskonstanten $k_i$ absorbiert werden.
Außerdem wird angenommen, dass die Konzentrationen der Spezies $[\mathrm{A}]$
und $[\mathrm{B}]$ konstant sind.

Das führt zu einem System von drei gekoppelten DGLs:
$$
  \begin{align}
    [X]'(t) &= k_1 c_A^0 [Y] - k_2 [X] [Y] + k_3 c_A^0 [X] - 2 k_4 [X]^2 \\
    [Y]'(t) &= -k_1 c_A^0 [Y] - k_2 [X] [Y] + k_5 c_B^0 [Z] \\
    [Z]'(t) &= k_3 c_A^0 [X] - k_5 c_B^0 [Z]
  \end{align}
    {{numeq}}{eq:oregonator}
$$

Ein wichtiger Unterschied zu Mn-Zerfall ist, dass wir hier ein System von
DGLs haben. Dankenswerterweise ist Gl. {{eqref: eq:euler_method}} für
DGL-Systeme genauso gültig wie für einzelne DGLs, wenn wir für $y(x)$
eine vektorwertige Funktion einsetzen, also in diesem Fall
$y(x) \equalhat ([X](t), [Y](t), [Z](t))^\intercal$.

Nach dem Importieren der benötigten Libraries
```python
{{#include ../codes/02-differential_equations/euler_bz.py:imports}}
```
können wir die Funktion `dydx` für das Oregonator-Modell z.B. so 
implementieren:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:dydx}}
```
Beachten Sie, dass der Datentyp für das Argument `y` ein `np.ndarray` ist.
Dann implementieren wir die Funktion `euler_step` und `euler_method` genau
so wie für den Mn-Zerfall, nur dass wir hier die Typ-Deklaraionen in 
den Signaturen angepasst werden müssen.
```python
{{#include ../codes/02-differential_equations/euler_bz.py:euler_step}}
```
```python
{{#include ../codes/02-differential_equations/euler_bz.py:euler_method}}
```

Ein weiterer Unterschied ist, dass die Variablen `y` in der Funktion
`euler_method` als ein 2D-Array initialisiert wird, wo die erste Dimension
die Anzahl der Komponenten und die zweite Dimension die Anzahl der Schritte
ist. Das Grid `x` bleibt ein 1D-Array, da die Zeit für alle Komponenten
gleich ist. 

Nun lösen wir das DGL-System in Gl. {{eqref: eq:oregonator}}:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:solve_ode_bad}}
```
Hier haben wir die Anfangsbedingungen `C0 = np.array([0.0, 0.001, 0.0])` 
definiert, welche bedeutet, dass in der Anfangsmischung nur die Spezies 
$\mathrm{Br^-}$ mit einer Konzentration von $0.001\ \mathrm{M}$ vorhanden ist, 
und die Spezies $\mathrm{HBrO_ 2}$ und $\mathrm{Ce^{4+}}$ liegen nicht vor.

Dieser Code-Block schafft es aber nicht, uns die richtige Lösung zu liefern.
Im Gegenteil erhalten wir mehrere Warnungen wie
```
/path/to/your/python_script.py:line: RuntimeWarning: overflow encountered in xxxxx
```
Sollten Sie das Ergebnis plotten und anschauen, werden Sie feststellen, dass 
die Komponenten in `y` betragsmäßig sehr groß werden, was zum
[arithmetischen Überläuf](https://de.wikipedia.org/wiki/Arithmetischer_Überlauf)
(*engl. overflow*) führt. 
Das ist ein Zeichen dafür, dass die Schrittweite $h$ zu groß ist und das
Euler-Verfahren instabil wird.

Um also die numerische Lösung zu erhalten, müssen wir die Schrittweite $h$
verkleinern. Tatächlich brauchen wir hier eine Schrittweite von `h = 0.0005`,
um eine stabile Lösung zu erhalten:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:solve_ode}}
```
Obwohl wir hier eine korrekte Lösung erhalten, haben wir
$200/0.0005 = 400.000$ Schritte benötigt. Wäre die Anfangsbedingung
so gewählt, dass $c_Y^0$ noch größer ist als $0.001\ \mathrm{M}$, dann
wäre die DGLs noch schwieriger zu lösen und wir müssten noch kleinere
Schrittweiten verwenden.

Nun plotten wir das Ergebnis:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:plot}}
```
Wir haben hier mit den Funktionen `ax.set_xlim` und `ax.set_ylim` die Achsen
eingeschränkt. Außerdem haben wir das Argument `loc='upper right'` an die 
Funktion `ax.legend` gegeben, um die Legende nach oben rechts zu verschieben. 
Das voreingestellte Argument ist `loc='best'`, was die "beste" Position für 
die Legende automatisch berechnet.

Das Diagramm sieht dann so aus:
![Euler-Verfahren für Oregonator-Modell](../assets/figures/02-differential_equations/euler_bz.svg)
Man erkennt hier der periodische, impulsartige Verlauf von $[\mathrm{Z}]$,
was in der Originalreaktion $[\mathrm{Ce^{4+}}]$ entspricht und in dem obigen
Video die Konzentration vom Fe(III) in Ferroin darstellt und dort als die 
Blaufärbung zu sehen ist.

Woher weiß man denn, ob die Schrittweite $h$ klein genug ist?
Eine Faustregel besagt, dass man die Lösung mit halbierter Schrittweite $h/2$
berechnen soll. Bleibt die Lösung gleich wie bei $h$, dann ist $h$ klein genug.

Unabhängig davon, ob $h = 0.0005$ klein genug ist oder nicht, können wir
uns darauf einigen, dass das Euler-Verfahren Schwierigkeiten mit Gl. 
{{eqref: eq:oregonator}} hat. Gibt es Methoden, die trotz größerer
Schrittweiten stabile Lösungen liefern können? 
Da das Euler-Verfahren nur den konstanten und linearen Term der 
Taylor-Entwicklung berücksichtigt (vgl. Gl. {{eqref: eq:ode_taylor_first_order}}),
kommt man unschwer auf die Idee, weitere Ordnungen mitzunehmen. Das führt
zu einer Familie von Methoden, die als *Runge-Kutta-Verfahren* bekannt sind.



[^field1986]: R. J. Field, H.-D. Försterling, *J. Phys. Chem.* **1986**, *90*, 5400&ndash;5407.
[^schneider1996]: F. W. Schneider, A. F. Münster, in *Nichelineare Dynamik in der Chemie*, Spektrum Akademischer Verlag, Heidelberg, **1996**, pp. 67&ndash;72.

---

### Übung

#### Aufgabe 2.1
<!--
{{#include ../psets/02.md:aufgabe_1}}
-->
