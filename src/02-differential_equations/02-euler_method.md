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

Erst durch Vergrößerung des Konzentrationsverlaufes können wir den Unterschied zwischen der analytischen
und der numerischen Lösung erkennen, was hier mit Hilfe des Befehls `ax.inset_axes` erreicht wurde
(die genaue Funktionsweise dieses Befehls ist an dieser Stelle unwichtig). Das Euler-Verfahren mit $h = 1$ liefert
demnach eine sehr gute Approximation.

#### Belousov-Zhabotinsky-Reaktion

Nun wollen wir eine Reaktion mit deutlich komplizierterer Kinetik betrachten. Die 
[Belousov-Zhabotinsky-Reaktion](https://de.wikipedia.org/wiki/Belousov-Zhabotinsky-Reaktion)
ist ein klassisches Beispiel einer oszillierenden Reaktion, wobei
ein Redox-System ($\mathrm{Ce^{3+}/Ce^{4+}}$ in der Originalreaktion) abwechselnd
in der oxidierten und der reduzierten Form vorliegt.
In dem folgenden Video können Sie die Reaktion beobachten:

<iframe width="750" height="422" src="https://www.youtube.com/embed/kw9wF-GNjqs?start=141&end=154" title="Everything about the Belousov Zhabotinsky reaction" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
<!--
<div onclick="this.nextElementSibling.style.display='block'; this.style.display='none'">
   <img src="../assets/figures/02-differential_equations/bz_thumbnail.png" style="cursor:pointer" />
</div>
<div style="display:none">
<iframe width="750" height="422" src="https://www.youtube.com/embed/kw9wF-GNjqs?start=141&end=154" title="Everything about the Belousov Zhabotinsky reaction" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowfullscreen></iframe>
</div>
-->

*Hier wird allerdings als Redox-System [Ferroin](https://de.wikipedia.org/wiki/Ferroin) 
verwendet, welches eine stärkere Farbänderung zeigt.*

Der detaillierte Mechanismus dieser Reaktion ist sehr
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
  \mathrm{B} + \mathrm{Z} &\xrightarrow{k_5} \mathrm{Y}\,,
\end{align}
$$
wobei
$\mathrm{A} = \mathrm{BrO_ 3^-}$, 
$\mathrm{B} = \mathrm{Brommalonsäure}$, 
$\mathrm{X} = \mathrm{HBrO_ 2}$, 
$\mathrm{Y} = \mathrm{Br^-}$, 
$\mathrm{Z} = \mathrm{2\ Ce^{4+}}$ 
und $\mathrm{P} = \mathrm{HOBr}$. 
Möglicherweise fällt Ihnen auf, dass die obigen Reaktionsgleichungen nicht ausbalanciert sind. 
Der Grund dafür ist, dass die Konzentrationen einiger Spezies entweder aufgrund
ihrer hohen Konzentration ($\mathrm{H^+}$, Malonsäure, etc.) oder ihrer schnellen
Reaktionen ($\mathrm{BrO_ 2^\cdot}$) als konstant angenommen werden können.
Insbesondere wird angenommen, dass die Konzentrationen $[\mathrm{A}]$
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

Ein wichtiger Unterschied zum Mn-Zerfall ist, dass wir hier ein System von
DGLs vorliegen haben. Glücklicherweise ist Gl. {{eqref: eq:euler_method}} für
DGL-Systeme genauso gültig wie für einzelne DGLs, sofern wir $y(x)$
als eine vektorwertige Funktion ansetzten, d.h.
$y(x) \equalhat ([X](t), [Y](t), [Z](t))^\intercal$.

Nachdem wir die benötigten Libraries importiert haben
```python
{{#include ../codes/02-differential_equations/euler_bz.py:imports}}
```
können wir die Funktion `dydx` für das Oregonator-Modell implementieren:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:dydx}}
```
Beachten Sie, dass der Datentyp des Arguments `y` ein `np.ndarray` ist.
Die Funktionen `euler_step` und `euler_method` sind analog zu den
Funktionen für den Mn-Zerfall implementiert, nur dass wir hier die 
Typ-Deklarationen in den Signaturen anpassen müssen.
```python
{{#include ../codes/02-differential_equations/euler_bz.py:euler_step}}
```
```python
{{#include ../codes/02-differential_equations/euler_bz.py:euler_method}}
```

Ein weiterer Unterschied ist, dass die Variable `y` in der Funktion
`euler_method` als ein 2D-Array initialisiert wird, wobei die erste Dimension
die Anzahl der Komponenten angibt (hier: `ndim = 3`) und die zweite Dimension 
die Anzahl der Schritte enthält. Das Grid `x` bleibt ein 1D-Array, 
da die Zeit für alle Komponenten gleichermaßen gilt. 

Nun lösen wir das DGL-System in Gl. {{eqref: eq:oregonator}} mit dem Euler-Verfahren:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:solve_ode_bad}}
```
Hier haben wir die Anfangsbedingungen `C0 = np.array([0.0, 0.001, 0.0])` 
gewählt, was bedeutet, dass zu Beginn nur die Spezies 
$\mathrm{Br^-}$ mit einer Konzentration von $0.001\ \mathrm{M}$ vorhanden ist. 
Die Spezies $\mathrm{HBrO_ 2}$ und $\mathrm{Ce^{4+}}$ liegen nicht vor.

Dieser Code-Block schafft es jedoch nicht, uns die richtige Lösung zu liefern.
Der Code wird zwar ohne Fehlermeldung ausgeführt, wir erhalten 
jedoch eine Reihe von Warnungen, wie z.B.
```
/path/to/your/python_script.py:line: RuntimeWarning: overflow encountered in xxxxx
```
Sollten Sie das Ergebnis plotten, werden Sie feststellen, dass 
die Komponenten in `y` betragsmäßig sehr groß werden, was zu einem
[arithmetischen Überlauf](https://de.wikipedia.org/wiki/Arithmetischer_Überlauf)
(*engl. overflow*) führt. 
Das kann ein Zeichen dafür sein, dass wir die Schrittweite $h$ zu groß gewählt haben 
und das Euler-Verfahren instabil wird.

In diesem Fall müssen wir die Schrittweite $h$ tatsächlich auf `h = 0.0005` verkleinern,
um eine stabile numerische Lösung zu erhalten:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:solve_ode}}
```
Es sei angemerkt, dass sich dadurch natürlich auch die Anzahl der Schritte erhöht;
hier waren es $200/0.0005 = 400.000$ Schritte. Für Anfangsbedingungen
$c_Y^0 > 0.001$ wären die DGLs sogar noch schwieriger zu lösen und wir müssten eine noch 
kleinere Schrittweite verwenden.

Nun plotten wir das Ergebnis:
```python
{{#include ../codes/02-differential_equations/euler_bz.py:plot}}
```
Wir haben hier mit den Funktionen `ax.set_xlim` und `ax.set_ylim` die Achsen
eingeschränkt. Außerdem haben wir das Argument `loc='upper right'` an die 
Funktion `ax.legend` übergeben, um die Legende nach oben rechts zu verschieben. 
Das voreingestellte Argument ist `loc='best'`, was die "beste" Position für 
die Legende automatisch bestimmt.

Wir erhalten das folgende Diagramm:
![Euler-Verfahren für Oregonator-Modell](../assets/figures/02-differential_equations/euler_bz.svg)
Man erkennt hier den periodischen, impulsartigen Verlauf von $[\mathrm{Z}]$,
was in dem obigen Video der Konzentration von Fe(III) entspricht ($[\mathrm{Ce^{4+}}]$ in der Originalreaktion) 
und dort als Blaufärbung zu sehen ist.

Woher weiß man nun, ob die Schrittweite $h$ klein genug ist?
Eine Faustregel besagt, dass man bei erfolgter Rechnung mit gegebenem $h$ die Rechnung mit halbierter Schrittweite $h/2$
erneut durchführen soll. Bleibt das Ergebnis gleich wie mit $h$, dann ist $h$ klein genug.

Unabhängig davon, ob $h = 0.0005$ klein genug ist oder nicht, konnten wir erkennen, 
dass das Euler-Verfahren Schwierigkeiten hat das Gleichungssystem {{eqref: eq:oregonator}} 
zu lösen. Wir werden uns im folgenden Abschnitt mit Methoden befassen, die trotz größerer
Schrittweiten stabile Lösungen liefern können. Dabei sei daran erinnert, dass das Euler-Verfahren 
nur den konstanten und linearen Term der Taylor-Entwicklung berücksichtigt (vgl. Gl. {{eqref: eq:ode_taylor_first_order}}).
Es liegt daher nahe, auch höhere Ordungen zu berücksichitgen, was zu einer Familie von Methoden führt, die als
*Runge-Kutta-Verfahren* bekannt sind.

[^field1986]: R. J. Field, H.-D. Försterling, *J. Phys. Chem.* **1986**, *90*, 5400&ndash;5407.
[^schneider1996]: F. W. Schneider, A. F. Münster, in *Nichelineare Dynamik in der Chemie*, Spektrum Akademischer Verlag, Heidelberg, **1996**, pp. 67&ndash;72.

---

### Übung

#### Aufgabe 2.1
<!--
{{#include ../psets/02.md:aufgabe_1}}
-->
