## The Euler Method

Let's dive into the world of initial value problems (IVPs). We start with a problem of the form:
$$
\begin{cases}
  y'(x) = f(x, y(x))\\ 
  y(x_0) = y_0
\end{cases}\,,
{{numeq}}{eq:ivp_first_order}
$$
where $f(x, y(x))$ is a given function of the independent variable $x$ 
and the unknown function $y(x)$. Equation {{eqref: eq:ivp_first_order}} represents a
first-order IVP, consisting of a first-order differential equation (as in 
Equation {{eqref: eq:ode_first_order}}) and an initial condition.

```admonish info title="Example"
Let's illustrate Equation {{eqref: eq:ivp_first_order}} with a concrete example:
$$
\begin{cases}
  y'(x) = x * y(x)\\ 
  y(0) = 2
\end{cases}\,,
$$

Try solving this IVP analytically using separation of variables and integration!
```

While analytical solutions are elegant, computers can't easily find particular solutions to IVPs analytically. Therefore, in the following sections, we'll explore how to find the particular solution $y(x)$ numerically.

Following the pattern of many numerical methods, our first step involves discretizing the function $y(x)$. This means selecting a set of points $x_i$ and
examining the function only at these points rather than across its entire domain. The most straightforward choice of points is
a uniform grid, where a starting point $x_0$ and a step size $h$ define subsequent
points as $x_n = x_0 + n h$ for $n = 0, 1, 2, \ldots$. The goal in the following sections will be to calculate, or approximate, the function values of
$y(x)$ at these points $x_1, x_2, \ldots$.

### Theoretical Foundations

The function value of $y$ at the point $x_{n+1} = x_n + h$ can (under certain
conditions on $y(x)$) be approximated using a Taylor expansion around the point
$x_n$:
$$
  y(x_{n+1}) = y(x_n) 
  + h y'(x_n) 
  + \frac{h^2}{2} y''(x_n) 
  + \frac{h^3}{6} y'''(x_n)
  + \cdots\,.
$$

```admonish info title="For the Mathematically Curious" collapsible=true
The condition is that $y(x)$ must be analytic at the point $x_n$, meaning there exists
a power series
$$
  \sum_{k=0}^\infty a_k (x - x_n)^k
$$
that converges for all $\|x - x_n\| < R$, where $R \in \mathbb{R}^+$ 
is the radius of convergence of the power series. Additionally, the radius of convergence
must be larger than the step size $h$, i.e., $h < R$.
```

We can approximate the function $y(x)$ using its first-order Taylor polynomial:
$$
  y(x_{n+1}) = y(x_n) + h y'(x_n) + \mathcal{O}(h^2)
$$
where the error $\mathcal{O}(h^2)$ is proportional to $h^2$. Assuming the error term is negligible, we can write:
$$
  y(x_{n+1}) \approx y(x_n) + h y'(x_n)\,.
  {{numeq}}{eq:ode_taylor_first_order}
$$
Equation {{eqref: eq:ode_taylor_first_order}} shows that knowledge of the
function value $y$ and its derivative $y'$ at point $x_n$ enables
calculation of the function value at the next point $x_{n+1}$.
Using the differential equation from Equation {{eqref: eq:ivp_first_order}}, the
derivative $y'(x_n)$ can be replaced with $f(x_n, y(x_n))$, leading to
$$
  y(x_{n+1}) = y(x_n) + h f(x_n, y(x_n))
$$
To emphasize the discrete nature of the Euler method's output, let's denote $y_n = y(x_n)$ and $y_{n+1} = y(x_{n+1})$, resulting in
$$
  y_{n+1} = y_n + h f(x_n, y_n)
  {{numeq}}{eq:euler_method}
$$
Equation {{eqref: eq:euler_method}} defines the *explicit Euler method* 
for numerically solving a first-order IVP. The initial condition
provides $y(x_0)$, allowing step-by-step calculation of all subsequent
$y(x_{n})$ using Equation {{eqref: eq:euler_method}}.

### Implementation

Let's bring these concepts to life by implementing the Euler method and applying it to some exciting examples from chemistry!

#### Manganese Decay

For our first example, we revisit the manganese decay example. Here we have a 
first-order reaction described by the differential equation
$$
  \dot{c}(t) = -k c(t)
  {{numeq}}{eq:mn_decay}
$$
where $c(t)$ represents the concentration of the manganese complex at time $t$. Note that we often write $\dot{c}(t)$ instead of $c'(t)$ to emphasize that we're taking the derivative with respect to time.

Starting with our required libraries
```python
{{#include ../codes/02-differential_equations/euler_mn.py:imports}}
```
the `dydx` function calculates $f(x_n, y(x_n))$:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:dydx}}
```
This function takes the arguments `x` ($t_n$) and `y` ($c(t_n)$) and
returns the derivative $\dot{c}(t_n) = -k c(t_n)$. Here we've used the fitted
rate constant $k$ from Section 
[1.4](../01-regression/04-nonlinear_regression.md#reaktionskinetik).
Following Equation {{eqref: eq:euler_method}}, the `euler_step` function calculates the function value $y(x_{n+1})$:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:euler_step}}
```

The complete Euler method implementation follows:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:euler_method}}
```

This function takes the initial conditions `x0` and `y0`, the
step size `h`, the derivative function `dydx`, and the number of steps
`n`. We first create the grid `x` using the 
[`np.arange`](https://numpy.org/doc/stable/reference/generated/numpy.arange.html)
function and initialize the zero array `y` to store our solution. Then we 
set the first entry `y[0]` of this array to
the initial value `y0`. Next, we use a for loop
over the number of steps, calling the `euler_step` function in each iteration
to calculate the function value at the next point, 
storing it in `y[i + 1]`. Finally, we return both the grid `x` and the solution `y`.

Let's apply the Euler method to Equation {{eqref: eq:mn_decay}}:
```python
{{#include ../codes/02-differential_equations/euler_mn.py:solve_ode}}
```
We set the initial conditions `C0 = 1.0` and `T0 = 0.0`, the step size
`h = 1.0`, and the maximum time `MAXTIME = 900.0`. The number of steps
`nsteps` is calculated by `int(MAXTIME / h)`. The `int` function
rounds down the result of the division and converts it to an integer. Then we
call the `euler_method` function and store the result in `x` and `y`.

Finally, we can compare our numerical result with the analytical solution.
Note that the analytical solution of Equation 
{{eqref: eq:mn_decay}} is easily calculated as $c(t) = c_0 \exp(-k t)$.
```python
{{#include ../codes/02-differential_equations/euler_mn.py:plot}}
```
This produces the following plot:
![Euler method for Mn decay](../assets/figures/02-differential_equations/euler_mn.svg)

Only by zooming in on the concentration profile can we see the difference between the analytical
and numerical solutions, which we achieved here using the `ax.inset_axes` command
(the exact workings of this command are not important at this point). The Euler method with $h = 1$ thus provides
a very good approximation.

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

**Self-Study Questions**

1. Show that the solution of Equation {{eq:mn_decay}} is given by $c(t) = c_0 \exp(-k t)$.

