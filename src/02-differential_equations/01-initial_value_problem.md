## Anfangswertproblem

Eine Differentialgleichung ist eine Gleichung, die eine unbekannte Funktion 
und ihre Ableitungen enthält. Eine allgemeine Form solcher Gleichungen ist
kompliziert, da die Funktion und ihre (partiellen) Ableitungen von 
verschiedester Ordnung und zu unterschiedlichen Potenzen auftreten können.
In diesem Abschnitt werden wir uns auf die 
*gewöhnlichen Differentialgleichungen* (GDGLs) beschränken, bei denen die
Funktion nur von einer Variablen abhängt. Außerdem werden wir nur 
*lineare* Differentialgleichungen betrachten, bei denen die Funktion und
ihre Ableitungen nur in der ersten Potenz auftreten. 
Die allgemeine Form einer linearen gewöhnlichen Differentialgleichung
lautet:
$$
  a_n(x) y^{(n)}(x) + a_{n-1}(x) y^{(n-1)}(x) + \ldots + a_1(x) y'(x) + a_0(x) y(x) = b_0(x) u(x)\,,
  {{numeq}}{eq:lode_general}
$$
wobei $x$ die unabhängige Variable, $y(x)$ die gesuchte Funktion und
$a_i(x)$, $b_0(x)$ und $u(x)$ gegebene Funktionen sind. Das größte $n$
mit $a_n(x) \neq 0$ ist der *Grad* der Differentialgleichung. Eine 
lineare Differentialgleichung erster Ordnung hat demnach die Form
$$
  a_1(x) y'(x) + a_0(x) y(x) = b_0(x) u(x)\,.
  {{numeq}}{eq:lode_first_order}
$$

Eine allgemeine DGL erster Ordnung kann in der Form
$$
  y'(x) = f(x, y(x))
  {{numeq}}{eq:ode_first_order}
$$
geschrieben werden, wobei $f(x, y(x))$ eine gegebene Funktion von der
unabhängigen Variablen $x$ und der gesuchten Funktion $y(x)$ ist.
Es ist offensichtlich, dass Gl. {{eqref: eq:lode_first_order}} ein
spezieller Fall von Gl. {{eqref: eq:ode_first_order}} ist, mit
$f(x, y(x)) = -\frac{a_0(x)}{a_1(x)} y(x) + \frac{b_0(x)}{a_1(x)} u(x)$.

Da bei der Bildung der Ableitung der konstante Term verschwindet,
ist die (allgemeine) Lösung einer DGL $y(x)$ nicht eine eindeutige Funktion, 
sondern eine Funktionsschar, die durch Integrationskonstanten parametrisiert 
ist. Um eine eindeutige Lösung zu erhalten, müssen diese Konstanten bestimmt 
werden, wie z.B. durch Anfangsbedingungen.

Ein *Anfangswertproblem* (AWP) ist ein Problem, bei dem zusätzlich zu der
Differentialgleichung eine Anfangsbedingung gegeben ist. Für eine
lineare GDL erster Ordnung ist die Anfangsbedingung durch den Funktionswert 
$y(x_0) = y_0$ gegeben, wobei $x_0$ ein gegebener Punkt und $y_0$ ein
gegebener Funktionswert sind. Obwohl das Anfangsbedingung heißt, muss der
gegebene Punkt nicht notwendigerweise der Anfangspunkt sein. Die Lösung
eines AWP ist dann eine eindeutige Funktion, die auch als *spezielle Lösung* 
der DGL bezeichnet wird.

Für eine DGL $n$-ter Ordnung benötigen wir $n$ Anfangsbedingungen, um
eine spezielle Lösung zu erhalten. Diese sind gegeben durch
$$
  y(x_0) = y_0\,,\quad y'(x_0) = y'_0\,,\quad \ldots\,,\quad y^{(n-1)}(x_0) = y^{(n-1)}_0\,.
$$

Es gibt noch weitere Möglichkeiten, eine spezielle Lösung zu erhalten,
z.B. durch die Angabe von Randbedingungen, was Sie in der Übung kennenlernen
werden.

