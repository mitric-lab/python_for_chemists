## Anfangswertproblem

Eine Differentialgleichung ist eine Gleichung, die eine unbekannte Funktion 
und ihre Ableitungen enthält. Eine allgemeine Form solcher Gleichungen ist
kompliziert, da die Funktion und ihre (partiellen) Ableitungen von 
verschiedester Ordnung und zu unterschiedlichen Potenzen auftreten können.
In diesem Abschnitt werden wir uns auf die 
*gewöhnlichen Differentialgleichungen* (GDGLs, engl. *ODEs*) beschränken, bei denen die
Funktion nur von einer Variablen abhängt. Außerdem werden wir nur 
*lineare* Differentialgleichungen betrachten, bei denen die Funktion und
ihre Ableitungen nur in der ersten Potenz auftreten (d.h. ersten *Grades* sind). 
Die allgemeine Form einer linearen gewöhnlichen Differentialgleichung
lautet:
$$
  a_n(x) y^{(n)}(x) + a_{n-1}(x) y^{(n-1)}(x) + \ldots + a_1(x) y'(x) + a_0(x) y(x) = b_0(x) u(x)\,,
  {{numeq}}{eq:lode_general}
$$
wobei $x$ die unabhängige Variable, $y(x)$ die gesuchte Funktion, $y^{(i)}(x)$ 
die $i$-te Ableitung dieser Funktion ist und
$a_i(x)$, $b_0(x)$ und $u(x)$ gegebene Funktionen von $x$ sind. Das größte $i$,
für das $a_i(x) \neq 0$, ist die *Ordnung* der Differentialgleichung. Eine 
lineare Differentialgleichung erster Ordnung hat demnach die Form
$$
  a_1(x) y'(x) + a_0(x) y(x) = b_0(x) u(x)\,.
  {{numeq}}{eq:lode_first_order}
$$

Eine allgemeine DGL erster Ordnung kann in der expliziten Form
$$
  y'(x) = f(x, y(x))
  {{numeq}}{eq:ode_first_order}
$$
geschrieben werden, wobei $f(x, y(x))$ eine gegebene Funktion von der
unabhängigen Variablen $x$ und der gesuchten Funktion $y(x)$ ist.
Es ist offensichtlich, dass Gl. {{eqref: eq:ode_first_order}} ein
spezieller Fall von Gl. {{eqref: eq:lode_first_order}} ist, mit
$f(x, y(x)) = -\frac{a_0(x)}{a_1(x)} y(x) + \frac{b_0(x)}{a_1(x)} u(x)$.

Da bei der Bildung der Ableitung der konstante Term verschwindet,
ist die (allgemeine) Lösung $y(x)$ einer DGL nicht eine eindeutige Funktion, 
sondern eine Funktionsschar, die durch Integrationskonstanten parametrisiert 
ist. Um eine eindeutige Lösung zu erhalten, müssen diese Konstanten bestimmt 
werden, wie z.B. durch die Wahl von Anfangsbedingungen.

Ein *Anfangswertproblem* (AWP) ist dann ein Problem, bei dem zusätzlich zu der
Differentialgleichung eine oder mehrere Anfangsbedingungen gegeben sind. Für eine
lineare GDL erster Ordnung ist die Anfangsbedingung durch den Funktionswert 
$y(x_0) = y_0$ gegeben, wobei $x_0$ ein gegebener Punkt und $y_0$ ein
gegebener Funktionswert sind. Obwohl wie dies als Anfangsbedingung bezeichnen, muss der
gegebene Punkt nicht notwendigerweise der Anfangspunkt, z.B. Zeitpunkt, sein. Die Lösung
eines AWP ist dann eine eindeutige Funktion, die auch als *spezielle Lösung* 
der DGL bezeichnet wird.

Für eine DGL $n$-ter Ordnung benötigen wir $n$ Anfangsbedingungen, um
eine spezielle Lösung zu erhalten. Diese sind gegeben durch
$$
  y(x_0) = y_0\,,\quad y'(x_0) = y_1\,,\quad \ldots\,,\quad y^{(n-1)}(x_0) = y_{n-1}\,.
$$

Es gibt noch weitere Möglichkeiten, eine spezielle Lösung zu erhalten,
z.B. durch die Angabe von Randbedingungen, was Sie in der Übung kennenlernen
werden.

