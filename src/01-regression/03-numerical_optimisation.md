## Numerische Optimierung

Die numerishe Optimierung bietet uns die Möglichkeit, komplexe Funktionen,
wie die Verlustfunktionen der kleinsten Quadrate
(Gl. {{eqref: eq:least_squares_loss}}), aber auch Funktionen im anderen
Kontext, wie z.B. die Energie eines Moleküls, zu minimieren oder maximieren.
Da die Maximierung eienr Funktion $f$ äquivalent zur Minimierung von $-f$ 
ist, werden wir im Folgenden nur noch vom Minimieren sprechen.
Die Fähigkeit, Optimierungsprobleme anzugehen, 
für die keine geschlossenen Lösungsformeln existieren, oder die Auswertung
des analytischen Ausdrucks zu aufwendig ist, 
erweitert signifikant unser Werkzeugset in der Datenanalyse. Insbesondere 
erlaubt sie es, Lösungen für Modelle zu finden, die durch 
Nichtlinearitäten, hohe Dimensionalitäten oder ungewöhnliche 
Datenverteilungen charakterisiert sind.

### Theoretische Grundlagen

Ein besonders zugänglicher und grundlegender Ansatz für die numerische
Optimierung ist das Verfahren das 
[Gradientenverfahren](https://de.wikipedia.org/wiki/Gradientenverfahren)
(engl. *Gradient Descent*). Das ist ein iteratives Verfahren, welches
von einem geschätzten Startpunkt $x^0$ ausgeht und in jedem Schritt die
Richtung des steilsten Abstiegs der Funktion folgt. Mathematisch
formuliert heißt es
$$
  x^{k+1} = x^k - \alpha \nabla f(x^k)
  {{numeq}}{eq:gradient_descent}
$$
wobei $x^k$ der Schätzwert des Minimums im $k$-ten Schritt ist. Das 
hochgestellte $k$ hat hier nichts mit Potenzierung zu tun, sondern ist 
lediglich eine Notation,
um die Verwechselung mit der $i$-ten Komponente, die hier tiefgestellt wird,
zu vermeiden. Der Gradient der *Objektivfunktion* $f$ bei $x^k$ kann dann als
$\nabla f(x^k) = \left(\frac{\partial f}{\partial x^k_1}, \ldots, \frac{\partial f}{\partial x^k_n}\right)^{\intercal}$ 
notiert werden. 
Die Propotionialitätskonstante $\alpha$ wird als *Schrittweite* oder
*learning rate* bezeichnet. Die Schrittweite ist ein wichtiger Parameter
des Verfahrens, da sie die Konvergenzgeschwindigkeit und die Stabilität
des Verfahrens beeinflusst. Ein zu kleiner Wert für $\alpha$ kann dazu
führen, dass das Verfahren sehr langsam konvergiert, während ein zu
großer Wert dazu führen kann, dass das Verfahren divergiert.

Damit wir das Gradientenverfahren nutzen können, müssen wir Zugang zum 
Gradient der Objektivfunktion haben. Liegt dieser nicht analytisch vor, 
muss eine numerische Approximation verwendet werden. Ein einfacher Ansatz 
liefert die 
