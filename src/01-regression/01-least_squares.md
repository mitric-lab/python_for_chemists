## Methode der kleinsten Quadrate

Das Ziel aller Regressionsprobleme besteht darin, ein Modell zu finden, 
welches am besten zu den Daten passt und hoffentlich zur Vorhersage neuer 
Datenpunkte verwendet werden kann. Einfachheitshalber gehen wir davon aus, 
dass unsere Daten nur eine unabhängige Variable ($x$) und eine abhängige
Variable ($y$) haben. Eine Verallgemeinerung auf mehrere unabhängige
Variablen ist simpel und wird in einem späteren Kapitel diskutiert.
Darüber hinaus gehen wir davon aus, dass die Daten reellwertig sind.
Für insgesamt $N$ Datenpunkte schreiben wir
$$
  \begin{align}
    x &= (x_1, x_2, \ldots, x_N) \in \mathbb{R}^N \\
    y &= (y_1, y_2, \ldots, y_N) \in \mathbb{R}^N \,,
  \end{align}
$$ 
d.h. $x$ und $y$ sind $N$-dimensionale Vektoren mit reellwertigen
Komponenten $x_i$ und $y_i$.

Wir können die Beziehung zwischen $x_i$ und $y_i$ als eine Funktion $f$
betrachten. Unter idealisierten Bedingungen haben wir
$$
  y_i = f(x_i) + \epsilon_i \,, {{numeq}}{eq:regression_true}
$$
wobei die Funktion $f$ die *wahre* Beziehung zwischen den Datenpunkten
$x_i$ und $y_i$ beschreibt, während $\epsilon_i$ den zufälligen
statistischen Fehler darstellt.

In der realen Welt ist die Beziehung $f$ jedoch oft unbekannt. Daher
müssen wir ein Modell verwenden, um sie abzuschätzen. In diesem Fall
schreiben wir
$$
  y_i = \hat{f}(\beta; x_i) + \hat{\epsilon}_ i\,, {{numeq}}{eq:regression_estimated}
$$
wobei $\hat{f}(\beta; x_i)$ der Abschätzer des Modells von $f$ ist,
parametrisiert durch $\beta \in \mathbb{R}^n$ mit $n$ Parametern im
Modell. Der Fehlerterm $\hat{\epsilon}_i$ ist nun eine Kombination aus
dem statistischen Fehler und dem unmodellierten Teil der echten
Beziehung. Das Ziel der Regressionsanalyse besteht darin, die
Parameter $\beta$ zu finden, sodass das Modell die Daten am besten
wiedergibt.

Aber was bedeutet "am besten"? Ein üblicher Ansatz ist die
*Methode der kleinsten Quadrate*, die darauf abzielt, die Summe der
quadratischen Fehler $\hat{\epsilon}_i$ zu minimieren, d.h. die
Verlustfunktion der kleinsten Quadrate
$$
  L(\beta; x, y) 
    = \sum_{i=1}^N\, \hat{\epsilon}_i^2
    = \sum_{i=1}^N\, (y_i - \hat{f}(\beta; x_i))^2 \,. {{numeq}}{eq:least_squares_loss}
$$
zu minimieren. Mathematisch ausgedrückt wollen wir die Parameter
$\beta^* $ durch Lösen des Optimierungsproblems
$$
  \beta^* 
    = \argmin{\beta\in\mathbb{R}^n} L(\beta; x, y)
    = \argmin{\beta\in\mathbb{R}^n} \sum_{i=1}^N\, (y_i - \hat{f}(\beta; x_i))^2 \,. {{numeq}}{eq:least_squares_opt}
$$
finden. Wenn wir den Vektor der vorhergesagten Werte $\hat{y}$ als
$$
  \hat{y}_ i = \hat{f}(\beta; x_i)\,,
$$
schreiben, können wir Gl. {{eqref: eq:least_squares_opt}} als
$$
  \beta^{* } 
    = \argmin{\beta\in\mathbb{R}^n} \sum_{i=1}^N\, (y_ i - \hat{y}_ i)^2 
    = \argmin{\beta\in\mathbb{R}^n} \|y - \hat{y}\|_ 2^2 \,, {{numeq}}{eq:least_squares_opt2}
$$
formulieren, wobei $\| v \|_2$ die euklidische Norm oder die $\ell_2$-Norm
eines Vektors $v$ bezeichnet, definiert als
$$
  \| v \|_2 = \sqrt{\sum_{i=1}^N\, v_i^2} \,.
$$

Die Methode der kleinsten Quadrate ist eine beliebte Wahl für die 
Regressionsanalyse, da sie geschlossene Lösungen für einige einfache, aber 
wichtige Modelle bietet und mit Methoden der numerischen linearen 
Algebra für komplexere Modelle effizient gelöst werden kann. Es weist jedoch 
einige Nachteile auf, z.B. die Empfindlichkeit gegenüber Ausreißern und die 
Anfälligkeit für Überanpassungen. Eine Alternative zur Methode der kleinsten 
Quadrate ist die Methode der kleinsten absoluten Abweichungen, welche die
$\ell_1$-Norm anstelle der $\ell_2$-Norm verwendet. Diese ist für einen
Vektor $v$ definiert als
$$
  \| v \|_1 = \sum_{i=1}^N\, |v_i| \,.
$$

Das Optimierungsproblem wird dann formuliert als
$$
  \beta^* 
    = \argmin{\beta\in\mathbb{R}^n} \|y - \hat{y}\|_1 \,. {{numeq}}{eq:least_absolute_deviations_opt}
$$

Es gibt viele andere Verlustfunktionen, die für die Regression verwendet
werden können, z.B. 
die [Huber-Verlustfunktion](https://en.wikipedia.org/wiki/Huber_loss),
die log-cosh-Verlustfunktion,
die [Quantilverlustfunktion](https://de.wikipedia.org/wiki/Quantilsregression#Optimierungsproblem)
usw. 
Obwohl die Methode der kleinsten Quadrate in vielen Fällen gut funktioniert,
ist es wichtig, sich der Grenzen der Methode bewusst zu sein und
gegebenenfalls andere Verlustfunktionen heranzuziehen.