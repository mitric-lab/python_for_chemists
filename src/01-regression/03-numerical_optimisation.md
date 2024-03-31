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
(engl. *Gradient Descent*). Das ist ein *iteratives Verfahren*, welches
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
*learning rate* bezeichnet. Das Verfahren wird solange wiederholt, bis
eine oder mehrere Abbruchbedingungen erfüllt sind. Typische 
Abbruchbedingungen für iterative Optimierungsverfahren sind:
- Die Änderung des Funktionswertes ist kleiner als ein Schwellwert
- Die Änderung des Schätzwertes ist kleiner als ein Schwellwert
- Eine maximale Anzahl an Iterationen ist erreicht
- Die Norm des Gradienten ist kleiner als ein Schwellwert.

Die Schrittweite $\alpha$ ist der einzige und zugleich ein wichtiger 
Parameter des Gradientenverfahrens, da sie die Konvergenzgeschwindigkeit 
und die Stabilität des Verfahrens beeinflusst. 
Ein zu kleiner Wert für $\alpha$ kann dazu
führen, dass das Verfahren sehr langsam konvergiert, während ein zu
großer Wert dazu führen kann, dass das Verfahren divergiert.

Damit wir das Gradientenverfahren nutzen können, müssen wir Zugang zum 
Gradienten der Objektivfunktion haben. Liegt dieser nicht analytisch vor, 
muss eine numerische Approximation verwendet werden. Ein einfacher Ansatz 
liefert die 
[Finite Differenz](https://en.wikipedia.org/wiki/Finite_difference)
. Hierbei wird die tangente der partiellen Ableitung durch die
Sekante ersetzt, was zu einer Approximation der Form
$$
  \frac{\partial f}{\partial x^k_i} \approx 
    \frac{f(x^k + h \hat{e}_i) - f(x^k - h \hat{e}_i)}{2h}
  {{numeq}}{eq:finite_difference_symmetric}
$$
führt. Hierbei ist $\hat{e}_i$ der $i$-te Einheitsvektor und $h$ ein
kleiner Wert, der die Schrittweite der Approximation bestimmt. Genauer 
genommen stellt Gl. {{eqref: eq:finite_difference_symmetric}} die 
*zentrale finite Differenz 2. Ordnung* dar. Es gibt auch *einseitige*
Approximationen und Approximationen höherer Ordnung, die wir hier aber 
nicht weiter betrachten.

### Implementierung

#### Finite Differenz

Als erstes implementieren wir die finite Differenz. Da wir mehrmals 
Ableitungen berechnen werden müssen, ist es sinnvoll, eine *Funktion* zu
implementieren. Funktionen im Programmierkontext sind ähnlich zu
mathematischen Funktionen, die eine Eingabe in eine Ausgabe umwandeln, 
können aber auch einiges mehr. 

Wir importieren zuerst `numpy` und die Objekt `Callable` und `Any`
aus dem Modul `typing`:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:import_numpy}}
```
Dann definieren die Funktion `finite_difference`, die den Gradienten einer
Funktion `func` an der Stelle `x0` berechnet. 
```python
{{#include ../codes/01-regression/numerical_optimisation.py:finite_difference}}
```
Der erste Block der Funktion ist die *Signatur*, die angibt, welche
Argumente die Funktion erwartet und welchen Typ die Rückgabe hat. Hierfür
wurde die Objekte `Callable` und `Any` benutzt. Obwohl die explizite Angabe 
der Datentypen in Python optional und immer noch selten ist, ist es nach
Meinung der Autoren eine gute Praxis. Näheres dazu finden Sie in dem 
folgenden
~~~admonish info title="Infokasten" collapsible=true
Eine Analogie zur Signatur in der Mathematik ist die Definitionsmenge und
der Wertebereich einer Funktion, die einschränken, welche Argumente
akzeptiert werden und welche Werte zurückgegeben werden. So funktioniert
auch die Signatur in Computerprogrammen. 

In statisch typisierten (engl. *statically typed*) Programmiersprachen, 
wie C++, Java oder Rust, müssen Funktionen explizite Signaturen haben.
Diese Angabe wird vom Compiler verwendet, um sicherzustellen, dass die
Funktion korrekt verwendet wird und ggf. Optimierungen durchzuführen.

Als Beispiel dient hier eine naive Implementierung der Funktion `powi`
in Rust, die die ganzzahlige Potenz einer Gleitkommazahl berechnet:
```rust,no_run,no_playground
fn powi(x: f64, n: i32) -> f64 {
    let mut result: f64 = 1.0;
    for _ in 0..n {
        result *= x;
    }
    result
}
```
Während die Details dieser Funktion uns nicht interessieren, ist die
erste Zeile des Codeblocks `fn powi(x: f64, n: i32) -> f64` wichtig,
da sie die Signatur der Funktion enthält. Man erkennt, dass die Funktion
`powi` zwei Argumente vom Typ `f64` und `i32` erwartet und einen Wert vom
Typ `f64` zurückgibt.

Eine übliche Implementierung dieser Funktion in Python könnte so aussehen:
```python
def powi(x, n):
    result = 1.0
    for _ in range(n):
        result *= x
    return result
```
Hier gibt es keine explizite Signatur, da Python eine dynamisch typisierte
(engl. *dynamically typed*) Programmiersprache ist. An dieser Funktion erkennt
man auf dem ersten Blick nicht, ob sie für Ganz- oder Gleitkommazahlen
gedacht ist. Um diese Funktion also korrekt zu verwenden, muss man sich mit 
der genauen Implementierung auseinandersetzen, während man durch die Signatur
sofort weiß, welche Argumente erwartet werden.

Deshalb empfehlen die Autoren, obwohl der Python-Interpreter keine
explizite Typisierung erfordert, die Verwendung von Typen in der Signatur
von Funktionen, um die Lesbarkeit und Wartbarkeit des Codes zu verbessern.
~~~

Die Signatur der Funktion `finite_difference` besagt, dass sie eine
Funktion `func` erwartet, die als Argument einen np.ndarray und 
irgendwas (`Any`) akzeptiert und einen Float zurückgibt.
Die weiteren Argumente sind der Punkt `x0` vom Typ np.ndarray, an dem
der Gradient berechnet werden soll, die Schrittweite `h` vom Typ Float
und ein *Tuple* (engl. *tuple*) `args`, das zusätzliche Argumente für die 
Funktion `func` liefern kann. Außerdem haben wir in der Signatur 
voreingestellte (engl. *default*) Werte für `h` und `args` definiert. 
Das bedeutet, dass wenn die Funktion `finite_difference` ohne die 
Argumente `h` und `args` aufgerufen wird, die Werte `h = 1e-5` und `args = ()`
verwendet werden.

In der Funktion wird zuerst die Dimension des Punktes `x0` bestimmt und
in den Integer `n` gespeichert. Der Gradient einer reellwertigen
Funktion mit $n$ Variablen ist ein Vektor der Länge $n$. Daher wird
die Variable `grad` als ein np.ndarray der Länge `n` mit Nullen
durch die Funktion 
[`np.zeros`](https://numpy.org/doc/stable/reference/generated/numpy.zeros.html)
initialisiert. 

Weil wir Gl. {{eqref: eq:finite_difference_symmetric}} an allen $n$ 
Komponenten des Punktes $x^k$ anwenden müssen, ist die Nutzung einer
*Schleife* (engl. *loop*) sinnvoll. Hier benutzen wir eine 
[*for*-Schleife](https://docs.python.org/3/tutorial/controlflow.html#for-statements),
die über alle Indizes `i` von `0` bis `n - 1` iteriert. Beachten Sie, dass
Python-Indizes bei `0` beginnen und die built-in Funktion 
[`range`](https://docs.python.org/3/library/functions.html#func-range)
den Endwert `n` nicht einschließt.

In jeder Iteration der Schleife definieren wir zuerst den Einheitsvektor
$\hat{e}_ i$, indem wir zuerst ein Null-Array der Länge `n` mit 
`np.zeros(n)` erstellen und dann die $i$-te Komponente auf `1` setzen.
Danach können wir den `i`-ten Eintrag des Gradientenarrays `grad` nach
Gl. {{eqref: eq:finite_difference_symmetric}} berechnen. Der einzige
Unterschied zwischen dem Code und der Gleichung ist, dass wir den
Zusatzargument `args` an die Funktion `func` übergeben. Hier 
steht ein Stern `*` vor dem Argument `args`. Das ist die Verwendung
von `*` als einen unären Operator, der ein Objekt in seine Bestandteile
entpackt
(engl. [*unpacking*](https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists)).
Das bedeutet, dass die Funktion `func` nach dem ersten Argument 
die weitere Argumente im Tuple `args` **einzelnd** akzeptiert.

Nach allen Iterationen geben wir den Gradienten gespeichert in `grad` 
zurück. 

#### Objektivfunktion

Als nächstes implementieren wir die Objektivfunktion, deren Wert wir
minimieren wollen. Gemäß Gl. {{eqref: eq:least_squares_loss_linear}} 
können wir die Funktion `objective_function` so implementieren:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:objective_function}}
```

Die Signatur der Funktion `objective_function` folgt die Typdefinition
vom Argument `func` in der Funktion `finite_difference` 
(`Callable[[np.ndarray, Any], float]`):
- Das Argument `beta` ist vom Typ np.ndarray,
- das Argument `args` ist vom Typ np.ndarray, also auch `Any` und
- die Funktion gibt einen Float zurück.

Hier sehen wir wieder die Verwendung von `*` vor dem Argument `args`. 
In diesem Fall heißt das, dass die Funktion `objective_function` beliebig
viele weitere Argumente nach `beta` akzeptiert.

Die Funktion `objective_function` definiert zuerst die Arrays
`concenctrations` und `absorbances` aus dem Argument `args`
und gibt anschließend den Wert der Verlustfunktion der kleinsten
Quadrate nach Gl. {{eqref: eq:least_squares_loss_linear}} zurück.
Hier gibt es wieder keine große Unterschiede zwischen der mathematischen
Formulierung und der programmatischen Implementierung.

Mit Hilfe der Funktionen `finite_difference` können wir die Funktion
`objective_function_gradient` implementieren, die den Gradienten der
Objektivfunktion berechnet:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:objective_function_gradient}}
```
Die Signatur dieser Funktion ist ähnlich zu der der Funktion 
`objective_function`, nur der Rückgabetyp ist ein np.ndarray, da die 
Objektivfunktion einen $n$-dimensionalen Vektor als Argument hat.
Nach der Definition der Arrays `concenctrations` und `absorbances`
berechnen wir den Gradienten der Objektivfunktion einfach durch
Aufrufen der Funktion `finite_difference`. Der Rückgabewert wird in der 
Variable `grad` gespeichert und zurückgegeben.

#### Gradientenverfahren

Anschließend implementieren wir das Gradientenverfahren. 
Gl. {{eqref: eq:gradient_descent}} sagt uns, dass wir als Argumente
den Gradienten der Objektivfunktion `func_grad`, den Startpunkt `x0`
und die Schrittweite `alpha` benötigen. Als Abbruchbedingung verwenden 
wir eine Kombination aus der maximalen Anzahl an Iterationen und der Norm 
des Gradienten. Das führt zu den zusätzlichen Argumenten `max_iter` und 
`max_norm`. Außerdem brauchen wir noch das Argument `args`, das die
weiteren Argumente für die Funktion `func_grad` enthält. Das führt zu
der folgenden Implementierung:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:gradient_descent}}
```
Als Rückgabewert haben wir hier nicht nur das Optimum `x`, sondern auch
die Anzahl der Iterationen `niter`. 

Zuerst setzen wir die Variable `x` auf den Startpunkt `x0`. Danach verwenden 
wir eine `for`-Schleife, die die Variable `niter` von `0` bis `max_iter - 1`
iteriert. In jeder Iteration berechnen wir den Gradienten `grad` durch 
Aufrufen der Funktion `func_grad` mit den Argumenten `x` und `args`.
Anschließend verwenden wir Gl. {{eqref: eq:gradient_descent}} um die Variable
`x` zu aktualisieren. Danach berechnen wir die Norm des Gradienten 
mit der Funktion 
[`np.linalg.norm`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html)
und prüfen, ob sie kleiner als `max_norm` ist. Wenn ja, brechen wir die
Schleife mit dem `break`-Befehl ab.

Nach allen Iterationen überprüfen wir, ob die Variable `niter` die maximale
Anzahl an Iterationen erreicht hat. Wenn ja, geben wir eine Warnung aus, da
die Güte des Ergebnisses fragwürdig ist. Am Ende geben wir das Optimum `x`
so wie die Anzahl der Iterationen `niter` zurück.

### Anwendung

Nun können wir den implementierten Algorithmus auf die Daten aus Kapitel
[1.2](02-linear_regression.md) anwenden. Dazu definieren wir als erstes 
wieder die Arrays `concentrations` und `absorbances`:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:data_array}}
```
Dann definieren wir den Startpunkt `x0` und rufen die Funktion
`gradient_descent` auf:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:gradient_descent_call}}
```

Die optimale Parameter `beta0` und `beta1` sind identisch wie die analytische
Lösung. Auf dem Computer des Autors wurden 2507 Iterationen benötigt, um
die Abbruchbedingung zu erfüllen. Die genaue Anzahl der Iterationen kann
je nach Hardware leicht variieren. Da dieses Optimierungsproblem sehr
"einfach" ist, kann die Schrittweite `alpha` größer gewählt werden, ohne
dass das Verfahren divergiert.

```admonish tip title="Tipp"
Spielen Sie gerne mit den Argumenten `alpha` herum und beobachten Sie, wie
sich die Anzahl der Iterationen und die Güte des Ergebnisses verändern.
```

Da Optimierung ein sehr allgemeines Problem ist, existieren viele
Implementierungen von verschiedensten Algorithmen in Bibliotheken wie
[`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html).
Wir wollen die Funktion `sci.optimize.minimize` verwenden, um die optimalen
Parameter $\beta$ zu finden:
```python
{{#include ../codes/01-regression/numerical_optimisation.py:scipy_minimize}}
```
Beim Aufruden der Funktion `minimize` müssen wir nur die Objektivfunktion
`objective_function` und den Startpunkt `beta_guess` angeben. Um die 
Berechnung der numerischen Gradienten kümmert sich die `minimize`-Funktion
selbst.

~~~admonish note title="Anmerkung zur Funktion `minimize`"
Die Funktion `minimize` akzeptiert auch den Argument `jac` 
([Jacobi-Matrix](https://de.wikipedia.org/wiki/Jacobi-Matrix)).
Das ist eine Funktion,
die den Gradienten der Objektivfunktion berechnet. Sollte man den analytischen
Gradienten zur Verfügung haben, der nicht sehr zeitaufwendig zu berechnen
ist, sollte man ihn als Argument `jac` übergeben. Das kann den 
Optimierungsprozess beschleunigen.
~~~


Mit dem Argument `method='CG'` haben wir die 
*Methode des nichtlinearen Konjugierten Gradienten*
(engl. [<i>nonlinear **C**onjugate **G**radient method</i>](https://en.wikipedia.org/wiki/Nonlinear_conjugate_gradient_method))
ausgewählt. Diese Methode ist ein Verbesseung des Gradientenverfahrens und 
erreicht das Minimum in nur 3 Iterationen auf dem Computer des Autors.

Bei `minimize` gibt es eine Reihe von weiteren Minimierungsmethoden, die
verwendet werden können. Eine Übersicht finden Sie in der 
[Dokumentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) dieser Funktion. Zwei wichtige Methoden davon sind:
- `method='Nelder-Mead'`: 
    Das [Nelder-Mead-Verfahren](https://de.wikipedia.org/wiki/Downhill-Simplex-Verfahren)
    ist eine heutristische Methode, die ohne die Berechnung des Gradienten
    auskommt. Deshalb ist sie besonders nützlich, wenn die Berechnung des
    Gradienten sehr aufwendig ist oder der Gradient durch Rauschen stark
    fluktuiert. Deshalb ist sie besonders geeignet für Regressionen mit
    experimentellen Daten.
- `method='BFGS'`: 
    Das [Broyden-Fletcher-Goldfarb-Shanno-Verfahren](https://de.wikipedia.org/wiki/BFGS-Verfahren)
    ist eine Methode, die den Gradienten benutzt, um die Hesse-Matrix der
    Objektivfunktion zu approximieren. Deshalb zeigt sie eine sehr schnelle
    Konvergenz in der Nähe des Minimums. In der Praxis benögtigt sie weniger 
    Iterationen als andere Optimierungsmethoden und ist deshalb eine der
    am häufigsten verwendeten Methoden.

Im Kontext der Regression erlaubt uns die numerische Optimierung 
einerseits den Einsatz komplizierterer Modelle, die keine analytische
Lösung haben, und andererseits die Verwendung sophistizierterer
Verlustfunktionen, wie z.B. die der Methode der kleinsten absoluten 
Abweichungen (vgl. Gl. {{eqref: eq:least_absolute_deviations_opt}}).
oder die der kleinsten Quadrate mit Regularisierung, was Sie in der Übung
kennenlernen werden.

~~~admonish note title="Funktion `scipy.optimize.curve_fit`"
Es gibt auch die Funktion 
[`scipy.optimize.curve_fit`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html),
die eine (nichtlineare) Regression durchführt. 
Allerdings ist sie nicht so flexibel wie die Funktion
`minimize`, da sie nur über wenige Optimierungsmethoden verfügt und die
Objektivfunktion als die Verlustfunktion der kleinsten Quadrate festlegt.
Deshalb ist die Funktion für einfachere Regressionen geeignet und benötigt
weniger Code. Bei komplizierteren Modellen ist es jedoch sinnvoller, die
Funktion `minimize` zu verwenden.
~~~

---

### Übung

#### Aufgabe 1.2
<!--
{{#include ../psets/01.md:aufgabe_2}}
-->

