## Eigenwertzerlegung

Für die Definition der Eigenwertzerlegung (engl. *eigenvalue decomposition*, EVD)
benötigen wir zunächst die 
Definition eines Vektorraums. Obwohl dieser viel allgemeiner definiert wird,
beschränken wir uns hier auf $\mathbb{C}^n$, also den $n$-dimensionalen
Vektorraum über die komplexen Zahlen. 
Die Vektoren in $\mathbb{C}^n$ können in der Standardbasis
als Spaltenvektoren dargestellt werden. Die Notation $\vec{v} \in \mathbb{C}^n$
bezeichnet einen solchen Vektor:
$$
  \vec{v} = \begin{pmatrix} v_1 \\ v_2 \\ \vdots \\ v_n \end{pmatrix}
$$
mit den Komponenten $v_i \in \mathbb{C}$. Ein komplex konjugierter 
Zeilenvektor wird als der hermitesch Konjugierte eines Spaltenvektors notiert:
$$
  \vec{v}^\dag = \begin{pmatrix} v_1^* & v_2^* & \cdots & v_n^* \end{pmatrix}
$$
mit den komplex konjugierten Komponenten $v_i^*$.
Das Standard-Skalarprodukt zweier Vektoren $\vec{v}, \vec{w} \in \mathbb{C}^n$ 
kann in dieser Notation als
$$
  \langle \vec{v}, \vec{w} \rangle = \vec{v}^\dag \vec{w} = \sum_{i=1}^n v_i^* w_i
$$
geschrieben werden, wobei die Regel für Matrixmultiplikation gilt. 

Sei $\bm{A} \in \C{n}{n}$ nun eine quadratische Matrix. Ein Vektor 
$\vec{v}_i \in \mathbb{C}^n$ heißt Eigenvektor von $\bm{A}$, wenn
$$
  \bm{A} \vec{v}_i = \lambda_i \vec{v}_i
$$
für ein $\lambda \in \mathbb{C}$ gilt. Die Zahl $\lambda_i$ wird als Eigenwert 
von $\bm{A}$ zum Eigenvektor $\vec{v}_i$ bezeichnet.

### Theoretische Grundlagen

Definieren wir nun die Diagonalmatrix der Eigenwerte $\bm{\Lambda}$ als
$$
  \bm{\Lambda} = \begin{pmatrix}
    \lambda_1 & 0 & \cdots & 0 \\
    0 & \lambda_2 & \cdots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    0 & 0 & \cdots & \lambda_n
  \end{pmatrix}
$$
und die Matrix der Eigenvektoren $\bm{Q}$ als
$$
  \bm{Q} = \begin{pmatrix}
    \vert & \vert &  & \vert \\
    \vec{v}_1 & \vec{v}_2 & \cdots & \vec{v}_n \\
    \vert & \vert &  & \vert
  \end{pmatrix}\,,
$$
dann kann die Eigenwertzerlegung von $\bm{A}$ als
$$
  \bm{A} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}
$$
geschrieben werden. 

#### Eigenwertzerlegung für normale Matrizen

Obwohl die Eigenwertzerlegung für alle diagonalisierbaren
Matrizen existiert, wollen wir uns hier auf eine Untermenge solcher Matrizen
beschränken, die sogenannten 
[*Normalen Matrizen*](https://de.wikipedia.org/wiki/Normale_Matrix). Eine 
Matrix $\bm{A}\in\C{n}{n}$ heißt normal, wenn sie mit ihrer Transponierten
kommutiert, also wenn $\bm{A}^\dag \bm{A} = \bm{A} \bm{A}^\dag$ gilt.
Solche Matrizen sind unitärdagonalisierbar, d.h. die Matrix $\bm{Q}$ ist
orthogonal. Damit gilt
$$
  \bm{Q}^{-1} = \bm{Q}^\dag = \begin{pmatrix}
    \text{---}\ \vec{v}_1^\dag\ \text{---} \\
    \text{---}\ \vec{v}_2^\dag\ \text{---} \\
    \vdots \\
    \text{---}\ \vec{v}_n^\dag\ \text{---}
  \end{pmatrix}\,.
$$

Für eine normale Matrix $\bm{A}$ kann die Eigenwertzerlegung also als
$$
  \begin{align}
    \bm{A} &= \bm{Q} \bm{\Lambda} \bm{Q}^\dag \\
    &= \begin{pmatrix}
      \vert & \vert &  & \vert \\
      \vec{v}_1 & \vec{v}_2 & \cdots & \vec{v}_n \\
      \vert & \vert &  & \vert
    \end{pmatrix}
    \begin{pmatrix}
      \lambda_1 & 0 & \cdots & 0 \\
      0 & \lambda_2 & \cdots & 0 \\
      \vdots & \vdots & \ddots & \vdots \\
      0 & 0 & \cdots & \lambda_n
    \end{pmatrix}
    \begin{pmatrix}
      \text{---}\ \vec{v}_1^\dag\ \text{---} \\
      \text{---}\ \vec{v}_2^\dag\ \text{---} \\
      \vdots \\
      \text{---}\ \vec{v}_n^\dag\ \text{---}
    \end{pmatrix} \\
    &= \lambda_1 \vec{v}_1 \vec{v}_1^\dag + \lambda_2 \vec{v}_2 \vec{v}_2^\dag + \cdots + \lambda_n \vec{v}_n \vec{v}_n^\dag \\
    &= \sum_{i=1}^n \lambda_i \vec{v}_i \vec{v}_i^\dag 
    {{numeq}}{eq:eigenvalue_decomposition_normal}
  \end{align}
$$

geschrieben werden. Das Produkt eines Spaltenvektors mit einem Zeilenvektor, wie z.B.
$\vec{v}_i \vec{v}_i^\dag$ in der obigen Gleichung, ist dabei als das
[*dyadische Produkt*](https://de.wikipedia.org/wiki/Dyadisches_Produkt) 
bekannt und produziert eine (Rang-1-)Matrix.

Um die Eigenwertzerlegung durchführen zu können, müssen wir demnach zunächst die 
Eigenwerte und Eigenvektoren der Matrix bestimmen. Während sie für kleine
Matrizen noch analytisch bestimmt werden können, sind numerische Verfahren
für größere Matrizen notwendig. Wir werden uns im Folgenden mit zwei
solcher Verfahren beschäftigen.

#### Eigenwertalgorithmen

Der wohl einfachste Algorithmus zur Bestimmung der Eigenwerte einer Matrix
ist die [*Potenzmethode*](https://de.wikipedia.org/wiki/Potenzmethode).
Dieser Algorithmus findet den Eigenvektor der Matrix $\bm{A}\in\C{n}{n}$ zum
Eigenwert mit dem größten Betrag, sofern $\bm{A}$ diagonalisierbar ist.

Wir nehmen dazu an, dass $\bm{A}\in\C{n}{n}$ eine Matrix mit (unbekannten) Eigenwerten
$$
  |\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|
$$
ist. Dann liefert die Iteration
$$
  \vec{x}^{k+1} = \frac{\bm{A} \vec{x}^k}{\|\bm{A} \vec{x}^k\|}\quad k=0,1,2,\ldots
  {{numeq}}{eq:power_iteration}
$$
für einen beliebigen Startvektor $\vec{x}^0 \in \mathbb{C}^n$ eine Folge
von Vektoren $(\vec{x}^k)_{k\in\mathbb{N}}$, die gegen den Eigenvektor
zum Eigenwert mit dem größten Betrag, $\lambda_1$, konvergiert.

Den zugehörigen Eigenwert kann dann durch
$$
  \lambda_1 = \frac{x^{k\dag} \bm{A} \vec{x}^k}{x^{k\dag} \vec{x}^k}
  {{numeq}}{eq:power_iteration_eigenvalue}
$$
bestimmt werden.

Für matheamtisch Interessierte ist hier ein Beweis der Konvergenz der Potenzmethode 
zu finden.

```admonish proof title="Beweis der Konvergenz der Potenzmethode" collapsible=true
Eine Matrix ist diagonalisierbar, genau denn wenn ihre Eigenvektoren eine Basis 
bilden. Deshalb kann man den beliebigen Startvektor $\vec{x}^0$ als 
Linearkombination der Eigenvektoren schreiben:
$$
  \vec{x}^0 = \sum_{i=1}^n a_i \vec{v}_i
$$
mit Entwicklungskoeffizienten $a_i\in\mathbb{C}$ und Eigenvektoren $\vec{v}_i$
zu den Eigenwerten $\lambda_i$. Dann gilt
$$
  \begin{align}
    \vec{x}^k &= \bm{A} \vec{x}^{k-1} = \bm{A}^2 \vec{x}^{k-2} = \cdots = \bm{A}^k \vec{x}^0 \\
    &= \sum_{i=1}^n a_i \bm{A}^k \vec{v}_i \\
    &= \sum_{i=1}^n a_i \lambda_i^k \vec{v}_i \\
    &= \lambda_1^k a_1 \vec{v}_1 + \sum_{i=2}^n \lambda_i^k a_i \vec{v}_i \\
    &= \lambda_1^k \left[ 
      a_1 \vec{v}_1 + \sum_{i=2}^n \left(\frac{\lambda_i}{\lambda_1}\right)^k a_i \vec{v}_i
    \right]\,.
  \end{align}
$$
Da $|\lambda_i| < |\lambda_1|$ für $i>1$, geht der zweite Term in der Klammer für 
$k\to\infty$ gegen Null, sodass
$$
  \lim_{k\to\infty} \vec{x}^k = \lim_{k\to\infty} \lambda_1^k a_1 \vec{v}_1 \propto \vec{v}_1\,.
$$
Mit einer endlichen Anzahl an Iterationen $k$ kann $\vec{x}^k$ also als eine
Approximation des Eigenvektors $\vec{v}_1$ angesehen werden, zumal Eigenvektoren
beliebig skaliert werden können.

Für die numerische Stabilität ist es sinnvoll, den Vektor $\vec{x}^k$ in jedem Schritt 
zu normieren, was zu der Iteration
$$
  \vec{x}^{k+1} = \frac{\bm{A} \vec{x}^k}{\|\bm{A} \vec{x}^k\|}
$$
führt.
```

Eine Erweiterung der Potenzmethode ist die sog. *inverse Iteration*, die
zur Bestimmung eines Eigenwertes nahe einer vorgegeben aber beliebigen Zahl
$\mu$ sowie dessen Eigenvektor verwendet werden kann.

In den meisten Fällen sind wir allerdings an allen Eigenvektoren und Eigenwerten 
einer Matrix interessiert. Ein häufig verwendetes Verfahren, welches alle Eigenpaare 
einer diagonalisierbaren Matrix bestimmen kann, ist der 
[*QR-Algorithmus*](https://de.wikipedia.org/wiki/QR-Algorithmus). 

Der Schlüsselschritt des QR-Algorithmus ist die Zerlegung der Matrix $\bm{A}$
in eine unitäre Matrix $\bm{Q}$ und eine obere Dreiecksmatrix $\bm{R}$, also
$\bm{A} = \bm{Q} \bm{R}$. Diese Zerlegung wird als
[*QR-Zerlegung*](https://de.wikipedia.org/wiki/QR-Zerlegung) bezeichnet. Wir 
werden hier nicht auf die Details der QR-Zerlegung eingehen, sondern annehmen,
dass wir sie für jede quadratische Matrix $\bm{A}$ berechnen können.

Beim QR-Algorithmus wird die Ausgangsmatrix als $\bm{A}_0 = \bm{A}$ gesetzt und 
ihre QR-Zerlegung $\bm{A}_0 = \bm{Q}_0 \bm{R}_0$ berechnet. Dann wird die
Iteration
$$
  \bm{A}_{k+1} = \bm{R}_k \bm{Q}_k
  {{numeq}}{eq:qr_algorithm}
$$
durchgeführt; die Matrix wird im nächsten Schritt $\bm{A}_{k+1}$ also durch
die Multiplikation der Faktoren der QR-Zerlegung in umgekehrter Reihenfolge
berechnet. Es gilt
$$
  \begin{align}
    \bm{A}_{k+1} &= \bm{R}_k \bm{Q}_k = \textcolor{blue}{\identity} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \textcolor{blue}{\bm{Q}_k^\dag \bm{Q}_k} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{darkcyan}{\bm{A}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\identity} \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\bm{Q}_{k-1}^\dag \bm{Q}_{k-1}} \textcolor{darkcyan}{\bm{R}_{k-1} \bm{Q}_{k-1}} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{orange}{\bm{Q}_{k-1}^\dag} \textcolor{lightgreen}{\bm{A}_{k-1}} \textcolor{darkcyan}{\bm{Q}_{k-1}} \bm{Q}_k \\
    &  \vdots \\
    &= \bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{A}_0} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k \\
    &= \underbrace{\bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{Q}}}_{:=\bm{P}_k} 
      \textcolor{red}{\bm{\Lambda}}
      \underbrace{\textcolor{red}{\bm{Q}^{-1}} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k}_{=\bm{P}_k^{-1}} \\
    &= \bm{P}_k \bm{\Lambda} \bm{P}_k^{-1} \,,
  \end{align}
$$
wobei $\bm{A}_0 = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$ die exakte aber unbekannte 
Eigenwertzerlegung der Ausgangsmatrix ist. Insgesamt gilt also, dass die Matrix $\bm{A}_k$ 
die gleichen Eigenwerte wie die Ausgangsmatrix $\bm{A}$ hat.

Es kann gezeigt werden, dass die Folge der Matrizen 
$(\bm{A}_k)_{k\in\mathbb{N}_0}$ gegen eine obere Dreiecksmatrix konvergiert,
deren Eigenwerte auf der Diagonalen stehen. Diese Eigenwerte sind dann
gleichzeitig die Eigenwerte der Ausgangsmatrix $\bm{A}$. Ist die Ausgangsmatrix
$\bm{A}$ zudem normal, so konvergiert die Folge der Matrizen gegen eine 
Diagonalmatrix, da alle normalen Dreiecksmatrizen diagonal sind.
Da die Matrix der Eigenvektoren $\bm{Q}$ einer Diagonalmatrix die 
Einheitsmatrix $\identity$ ist, also $\bm{P}_k \approx \identity$, folgt aus
der obigen Definition
$$
  \bm{P}_k = \bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \bm{Q} \approx \identity\,,
$$
dass 
$$
  \bm{Q} \approx \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k = \prod_{i=0}^k \bm{Q}_i\,, 
  {{numeq}}{eq:qr_algorithm_eigenvectors}
$$
womit wir die Eigenvektoren der normalen Matrix $\bm{A}$ approximieren können.

Der QR-Algorithmus lässt sich auf verschiedene Weisen beschleunigen, die
wir hier nicht weiter diskutieren wollen. 

### Implementierung

Nun wollen wir die Eigenwertzerlegung mit der Potenzmethode und dem 
QR-Algorithmus implementieren. Zum Testen wollen wir zunächst zufällige reelle 
Matrizen wählen. Damit diese Matrix auf jeden Fall diagonalisierbar ist, 
werden wir sie mit ihrer Transponierten addieren,
um eine symmetrische Matrix zu erhalten. 

#### Potenzmethode
Nach dem Importieren von NumPy
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:import}}
```
implementieren wir die Potenzmethode nach Gl. {{eqref: eq:power_iteration}}
als eine Funktion:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:power_iteration_function}}
```
Diese Funktion akzeptiert die Argumente `a` (Matrix $\bm{A}$), 
`eps` (Fehlergrenze für die Konvergenz) und `maxiter` (maximale Anzahl an
Iterationen). Sie gibt den betragsmäßig größten Eigenwert und den zugehörigen
Eigenvektor zurück.

Als erstes wird die Dimension der Matrix mit `a.shape[0]` bestimmt und in der
Variable `n` gespeichert. Dann haben wir die Funktion
[`np.random.rand`](https://numpy.org/doc/stable/reference/random/generated/numpy.random.rand.html)
verwendet, um einen zufälligen Startvektor zu generieren, welcher 
anschließend durch deine Norm geteilt und damit normiert wird. Dazu wurde die
Funktion
[`np.linalg.norm`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html)
verwendet, um die (euklidische) Norm eines Vektors zu berechnen.

Eine mögliche Abbruchbedingung stellt dar, wenn die Norm der Differenz zwischen den
Vektoren in zwei aufeinanderfolgenden Iterationen kleiner als eine gegebene
Fehlergrenze `eps` ist. Um dies zu implmentieren, haben wir die Variable `x_prev`
definiert, die den Vektor in der vorherigen Iteration speichert. 

Danach beginnt die Iteration mit einer `for`-Schleife, die bis zur maximalen
Anzahl an Iterationen `maxiter` läuft. In jeder Iteration wird zuerst der
wert von `x` in `x_prev` gespeichert. Dann wird der neue Vektor `x` durch
`np.dot(a, x)` berechnet und anschließend normiert. Danach wird überprüft, 
ob die Norm der Differenz zwischen `x` und `x_prev` kleiner als `eps` ist.
Ist dies der Fall, wird die Schleife abgebrochen.

Nimmt die Laufvariable `i` den Wert `maxiter - 1` an, so ist die Potenzmethode
möglicherweise nicht konvergiert und es wird eine Warnung ausgeschrieben.
Anschließend wird der Eigenwert mit Gl. {{eqref: eq:power_iteration_eigenvalue}}
berechnet. Am Ende werden die Ergebnisse zurückgegeben.

```admonish tip title="Tipp zum Programmierstil"
Manchmal kommt es vor, dass der gewünschte Variablenname bereits ein
reservierter Befehl in Python ist. In diesem Fall kann entweder ein anderer
Name gewählt werden oder, gemäß der Konvention, ein Unterstrich an den Namen
angehängt werden. Hier wurde z.B. die Variable `lambda_` verwendet, um den
Eigenwert zu speichern.
```

Im Folgenden erzeugen wir eine symmetrisierte Zufallsmatrix `a_mat` und rufen
die Funktion `power_iteration` mit `eps=1e-8` auf:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:test_power_iteration}}
```
Um unsere Ergebnisse zu verifizieren, verwenden wir die Funktion
[`np.linalg.eigh`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html),
welche die Eigenwerte und Eigenvektoren einer symmetrischen Matrix berechnet.
Dabei werden die Eigenwerte in aufsteigender Reihenfolge in dem 1D-Array `eigvals` und die 
normierten Eigenvektoren als Spalten in dem 2D-Array `eigvecs` zurückgegeben. 
Aus diesem Grund nehmen wir den letzten Eintrag der Eigenwerte und den zugehörigen
Eigenvektor, um sie mit den Ergebnissen der Potenzmethode zu vergleichen:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:verify_results}}
```
Während wir die Eigenwerte direkt mit `np.isclose` vergleichen können, 
müssen wir beim Eigenvektor berücksichtigen, dass dieser nur 
bis auf eine Skalierung eindeutig sind. Da unsere Eigenvektoren normiert
sind, können wir die Projektion der Eigenvektoren aufeinander
mit `np.dot` berechnen und der Absolutwert der Projektion sollte 
nahe bei 1 sein.

Ein wiederholtes Ausführen unseres Codes würde unsere Implementierung der Potenzmethode
an verschiedenen Beispielen testen.

#### QR-Algorithmus
Wir implementieren nun den QR-Algorithmus nach 
Gl. {{eqref: eq:qr_algorithm}} und Gl. {{eqref: eq:qr_algorithm_eigenvectors}}.
Nach dem Importieren von NumPy definieren wir die Funktion `qr_algorithm`,
welche die gleiche Signatur wie die Funktion `power_iteration` hat:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:qr_algorithm_function}}
```
Nach der Bestimmung der Dimension wird die Ausgangsmatrix `a` in `a_k`
kopiert. Es ist hier wichtig, dass eine explizite Kopie von `a` gemacht wird,
da der QR-Algorithmus die Matrix `a` ansonsten verändern würde. Um die 
zukünftigen $\bm{Q}_k$ sequentiell multiplizieren zu können, wird die
Variable `eigvecs` als eine Einheitsmatrix initialisiert.

Innerhalb der `for`-Schleife, die bis `maxiter` läuft, wird die zuerst die QR-Zerlegung von
`a_k` berechnet. Dazu haben wir die Funktion
[`np.linalg.qr`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.qr.html)
verwendet. Danach überschreiben wir `a_k` mit dem Produkt von `r_k` und `q_k` 
und multiplizieren `eigvecs` mit `q_k`.

Da `a_k` gegen eine obere Dreiecksmatrix (bzw. in diesem Fall Diagonalmatrix)
anstrebt, können wir als Abbruchbedingung den größten Betrag der Elemente
im strikten unteren Dreieck der Matrix verwenden. Hierzu wird die Funktion
[`np.tril`](https://numpy.org/doc/stable/reference/generated/numpy.tril.html)
mit dem Argument `k=-1` verwendet, um das untere Dreieck der Matrix
zu erhalten. Ist dieser Betrag kleiner als die Fehlergrenze `eps`, wird die
Schleife abgebrochen. Anschließend wird, wie bei der Implementierung der 
Potenzmethode, überprüft, ob der Algorithmus konvergiert ist und ggf. eine 
Warnung ausgegeben.

Die Eigenwerte werden aus der Diagonale der Matrix `a_k` mit `np.diag`
extrahiert. Am Ende wird die Reihenfolge der Indizes der Eigenwerte mit
`np.argsort` berechnet und die Eigenwerte, sowie die Eigenvektoren entsprechend
sortiert, bevor sie zurückgegeben werden.

Auch hier testen wir die Implementierung mit einer symmetrischen Zufallsmatrix
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:test_qr_algorithm}}
```
und vergleichen unsere Ergebnisse mit denen von `np.linalg.eigh`:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:verify_results}}
```
Die paarweise Projektion unserer Eigenvektoren auf alle Eigenvektoren von
`np.linalg.eigh` kann in diesem Fall elegant mit der Matrixmultiplikation
`eigvecs.T @ eigvecs_ref` durchgeführt werden. Im Falle eines erfolgreich
konvergierten QR-Algorithmus
sollte dieses Produkt einer Diagonalmatrix entsprechen, deren Diagonalelemente nahe
bei 1 oder -1 liegen. Deshalb nehmen wir den Absolutbetrag der Elemente aus 
der Produktmatrix und vergleichen sie mit der Einheitsmatrix, die wir mit
[`np.eye`](https://numpy.org/doc/stable/reference/generated/numpy.eye.html)
erzeugen.

Die Eigenwertzerlegung liefert uns viele essentielle Informationen über
die korrespondierende quadratische Matrix. Können wir etwas ähnliches 
auch für nicht-quadratische Matrizen finden? 
Wie Sie im nächsten Abschnitt sehen werden, ist die Antwort ist ein klares Ja. 
Wir werden mit der Singulärwertzerlegung eine Methode kennenlernen, die als die Verallgemeinerung der
Eigenwertzerlegung für allgemeine Matrizen betrachtet werden kann.

---

### Übung

#### Aufgabe 1: Hückel-Theorie

{{#include ../psets/04.md:aufgabe_1}}

