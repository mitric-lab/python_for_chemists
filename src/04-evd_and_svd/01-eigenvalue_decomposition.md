## Eigenwertzerlegung

Für die Definition der Eigenwertzerlegung benötigen wir zunächst die 
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
Zeilenvektor wird als das hermitesch Konjugierte eines Spaltenvektors notiert:
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

Sei $\bm{A} \in \C{n}{n}$ eine quadratische Matrix. Ein Vektor 
$\vec{v}_i \in \mathbb{C}^n$ heißt Eigenvektor von $\bm{A}$, wenn
$$
  \bm{A} \vec{v}_i = \lambda_i \vec{v}_i
$$
für ein $\lambda \in \mathbb{C}$ gilt. Die Zahl $\lambda_i$ heißt in diesem
Fall Eigenwert von $\bm{A}$ zum Eigenvektor $\vec{v}_i$.

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
kann die Eigenwertzerlegung von $\bm{A}$ als
$$
  \bm{A} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}
$$
geschrieben werden. 

#### Eigenwertzerlegung für normale Matrizen

Obwohl die Eigenwertzerlegung für alle diagonalisierbaren
Matrizen existiert, wollen wir uns hier auf eine Untermenge solcher Matrizen
beschränken, die sogenannten 
[*Normalen Matrizen*](https://de.wikipedia.org/wiki/Normale_Matrix). Eine 
Matrix $\bm{A}\in\C{n}{n}$ heißt normal, wenn sie mit ihrem Transponierten
kommutiert, also wenn $\bm{A}^\dag \bm{A} = \bm{A} \bm{A}^\dag$ gilt.
Solche Matrizen sind unitärdagonalisierbar, d.h. die Matrix $\bm{Q}$ ist
orthogonal. Damit gilt aber
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
  \end{align}
$$

Das Produkt eines Spaltenvektors mit einem Zeilenvektor, wie z.B.
$\vec{v}_i \vec{v}_i^\dag$ in der obigen Gleichung, ist als das
[*dyadische Produkt*](https://de.wikipedia.org/wiki/Dyadisches_Produkt) 
bekannt. Dieses Produkt produziert eine (Rang-1-)Matrix.

Um die Eigenwertzerlegung durchführen zu können, müssen wir aber die 
Eigenwerte und Eigenvektoren der Matrix bestimmen. Während sie für kleine
Matrizen noch analytisch bestimmt werden können, sind numerische Verfahren
für größere Matrizen notwendig. Wir werden uns im Folgenden mit zwei
solchen Verfahren beschäftigen.

#### Eigenwertalgorithmen

Der wohl einfachste Algorithmus zur Bestimmung der Eigenwerte einer Matrix
ist die [*Potenzmethode*](https://de.wikipedia.org/wiki/Potenzmethode).
Dieser Algorithmus findet den Eigenvektor der Matrix $\bm{A}\in\C{n}{n}$ zum
Eigenwert mit dem größten Betrag, wenn $\bm{A}$ diagonalisierbar ist.

Wir nehmen nun an, dass $\bm{A}\in\C{n}{n}$ eine Matrix mit Eigenwerten
$$
  |\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|
$$
ist. Dann liefert die Iteration
$$
  \vec{x}^{k+1} = \frac{\bm{A} \vec{x}^k}{\|\bm{A} \vec{x}^k\|}\quad k=0,1,2,\ldots\,.
  {{numeq}}{eq:power_iteration}
$$
für einen beliebigen Startvektor $\vec{x}^0 \in \mathbb{C}^n$ eine Folge
von Vektoren $(\vec{x}^k)_{k\in\mathbb{N}}$, die gegen den Eigenvektor
zum Eigenwert mit dem größten Betrag, $\lambda_1$, konvergiert.

Das kann man am
```admonish proof title="Beweis der Konvergenz der Potenzmethode" collapsible=true
Eine Matrix ist diagonalisierbar, genau wenn ihre Eigenvektoren eine Basis 
bildet. Deshalb kann man den beliebigen Startvektor $\vec{x}^0$ als 
Linearkombination der Eigenvektoren schreiben:
$$
  \vec{x}^0 = \sum_{i=1}^n a_i \vec{v}_i\,.
$$
mit Entwicklungskoeffizienten $a_i\in\mathbb{C}$ und Eigenvektoren $\vec{v}_i$
zu den Eigenwerten $\lambda_i$. Dann gilt
$$
  \begin{align}
    \vec{x}^k &= \bm{A} \vec{x}^{k-1} = \cdots = \bm{A}^k \vec{x}^0 \\
    &= \sum_{i=1}^n a_i \lambda_i^k \vec{v}_i \\
    &= \lambda_1^k a_1 \vec{v}_1 + \sum_{i=2}^n \left(\frac{\lambda_i}{\lambda_1}\right)^k a_i \vec{v}_i \\
    &= \lambda_1^k \left[ 
      a_1 \vec{v}_1 + \sum_{i=2}^n \left(\frac{\lambda_i}{\lambda_1}\right)^k a_i \vec{v}_i
    \right]\,.
  \end{align}
$$
Da $|\lambda_i| < |\lambda_1|$ für $i>1$ gilt, gilt
$$
  \lim_{k\to\infty} \left(\frac{\lambda_i}{\lambda_1}\right)^k = 0\,,
$$
also geht der zweite Term in der Klammer gegen Null. Das führt zu
$$
  \lim_{k\to\infty} \vec{x}^k = \lim_{k\to\infty} \lambda_1^k a_1 \vec{v}_1 \propto \vec{v}_1\,.
$$
Mit einer endlichen Anzahl an Iterationen $k$ kann also $\vec{x}^k$ als eine
Approximation des Eigenvektors $\vec{v}_1$ angesehen werden, da Eigenvektoren
beliebig skaliert werden können.

Für die numerische Stabilität ist es sinnvoll, in jedem Schritt den Vektor
$\vec{x}^k$ zu normieren, was zu der Iteration
$$
  \vec{x}^{k+1} = \frac{\bm{A} \vec{x}^k}{\|\bm{A} \vec{x}^k\|}
$$
führt.
```
sehen. Den zugehörigen Eigenwert kann dann durch
$$
  \lambda_1 = \frac{x^{k\dag} \bm{A} \vec{x}^k}{x^{k\dag} \vec{x}^k}
  {{numeq}}{eq:power_iteration_eigenvalue}
$$
bestimmt werden.

Eine Erweiterung der Potenzmethode ist die sog. *inverse Iteration*, die
zur Bestimmung eines Eigenwertes nahe einer vorgegeben aber beliebigen Zahl
$\mu$ sowie seinen Eigenvektor verwendet werden kann, Diese Methode werden Sie
in der Übung näher kennenlernen.

Wir wollen uns aber nicht zufrieden geben mit Methoden, die nur ein Eigenpaar
bestimmen kann. Ein häufig verwendetes Verfahren, welches alle Eigenpaare einer
diagonalisierbaren Matrix bestimmen kann, ist der
[*QR-Algorithmus*](https://de.wikipedia.org/wiki/QR-Algorithmus). 

Der Schlüsselschritt des QR-Algorithmus ist die Zerlegung der Matrix $\bm{A}$
in eine unitäre Matrix $\bm{Q}$ und eine obere Dreiecksmatrix $\bm{R}$, also
$\bm{A} = \bm{Q} \bm{R}$. Diese Zerlegung wird als
[*QR-Zerlegung*](https://de.wikipedia.org/wiki/QR-Zerlegung) bezeichnet. Wir 
werden hier nicht auf die Details der QR-Zerlegung eingehen, sondern annehmen,
dass wir sie für jede quadratische Matrix $\bm{A}$ berechnen können.

Beim QR-Algorithmus wird die Ausgangsmatrix $\bm{A}_0 = \bm{A}$ gesetzt und 
ihre QR-Zerlegung $\bm{A}_0 = \bm{Q}_0 \bm{R}_0$ berechnet. Dann wird die
Iteration
$$
  \bm{A}_{k+1} = \bm{R}_k \bm{Q}_k
  {{numeq}}{eq:qr_algorithm}
$$
durchgeführt, also wird die Matrix im nächsten Schritt $\bm{A}_{k+1}$ durch
die Multiplikation der Faktoren der QR-Zerlegung in umgekehrter Reihenfolge
berechnet. Es gilt
$$
  \begin{align}
    \bm{A}_{k+1} &= \bm{R}_k \bm{Q}_k = \textcolor{blue}{\identity} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \textcolor{blue}{\bm{Q}_k^\dag \bm{Q}_k} \textcolor{green}{\bm{R}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \textcolor{darkcyan}{\bm{A}_k} \bm{Q}_k \\
    &= \bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{A}_0} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k \\
    &= \underbrace{\bm{Q}_k^\dag \bm{Q}_{k-1}^\dag \cdots \bm{Q}_0^\dag \textcolor{red}{\bm{Q}}}_{:=\bm{P}_k} 
      \textcolor{red}{\bm{\Lambda}}
      \underbrace{\textcolor{red}{\bm{Q}^{-1}} \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k}_{=\bm{P}_k^{-1}} \\
    &= \bm{P}_k \bm{\Lambda} \bm{P}_k^{-1} \,,
  \end{align}
$$
wobei $\bm{A}_0 = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$ die Eigenwertzerlegung
der Ausgangsmatrix ist. Insgesamt gilt also, dass die Matrix $\bm{A}_k$ 
die gleichen Eigenwerte wie die Ausgangsmatrix $\bm{A}$ hat.

Es kann gezeigt werden, dass die Folge der Matrizen 
$(\bm{A}_k)_{k\in\mathbb{N}_0}$ gegen eine obere Dreiecksmatrix konvergiert,
deren Eigenwerte auf der Diagonalen stehen. Diese Eigenwerte sind dann
die Eigenwerte der Ausgangsmatrix $\bm{A}$. Ist zudem die Ausgangsmatrix
$\bm{A}$ normal, so konvergiert die Folge der Matrizen gegen eine 
Diagonalmatrix, da alle normalen Dreiecksmatrizen diagonal sind.
Da die Matrix der Eigenvektoren $\bm{Q}$ einer Diagonalmatrix die 
Einheitsmatrix $\identity$ ist, also $\bm{P}_k \approx \identity$ gilt, ist
$$
  \bm{Q} \approx \bm{Q}_0 \bm{Q}_1 \cdots \bm{Q}_{k-1} \bm{Q}_k = \prod_{i=0}^k \bm{Q}_i\,, 
  {{numeq}}{eq:qr_algorithm_eigenvectors}
$$
womit wir die Eigenvektoren der normalen Matrix $\bm{A}$ approximieren können.

Der QR-Algorithmus lässt sich auf verschiedene Weisen beschleunigen, die
wir hier nicht weiter diskutieren wollen. 

### Implementierung

Nun implementieren wir die Eigenwertzerlegung mit der Potenzmethode und dem 
QR-Algorithmus implementieren. Als Testfall wollen wir zufällige reelle 
Matrizen wählen. Damit diese Matrix auf jeden Fall diagonalisierbar ist, 
werden wir die zufällig generierte Matrix mit ihrer Transponierten addieren,
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

Als erstes wird die Dimension der Matrix mit `a.shape[0]` bestimmt und in die
Variable `n` gespeichert. Dann haben wir die Funktion
[`np.random.rand`](https://numpy.org/doc/stable/reference/random/generated/numpy.random.rand.html)
verwendet, um einen zufälligen Startvektor zu generieren, welcher 
anschließend durch ihre Norm geteilt und damit normiert wird. Hier wurde die
Funktion
[`np.linalg.norm`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.norm.html)
verwendet, um die (euklidische) Norm eines Vektors zu berechnen.

Eine mögliche Abbruchbedingung ist, wenn die Differenz zwischen den
Vektoren in zwei aufeinanderfolgenden Iterationen kleiner als eine gegebene
Fehlergrenze `eps` ist. Um es zu implmentieren, haben wir die Variable `x_prev`
definiert, die den Vektor in der vorherigen Iteration speichert. 

Danach beginnt die Iteration mit einer `for`-Schleife, die bis zur maximalen
Anzahl an Iterationen `maxiter` läuft. In jeder Iteration wird zuerst der
wert von `x` in `x_prev` gespeichert. Dann wird der neue Vektor `x` durch
`np.dot(a, x)` berechnet und anschließend normiert. Danach wird überprüft, 
ob die Norm der Differenz zwischen `x` und `x_prev` kleiner als `eps` ist.
Ist es der Fall, wird die Schleife abgebrochen.

Nimmt die Laufvariable `i` den Wert `maxiter - 1` an, so hat die Potenzmethode
möglicherweise nicht konvergiert und es wird eine Warnung geprintet.
Anschließend wird der Eigenwert mit Gl. {{eqref: eq:power_iteration_eigenvalue}}
berechnet. Am Ende werden die Ergebnisse zurückgegeben.

```admonish tip title="Tipp zum Programmierstil"
Manchmal kann es vorkommen, dass die gewünschte Variablename bereits ein
reservierter Befehl in Python ist. In diesem Fall kann entweder einen anderen
Namen gewählt werden oder gemäß der Konvention ein Unterstrich an den Namen
angehängt werden. Hier wurde z.B. die Variable `lambda_` verwendet, um den
Eigenwert zu speichern.
```

Anschließend erzeugen wir eine symmetrisierte Zufallsmatrix `a_mat` und rufen
die Funktion `power_iteration` auf mit `eps=1e-8`:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:test_power_iteration}}
```
Um unsere Ergebnisse zu verifizieren, verwenden wir die Funktion
[`np.linalg.eigh`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.eigh.html),
welche die Eigenwerte und Eigenvektoren einer symmetrischen Matrix berechnet.
Dabei werden die Eigenwerte in aufsteigender Reihenfolge und die normierten
Eigenvektoren zurückgegeben. 
Deshalb nehmen wir den letzten Eintrag der Eigenwerte und den zugehörigen
Eigenvektor, um sie mit den Ergebnissen der Potenzmethode zu vergleichen:
```python
{{#include ../codes/04-evd_and_svd/power_iteration.py:verify_results}}
```
Während wir die Eigenwerte direkt mit `np.isclose` vergleichen können, 
müssen wir beim Eigenvektor berücksichtigen, dass die Eigenvektoren nur 
bis auf eine Skalierung eindeutig sind. Weil unsere Eigenvektoren normiert
sind, können wir die Projektion der Eigenvektoren aufeinander
mit `np.dot` berechnen und der Absolutwert der Projektion sollte 
nahe bei 1 sein.

Wiederholtes Ausführen des Codes wird unsere Implementierung der Potenzmethode
an verschiedenen Beispielen testen.

#### QR-Algorithmus
Wir implementieren nun den QR-Algorithmus nach 
Gl. {{eqref: eq:qr_algorithm}} und Gl. {{eqref: eq:qr_algorithm_eigenvectors}}.
Nach dem Importieren von NumPy schreiben wir die Funktion `qr_algorithm`,
die die gleiche Signatur wie die Funktion `power_iteration` hat:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:qr_algorithm_function}}
```
Nach der Bestimmung der Dimension wird die Ausgangsmatrix `a` in `a_k`
kopiert. Es ist hier wichtig, dass eine explizite Kopie von `a` gemacht wird,
da der QR-Algorithmus die Matrix `a` sonst verändern würde. Um die 
zukünftigen $\bm{Q}_k$'s sequentiell multiplizieren zu können, wird die
Variable `eigvecs` als eine Einheitsmatrix initialisiert.

In der `for`-Schleife bis `maxiter` wird die zuerst die QR-Zerlegung von
`a_k` berechnet. Dazu haben wir die Funktion
[`np.linalg.qr`](https://numpy.org/doc/stable/reference/generated/numpy.linalg.qr.html)
verwendet. Danach überschreiben wir `a_k` mit dem Produkt von `r_k` und `q_k` 
und multiplizieren `eigvecs` mit `q_k`.

Da `a_k` gegen eine obere Dreiecksmatrix (bzw. Diagonalmatrix in diesem Fall)
anstrebt, können wir als Abbruchbedingung den größten Betrag der Elemente
im strikten unteren Dreieck der Matrix verwenden. Hierzu wird die Funktion
[`np.tril`](https://numpy.org/doc/stable/reference/generated/numpy.tril.html)
mit dem Argument `k=-1` verwendet, um das strikte untere Dreieck der Matrix
zu erhalten. Ist dieser Betrag kleiner als die Fehlergrenze `eps`, wird die
Schleife abgebrochen. Anschließend wird wie bei der Implementierung der 
Potenzmethode überprüft, ob der Algorithmus konvergiert hat ggf. eine 
Warnung geprintet.

Die Eigenwerte werden aus der Diagonale der Matrix `a_k` mit `np.diag`
extrahiert. Am Ende werden die Indizes der Sortierung der Eigenwerte mit
`np.argsort` berechnet und die Eigenwerte sowie die Eigenvektoren entsprechend
sortiert, bevor sie zurückgegeben werden.

Auch hier testen wir die Implementierung mit einer symmetrischen Zufallsmatrix:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:test_qr_algorithm}}
```
und vergleichen unsere Ergebnisse mit denen von `np.linalg.eigh`:
```python
{{#include ../codes/04-evd_and_svd/qr_algorithm.py:verify_results}}
```
Die paarweise Projektion unserer Eigenvektoren auf die Eigenvektoren von
`np.linalg.eigh` kann in diesem Fall elegant mit der Matrixmultiplikation
`eigvecs.T @ eigvecs_ref` berechnet werden. Bei korrekter Implementierung
sollte dieses Produkt eine Diagonalmatrix sein, deren Diagonalelemente nahe
bei 1 oder -1 liegen. Deshalb nehmen wir den Absolutbetrag der Elemente aus 
der Produktmatrix und vergleichen sie mit der Einheitsmatrix, die mit
[`np.eye`](https://numpy.org/doc/stable/reference/generated/numpy.eye.html)
erzeugt wird.

Eigenwerte und Eigenvektoren enthalten viele essentielle Informationen über
die korrespondierende quadratische Matrix. Können wir auch etwas ähnliches 
finden, die für nicht-quadratische Matrizen gelten? 
Die Antwort ist ein klares Ja, und wir werden uns im nächsten Abschnitt mit
der Singulärwertzerlegung beschäftigen, die als die Verallgemeinerung der
Eigenwertzerlegung für allgemeine Matrizen betrachtet werden kann.

