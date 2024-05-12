## Hauptkoordinatenanalyse

Wie im letzten Kapitel beschrieben, ist die *Hauptkoordinatenanalyse* 
(engl. *principal Coordinate Analysis*, PCoA) ein Verfahren, welches 
die Hauptkomponenten des Datensatzes aus den Abständen zwischen den
Datenpunkten bestimmt.

### Theoretische Grundlagen

Damit wir am Ende die PCA durchführen können, benötigen wir eine 
Koordinatendarstellung der Datenpunkte $\widetilde{\bm{X}}$. Nehmen wir zuerst
an, dass wir eine solche Koordinatendarstellung durch eine magische Kraft
erhalten haben. Dann können wir die an der zentrierten Datenmatrix $\bm{X}$
die SVD durchführen ($\bm{X} = \bm{U} \bm{\Sigma} \bm{V}^\dag$) und erhalten 
die Projektion auf die Hauptkomponenten durch $\bm{U} \bm{\Sigma}$.

Betrachten wir nun die *Gram-Matrix* $\bm{G}$ der zentrierten Datenmatrix 
$\bm{X}$, die durch $\bm{G} = \bm{X} \bm{X}^\dag$ gegeben ist. Setzen wir
die SVD von $\bm{X}$ ein, so erhalten wir
$$
  \bm{G} = \bm{X} \bm{X}^\dag 
  = \bm{U} \bm{\Sigma} \bm{V}^\dag \bm{V} \bm{\Sigma}^\dag \bm{U}^\dag 
  = \bm{U} \underbrace{\bm{\Sigma} \bm{\Sigma}^\dag}_ {:=\bm{\Lambda}} \bm{U}^\dag\,,
$$
ganz nach Gl. {{eqref: eq:svd_and_evd}}. Die Projektion der Datenpunkte 
auf die Hauptkomponenten kann also aus der Eigenwertzerlegung der Gram-Matrix
berechnet werden als $\bm{U} \bm{\Lambda}^{1/2}$.

Setzen wir nun Gl. {{eqref: eq:centre_data_matrix}} in die Definition der
zentrierten Gram-Matrix ein, so erhalten wir
$$
  \begin{align}
    \bm{G} &= \bm{X} \bm{X}^\dag 
    = (\identity_N - \frac{1}{N} \mathbf{1}_ N) \widetilde{\bm{X}} \widetilde{\bm{X}}^\dag (\identity_N - \frac{1}{N} \mathbf{1}_ N) \\
    &= \widetilde{\bm{X}} \widetilde{\bm{X}}^\dag 
      - \frac{\mathbf{1}_ N}{N} \mathbf{1}_ N \widetilde{\bm{X}} \widetilde{\bm{X}}^\dag 
      - \widetilde{\bm{X}} \widetilde{\bm{X}}^\dag \frac{\mathbf{1}_ N}{N}
      + \frac{\mathbf{1}_ N}{N} \widetilde{\bm{X}} \widetilde{\bm{X}}^\dag \frac{\mathbf{1}_ N}{N} \\
    &= \widetilde{\bm{G}} - \frac{\mathbf{1}_ N}{N} \widetilde{\bm{G}} - \widetilde{\bm{G}} \frac{\mathbf{1}_ N}{N} + \frac{\mathbf{1}_ N}{N} \widetilde{\bm{G}} \frac{\mathbf{1}_ N}{N}\,,
  \end{align}
$$
wo wir die Gram-Matrix der unzentrierten Datenmatrix als
$\widetilde{\bm{G}} = \widetilde{\bm{X}} \widetilde{\bm{X}}^\dag$ definiert 
haben. Dieser Prozess wird als *Double Centering* bezeichnet, da wir von der
Matrix $\widetilde{\bm{G}}$ sowohl die Zeilen- als auch die Spaltenmittelwerte
durch $\frac{\mathbf{1}_ N}{N} \widetilde{\bm{G}}$ und 
$\widetilde{\bm{G}} \frac{\mathbf{1}_ N}{N}$ abziehen, und dann den Mittelwert
der gesamten Matrix, der doppelt abgezogen wurde, durch
$\frac{\mathbf{1}_ N}{N} \widetilde{\bm{G}} \frac{\mathbf{1}_ N}{N}$ wieder
hinzufügen. Die Matrix $\bm{G}$ hat also 0 als Spalten- und Zeilenmittelwert.

Das ist zwar schön und gut, dass wir aus den unverarbeiteten Datenkoordinaten
die Hauptkomponenten berechnen können, aber wie müssen erst (ohne Magie)
die Koordinaten erhalten. Wir betrachten nun das was wir haben: die
Abstände zwischen den Datenpunktpaaren. Diese lassen sich in einer
symmetrischen $N \times N$-Matrix $\bm{D}$ mit den Elementen $d_{ij}$ 
speichern, wobei $d_{ij}$ den Abstand zwischen den Datenpunkten $i$ und $j$
angibt. Des Weiteren nehmen wir an, dass die euklidischen Abstände
zwischen den Datenpunkten gegeben sind. Damit gilt
$$
  \begin{align}
  d_{ij}^2 = \|\vec{x}_ i - \vec{x}_ j\|_ 2^2
    &= \|\vec{x}_ i - \vec{\mu}\|_ 2^2 + \|\vec{x}_ j - \vec{\mu}\|_ 2^2 
      - 2 \langle \vec{x}_ i - \vec{\mu}, \vec{x}_ j - \vec{\mu} \rangle \\
    &= \|\vec{x}_ i - \vec{\mu}\|_ 2^2 + \|\vec{x}_ j - \vec{\mu}\|_ 2^2
      -2 G_{ij}\,,
  \end{align}
$$
wobei wir den Kosinussatz verwendet haben. Das Skalarprodukt der zentrierten
Datenpunkte $\vec{x}_ i - \vec{\mu}$ und $\vec{x}_ j - \vec{\mu}$ ist
gerade das Element $G_{ij}$ der zentrierten Gram-Matrix $\bm{G}$.

Wenn wir die Matrix $\bm{D}^{(2)}$ der quadrierten Abstände durch 
$D^{(2)}_{ij} = d_{ij}^2$ definieren, so unterscheidet sich 
$-\frac{1}{2} \bm{D}^{(2)}$ von der zentrierten Gram-Matrix $\bm{G}$ nur
durch einen Zeilen- und einen Spaltenmittelwert. Führt man das Double Centering
auf $\bm{D}^{(2)}$ durch, so erhält man die Matrix $\bm{G}$:
$$
  \bm{G} = -\frac{1}{2}\left(\identity_N - \frac{1}{N} \mathbf{1}_ N\right) \bm{D}^{(2)} \left(\identity_N - \frac{1}{N} \mathbf{1}_ N\right)\,.
$$
Das ist genau die "magische Kraft", die wir benötigen, um die Abstände
in Koordinaten umzuwandeln. Es ergibt sich der folgende Algorithmus:
1. Berechne die Matrix $\bm{D}^{(2)}$ der quadrierten Abstände.
2. Führe das Double Centering auf $\bm{D}^{(2)}$ durch.
3. Berechne die Eigenwertzerlegung von $\bm{G}$ als 
   $\bm{G} = \bm{U} \bm{\Lambda} \bm{U}^\dag$.
4. Berechne projizierten Koordinaten auf die Hauptkomponenten mit
   $\bm{U} \bm{\Lambda}^{1/2}$.
Dieser Algorithmus wird als *Principal Coordinate Analysis* (PCoA) bezeichnet.

Damit ist PCoA äquivalent zur PCA, wenn der Abstand zwischen den Datenpunkten
euklidisch ist. Verwendet man aber eine andere Abstandsmetrik, so liefert 
die PCoA andere Projektionen der Datenpunkte als die PCA. In diesem Fall ist 
die erhaltene Projektion oft eine gute Approximation der opti

Es sei noch angemerkt, dass die PCoA zu einer Familie von Verfahren gehört,
die als
```admonish info title="Multidimensionale Skalierung" collapsible=true
Die
[*Multidimensionale Skalierung*](https://de.wikipedia.org/wiki/Multidimensionale_Skalierung)
(engl. *Multidimensional Scaling*, MDS) versucht,
eine Koordinatendarstellung von Datenpunkten in 
$k$ Dimensionen zu finden, so dass die Abstände zwischen den Datenpunkten
möglichst gut erhalten bleiben. Sei also der Abstand zwischen dem $i$-ten
und $j$-ten Datenpunkt in den ursprünglichen Koordinaten durch $d_{ij}$
und ihre Koordinaten im $k$-dimensionalen Raum durch $\vec{x}_ i$ und
$\vec{x}_ j$ gegeben. Dann wird der *Stress*-Wert
$$
  \text{Stress}(\vec{x}_ 1, \ldots, \vec{x}_ n) = \sqrt{
    \sum_{i\neq j} (d_{ij} - \|\vec{x}_ i - \vec{x}_ j\|)^2
  }
$$
durch die MDS minimiert. Hier wurde die genaue Abstandsmetrik für die 
Berechnung von $d_{ij}$ und die Norm für die Berechnung der Distanz
zwischen den Koordinaten $\vec{x}_ i$ und $\vec{x}_ j$ nicht angegeben,
da die MDS für beliebige Metriken und Normen definiert werden kann. 

Streng genommen ist die PCoA keine MDS, auch wenn die in diesem Kontext als
*Klassische MDS* (CMDS) bezeichnet wird. Es liegt daran, dass die PCoA
den *Strain*-Wert
$$
  \text{Strain}(\vec{x}_ 1, \ldots, \vec{x}_ n) = \sqrt{
     \frac{\sum_{i\neq j} (G_{ij} - \langle \vec{x}_ i, \vec{x}_ j \rangle)^2}
          {\sum_{i\neq j} G_{ij}^2}
  }
$$
minimiert. Diese Funktion ist nicht äuquivalent zum *Stress*-Wert der MDS.
Aber weil die Idee des Strain-Werts sehr ähnlich zum Stress-Wert ist, wird
die PCoA oder die CMDS oft als eine variante der MDS betrachtet.
```
bekannt sind.

### Implementierung

Für die Implementierung der PCoA nehmen wir ein Beispieldatensatz ohne
Koordinaten, die auf dem ersten Blick sinnvoll erscheinen: der 
GDB-9 Datensatz, bestehend aus 133885 kleinen organischen Molekülen 
aus den Atomen H, C, N, O, F bis zu einer Größe von 9 schweren Atomen.
Damit die Anzahl der Datenpunkte nicht zu groß wird, wählen wir nur
eine Untermenge dieses DAtensatztes aus: alle Moleküle aus den
Atomen H, C, N, O mit maximal 5 schweren Atomen. Das liefert uns einen
Datensatz mit 177 Molekülen, der 
<a href="../codes/04-evd_and_svd/gdb9_subset_5.sdf" download>hier</a>
heruntergeladen werden kann.

WIP
