## Unüberwachtes Lernen

Bisher haben wir stets Daten betrachtet, die aus Datenpaaren $(\vec{x}_i, y_i)$ für 
$i = 1, \ldots, N$ bestanden, wobei $\vec{x}_i \in \mathbb{R}^n$ ein Featurevektor 
und $y_i \in \mathbb{R}$ ein quantitatives (Regression) oder qualitatives (Klassifikation)
Label oder Target war. Beim so genannten *unüberwachten Lernen* hingegen fehlt diese 
zweite Information und wir haben lediglich Datenpunkte $\{(\vec{x}_i)\}_{i=1, \dots, N}$
zur Verfügung. Dabei interessieren uns insbesondere die *Dimensionsreduktion* oder das 
*Clustering* von Datenpunkten.

Unter Dimensionsreduktion versehen wir die Transformation der Daten in einen Raum
geringerer Dimension mit möglichst geringem Informationsverlust. Dazu könnten wir 
uns z.B. am Beispiel des Wine Quality Datensatzes fragen, ob wir zur Beschreibung der 
Daten nicht nur die 11 Features, sondern auch weniger verwenden können. Zudem sind 
manche der Features möglicherweise stark korreliert und somit überflüssig. Mit der 
Hauptkomponentenanalyse (PCA) haben Sie in Kapitel 
[(4.3)](../04-evd_and_svd/03-principal_component_analysis.md) bereits die 
wichtigste Methode zur Dimensionsreduktion kennengelernt, welche die Daten 
durch die Hauptkomponenten der größten Varianz beschreibt. Mit der 
Hauptkoordinatenanalyse (PCoA) aus Kapitel 
[(4.4)](../04-evd_and_svd/04-principal_coordinate_analysis.md) kennen Sie zudem 
eine Methode zur Dimensionsreduktion, die auf Ähnlichkeiten zwischen den Datenpunkten
basiert. Wir werden daher direkt zum Clustering übergehen.

```admonish note title="Hinweis"
Häufig werden Methoden des überwachten und unüberwachten Lernens kombiniert oder 
nacheinander angewendet. So kann z.B. ein hochdimensionaler Datensatz zunächst
mit PCA auf wenige Dimensionen reduziert werden, bevor ein Klassifikator darauf
trainiert wird oder Clustering durchgeführt wird.
```

Unter Clustering verstehen wir die Gruppierung von Datenpunkten in Cluster, wobei
die Datenpunkte innerhalb eines Clusters möglichst ähnlich und zwischen den Clustern
möglichst verschieden sein sollen. Dazu betrachten wir im Folgenden den so genannten
*k*-Means Algorithmus, der auch als *Lloyd's Algorithmus* bekannt ist.

### $k$-Means Clustering

Gegeben seien wie zuvor $N$ Datenpunkte $\mathcal{D} := \{\vec{x}_i\}_{i=1, \dots, N}$, 
welche wir in eine vorgegebene Anzahl $K$ von Clustern 
gruppieren möchten. Wir beginnen zunächst mit einer Wunschliste der Eigenschaften, die 
wir von einem Clustering-Algorithmus erwarten:

1. Eine allgemeine **Zuweisungsregel**, die jedem Datenpunkt einen Cluster zuordnet, d.h.
   $\vec{x}_i \mapsto k \in \{1, \ldots, K\}$ für $i = 1, \ldots, N$.
2. Eine **Rekonstruktionsregel**, die für jedes Cluster $k \in \{1, \ldots, K\}$ ein 
   repräsentatives Element $\vec{m}_k$ bestimmt, d.h. $k \mapsto \vec{m}_k \in \mathbb{R}^n$
   für $k \in \{1, \ldots, K\}$.

Dabei bezeichnen wir $\vec{m}_k$ auch als Mittelwert des Clusters $k$. Um den $k$-Means
Algorithmus zu formulieren, fixieren wir zunächst die Anzahl der Cluster $K$ und definieren 
zwei Größen, $\mathbf{C} := (C_1, \ldots, C_K)$ und 
$\mathbf{M} := (\vec{m}_1, \ldots, \vec{m}_K)$. Die Cluster-variable $\mathbf{C}$ enthält 
die Teilmengen $C_k \subseteq \mathcal{D}$ der Datenpunkte, die dem Cluster $k$ zugeordnet
sind, während $\mathbf{M}$ die Mittelwerte der Cluster enthält. Die Vereinigung der Cluster
muss dabei die gesamte Datenmenge $\mathcal{D}$ ergeben, d.h. $\bigcup_{k=1}^K C_k = \mathcal{D}$
und $C_i \cap C_j = \emptyset$ für $i \neq j$, d.h. ein Datenpunkt kann nicht gleichzeitig
mehreren Clustern zugeordnet sein. Der $k$-Means Algorithmus ist ein iteratives Verfahren,
welches die Cluster-Variable $\mathbf{C}$ und die Mittelwerte $\mathbf{M}$ abwechselnd 
akualisiert. Für ein initiales Clustering $\mathbf{C}$ wird dabei zunächst der Mittelwert jedes
Clusters als der Mittelwert der Datenpunkte in diesem Cluster berechnet:

$$
\vec{m}_k \leftarrow \frac{1}{|C_k|} \sum_{\vec{x}_i \in C_k} \vec{x}_i\,,
$$

wobei $|C_k|$ die Anzahl der Datenpunkte im Cluster $k$ bezeichnet. Dies entspricht der 
Rekonstruktionsregel. Anschließend werden 
die berechneten Mittelwerte $\mathbf{m}$ festgehalten und die $k$-te Gruppe als diejenige
Menge von Datenpunkten definiert, die dem Mittelwert $\vec{m}_k$ näher ist als jedem anderen
Mittelwert $\vec{m}_j$ für $j \neq k$. Formal ausgedrückt bedeutet das:

$$
C_k \leftarrow \{\vec{x}_i \in \mathcal{D} \mid \|\vec{x}_i - \vec{m}_k\| \leq \|\vec{x}_i - \vec{m}_j\| \text{ für alle } j \neq k\}\,,
$$

was der Zuweisungsregel entspricht. Diese beiden
Schritte werden dann iterativ für eine vorgegebene Anzahl von Iterationen wiederholt, 
oder bis sich die Cluster nicht mehr ändern. Für Daten in $\mathbb{R}^2$ ist in der 
folgenden Abbildung ein Beispiel für den $k$-Means Algorithmus dargestellt:

<figure>
    <center>
    <img src="../assets/figures/05-machine_learning/k_means.pdf"
         alt="k_means"
         width="600"\>
    <figcaption>Illustration des $k$-Means Algorithmus für $K = 2$ Cluster. Quelle: Christopher M. Bishop, Pattern Recognition and Machine Learning.</figcaption>
    </center>
</figure>

Eine andere Blickweise auf die durch den $k$-Means Algorithmus zugeteilten Cluster
stellen übrigens die so genannten *Voronoi-Zellen* dar, die als 

$$
V_k := \{\vec{x} \in \mathbb{R}^n \mid \|\vec{x} - \vec{m}_k\| \leq \|\vec{x} - \vec{m}_j\| \text{ für alle } j \neq k\}
$$

definiert sind und in Voronoi-Diagrammen dargestellt werden können:

<figure>
    <center>
    <img src="../assets/figures/05-machine_learning/voronoi.pdf"
         alt="voronoi"
         width="300"\>
    <figcaption>Voronoi-Zellen, wobei die Punkte die Mittelwerte der Cluster darstellen und die Linien die Grenzen der Zellen. Quelle: Kevin P. Murphy, Machine Learning - A Probabilistic Perspective.</figcaption>
    </center>
</figure>
