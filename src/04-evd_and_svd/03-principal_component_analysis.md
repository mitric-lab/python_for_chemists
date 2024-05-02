## Hauptkomponentenanalyse

Die Singulärwertzerlegung einer Matrix $\bm{A}$ in der Form
$$
  \bm{A} = \sum_{i=1}^p \sigma_i \vec{u}_ i \vec{v}_ i^\dag
$$
mit Singulärvektoren $\vec{u}_ i$ und $\vec{v}_ i$ sowie Singulärwerten 
$\sigma_i$ 
(vgl. Gl. {{eqref: eq:singular_value_decomposition}})
suggestiert eine Approximation der Matrix $\bm{A}$ durch Abschneiden der Summe
bis zum $k$-ten Summanden:
$$
  \bm{A_k} = \sum_{i=1}^k \sigma_i \vec{u}_ i \vec{v}_ i^\dag\,.
$$
mit $k \leq p$ und absteigend sortierten Singulärwerten 
$\sigma_1 \geq \sigma_2 \geq \ldots \geq \sigma_p$.
Man spricht hier von einer Rang-$k$-Approximation der Matrix $\bm{A}$,
da $\bm{A_k}$ aus der Summe von $k$ (linear unabhängigen) Rang-1-Matrizen besteht.
Wie gut ist denn diese Approximation? Eine Antwort auf diese Frage liefert das
```admonish note title="Eckart-Young-Theorem"
Sei $\bm{A} \in \C{m}{n}$ eine beliebige Matrix mit der SVD
$\bm{A} = \sum_{i=1}^p \sigma_i \vec{u}_ i \vec{v}_ i^\dag$, wobei 
$p = \min(m,n)$ und $\sigma_1 \geq \sigma_2 \geq \ldots \geq \sigma_p$.
Dann ist die Matrix $\bm{A_k} := \sum_{i=1}^k \sigma_i \vec{u}_ i \vec{v}_ i^\dag$
die **beste** Rang-$k$-Approximation von $\bm{A}$ im Sinne der 
Frobenius-Norm $\|\cdot\|_F$, also
$$
  \|\bm{A} - \bm{A_k}\|_F = \min_{\substack{\bm{B} \in \C{m}{n} \\ \text{rank}(\bm{B}) \leq k}} \|\bm{A} - \bm{B}\|_F\,.
$$

Die Frobenius-Norm einer Matrix $\bm{A} \in \C{m}{n}$ ist definiert als
$$
  \|\bm{A}\|_F = \sqrt{\sum_{i=1}^m \sum_{j=1}^n |a_{ij}|^2}\,.
$$

Der Beweis dieses Theorems erfordert einige Kenntnisse der linearen Algebra
und wird sehr schnell sehr technisch, wenn man alles "zu Fuß" zeigt. Deshalb
wird hier auf den Beweis verzichtet, die Interessierten können z.B.
[hier](https://en.wikipedia.org/wiki/Low-rank_approximation#Proof_of_Eckart–Young–Mirsky_theorem_(for_Frobenius_norm)
einen Beweis finden.
```


```admonish note title="Eckart-Young-Mirsky-Theorem" collapsible=true
Leon Mirsky konnte die obige Approximationseigenschaft auf beliebigt
*unitär invariante Norm* erweitern. 
Eine Norm $\|\cdot\|$ ist *unitär invariant*, wenn für beliebige unitäre
Matrizen $\bm{U}\in \C{m}{m}$ und $\bm{V}\in \C{n}{n}$ die Bedingung
$$
  \|\bm{UAV}\| = \|\bm{A}\|
$$
für alle Matrizen $\bm{A}\in \C{m}{n}$ erfüllt ist.

Einige gebräuchliche unitär invariante Normen seien hier für $\bm{A}\in \C{m}{n}$
mit $p = \min(m,n)$ aufgeführt:
- Die [Frobenius-Norm](https://de.wikipedia.org/wiki/Frobeniusnorm)
  $$
    \|\bm{A}\|_F = \sqrt{\sum_{i=1}^m \sum_{j=1}^n |a_{ij}|^2}
    = \sqrt{\sum_{i=1}^{p} \sigma_i^2}\,,
  $$
- Die [Spektralnorm](https://de.wikipedia.org/wiki/Spektralnorm)
  $$
    \|\bm{A}\|_2 = \max_{\|\vec{x}\|_2 = 1} \|\bm{A}\vec{x}\|_2
    = \sigma_1\,,
  $$
  und
- Die [Ky-Fan-Norm](https://en.wikipedia.org/wiki/Matrix_norm#Ky-Fan_norm)
oder auch Spurnorm
  $$
    \|\bm{A}\|_{*} = \sum_{i=1}^{p} \sigma_i\,.
  $$
Alle drei Normen sind Spezialfälle der 
[Schatten-Normen](https://de.wikipedia.org/wiki/Matrixnorm#Schatten-Normen).
```

Das [Eckart-Young-Theorem](https://en.wikipedia.org/wiki/Low-rank_approximation#Basic_low-rank_approximation_problem),
welches zuerst von Erhard Schmidt bewiesen wurde, besagt also, dass
die SVD nicht nur eine gute Approximation der Matrix $\bm{A}$ liefert,
sondern sogar die **beste** Approximation bis zum Rang $k$ bezüglich der 
Frobenius-Norm.

Die SVD liefert uns aber noch mehr als nur eine Approximation der Matrix 
$\bm{A}$. Da die Approximation sukzessive durch Abschneiden der Singulärwerte
erfolgt, kann die Rang-1-Matrix korrespondierend zum Größten Singulärwert
$\sigma_1 \vec{u}_ 1 \vec{v}_ 1^\dag$ als die "wichtigste" Komponente, 
die *Hauptkomponente*, der Matrix $\bm{A}$ interpretiert werden.
Der Singulärwert $\sigma_1$ gibt dabei an, wie wichtig diese Hauptkomponente
für die Matrix $\bm{A}$ ist. Die weiteren Rang-1-Matrizen 
$\sigma_i \vec{u}_ i \vec{v}_ i^\dag$ mit $i > 1$ können als die
nächsten wichtigen Komponenten interpretiert werden, mit dem jeweiligen
Gewicht $\sigma_i$. Diese Interpretation führt uns zur 
*Hauptkomponentenanalyse* (engl. *Principal Component Analysis*, PCA),
welche in der Praxis häufig zur Dimensionsreduktion von Daten verwendet wird.

### Theoretische Grundlagen

#### Principal Component Analysis (PCA)

Wir betrachten eine Messung von $P$ Merkmalen (engl. *features*) an $N$ 
Proben (engl. *samples*), welche durch eine $N \times P$-Matrix 
$\widetilde{\bm{X}}$ repräsentiert wird. Die $i$-te Zeile der Matrix 
$\widetilde{\bm{X}}$ entspricht dabei den $P$ Merkmalen der $i$-ten Probe.

Der erste Schritt der PCA ist die Vorverarbeitung der Daten auf einer
der zwei möglichen Arten: Zentrierung oder Standardisierung. Hierzu
schreiben wir die Datenmatrix $\widetilde{\bm{X}}$ als
$$
  \widetilde{\bm{X}} =
    \begin{pmatrix}
      \text{---}\ \vec{x}_ 1\ \text{---} \\
      \text{---}\ \vec{x}_ 2\ \text{---} \\
      \vdots \\
      \text{---}\ \vec{x}_ N\ \text{---}
    \end{pmatrix}
$$
mit den Messdaten $\vec{x}_ i \in \mathbb{R}^P$. Jeder Datenpunkt ist also 
ein $P$-dimensionaler Vektor, wobei die Features die Basisvektoren darstellen. 
Nun definieren wir die Mittelwerte Messung $\vec{\mu}$ als
$$
  \vec{\mu} = \frac{1}{N} \sum_{i=1}^N \vec{x}_ i\,.
$$
Bei der **Zentrierung** subtrahieren wir von jedem Datenpunkt den Mittelwert:
$$
  \bm{X} = 
    \begin{pmatrix}
      \text{---}\ (\vec{x}_ 1 - \vec{\mu})\ \text{---} \\
      \text{---}\ (\vec{x}_ 2 - \vec{\mu})\ \text{---} \\
      \vdots \\
      \text{---}\ (\vec{x}_ N - \vec{\mu})\ \text{---}
    \end{pmatrix}\,,
$$ 
um die zentrierte Datenmatrix $\bm{X}$ zu erhalten.

Manchmal kommt es vor, dass die Features in unterschiedlichen 
Größenordnungen auftreten bzw. in unterschiedlichen Einheiten gemessen 
werden. In diesem Fall ist es sinnvoll, neben der Mittelwertsubtraktion
noch eine Normierung der Daten durchzuführen. Die Kombination dieser
zwie Schritte wird als **Standardisierung** bezeichnet.
Nach der Berechnung der Standardabweichung der $j$-ten Feature durch
$$
  \sigma_j = \sigma(\{x_{1j}, x_{2j}, \ldots, x_{Nj}\}) = \sqrt{\frac{1}{N} \sum_{i=1}^N (x_{ij} - \mu_j)^2}\,.
$$
kann die Standardisierte Datenmatrix $\bm{X}$ durch
$$
  \bm{X} = 
    \begin{pmatrix}
      \text{---}\ (\vec{x}_ 1 - \vec{\mu})\ \text{---} \\
      \text{---}\ (\vec{x}_ 2 - \vec{\mu})\ \text{---} \\
      \vdots \\
      \text{---}\ (\vec{x}_ N - \vec{\mu})\ \text{---}
    \end{pmatrix}
    \begin{pmatrix}
      \sigma_1^{-1} &  &  &  \\
       & \sigma_2^{-1} &  &  \\
       &  & \ddots &  \\
       &  &  & \sigma_P^{-1}
    \end{pmatrix}
$$
berechnet werden. In anderen Worten heißt das, dass wir zusätzlich zu
der zeilenweisen Mittelwertsubtraktion noch 
die $j$-te Spalte durch die Standardabweichung $\sigma_j$ teilen müssen.
Je nach Art der Daten soll entschieden werden, ob eine Zentrierung
oder Standardisierung durchgeführt werden soll.

Anschließend berechnen wir die Singulärwertzerlegung der standardisierten
Datenmatrix $\bm{X}$:
$$
  \bm{X} = \bm{U} \bm{\Sigma} \bm{V}^\dag
$$
mit $\bm{U} \in \C{N}{N}$, $\bm{V} \in \C{P}{P}$ und $\bm{\Sigma} \in \R{N}{P}$.
Die Rechts-Singulärvektoren $\vec{v}_ i$ sind die Hauptkomponenten,
die auch als *loadings* bezeichnet werden,
während die Singulärwerte $\sigma_i$'s die Wichtigkeit der jeweiligen
Hauptkomponente angeben. Gebräuchlich ist noch der Varianzanteil
(engl. *explained variance*) $\eta_i$ der $i$-ten Hauptkomponente,
die durch 
$$
  \eta_i = \frac{\sigma_i^2}{\sum_{j} \sigma_j^2}
$$
gegeben ist, wobei die Summe über alle Singulärwerte läuft.

Weil die Rechts-Singulärvektoren $\vec{v}_ i$ orthonormal sind, können
wir diese als eine Orthonormalbasis des $P$-dimensionalen Raums interpretieren,
die möglicherweise eine natürlichere Repräsentation der Datenpunkte 
darstellt. 

Die Projektion der Datenpunkte auf die Hauptkomponenten berechnet sich
durch das Matrixprodukt 
$$
  \bm{X} \bm{V} = \bm{U} \bm{\Sigma} \bm{V}^\dag \bm{V} = \bm{U} \bm{\Sigma}\,,
$$
also ist die Projektion der Datenpunkte auf die $i$-te Hauptkomponente
durch das Produkt des $i$-ten Links-Singulärvektors $\vec{u}_ i$ mit
dem $i$-ten Singulärwert $\sigma_i$ gegeben. Dies Projektion
$\vec{u}_ i \sigma_i$ wird auch als *score* der $i$-ten Hauptkomponente
bezeichnet.

Zwar ist die Transformation der Datenpunkte in die Basis der Hauptkomponenten
allein schon ein nützliches Werkzeug, aber die Hauptkomponentenanalyse
kann auch zur Dimensionsreduktion verwendet werden. Die Idee ist, dass
wir nur die ersten $k$ Hauptkomponenten behalten und die Datenpunkte in den
$k$-dimensionalen Raum der Hauptkomponenten projizieren. Dank des 
Eckart-Young-Theorems wissen wir, dass diese Projektion immer die beste
Beschreibung der ursprünglichen Datenpunkte in einem $k$-dimensionalen Raum
liefert. Also können wir PCA verwenden, um einen hochdimensionalen Datensatz
mit wenigen Hauptkomponenten zu approximieren, ohne dabei zu viel Information
zu verlieren.





Hat der Datensatz sehr viele Features (d.h. $P$ ist groß), so k

#### Principal Coordinate Analysis (PCoA)

WIP



