# Eigenwert- und Singulärwertzerlegung

Eigenwerte und Eigenvektoren einer Quadratischen Matrix sollten Ihnen bekannt 
sein. Hat man alle Eigenwerte und Eigenvektoren einer diagonalisierbaren Matrix
$\bm{A}$ gefunden, kann diese als $\bm{A} = \bm{Q} \bm{\Lambda} \bm{Q}^{-1}$ 
dargestellt werden, wobei $\bm{\Lambda}$ eine Diagonalmatrix der Eigenwerte 
und $\bm{Q}$ die Matrix der Eigenvektoren ist. Diese Darstellung wird als 
*Eigenwertzerlegung* bezeichnet.

Die *Singulärwertzerlegung* kann als eine Verallgemeinerung der 
Eigenwertzerlegung betrachtet werden, da sie für Matrizen beliebiger Dimension
definiert ist. Die Singulärwertzerlegung einer Matrix 
$\bm{A}\in \R{m}{n}$ ist gegeben durch $\bm{A} = \bm{U} \bm{\Sigma} \bm{V}^\intercal$, 
wobei $\bm{\Sigma}\in \R{m}{n}$ eine Diagonalmatrix der Singulärwerte
und $\bm{U}\in \R{m}{m}$ und $\bm{V}\in \R{n}{n}$ orthogonale Matrizen sind.
Die genauen Hintergründe werden wir noch in diesem Kapitel kennenlernen.

Nach der Meinung der Autoren ist die Singulärwertzerlegung **die** wichtigste
Matrixzerlegung in der (numerischen) linearen Algebra und umfasst unzählige
Anwendungen in der Naturwissenschaft, Technik und darüber hinaus. 
Gleichzeitig bildet sie die Basis für viele moderne Algorithmen.
In diesem Kapitel werden wir uns mit der
Eigenwert- und Singulärwertzerlegung beschäftigen und ihre Bedeutung für die
Quantenchemie, die Datenanalyse und das maschinelle Lernen diskutieren.

