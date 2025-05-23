# Klausurvorbereitung

Die folgenden Aufgaben sollen Ihnen dabei helfen, sich auf die Klausur vorzubereiten. 
Die Aufgaben sind so gewählt, dass sie den Prüfungsfragen ähneln und Themen aus 
den Übungen und Vorlesungen abdecken.

## Aufgabe 1: Code-Snippets

<!--- ANCHOR: aufgabe_1 --->
Wählen Sie für die folgenden Code-Snippets die erwartete Ausgabe aus.

| a) | b) | c) |
|----|----|----|
| <pre><code>x = [a for a in 'pkc']<br>y = len(x)<br>print(y ** 2)</code></pre>   | <pre><code>x = 99<br>y = x // -10<br>z = int(y) ** 2<br>print(z)</code></pre> | <pre><code>x = [[i] for i in range(4)]<br>x[-1].append(4)<br>y = np.array(x)<br>print(y)</code></pre> |
| &#9744; <code>3</code>     | &#9744; <code>98.01</code> | &#9744; <code>[[1], [2], [3, 4]]</code> |
| &#9744; <code>6</code>     | &#9744; <code>81</code>    | &#9744; <code>[[1], [2], [3], [4]]</code> |
| &#9744; <code>9</code>     | &#9744; <code>-81</code>   | &#9744; <code>[[0], [1], [2], [3, 4]]</code> |
| &#9744; <code>Error</code> | &#9744; <code>100</code>   | &#9744; <code>Error</code> |
<!--- ANCHOR_END: aufgabe_1 --->

## Aufgabe 2: Euler-Verfahren

<!--- ANCHOR: aufgabe_2 --->
In der folgenden Implementierung des Euler-Verfahrens haben sich **fünf Fehler** 
eingeschlichen. Finden Sie die Fehler und korrigieren Sie diese. Nehmen Sie dabei an, 
dass alle benötigten Module importiert wurden.

```python
{{#include ../codes/07-summary/exam_preparation.py:euler_errors}}
```
<!--- ANCHOR_END: aufgabe_2 --->

## Aufgabe 3: $k$-Nearest Neighbors

<!--- ANCHOR: aufgabe_3 --->
Der $k$-Nearest Neighbors (KNN) Algorithmus ist ein einfacher Algorithmus 
für Klassifikationsprobleme, der aber auch für 
Regressionsprobleme verwendet werden kann. Er ist ein 
parameterfreier Algorithmus, das bedeutet, dass er keine Trainingsphase hat. 

Für einen gegebenen Punkt $\vec{x}_i$ sagt der KNN-Algorithmus die Klasse $y_i$ 
voraus, indem die $k \in \mathbb{N}$ nächsten Nachbarn von $\vec{x}_i$ im 
Trainingsdatensatz basierend auf ihrer euklidischen Distanz zu $\vec{x}_i$ bestimmt werden.
Die Klasse $y_i$ von $\vec{x}_i$ ist dann diejenige Klasse, die unter den $k$ 
nächsten Nachbarn am häufigsten vorkommt. Dies ist im folgenden Bild für zwei Klassen 
veranschaulicht, wobei für $k = 3$ (innerer schwarzer Kreis) die Klasse $y_i$ 
des grünen Punktes als *rot* vorhergesagt wird.

![kNN](../assets/figures/07-summary/kNN.svg)

**(a)**

Welche Werte für $k$ sind mehr oder weniger sinnvoll, wenn die Anzahl der Klassen 
zwei ist? Begründen Sie Ihre Antwort.

**(b)**

Vervollständigen Sie die Methode `predict` der Klasse `kNN_Classifier` in der folgenden 
Implementierung.

```python
{{#include ../codes/07-summary/exam_preparation.py:knn_incomplete}}
```

Das Model soll wie folgt verwendet werden:

```python
{{#include ../codes/07-summary/exam_preparation.py:knn_example}}
```

<!-- 
Lösung:
```
{{#include ../codes/07-summary/exam_preparation.py:knn_complete}}
```
 -->

<!--- ANCHOR_END: aufgabe_3 --->
