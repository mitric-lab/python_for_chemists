# Maschinelles Lernen

**Maschinelles Lernen** (engl. *machine learning*, ML) wird häufig als ein Teilgebiet der **Künstlichen 
Intelligenz** (engl. *artificial intelligence*, AI) betrachtet, welches sich mit der Entwicklung 
von Algorithmen beschäftig, die es uns erlauben, aus Daten zu lernen und Vorhersagen zu treffen. 
Dabei wollen wir die grundlegenden Muster und Strukturen in den Daten erkennen, um diese für
zukünftige Vorhersagen zu nutzen.
Als Künstliche Intelligenz wiederum verstehen wir diejenigen Systeme, die in der Lage sind, Aufgaben 
zu erfüllen, die sonst *menschliche Intelligenz* erfordern würden. Dazu gehören beispielsweise das Erkennen 
von Sprache, das Verstehen von Texten oder das Erkennen von Objekten in Bildern. 

Die Anwendung von ML-Methoden ist in den letzten Jahren immer populärer geworden und hat in unzähligen
Bereichen Einzug gehalten, auch in der chemischen Forschung. Wir möchten jedoch zu Beginn anmerken, 
dass wir hier nur einen sehr groben Überblick über das Thema geben können. Maschinelles Lernen ist ein 
interdisziplinäres Forschungsgebiet, das Methoden aus der Statistik, der Informatik und der Mathematik 
vereint und daher sehr umfangreich ist. Vielmehr möchten wir Ihnen zeigen, wie Sie mit Ihren neu erworbenen
Python-Kenntnissen einfache ML-Modelle implementieren und anwenden können.

Um zu verstehen, was maschinelles Lernen von anderen, Ihnen bereits bekannten Algorithmen unterscheidet,
betrachten wir das folgende Schaubild:

<figure>
    <center>
    <img src="./assets/figures/05-machine_learning/ml_scheme.pdf"
         alt="Classical Algorithms vs. ML"
         width="400"\>
    <figcaption>Klassisches Programmieren vs. ML. Quelle: François Chollet, Deep Learning with Python.</figcaption>
    </center>
</figure>

Der herkömliche Weg einen Computer nützliche Aufgaben erfüllen zu lassen, besteht darin, ihm
eine Reihe von Regeln zu übergeben. Diese Regeln werden von einem Programmierer festgelegt und 
in Form eines Programms umgesetzt. Für einen gegebenen Input, soll der Computer diese Regeln dann
befolgen, um ein bestimmtes Resultat zu erzeugen.

Im Gegensatz dazu, *lernt* ein ML-Algorithmus aus Daten. Anstatt ihm Regeln vorzugeben, wird ihm
eine Menge von Daten gegeben, die aus den Inputs und erwarteten Outputs bestehen. Der Algorithmus lernt dann
aus diesen Daten, wie er den Input in den Output umwandeln kann. Das bedeutet, dass der Algorithmus
selbstständig Regeln aus den Daten extrahiert, anstatt dass diese ihm vorgegeben werden. Diese Regeln 
werden in Form von *Modellen* repräsentiert, die aus den Daten gelernt werden. Wir werden in den folgenden
Abschnitten einen Programmierstil kennenlernen, der es uns erlaubt, solche Modelle zu erstellen und anzuwenden.

Man unterscheidet im Allgemeinen zwischen zwei oder mehr Arten von ML-Verfahren, die sich in der Art und Weise
unterscheiden, wie sie aus den Daten lernen:

- **Überwachtes Lernen** (engl. *supervised learning*): Hierbei werden dem Algorithmus Daten gegeben, die aus
  Inputs und den dazugehörigen Outputs bestehen. Der Algorithmus lernt dann, wie er die Inputs in die Outputs
  umwandeln kann. Ein Beispiel für ein solches Problem ist die Vorhersage der Synthetisierbarkeit eines Moleküls
  basierend auf dessen Struktur und weiteren Eigenschaften.

- **Unüberwachtes Lernen** (engl. *unsupervised learning*): Hierbei werden dem Algorithmus nur die Inputs gegeben,
  ohne dass die dazugehörigen Outputs bekannt sind. Der Algorithmus lernt dann, wie er die Inputs in sinnvolle
  Gruppen einteilen kann. Ein Beispiel für ein solches Problem ist die Identifikation von Gruppen von Molekülen
  mit ähnlichen Eigenschaften.

Manchmal werden auch weitere Arten von Lernverfahren unterschieden, wie z.B. **verstärkendes Lernen** (engl. 
*reinforcement learning*) oder **semi-überwachtes Lernen** (engl. *semi-supervised learning*). Wir werden uns in 
diesem Kapitel jedoch auf die beiden oben genannten Arten beschränken.

```admonish tip title="Tipp"
Versuchen Sie, für die bereits in diesem Kurs behandelten Algorithmen zu überlegen, ob sie eher dem klassischen
Programmieransatz oder dem maschinellen Lernen (überwacht oder unüberwacht) entsprechen. Sie werden feststellen, 
dass einige Ansätze dem maschinellen Lernen zugeordnet werden können, sodass Sie bereits jetzt von sich
behaupten können, ein wenig über ML zu wissen. Herzlichen Glückwunsch!
```

