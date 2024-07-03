# Neuronale Netzwerke

In diesem letzten Kapitel werden wir uns mit *neuronalen Netzwerken* beschäftigen, einer 
speziellen Klasse von Algorithmen des maschinellen Lernens (ML), die in den letzten Jahren
besonders populär geworden sind. Dieses Teilgebiet des maschinellen Lernens wird auch als 
*Deep Learning* bezeichnet, also als ML mit *tiefen neuronalen Netzen*, wobei die Zusammenhänge 
zwischen Künstlicher Intelligenz (AI), ML und Deep Learning in der Abbildung unten
dargestellt sind. Deep Learning ist zum allergrößten Teil 
verantwortlich für den aktuellen Hype um AI und hat in vielen Anwendungen zu großen 
Fortschritten geführt (z.B. ChatGPT in der Sprachverarbeitung, StableDiffusion in der 
Bildgenerierung und AlphaFold in der Proteinstrukturvorhersage). Die Grundidee von neuronalen
Netzwerken wurde bereits in den 50er Jahren entwickelt, ist aber erst in den letzten Jahren durch
Fortschritte in der Hardware, Software und Datenverfügbarkeit populär geworden. Bevor wir 
uns jedoch den Details von neuronalen Netzwerken zuwenden, wollen eine kurze 
Motivation geben, warum wir uns überhaupt mit diesem Thema beschäftigen sollten.

<figure>
    <center>
    <img src="./assets/figures/06-neural_networks/deep_learning_scheme.pdf"
         alt="AI vs. ML vs. DL"
         width="400"\>
    <figcaption>Einordnung von AI, ML und Deep Learning. Quelle: François Chollet, 
    Deep Learning with Python.</figcaption>
    </center>
</figure>

Wir haben in den vorherigen Abschnitten bereits einige ML-Methoden kennengelernt und 
angewendet. Dazu gehören Methoden des überwachten Lernes, wie z.B. die lineare Regression, 
und des unüberwachten Lernens, wie z.B. PCA. Was all diese Methoden gemeinsam haben, ist, 
dass eine Reihe an Schritten durchgeführt werden müssen, bevor ein Modell
für die Vorhersage genutzt werden kann. Diese umfassen u.a. 

- die Verarbeitung der Daten (*Data Preprocessing*),
- die Aufarbeitung von fehlenden Daten (*Data Annotation*),
- das Extrahieren geeigneter Features (*Feature Extraction*),
- oder sogar die Generierung neuer Features (*Feature Engineering*),
- die Auswahl eines geeigneten Modells (*Model Selection*),
- das Training des Modells (*Model Training*),
- die Evaluation des Modells (*Model Evaluation*),
- die Analyse der Ergebnisse (*Error Analysis*),
- bis hin zur Anwendung des Modells auf neue Daten (*Deployment*).

Viele dieser Schritte haben wir in der Vergangenheit implizit durchgeführt ohne es 
vielleicht zu merken, oder waren nicht notwenig, da die Daten bereits in einem
geeigneten Format vorlagen. In der Praxis ist dies jedoch in den seltensten Fällen der
Fall und viel Arbeit geht in die Vorbereitung der Daten. Die korrekte manuelle Auswahl 
geeigneter Darstellungen ist jedoch nicht immer eindeutig (z.B. bei Bildern, Texten
oder Molekülen) oder geht mit erheblichen Informationsverlust einher (z.B. bei PCA).

Neuronale Netzwerke sind eine Klasse von Algorithmen, die es ermöglichen, viele dieser
Schritte zu automatisieren. Sie sind in der Lage, *selbstständig* geeignete Darstellungen
aus den Daten zu lernen, die für die Vorhersage notwendig sind. Dieser Prozess wird als
*Feature Learning* oder *Representation Learning* bezeichnet und ist ein wesentlicher 
Vorteil gegenüber anderen ML-Methoden. Ein Nachteil ist jedoch, dass diese selbstständig
gelernten Darstellungen nicht immer einfach zu interpretieren sind.

