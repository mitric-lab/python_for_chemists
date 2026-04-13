---
numbering:
  headings: false
---

(sec:what_is_programming)=

# Was ist Programmieren?

:::{caution} Under Construction
Diese Seite befindet sich noch im Aufbau. Bitte haben Sie etwas Geduld, bis sie vollständig ist.
:::

Programmieren bedeutet, einem Computer Anweisungen so aufzuschreiben, dass er sie Schritt für Schritt ausführen kann. Diese Anweisungen werden in einer Textdatei gespeichert und bilden zusammen ein Programm.

Im Gegensatz zu Menschen, deren Anweisungen eine gewisse Mehrdeutigkeit zulassen, verarbeitet ein Computer Informationen nach genau festgelegten Regeln. Er "versteht" nicht, was wir meinen, sondern führt nur genau die Befehle aus, die ihm gegeben werden.

Die zentralen Komponenten eines Computers bestehen aus Eingabe, Verarbeitung, Speicher und Ausgabe und sind in Abbildung \ref{fig:computer_components} dargestellt. Die Central Processing Unit (CPU) ist für die Verarbeitung verantwortlich. Sie besteht aus einer Steuereinheit (Control Unit), die den Ablauf eines Programms steuert, und einer arithmetisch-logischen Einheit (ALU), die die eigentlichen Rechenoperationen ausführt. Alle dafür benötigten Daten und Programme liegen im Speicher (Memory). Von dort werden sie von der Steuereinheit abgerufen und nach der Verarbeitung wieder abgelegt. Die Ergebnisse werden schließlich über die Ausgabe bereitgestellt, etwa auf dem Bildschirm oder als Datei.

:::{figure} ../../assets/figures/preface/computer_figma.pdf
:align: center
:::

Doch wie wird ein solches Programm tatsächlich ausgeführt? Abbildung \ref{fig:program_execution} zeigt zwei typische Wege. Bei kompilierten Sprachen (z. B. C++) wird der Quellcode vor der Ausführung über mehrere Zwischenschritte (Parser, Optimierung und Codegenerierung) vollständig in Maschinencode übersetzt. Diesen Prozess nennt man Kompilierung. Das Ergebnis ist eine eigenständige, ausführbare Datei, die anschließend direkt von der CPU geladen und ausgeführt werden kann. Python unterscheidet sich davon grundlegend: Der Quellcode wird nicht vorab in Maschinencode übersetzt, sondern zunächst in eine Zwischendarstellung (Bytecode) überführt. Der Bytecode ist jedoch nicht direkt ausführbar, sondern wird von einem Interpreter während der Ausführung schrittweise verarbeitet und in konkrete Operationen auf der CPU umgesetzt.

:::{figure} ../../assets/figures/preface/code_figma.pdf
:align: center
:::

Obwohl Python durch diese Art der Ausführung tendenziell langsamer ist als kompilierte Sprachen, ist ein wesentlicher Vorteil die große Anzahl verfügbarer Bibliotheken. Viele davon sind intern in effizienten Sprachen wie C oder C++ implementiert, wodurch sich auch rechenintensive Aufgaben mit Python performant umsetzen lassen. Allerdings sind diese Bibliotheken oft versionsabhängig. Unterschiedliche Versionen von Python oder einzelnen Paketen können zu unterschiedlichem Verhalten führen oder inkompatibel sein. Aus diesem Grund verwendet man sogenannte Python-Umgebungen: Sie verknüpfen eine bestimmte Python-Version mit genau definierten Paketversionen und erlauben es, mehrere solcher Konfigurationen parallel zu verwalten. Eine Python-Datei kann daher in unterschiedlichen Umgebungen ausgeführt werden. Je nach installierten Paketen und Versionen kann das Programm identische Ergebnisse liefern, leicht unterschiedliche Resultate erzeugen oder sogar fehlschlagen, etwa wenn eine benötigte Bibliothek nicht verfügbar ist.

:::{figure} ../../assets/figures/preface/python_figma.pdf
:align: center
:::