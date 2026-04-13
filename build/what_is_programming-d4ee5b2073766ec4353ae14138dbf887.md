---
numbering:
  headings: false
---

(sec:what_is_programming)=

# Was ist Programmieren?

Programmieren bedeutet, einem Computer Anweisungen so aufzuschreiben, dass er 
sie Schritt für Schritt ausführen kann. Diese Anweisungen werden in einer Textdatei 
gespeichert und bilden zusammen ein Programm.

Im Gegensatz zu Menschen, deren Anweisungen eine gewisse Mehrdeutigkeit zulassen, 
verarbeitet ein Computer Informationen nach genau festgelegten Regeln. 
Er "versteht" nicht, was wir meinen, sondern führt nur genau die Befehle aus, 
die ihm gegeben werden.

## Grundaufbau eines Computers

Die zentralen Komponenten eines Computers sind Eingabe, Verarbeitung, Speicher 
und Ausgabe. [Abb.&nbsp;%s](#fig:computer_components) zeigt dieses Grundprinzip. 
Die Verarbeitung übernimmt die _**C**entral **P**rocessing **U**nit_ (CPU). 
Sie besteht aus einer Steuereinheit (Control Unit), die den Ablauf eines Programms 
steuert, und einer arithmetisch-logischen Einheit (_**A**ritmetic **L**ogic **U**nit_, ALU), 
die die Rechenoperationen ausführt. Daten und Programme liegen im Speicher und 
werden von dort gelesen, verarbeitet und anschließend wieder abgelegt. Ergebnisse 
werden schließlich über die Ausgabe bereitgestellt, etwa auf dem Bildschirm oder 
in einer Datei.

:::{figure} ../../assets/figures/preface/computer_figma.svg
:align: center
:label: fig:computer_components
:width: 590.4 px
Grundlegende Komponenten eines Computers.
:::

## Ausführung von Programmen

Doch wie wird ein solches Programm tatsächlich ausgeführt? [Abb.&nbsp;%s](#fig:program_execution) 
zeigt zwei typische Wege. Vereinfacht gesagt wird bei kompilierten Sprachen 
wie C++ der Quellcode vor der Ausführung in Maschinencode übersetzt. Das Ergebnis 
ist eine ausführbare Datei, die anschließend vom Betriebssystem geladen und 
von der CPU ausgeführt werden kann. 
Python[^cpython] arbeitet anders: Der Quellcode wird zunächst in Bytecode überführt. 
Dieser Bytecode wird anschließend von einem Interpreter schrittweise ausgeführt 
und in konkrete Operationen auf der CPU umgesetzt. Diese Art der Ausführung ist 
ein Grund dafür, dass Python oft langsamer ist als kompilierte Sprachen.

[^cpython]: Streng genommen ist hier nur [*CPython*](wiki:CPython) gemeint, die verbreitetste 
Implementierung der Programmiersprache Python. Andere Implementierungen wie 
[*PyPy*](wiki:PyPy) oder [*Jython*](wiki:Jython) funktionieren teilweise anders.

:::{figure} ../../assets/figures/preface/code_figma.svg
:align: center
:label: fig:program_execution
:width: 534.6 px
Ausführung eines Programms in kompilierten und interpretierten Sprachen.
:::

Für Python ist damit jedoch noch nicht die ganze Ausführungssituation beschrieben. 
Neben dem Quellcode ist auch die Umgebung wichtig, in der ein Programm läuft.

## Python und seine Umgebung

Python ist nicht nur eine Programmiersprache, sondern auch Teil eines sehr großen 
Ökosystems aus Bibliotheken und Werkzeugen. Die vergleichsweise einfache Syntax 
und die hohe Ausdrucksstärke machen Python besonders angenehm zum Schreiben und 
Erweitern von Programmen. Das hat dazu beigetragen, dass für sehr viele Anwendungen 
bereits passende Bibliotheken existieren. Viele davon sind intern in effizienten 
Sprachen wie C oder C++ implementiert, sodass auch rechenintensive Aufgaben gut 
bearbeitet werden können.

In diesem großen Ökosystem laufen Python-Programme daher nicht isoliert, sondern immer in einer 
bestimmten Python-Umgebung. Eine solche Umgebung verknüpft eine konkrete Python-Version 
mit genau festgelegten Paketversionen. [Abb.&nbsp;%s](#fig:python_environment) 
verdeutlicht das: Dieselbe Python-Datei kann in verschiedenen Umgebungen unterschiedlich 
ausgeführt werden. Je nach installierten Versionen erhält man identische Ergebnisse, 
leicht andere Resultate oder eine Fehlermeldung, wenn ein benötigtes Paket fehlt oder 
inkompatibel ist. Deshalb ist die verwendete Umgebung ein wichtiger Bestandteil 
eines Programms.

:::{figure} ../../assets/figures/preface/python_figma.svg
:align: center
:label: fig:python_environment
:width: 600 px
Dieselbe Python-Datei kann je nach Umgebung unterschiedlich ausgeführt werden.
:::
