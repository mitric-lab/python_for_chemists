# Motivation

## Warum Python?

### Beliebtheit

Python ist eine der am weitesten verbreiteten Programmiersprachen und
besonders anfängerfreundlich.
Die folgende Abbildung zeigt die Beliebtheit einiger Programmiersprachen
gemäß des ***P**opularit**Y** of **P**rogramming **L**anguage*
([PYPL](https://pypl.github.io/PYPL.html)) Indexes (Stand: Jan. 2024).

![Popularity of Programming languages](figures/00_preface/popularity_pypl_202401.svg)
*Beliebtheit von Programmiersprachen. Die Beliebtheiten sind dem
[PYPL Index](https://pypl.github.io/PYPL.html) entnommen.*

### Weniger Code, weniger Arbeit

Python benötigt oft deutlich weniger Code als kompilierte Sprachen wie
Java oder C/C++ um die gleichen Algorithmen zu implementieren.
![LOC of PL](figures/00_preface/loc.png)
*Programmlänge, gemessen in der Anzahl der nicht-kommentierten
Quellcodezeilen* (**L**ines **O**f **C**ode (LOC))*.*

Darüber hinaus ist es oft möglich, die gleiche Aufgabe mit einer 
Skriptsprache, wie z.B. Python, in deutlich weniger Zeit zu erledigen,
als mit einer kompilierten Sprache.
![Hours of work to code](figures/00_preface/hours.png)
*Entwicklungszeit um eine bestimme Programmieraufgabe zu erledigen, 
gemessen in Stunden.*

### Aber Python ist langsam?!

Ein häufiges Argument gegen die Verwendung von Python ist, dass es sich
um eine vergleichsweise langsame Sprache handelt.
An dieser Aussage ist etwas dran.
Daher wollen wir uns einen der Gründe dafür ansehen und dann erklären,
warum (in den meisten Fällen) dies für uns kein Problem darstellt.

#### Python ist (auch) eine interpretierte Sprache
Python selbst ist ein C-Programm, das zuerst den Quellcode zum sog. 
[*Bytecode*](https://de.wikipedia.org/wiki/Bytecode) kompiliert
("übersetzt") und diesen dann interpretiert und ausführt.
Dies ist im Gegensatz zu kompilierten Sprachen wie C, C++, Rust, etc., 
bei denen der Quellcode in Maschinencode kompiliert wird. Der Compiler 
kann dabei viele Optimierungen am Code vornehmen, was zu einer kürzeren
Laufzeit führt.
Dieses Verhalten lässt sich an einem einfachen Beispiel zeigen: 
Eine naive Implementierung, die alle ungeraden Zahlen bis 100 Millionen
aufsummiert. Diese könnte wie folgt aussehen:

```python
s = 0
for i in range(100_000_000):
    if i % 2 == 1:
        s += i
```
Dieser Code benötigt auf dem Computer des Autors ca. 8 Sekunden.
Nun wird der gleiche Algorithmus in einer kompilierten Sprache 
(in diesem Fall *Rust*) implementiert, um den Einfluss des Compilers 
zu zeigen. 
Die Details dieses Codes und die Programmiersprache sind an dieser
Stelle nicht wichtig.

```rust,no_run,no_playground
let mut s: usize = 0;
for i in 0..100_000_000 {
    if i % 2 == 1 {
        s += i;
    }
}
```

Dieser Code hat tatsächlich überhaupt keine Laufzeit und wird sofort
ausgewertet. Der Compiler ist schlau genug zu verstehen, dass alles zur
Compile-Zeit berechnet werden kann und ersetzt einfach den Wert für die
Variable `s`. Dies macht nun deutlich, dass kompilierte Sprachen von
Methoden profitieren können, die interpretierte Sprachen einfach aufgrund
ihrer Herangehensweise nicht haben. Wir haben jedoch bereits gesehen,
dass kompilierte Sprachen in der Regel mehr Codezeilen und mehr Arbeit
erfordern. Darüber hinaus gibt es in der Regel viel mehr Konzepte in
kompilierten Sprachen zu lernen.

#### Python kann sehr performant sein

Während dieser Veranstaltung werden wir oft Python-Bibliotheken wie NumPy 
oder SciPy für mathematische Algorithmen und insbesondere lineare Algebra.
Diese Pakete bringen zwei große Vorteile. Einerseits ermöglichen sie die
sehr einfache Verwendung komplizierter Algorithmen und andererseits sind
diese Pakete in kompilierten Sprachen wie C oder Fortran geschrieben.
Auf diese Weise können wir von den Leistungsvorteilen profitieren,
ohne eine möglicherweise kompliziertere Sprache lernen zu müssen.

