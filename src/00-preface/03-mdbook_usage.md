## Bedienung dieser Webseite

Dieser Abschnitt gibt eine Einführung in die Bedienung dieses 
Vorlesungsskripts. Dieses ist in *Kapitel* organisiert. Jedes Kapitel ist eine
eigene Seite. Kapitel sind in einer Hierarchie von Unterkapiteln
verschachtelt. Typischerweise ist jedes (Unter-)Kapitel in eine Reihe von
*Überschriften* unterteilt.

### Navigation

Es gibt mehrere Methoden, um durch die Kapitel eines Buches zu navigieren.

Die **Seitenleiste** auf der linken Seite bietet eine Liste aller Kapitel.
Durch Klicken auf einen der Kapiteltitel wird die entsprechende Seite geladen.

Die Seitenleiste erscheint möglicherweise nicht automatisch, wenn das 
Fenster zu schmal ist, insbesondere auf mobilen Displays.
In diesem Fall kann das Menüsymbol (drei horizontale Balken) oben links auf
der Seite gedrückt werden, um die Seitenleiste zu öffnen und zu schließen.

Die **Pfeiltasten** links und rechts neben der Seite können verwendet werden, 
um zum vorherigen oder nächsten Kapitel zu navigieren.

Die **Pfeiltasten** auf der Tastatur können verwendet werden, um zum vorherigen 
oder nächsten Kapitel zu navigieren.

### Obere Menüleiste

Die Menüleiste oben auf der Seite bietet einige Symbole zur Interaktion mit den Notizen.

| Symbol | Beschreibung |
|--------|--------------|
| <i class="fa fa-bars"></i> | Öffnet und schließt die Seitenleiste. |
| <i class="fa fa-paint-brush"></i> | Öffnet eine Liste zur Auswahl eines anderen Farbthemas. |
| <i class="fa fa-search"></i> | Öffnet eine Suchleiste für die Suche im Skript. |
| <i class="fa fa-print"></i> | Fordert den Webbrowser auf, das Skript zu drucken. |
| <i class="icon-uw"></i> | Öffnet den WueCampus-Kursraum zu dieser Veranstaltung. |

Ein Klick auf die Menüleiste scrollt die Seite nach oben.

### Suche

Dieses Vorlesungsskript verfügt über ein eingebautes Suchsystem.
Durch Drücken des Suchsymbols (<i class="fa fa-search"></i>) in der 
Menüleiste oder Drücken der Taste `S` auf der Tastatur wird ein Eingabefeld 
zum Eingeben von Suchbegriffen geöffnet.
Das Eingeben von Begriffen zeigt passende Kapitel und Abschnitte in 
Echtzeit an.

Durch Klicken auf eines der Ergebnisse wird zu diesem Abschnitt gesprungen.
Die Pfeil-nach-oben- und Pfeil-nach-unten-Tasten können verwendet werden, um
die Ergebnisse zu navigieren, und Enter öffnet den markierten Abschnitt.

Nach dem Laden eines Suchergebnisses werden die passenden Suchbegriffe 
im Text hervorgehoben.
Durch Klicken auf ein hervorgehobenes Wort oder Drücken der `Esc`-Taste
werden die Hervorhebungen entfernt.

### Codeblöcke

Codeblöcke besitzen ein Kopieren-Symbol <i class="fa fa-copy"></i>, das den
Codeblock in die lokale Zwischenablage kopiert.

Hier ein Beispiel:
```python
print("Hello, World!")
```

Wir verwenden die `assert`-Anweisung häufig in Codeblöcken, um Ihnen den
Wert einer Variablen anzuzeigen. Da die Codeblöcke in diesem Dokument nicht 
interaktiv sind (Sie können sie nicht einfach in Ihrem Browser ausführen), 
ist es nicht möglich, den Wert der Variablen auf dem Bildschirm auszugeben. 
Deshalb stellen wir für Sie sicher, dass alle Codeblöcke in diesem 
Vorlesungsskript fehlerfrei laufen und wir auf diese Weise den Wert einer 
Variablen durch die Verwendung von `assert` darstellen können. Der folgende 
Codeblock zeigt beispielsweise, dass die Variable „a“ den Wert 2 hat:
```python
a = 2
assert a == 2
```

Falls die Bedingung `False` ergeben würde, würde die `assert`-Anweisung 
einen `AssertionError` auslösen.

