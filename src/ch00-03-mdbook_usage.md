## Bedienung dieser Webseite

**WARNING: Translate and adapt!**

Dieser Abschnitt gibt eine Einführung in die Bedienung dieses 
Vorlesungsskripts. Die Inhalte sind in *Kapitel* organisiert und 
jedes Kapitel stellt eine eigene Seite dar. Die Kapitel sind in einer
Hierarchie von *Unterkapiteln* organisiert. Typischerweise ist jedes
Kapitel in eine Reihe von *Überschriften* unterteilt.

### Navigation

Es gibt mehrere Möglichkeiten, um durch die Kapitel eines Buches zu
navigieren.

Die **Seitenleiste** auf der linken Seite bietet eine Liste aller
Kapitel. Wenn Sie auf einen der Kapiteltitel klicken, wird diese Seite
geladen.

Die Seitenleiste erscheint möglicherweise nicht automatisch, wenn das
Fenster zu schmal ist, insbesondere auf Mobilgeräten. In diesem
Fall kann das Menüsymbol (drei horizontale Balken) oben links auf der
Seite gedrückt werden, um die Seitenleiste zu öffnen und zu schließen.

Die **Pfeilschaltflächen** links und rechts mittig auf der Seite können 
verwendet werden, um zum vorherigen oder nächsten Kapitel zu navigieren.

Die **Pfeiltasten** links und rechts auf der Tastatur können verwendet
werden, um zum vorherigen oder nächsten Kapitel zu navigieren.

### Top-Menüleiste

Die Menüleiste oben auf der Seite bietet einige Symbole zur Interaktion
mit dem Skript.

| Icon | Beschreibung |
|------|-------------|
| <i class="fa fa-bars"></i> | Öffnet und schließt die Seitenleiste. |
| <i class="fa fa-paint-brush"></i> | Öffnet eine Dropdown-Liste, um ein anderes Farbschema auszuwählen. |
| <i class="fa fa-search"></i> | Öffnet eine Suchleiste zum Suchen im Buch. |
| <i class="fa fa-print"></i> | Fordert den Webbrowser auf, das Skript für den Druck zu formatieren. |

Das Tippen auf die Menüleiste scrollt die Seite nach oben.

### Suchen

Dieses Vorlesungsskript verfügt über ein integriertes Suchsystem.
Durch Drücken des Suchsymbols (<i class="fa fa-search"></i>) in der
Menüleiste oder Drücken der Taste `S` auf der Tastatur wird ein
Eingabefeld zum Eingeben von Suchbegriffen geöffnet.
Das Eingeben von Begriffen zeigt dann übereinstimmende Kapitel und
Abschnitte in Echtzeit an.

Durch Klicken auf eines der Ergebnisse können Sie zu diesem Abschnitt
wechseln. Die Pfeiltasten nach oben und unten können verwendet werden,
um die Ergebnisse zu durchsuchen, und die Eingabetaste öffnet den
hervorgehobenen Abschnitt.

Nach dem Laden eines Suchergebnisses werden die übereinstimmenden
Suchbegriffe im Text hervorgehoben. Durch Klicken auf ein hervorgehobenes
Wort oder Drücken der Taste `Esc` wird die Hervorhebung entfernt.

### Codeblöcke

Codeblöcke enthalten ein Kopiersymbol <i class="fa fa-copy"></i>, das den
Codeblock in die lokale Zwischenablage kopiert.

Hier ist ein Beispiel:

```python
print("Hello, World!")
```

Oft wird die Anweisung `assert` in Codebeispielen verwendet, um Ihnen
den Wert einer Variablen zu zeigen. Da die Codeblöcke in diesem Dokument
nicht interaktiv sind (Sie können sie nicht einfach in Ihrem Browser
ausführen), ist es nicht möglich, den Wert der Variablen auf dem Bildschirm
anzuzeigen. Daher stellen wir sicher, dass alle Codeblöcke in diesem
Skript ohne Fehler ausgeführt werden können, wenn es nicht anders 
angegeben ist. Auf diese Weise können wir den Wert einer Variablen durch
die Verwendung von `assert` darstellen. Das folgende Codebeispiel zeigt
beispielsweise, dass die Variable `a` den Wert 2 hat:

```python
a = 2
assert a == 2
```

Die `assert`-Anweisung überprüft, ob die Bedingung, die ihr übergeben
wird, `True` ist. Wenn die Bedingung `False` ist, wird eine
`AssertionError`-Ausnahme ausgelöst. 

