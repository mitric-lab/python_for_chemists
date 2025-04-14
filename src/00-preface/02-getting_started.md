## Erste Schritte

### Intallieren der Python-Distribution

Auf vielen Computern ist Python zwar bereits vorinstalliert, dennoch empfehlen wir, Python über die [Miniforge](https://github.com/conda-forge/miniforge) Distribution zu installieren. Sollten Sie bereits die [Anaconda](https://www.anaconda.com) Distribution nutzen, ist auch diese für den Kurs geeignet. Um Miniforge zu installieren, besuchen Sie die [offizielle Miniforge-Webseite](https://github.com/conda-forge/miniforge) und folgen Sie den Anweisungen im Abschnitt *Install* für Ihr Betriebssystem. Bei Fragen oder Schwierigkeiten während der Installation können Sie diese gerne in der ersten Vorlesung ansprechen.

### Python-Pakete

Während des Kurses verwenden wir verschiedene Python-Pakete, die einfach mit dem Paketmanager `mamba` installiert werden können. Dieser wird automatisch mit Miniforge installiert. Die wichtigsten Pakete sind in der folgenden Tabelle aufgelistet. Bitte installieren Sie diese, indem Sie die angegebenen Befehle in Ihrem Terminal ausführen.

| Paket | Befehl |
| ----- | ------------ |
| [NumPy](https://numpy.org) | `mamba install numpy` |
| [Matplotlib](https://matplotlib.org) | `mamba install -c conda-forge matplotlib` |
| [SciPy](https://scipy.org) | `mamba install -c conda-forge scipy` |
| [Jupyter](https://jupyter.org) | `mamba install jupyter` |

```admonish note title="Hinweis für Windows-Benutzer"
Unter Windows führen Sie die Befehle bitte in der Miniforge Prompt aus, die mit der Miniforge-Installation bereitgestellt wird.
```

#### Jupyter Notebook

Für Python-Einsteiger empfehlen wir [Jupyter Notebook](https://jupyter.org). Jupyter Notebooks ermöglichen es, Python-Code blockweise auszuführen und dabei Code, Text und Grafiken übersichtlich zu kombinieren. Um Ihre Installation zu testen und sich mit Jupyter Notebook vertraut zu machen, laden Sie die Datei *Crashkurs: Jupyter Notebooks und Python* herunter, die auf WueCampus bereitgestellt ist.

Öffnen Sie das Notebook, indem Sie folgenden Befehl in Ihrem Terminal eingeben:
```bash
jupyter notebook
```
Ihr Webbrowser öffnet sich anschließend automatisch. Navigieren Sie zum Speicherort des heruntergeladenen Notebooks, öffnen Sie es durch Doppelklicken und führen Sie die Code-Zellen aus. Wenn die Installation korrekt war, sollten dabei keine Fehler auftreten. Weitere hilfreiche Tutorials zu Jupyter Notebooks finden Sie beispielsweise [hier](https://www.dataquest.io/blog/jupyter-notebook-tutorial/).

### Empfehlungen für integrierte Entwicklungsumgebungen (IDEs) oder Editoren

Wenn Sie bereits Erfahrungen mit Python haben, empfehlen wir Ihnen eine integrierte Entwicklungsumgebung (IDE) oder einen spezialisierten Editor. Diese bieten Funktionen wie Syntax-Highlighting, automatische Code-Vervollständigung und Debugging-Tools, die Ihre Arbeit deutlich effizienter gestalten. Beliebte IDEs und Editoren sind beispielsweise [Spyder](https://www.spyder-ide.org), [Visual Studio Code](https://code.visualstudio.com) oder [PyCharm](https://www.jetbrains.com/de-de/pycharm/).