## Erste Schritte

### Python-Distribution

Wir empfehlen die Verwendung der
[Miniforge](https://github.com/conda-forge/miniforge) Distribution.
Falls Sie bereits die [Anaconda](https://www.anaconda.com) Distribution
installiert haben, können Sie auch diese für den Kurs verwenden.
Folgen Sie die untenstehenden Anweisungen, um Miniforge zu installieren.
Sie sollten eine Python-Version von 3.9 oder höher verwenden.

#### Installation
1. Laden Sie den Miniforge-Installer gemäß Ihrem Betriebssystem
   [hier](https://github.com/conda-forge/miniforge#miniforge3) herunter.
2. Führen Sie den Installer aus:
    - Mac OS & Linux:
        Öffnen Sie das Terminal, navigieren Sie zu Ihrem Downloads-Ordner und
        rufen Sie
        ```bash
        bash [Miniforge-Installer].sh
        ```
        auf, wobei `[Miniforge-Installer]` durch den Namen des Installers
        ersetzt werden sollte, z.B. `Miniforge3-Linux-x86_64.sh` oder
        `Miniforge3-MacOSX-arm64.sh`.
    - Windows:
        Doppelklicken Sie auf den Installer im Explorer.

### Python-Pakete

Während der Vorlesung werden wir mehrere Python-Pakete verwenden, die
bequem mit dem Paketmanager `mamba` installiert werden können, der mit
der Miniforge-Distribution installiert wird.
Einige der wichtigsten Pakete sind in der Tabelle unten aufgeführt.
Bitte installieren Sie diese Pakete, indem Sie die folgenden Befehle
in Ihrem Terminal ausführen.

```admonish note title="Hinweis für Windows-Benutzer"
Wenn Sie Windows verwenden, führen Sie die Befehle in der
Miniforge Prompt aus, die mit Ihrer Miniforge-Installation
geliefert wird.
```

| Paket | Befehl |
| ----- | ------------ |
| [NumPy](https://numpy.org) | `mamba install numpy` |
| [Matplotlib](https://matplotlib.org) | `mamba install -c conda-forge matplotlib` |
| [SciPy](https://scipy.org) | `mamba install -c conda-forge scipy` |
| [Jupyter](https://jupyter.org) | `mamba install jupyter` |


#### Jupyter Notebook

Wenn Sie ein Anfänger in Python sind, empfehlen wir die Verwendung von
[Jupyter Notebook](https://jupyter.org) zum Schreiben und Ausführen Ihres
Python-Codes. Es handelt sich um ein blockweise ausführbares Dokument, das
Code, Text und Grafiken enthalten kann. Um Ihre Python-Installation zu
testen und sich mit Jupyter Notebook vertraut zu machen, laden Sie bitte
die Datei <em>Crashkurs: Jupyter Notebooks und Python</em> 
herunter, die in WueCampus bereitgestellt ist. Sie können das Notebook öffnen
indem Sie
```bash
jupyter notebook
```
in Ihrem Terminal eingeben. Anschließend öffnet sich ein neuer Tab in 
Ihrem Webbrowser, wo Sie zum Ordner navigieren können, in dem sich das 
Jupyter Notebook befindet. Öffnen Sie das Notebook und führen Sie alle 
Code-Zellen aus. Wenn Ihre Installation erfolgreich war, sollten Sie 
keinen Fehler erhalten. Zur Bedienung des Notebooks finden Sie viele 
weitere Tutorials im Internet, wie z.B. 
[hier](https://www.dataquest.io/blog/jupyter-notebook-tutorial/).

### Empfehlungen für weitere integrierte Entwicklungsumgebungen (IDEs) oder Editoren

Wenn Sie bereits mit Python vertraut sind und eine erweiterte 
Entwicklungsumgebung nutzen möchten, sollten Sie eine IDE oder einen 
spezialisierten Editor verwenden. Obwohl es möglich ist, Python-Code mit 
einem herkömmlichen Texteditor zu schreiben, kann eine gute IDE Ihre 
Programmiererfahrung erheblich verbessern. Deshalb stellen wir Ihnen hier 
einige Empfehlungen für IDEs und spezialisierte Editoren neben Jupyter 
Notebook vor.

- [Spyder](https://www.spyder-ide.org)
    - voll ausgestattete Python-IDE mit Schwerpunkt auf wissenschaftlichr
      Entwicklung
    - leichtgewichtig und einfach zu konfigurieren

- [Visual Studio Code](https://code.visualstudio.com)
    - einer der meistgenutzten Editoren trotz seiner Leichtigkeit
    - fast alle Funktionen einer IDE, obwohl es offiziell keine ist

- [PyCharm](https://www.jetbrains.com/de-de/pycharm/)
    - kommerzielle (für Studierende kostenlos) IDE mit vielen Funktionen
    - für sehr große und komplexe Python-Projekte geeignet
    - möglicherweise zu umständlich für kleine Projekte wie die in dieser
      Veranstaltung

- Vim/[NeoVim](https://neovim.io)
    - Kommandozeilen-Editor, der auf fast allen Unix-ähnlichen Computern 
      vorinstalliert ist
    - möglicherweise sehr unhandlich am Anfang aufgrund der vielen 
      Tastenkombinationen und der mangelnden Anfängerfreundlichkeit im Vergleich zu 
      einer typischen IDE oder Jupyter Notebooks
    - extrem konfigurierbar
    - mit einigem Aufwand zu einem sehr umfangreichen und 
      komfortablen Editor anpassbar
    - immer noch einer der meistgenutzten Editoren in der 
      Softwareentwicklung