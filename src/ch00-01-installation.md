## Installation

### Python-Distribution

Wir empfehlen die Verwendung der 
[Miniforge](https://github.com/conda-forge/miniforge)-Distribution.
Falls Sie bereits [Anaconda](https://www.anaconda.com) installiert haben,
können Sie auch diese Distribution problemlos für diese Veranstaltung
verwenden.
Folgen Sie die untenstehenden Anweisungen, um Mamba zu installieren.
Sie sollten eine Python-Version von 3.9 oder höher verwenden.

#### Installation
1. Laden Sie den Miniforge-Installer für Ihr Betriebssystem 
   [hier](https://github.com/conda-forge/miniforge#miniforge3) herunter.
2. Führen Sie den Installer aus:
    - Unixoide Plattformen (Mac OS & Linux):
        Öffnen Sie das Terminal und führen Sie
        ```bash
        bash [Miniforge-Installer].sh
        ```
        auf, wobei `[Miniforge-Installer]` durch den Namen des Installers 
        ersetzt werden sollte, z.B. `Miniforge3-Linux-x86_64.sh` oder 
        `Miniforge3-MacOSX-arm64.sh`.
    - Windows:
        Doppelklicken Sie auf den Installer im Explorer.


### Empfehlungen für eine Integrierte Entwicklungsumgebung (IDE) oder Editor
Auch wenn es möglich ist, Python-Code mit einem herkömmlichen Texteditor
zu schreiben, kann eine integrierte Entwicklungsumgebung das Programmieren 
erheblich erleichtern. Daher zeigen wir Ihnen hier einige Empfehlungen für
IDEs und spezialisierte Editoren.

- [Jupyter Notebook](https://jupyter.org)
    - blockweise ausführbare Dokumente, die Code, Text 
      und Grafiken enthalten können
    - viele Tutorials im Internet verfügbar, z.B.
      [hier](https://www.dataquest.io/blog/jupyter-notebook-tutorial/)

- [Spyder](https://www.spyder-ide.org)
    - vollwertige Python-IDE mit dem Fokus auf wissenschaftliche 
      Entwicklung
    - leichte und übersichtliche Konfiguration

- [Visual Studio Code](https://code.visualstudio.com)
    - einer der meistgenutzten Editoren trotz der Leichtgewichtigkeit
    - fast alle Funktionen einer IDE obwohl es offiziell keine echte IDE ist

- [PyCharm](https://www.jetbrains.com/de-de/pycharm/)
    - kommerzielle (für Studierende kostenlose) IDE mit vielen Funktionen
    - geeignet für sehr große und komplexe Python-Projekte
    - möglicherweise zu umständlich für kleine Projekte wie die in der 
      Vorlesung

- Vim/[NeoVim](https://neovim.io)
    - Kommandozeilen-Editor, normalerweise auf allen unixoiden 
      Computern vorinstalliert
    - möglicherweise sehr umständlich am Anfang wegen der vielen 
      Tastaturkürzel und der fehlenden Anfängerfreundlichkeit
      wie bei einer typischen IDE oder Jupyter Notebooks.
    - extrem konfigurierbar, mit etwas Aufwand zu einem sehr 
      umfangreichen und komfortablen Editor anpassbar
    - immer noch einer der meistgenutzten Editoren


### Python-Module

Während der Vorlesung werden wir verschiedene Python-Module verwenden, 
die bequem über den Paketmanager `mamba` installiert werden können,
welcher mit der Miniforge-Distribution installiert wird.
Bitte installieren Sie diese Pakete, indem Sie die folgenden Befehle
mit Kommandozeile ausführen.

```admonish note title="Hinweis für Windows-Nutzer"
Falls Sie Windows verwenden, führen Sie die Befehle in der
Miniforge Prompt aus, die mit Ihrer Miniforge-Installation
mitgeliefert wurde.
```

| Paket | Mamba-Befehl |
| ----- | ------------ |
| [NumPy](https://numpy.org) | `mamba install numpy` |
| [Matplotlib](https://matplotlib.org) | `mamba install matplotlib` |
| [SciPy](https://scipy.org) | `mamba install scipy` |
| [SymPy](https://www.sympy.org/en/index.html) | `mamba install sympy` |
| [Jupyter](https://jupyter.org) | `mamba install notebook` |
| [pytest](https://docs.pytest.org/en/7.4.x/) | `mamba install pytest` |
| [ipytest](https://github.com/chmp/ipytest) | `mamba install ipytest` |
| [Hypothesis](https://hypothesis.readthedocs.io/en/latest/) | `mamba install hypothesis` |

