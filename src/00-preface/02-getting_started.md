## Getting Started

### Python Distribution

We highly recommend to use the 
[Miniforge](https://github.com/conda-forge/miniforge) distribution.
If you already have the installed the 
[Anaconda](https://www.anaconda.com) distribution,
you can also use it for this course. 
Follow the instructions below to install Miniforge. 
You should use a Python version of 3.9 or higher.

#### Installation
1. Download the Miniforge installer according to your operating system
   [here](https://github.com/conda-forge/miniforge#miniforge3).
2. Execute the installer:
    - Unix-like platforms (Mac OS & Linux):
        Open the terminal, navigate to your Downloads folder and call
        ```bash
        bash [Miniforge-Installer].sh
        ```
        where `[Miniforge-Installer]` should be replaced by the name of the 
        installer, e.g. `Miniforge3-Linux-x86_64.sh` or 
        `Miniforge3-MacOSX-arm64.sh`.
    - Windows:
        Double-click the installer on the file explorer.

### Python Packages

During the lecture, we will use several Python packages that can be
conveniently installed using the package manager `mamba`, which is
installed with the Miniforge distribution.
Some of the most important packages are listed in the table below.
Please install these packages by executing the following commands
in your terminal.

```admonish note title="Note for Windows Users"
If you are using Windows, execute the commands in 
Miniforge Prompt that comes with your Miniforge installation.
```

| Package | Command |
| ----- | ------------ |
| [NumPy](https://numpy.org) | `mamba install numpy` |
| [Matplotlib](https://matplotlib.org) | `mamba install -c conda-forge matplotlib` |
| [SciPy](https://scipy.org) | `mamba install -c conda-forge scipy` |
| [Jupyter](https://jupyter.org) | `mamba install jupyter` |


#### Jupyter Notebook

If you are a beginner in Python, we recommend using 
[Jupyter Notebook](https://jupyter.org) for writing and executing your 
Python code. It is a blockwise executable document that can contain code, 
text, and graphics. To test your Python installation and to make yourself 
familiar with Jupyter Notebook, please download and open the Juptyer 
Notebook that is provided in WueCampus by calling

```bash
jupyter notebook
```

in your terminal. This will open a new tab in your web browser where you
can navigate to the folder where the Jupyter Notebook is located and open
it. If your installation was successful, you should be able to execute
all code cells in the notebook without any errors. You can can find many
more tutorials on the internet, e.g. 
[here](https://www.dataquest.io/blog/jupyter-notebook-tutorial/).


### Recommendations for Further Integrated Development Environments (IDEs) or Editors

If you are already familiar with Python and want to use a more
sophisticated development environment, you might want to use an IDE or a
specialised editor. Although it is possible to write Python code using a 
conventional text editor, a good IDE can boost your programming experience 
greatly. Therefore, we shall show you some recommendations of IDEs and 
specialised editors here besides Jupyter Notebook.

- [Spyder](https://www.spyder-ide.org)
    - full-featured Python IDE with a focus on scientific development
    - lightweight and easy to configure

- [Visual Studio Code](https://code.visualstudio.com)
    - one of the most used editors despite its lightweightness
    - almost all features of an IDE although it is officially not one

- [PyCharm](https://www.jetbrains.com/de-de/pycharm/)
    - commercial (free for students) IDE with lots of functionalities
    - most suitable for very large and complex Python projects
    - possibly too cumbersome for small projects like the ones in this 
      lecture

- Vim/[NeoVim](https://neovim.io)
    - command-line editor pre-installed on almost all Unix-like
      computers
    - possibly very unwieldy at the beginning due to the many 
      keyboard shortcuts and the lack of beginner-friendliness like 
      a typical IDE or Jupyter Notebooks
    - extremely configurable
    - adaptable to a very extensive and 
      comfortable editor with some effort
    - still one of the most used editors

