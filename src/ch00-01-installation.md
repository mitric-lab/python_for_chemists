# Installation

**WARNING: Translate and adapt!**

### Python Distribution

We highly recommend to use the [Anaconda](https://www.anaconda.com) or 
the [Mamba](https://github.com/mamba-org/mamba) Python distribution. Between 
these two options, Mamba is the top choice since it is much faster than 
Anaconda. 
Follow one of the following instructions to install Anaconda or Mamba. 
Do **NOT** install both distributions! 
Your Python version for this class should be 3.5 or higher.

#### Anaconda
Follow the installation guide for your operating system 
[here](https://docs.anaconda.com/anaconda/install/index.html).

#### Mamba
1. Download Mambaforge installer acoording to your operating system 
[here](https://github.com/conda-forge/miniforge#mambaforge).
2. Execute the installer:
  - Unix-like platforms (Mac OS & Linux):

    Open the terminal and call
    ```bash
    bash [Mambaforge-installer].sh
    ```
    where `[Mambaforge-installer]` should be replaced by your installer, 
    e.g. `Mambaforge-Linux-x86_64` or `Mambaforge-MacOSX-x86_64`.
  - Windows:

    Double click the installer on the file browser.

### Integrated Development Environment (IDE) or Editor recommendations

Although it is possible to write Python code using a conventional 
text editor, an integrated development environment can boost your 
programming experience greatly. Therefore, we shall show you some 
recommendations of IDEs and specialised editors here.

- [Jupyter Notebook](https://jupyter.org)
    - will be automatically available after you have installed 
      Anaconda/Miniconda.
    - if you are not yet familiar with Jupyter notebooks, you can find 
      many tutorials in the web, e.g. [Tutorial](https://www.dataquest.io/blog/jupyter-notebook-tutorial/).

- [Spyder](https://www.spyder-ide.org)
    - is a full Python IDE with the focus on scientific development.
    - will be automatically installed if you install the full Anaconda
      package.

- [Visual Studio Code](https://code.visualstudio.com)
    - is a more lightweight text editor for coding purposes and one 
      of the most used ones. 
    - Although it is officially not a real IDE, it has almost all
      features of an IDE. 

- [PyCharm](https://www.jetbrains.com/de-de/pycharm/)
    - commercial (free for students) IDE with a lot of functionalities. 
    - is usually more suitable for very large and complex Python projects
      and can be too cumbersome for small projects like the ones used here
      in the lecture.

- Vim/[NeoVim](https://neovim.io)
    - command-line text editor which is usually pre installed on all 
      unix-like computers.
    - can be very unwieldy at the beginning,
      since Vim has a very large number of keyboard shortcuts and is not
      as beginner friendly as a typical IDE or Jupyter notebooks. 
    - is extremely configurable and can be adapted with a little effort to
      a very extensive and comfortable editor.
    - still one of the most used editors.

### Dependencies

During the lecture we will use different Python libraries, which can be 
conveniently installed via a package manager. The Anaconda distribution 
comes with the package manager conda while the Mamba distribution has mamba. 
Please install these packages by executing the following commands in your 
terminal. 

```admonish note
If you are using Windows, execute the commands in Anaconda Prompt 
or Miniforge Prompt that comes with your Anaconda or Mambaforge installation.
```

| Library | conda | mamba |
| ------- | ----- | ----- |
| [NumPy](https://numpy.org) | `conda install numpy` | `mamba install numpy` |
| [Matplotlib](https://matplotlib.org) | `conda install -c conda-forge matplotlib` | `mamba install matplotlib` |
| [SciPy](https://scipy.org) | `conda install scipy` | `mamba install scipy` |
| [SymPy](https://www.sympy.org/en/index.html) | `conda install sympy` | `mamba install sympy` |
| [Jupyter](https://jupyter.org) | Preinstalled with Anaconda | `mamba install notebook` |
| [pytest](https://docs.pytest.org/en/7.2.x/) | `conda install pytest` | `mamba install pytest` |
| [ipytest](https://github.com/chmp/ipytest) | `conda install -c conda-forge ipytest` | `mamba install ipytest` |
| [Hypothesis](https://hypothesis.readthedocs.io/en/latest/) | `conda install -c conda-forge hypothesis` | `mamba install hypothesis` |


