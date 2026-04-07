# Programmierkurs für Chemiker

[![Deploy](https://github.com/mitric-lab/python_for_chemists_ss24/workflows/Deploy/badge.svg)](https://github.com/mitric-lab/python_for_chemists_ss24/actions/workflows/deploy.yml)

Repository for the bachelor course "Programmierkurs für Chemiker"
in summer term 2026 at the University of Würzburg.

## Development

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then run:

```bash
uv run jupyter book start
```

The site will be available at [http://localhost:3000](http://localhost:3000).

Before building the student-facing book, regenerate the student notebooks from
the instructor solution notebooks:

```bash
uv run python scripts/generate_student_notebooks.py
```

The generator looks for notebooks named `instructors/*_solution.ipynb` and
creates the corresponding student notebook in `chapters/problem_sets/`. Inside
code cells, content between `# SOLUTION-START` and `# SOLUTION-END` is replaced
with `# your code here`, preserving the original indentation and roughly the
same number of lines.

A typical local workflow is:

```bash
uv run python scripts/generate_student_notebooks.py
uv run jupyter book start
```

The GitHub Pages deploy workflow runs the same generator automatically before
building the HTML site, so the published book always uses the generated student
notebooks.
