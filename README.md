# Programmierkurs für Chemiker

[![Deploy](https://github.com/mitric-lab/python_for_chemists_ss24/workflows/Deploy/badge.svg)](https://github.com/mitric-lab/python_for_chemists_ss24/actions/workflows/deploy.yml)

Repository for the bachelor course "Programmierkurs für Chemiker"
in summer term 2026 at the University of Würzburg.

## Development

Install [uv](https://docs.astral.sh/uv/getting-started/installation/). To preview
the book locally after generating the student notebooks, run:

```bash
uv run python scripts/generate_student_notebooks.py
uv run jupyter book start
```

The site will be available at [http://localhost:3000](http://localhost:3000).

## Authoring Notes

Use colon fences for admonitions:

````md
:::{note}
Admonition content.
:::
````

Use MyST math directives only for numbered equations. Numbered equations must
have both `:enumerated: true` and a `:label:`:

````md
```{math}
:enumerated: true
:label: eq:example

E = mc^2
```
````

Write unlabelled display equations with double dollar signs instead of a math
directive:

```md
$$
E = mc^2
$$
```

## Notebook Generation

Student problem-set notebooks are generated from instructor solution notebooks.
Run the generator before building the book:

```bash
uv run python scripts/generate_student_notebooks.py
```

The generator looks for notebooks named `instructors/pset_X_solution.ipynb` and
creates two student-facing notebooks:

- A **rich** notebook for rendering the book page:
  - `chapters/problem_sets/pset_X_rich.ipynb`
- A **clean** notebook for students to download:
  - `assets/downloads/problem_sets/pset_X.ipynb`

For both generated notebooks, code-cell content between `# SOLUTION-START` and
`# SOLUTION-END` is replaced with `# your code here`, preserving the original
indentation and roughly the same number of lines. Code-cell content between
`# INSTRUCTOR-ONLY-START` and `# INSTRUCTOR-ONLY-END` is removed without adding a
student placeholder. Code-cell outputs are removed and execution counts are reset.

Markdown cells are normalised for students: headings starting with
`# Solution of ...` become regular problem-set headings, and solution-specific
section labels such as `(sec:pset_1_solution)=` are rewritten to the student label.
In the first markdown cell, the clean notebook also removes page frontmatter and
problem-set section labels.

The rich notebook gets standardised page frontmatter, a `(sec:pset_X)=` section
label, and a download link to the clean notebook. This keeps the rendered problem
set page and the downloadable notebook in sync.

Expected outputs can be shown on the rendered pages by tagging an instructor code cell
with the `expected-output` cell tag. The generator reads the original instructor
cell output and injects a small expected-output markdown block after the generated
code cell in the rich notebook only. Stream outputs and `text/plain` outputs are
embedded as HTML. PNG, JPEG, and GIF outputs are written to
`assets/figures/problem_sets/expected_outputs/` and referenced from the injected
HTML. Existing expected-output images for that problem set are cleaned up before
new ones are written.

The GitHub Pages deploy workflow runs the same generator automatically before
building the HTML site, so the published book always uses the generated student
notebooks.
