#!/usr/bin/env python3
"""Generate student notebooks from instructor solution notebooks.

For each instructor notebook in `instructors/*_solution.ipynb`, we generate:
1. A "clean" notebook (for students to download) at `assets/downloads/problem_sets/pset_X.ipynb`
2. A "rich" notebook (for the website) at `chapters/problem_sets/pset_X_rich.ipynb`

Rich notebooks may include injected "expected output" blocks based on a cell tag
(`expected-output`) in the instructor notebook.
"""

import base64
import html
import json
import re
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INSTRUCTORS_DIR = ROOT / "instructors"
EXPECTED_OUTPUT_TAG = "expected-output"
EXPECTED_OUTPUT_IMAGE_DIR = ROOT / "assets/figures/problem_sets/expected_outputs"
CLEAN_NOTEBOOK_DIR = ROOT / "assets/downloads/problem_sets"
INSTRUCTOR_ONLY_START = "# INSTRUCTOR-ONLY-START"
INSTRUCTOR_ONLY_END = "# INSTRUCTOR-ONLY-END"

# TODO (future): Convert MyST equation directives to plain markdown/LaTeX for clean notebooks.

_SEC_LABEL_RE = re.compile(r"^\(sec:pset_\d+\)=\s*$")
_PSET_INDEX_RE = re.compile(r"pset_(\d+)_solution\.ipynb$")
_YAML_DELIM = "---"
_IMAGE_MIME_TYPES = ("image/png", "image/jpeg", "image/gif")


def load_notebook(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_notebook(path: Path, notebook: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(notebook, ensure_ascii=False, indent=1) + "\n", encoding="utf-8")


def problem_set_index(solution_path: Path) -> str:
    """Extract the problem-set index from a solution notebook filename."""
    match = _PSET_INDEX_RE.search(solution_path.name)
    if not match:
        raise ValueError(f"Could not extract problem set index from {solution_path.name}")
    return match.group(1)


def clean_target_path_for(solution_path: Path) -> Path:
    """Map an instructor solution notebook to its clean student notebook path."""
    if not solution_path.name.endswith("_solution.ipynb"):
        raise ValueError(f"Unsupported solution notebook name: {solution_path.name}")
    return CLEAN_NOTEBOOK_DIR / solution_path.name.replace("_solution.ipynb", ".ipynb")


def rich_target_path_for(solution_path: Path) -> Path:
    """Map an instructor solution notebook to its rich student notebook path."""
    if not solution_path.name.endswith("_solution.ipynb"):
        raise ValueError(f"Unsupported solution notebook name: {solution_path.name}")
    return ROOT / "chapters/problem_sets" / solution_path.name.replace(
        "_solution.ipynb", "_rich.ipynb"
    )


def strip_solution_blocks(source: list[str]) -> list[str]:
    """Replace solution regions with a placeholder and blank lines.

    Solution regions are delimited by:
    - `# SOLUTION-START`
    - `# SOLUTION-END`
    """
    result: list[str] = []
    i = 0

    while i < len(source):
        if source[i].strip() == "# SOLUTION-START":
            start = i
            indent = source[i][: len(source[i]) - len(source[i].lstrip())]
            i += 1
            while i < len(source) and source[i].strip() != "# SOLUTION-END":
                i += 1

            if i >= len(source):
                result.extend(source[start:])
                break

            end = i
            block_line_count = end - start + 1
            result.append(f"{indent}# your code here\n")
            result.extend("\n" for _ in range(max(0, block_line_count - 1)))
            i += 1
            continue

        result.append(source[i])
        i += 1

    return result


def strip_instructor_only_blocks(source: list[str]) -> list[str]:
    """Remove instructor-only regions without inserting student placeholders.

    Instructor-only regions are delimited by:
    - `# INSTRUCTOR-ONLY-START`
    - `# INSTRUCTOR-ONLY-END`
    """
    result: list[str] = []
    skipping = False

    for line in source:
        stripped = line.strip()
        if stripped == INSTRUCTOR_ONLY_START:
            skipping = True
            continue
        if stripped == INSTRUCTOR_ONLY_END:
            skipping = False
            continue
        if not skipping:
            result.append(line)

    return result


def normalize_student_markdown(source: list[str]) -> list[str]:
    """Adjust solution-specific labels and titles for the student notebook."""
    normalized: list[str] = []
    for line in source:
        line = line.replace("_solution)=", ")=")
        if line.startswith("# Solution of "):
            line = line.replace("# Solution of ", "# ", 1)
        normalized.append(line)
    return normalized


def strip_frontmatter_and_labels(markdown: str) -> str:
    """Strip YAML frontmatter and section labels like (sec:pset_1)= from markdown."""
    lines = markdown.splitlines(keepends=True)
    if not lines:
        return markdown

    # Remove a leading YAML frontmatter block.
    if lines[0].strip() == _YAML_DELIM:
        end = None
        for i in range(1, len(lines)):
            if lines[i].strip() == _YAML_DELIM:
                end = i
                break
        if end is not None:
            lines = lines[end + 1 :]
            if lines and lines[0].strip() == "":
                lines = lines[1:]

    # Remove MyST section labels like (sec:pset_1)=
    filtered: list[str] = []
    for line in lines:
        if _SEC_LABEL_RE.match(line.strip()):
            continue
        filtered.append(line)

    while filtered and filtered[0].strip() == "":
        filtered = filtered[1:]

    return "".join(filtered)


def build_rich_notebook_header(markdown: str, *, pset_index: str) -> str:
    """Build the standardized rich-notebook header (frontmatter + section label) for a problem set."""
    # These pages render under /chapters/problem-sets/<slug>, so go up two levels to reach /assets/.
    # This keeps downloads working even when the site is served under a non-root BASE_URL.
    url = f"../../assets/downloads/problem_sets/pset_{pset_index}.ipynb"
    filename = f"pset_{pset_index}.ipynb"

    body = strip_frontmatter_and_labels(markdown)
    if body and not body.endswith("\n"):
        body += "\n"

    header = (
        "---\n"
        "numbering:\n"
        "  headings: false\n"
        "downloads:\n"
        "  -\n"
        f"    title: {filename}\n"
        f"    url: {url}\n"
        f"    filename: {filename}\n"
        "---\n"
        "\n"
        f"(sec:pset_{pset_index})=\n"
        "\n"
    )

    return header + body


def cell_has_tag(cell: dict, tag: str) -> bool:
    """Return whether a notebook cell has the given metadata tag."""
    tags = cell.get("metadata", {}).get("tags", [])
    return tag in tags


def markdown_cell(source: str) -> dict:
    """Create a markdown cell with the given source."""
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": source.splitlines(keepends=True),
    }


def _as_text(value: object) -> str:
    if isinstance(value, list):
        return "".join(str(v) for v in value)
    if value is None:
        return ""
    return str(value)


def _escape_output_text(text: str) -> str:
    # Preserve output whitespace in regular HTML without using <pre>, which
    # some themes treat as a copyable code block.
    return html.escape(text.rstrip()).replace(" ", "&nbsp;").replace("\n", "<br/>")


def _output_has_image(output: dict) -> bool:
    data = output.get("data", {})
    if not isinstance(data, dict):
        return False
    return any(bool(data.get(mime_type)) for mime_type in _IMAGE_MIME_TYPES)


def output_to_html(
    output: dict,
    *,
    image_dir: Path,
    image_prefix: str,
    image_counter: int,
) -> str | None:
    """Convert a notebook output object into expected-output HTML."""
    output_type = output.get("output_type")

    if output_type == "stream":
        escaped = _escape_output_text(_as_text(output.get("text", "")))
        return (
            '<div class="expected-output-container">'
            '<div class="myst-jp-stream-output expected-output">'
            f'<div class="expected-output-text">{escaped}</div>'
            "</div>"
            "</div>"
        )

    data = output.get("data", {})
    if not isinstance(data, dict):
        data = {}

    for mime_type in _IMAGE_MIME_TYPES:
        payload = data.get(mime_type)
        if not payload:
            continue
        payload = _as_text(payload)
        extension = mime_type.split("/", 1)[1]
        image_path = image_dir / f"{image_prefix}_{image_counter:02d}.{extension}"
        image_path.write_bytes(base64.b64decode(payload))
        relative_image_path = image_path.relative_to(ROOT).as_posix()
        return (
            '<div class="expected-output-container expected-output-image">'
            '<div class="expected-output-image-frame">'
            '<div class="expected-output-label">Expected Output</div>'
            f'<img src="/{relative_image_path}" alt="Expected output" />'
            "</div>"
            "</div>"
        )

    text_plain = data.get("text/plain")
    text_plain = _as_text(text_plain)
    if text_plain.strip():
        escaped = _escape_output_text(text_plain)
        return (
            '<div class="expected-output-container">'
            '<div class="myst-jp-stream-output expected-output">'
            f'<div class="expected-output-text">{escaped}</div>'
            "</div>"
            "</div>"
        )

    return None


def expected_output_cell(
    cell: dict,
    *,
    image_dir: Path,
    image_prefix: str,
) -> dict | None:
    """Build a markdown cell that shows saved outputs as expected output."""
    rendered_outputs: list[str] = []
    image_counter = 1
    for output in cell.get("outputs", []):
        rendered = output_to_html(
            output,
            image_dir=image_dir,
            image_prefix=image_prefix,
            image_counter=image_counter,
        )
        if rendered:
            rendered_outputs.append(rendered)
            if _output_has_image(output):
                image_counter += 1

    if not rendered_outputs:
        return None

    return markdown_cell("\n\n".join(rendered_outputs))


def cleanup_expected_output_images(pset_index: str) -> None:
    """Remove previously materialized expected-output images for a given problem set."""
    EXPECTED_OUTPUT_IMAGE_DIR.mkdir(parents=True, exist_ok=True)
    for path in EXPECTED_OUTPUT_IMAGE_DIR.glob(f"{pset_index}_*"):
        if path.suffix.lower() in {".png", ".jpeg", ".jpg", ".gif"}:
            path.unlink()


def generate_student_notebook(solution_path: Path) -> tuple[Path, Path]:
    """Generate clean + rich student notebooks from an instructor notebook."""
    notebook = load_notebook(solution_path)
    clean = deepcopy(notebook)
    rich = deepcopy(notebook)

    clean_cells: list[dict] = []
    rich_cells: list[dict] = []

    clean_path = clean_target_path_for(solution_path)
    rich_path = rich_target_path_for(solution_path)
    pset_index = problem_set_index(solution_path)

    # Rich notebooks: expected-output images are materialized here.
    cleanup_expected_output_images(pset_index)
    expected_output_block_index = 1

    for cell_idx, cell in enumerate(notebook.get("cells", [])):
        clean_cell = deepcopy(cell)
        rich_cell = deepcopy(cell)

        if cell.get("cell_type") == "markdown":
            clean_cell["source"] = normalize_student_markdown(clean_cell.get("source", []))
            rich_cell["source"] = normalize_student_markdown(rich_cell.get("source", []))

            # First cell: strip page metadata/labels in clean, add downloads in rich.
            if cell_idx == 0:
                clean_src = strip_frontmatter_and_labels("".join(clean_cell.get("source", [])))
                clean_cell["source"] = clean_src.splitlines(keepends=True)

                rich_src = build_rich_notebook_header("".join(rich_cell.get("source", [])), pset_index=pset_index)
                rich_cell["source"] = rich_src.splitlines(keepends=True)

            clean_cells.append(clean_cell)
            rich_cells.append(rich_cell)
            continue

        if cell.get("cell_type") != "code":
            clean_cells.append(clean_cell)
            rich_cells.append(rich_cell)
            continue

        # Code cell: strip solutions and clear outputs in both.
        clean_cell["source"] = strip_solution_blocks(
            strip_instructor_only_blocks(clean_cell.get("source", []))
        )
        clean_cell["outputs"] = []
        clean_cell["execution_count"] = None
        clean_cells.append(clean_cell)

        rich_cell["source"] = strip_solution_blocks(
            strip_instructor_only_blocks(rich_cell.get("source", []))
        )
        rich_cell["outputs"] = []
        rich_cell["execution_count"] = None
        rich_cells.append(rich_cell)

        # Rich-only: expected output blocks (read from the instructor notebook output).
        if cell_has_tag(cell, EXPECTED_OUTPUT_TAG):
            expected_output = expected_output_cell(
                cell,
                image_dir=EXPECTED_OUTPUT_IMAGE_DIR,
                image_prefix=f"{pset_index}_expected_output_{expected_output_block_index:02d}",
            )
            if expected_output:
                rich_cells.append(expected_output)
                expected_output_block_index += 1

    clean["cells"] = clean_cells
    rich["cells"] = rich_cells

    write_notebook(clean_path, clean)
    write_notebook(rich_path, rich)

    return clean_path, rich_path


def main() -> None:
    solution_paths = sorted(INSTRUCTORS_DIR.glob("*_solution.ipynb"))
    if not solution_paths:
        print("No instructor solution notebooks found.")
        return

    for solution_path in solution_paths:
        clean_path, rich_path = generate_student_notebook(solution_path)
        print(f"Generated {clean_path.relative_to(ROOT)}")
        print(f"Generated {rich_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
