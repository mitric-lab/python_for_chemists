#!/usr/bin/env python3
"""Generate student notebooks from instructor solution notebooks."""

from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
INSTRUCTORS_DIR = ROOT / "instructors"


def strip_solution_blocks(source: list[str]) -> list[str]:
    """Replace solution regions with a placeholder and blank lines."""
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


def target_path_for(solution_path: Path) -> Path:
    """Map an instructor solution notebook to its student notebook path."""
    if solution_path.name.endswith("_solution.ipynb"):
        return ROOT / "chapters/problem_sets" / solution_path.name.replace(
            "_solution.ipynb", ".ipynb"
        )
    raise ValueError(f"Unsupported solution notebook name: {solution_path.name}")


def normalize_student_markdown(source: list[str]) -> list[str]:
    """Adjust solution-specific labels and titles for the student notebook."""
    normalized: list[str] = []
    for line in source:
        line = line.replace("_solution)=", ")=")
        if line.startswith("# Solution of "):
            line = line.replace("# Solution of ", "# ", 1)
        normalized.append(line)
    return normalized


def generate_student_notebook(solution_path: Path) -> Path:
    """Generate one student notebook from an instructor notebook."""
    notebook = json.loads(solution_path.read_text(encoding="utf-8"))
    student = deepcopy(notebook)

    for cell in student.get("cells", []):
        if cell.get("cell_type") == "markdown":
            cell["source"] = normalize_student_markdown(cell.get("source", []))
            continue
        if cell.get("cell_type") != "code":
            continue
        cell["source"] = strip_solution_blocks(cell.get("source", []))
        cell["outputs"] = []
        cell["execution_count"] = None

    target_path = target_path_for(solution_path)
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target_path.write_text(
        json.dumps(student, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8",
    )
    return target_path


def main() -> None:
    solution_paths = sorted(INSTRUCTORS_DIR.glob("*_solution.ipynb"))
    if not solution_paths:
        print("No instructor solution notebooks found.")
        return

    for solution_path in solution_paths:
        target_path = generate_student_notebook(solution_path)
        print(f"Generated {target_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
