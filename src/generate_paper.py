"""Regenerates the downloadable PDF and Word (.docx) copies of
``docs/full_methodology_paper.md`` via `pandoc <https://pandoc.org>`_ (a
system dependency, not a pip package -- must already be installed and on
PATH). PDF generation additionally requires a LaTeX engine; this targets
``xelatex`` for its better Unicode support (em dashes, arrows) over the
default ``pdflatex``.

The app reads the generated PDF's bytes directly (`app/streamlit_app.py`'s
`_load_file_bytes`) for both the download button and the "Open PDF in New
Tab" button (which decodes them into a Blob client-side via injected
JavaScript) -- no server route or mirrored copy needed.

Run with: python src/generate_paper.py
"""
from __future__ import annotations

import shutil
import subprocess
import sys
import os

DOCS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs")
SOURCE_MD = os.path.join(DOCS_DIR, "full_methodology_paper.md")
OUTPUT_PDF = os.path.join(DOCS_DIR, "full_methodology_paper.pdf")
OUTPUT_DOCX = os.path.join(DOCS_DIR, "full_methodology_paper.docx")


def _run_pandoc(args: list[str]) -> None:
    subprocess.run(["pandoc", *args], check=True)


def generate() -> None:
    if shutil.which("pandoc") is None:
        print(
            "pandoc not found on PATH -- install it from https://pandoc.org/installing.html "
            "to regenerate the PDF/Word versions of the paper. The markdown source at "
            f"{SOURCE_MD} is unaffected and is what the app's in-app summary view uses.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Generating {OUTPUT_PDF} ...")
    _run_pandoc(
        [
            SOURCE_MD,
            "-o",
            OUTPUT_PDF,
            "--pdf-engine=xelatex",
            "--toc",
            "--toc-depth=2",
            "-V",
            "geometry:margin=1in",
            "-V",
            "colorlinks=true",
            "-V",
            "mainfont=Georgia",
        ]
    )

    print(f"Generating {OUTPUT_DOCX} ...")
    _run_pandoc([SOURCE_MD, "-o", OUTPUT_DOCX, "--toc", "--toc-depth=2"])

    print("Done.")


if __name__ == "__main__":
    generate()
