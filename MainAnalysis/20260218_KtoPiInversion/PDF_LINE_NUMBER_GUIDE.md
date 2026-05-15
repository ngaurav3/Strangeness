# PDF Line Number Guide

Date: 2026-03-13

## Scope

This guide defines how agents should interpret requests like:

- "What is at line 102?"
- "Please check L102."
- "Look at line 250 in the note."

for the two main documents in this repository:

- analysis note PDF: `overleaf_repo_git/main.pdf`
- paper draft PDF: `paper/overleaf_paper_draft/main.pdf`

## Default Rule

Unless the user explicitly says "source line", ".tex line", or names a `.tex`
file, line numbers must be interpreted as line numbers in the compiled PDF, not
line numbers in the LaTeX source.

Examples:

- "line 102 in the analysis note" means PDF line `102` in `overleaf_repo_git/main.pdf`
- "L102 in the paper" means PDF line `102` in `paper/overleaf_paper_draft/main.pdf`
- "main.tex:102" means source line `102` in the LaTeX file

## How To Read a PDF Line Number

Use `pdftotext -layout` so the left-margin line numbers are preserved.

Analysis note example:

```bash
pdftotext -layout overleaf_repo_git/main.pdf - | rg -n '^ *102 '
```

Paper example:

```bash
pdftotext -layout paper/overleaf_paper_draft/main.pdf - | rg -n '^ *102 '
```

If you need the surrounding context:

```bash
pdftotext -layout overleaf_repo_git/main.pdf - | sed -n '132,140p'
```

If you already know the page and want a page-local extraction:

```bash
pdftotext -f 4 -l 4 -layout overleaf_repo_git/main.pdf -
pdftotext -f 4 -l 4 -layout paper/overleaf_paper_draft/main.pdf -
```

## How To Answer

When the user asks for a line number in the note or paper:

1. Quote the text from the compiled PDF line number.
2. State the page number when possible.
3. If helpful, note whether the sentence continues onto the next PDF line.
4. Do not silently switch to source-line numbering.

Good answer pattern:

- "In the compiled analysis-note PDF, line 102 on page 4 reads: ..."

## Important Distinction

PDF line numbers and source line numbers are different and must not be mixed.

- PDF line number: what the user sees in the rendered PDF margin
- source line number: what appears in `main.tex` or another source file

If there is any ambiguity, ask whether the user means:

- compiled PDF line number
- LaTeX source line number

But if the user only says "line 102" for the analysis note or paper, default to
the compiled PDF line number.
