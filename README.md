## Control Theory Lab Report with PyLaTeX

This project demonstrates the use of the `pylatex` library to create a report for a control theory lab in college. The report analyzes various basic transfer functions using step, impulse response plots and pole-zero, Nyquist, Bode diagrams.

The report is written in Russian language. Adjustments may be required for other languages.

**Requirements:**

* `pylatex`
* `matplotlib`
* `sympy`
* `numpy`
* `pandas`

In addition, you should have some LaTeX distribution installed on your machine. I preferred `BasicTeX` to keep a size of the distribution small. I installed additional packages via `tlmgr` CLI in the course of work.

**PyLaTeX Features used:**

* `Section` and `Subsection`
* `Math` 
* `Figure`
* `Marker`, `Label`, and `Ref` for figure referencing

**Additional Features:**

* Optional cover page (Concatenate an existing cover page in PDF format with the rendered document)
* Custom figure captions
* Custom title styles

**Project Pipeline:**

1. [`step_0_create_figures.py`](python/step_0_create_figures.py): Performs calculations with SymPy and NumPy, generates plots using Matplotlib, and saves them as individual files
2. [`step_1_create_doc.py`](python/step_1_create_doc.py): Composes the report using PyLaTeX and exports it to a PDF named [`report.pdf`](report.pdf)

Also, project has:
* [`utils.py`](python/utils.py): utility functions and set-ups that could be reused 
* [`conclusions.csv`](misc/conclusions.csv): holds commentaries for analysis of each transfer function

**Encountered Issues and Solutions:**

* **Messed up styles with PDFLaTeX and babel:** Even using `babel` package for russian and english, `pdflatex` was rendering file badly: messed up styles, no bold, italic fonts, etc. I switched to `xelatex` renderer and it helped.
* **Misplaced text in rendered document:** Initially, multiple figures were placed using the `[h!]` option, which caused text to sticks only to the first figure. In general, text was filling empty spaces between figures. To resolve this, the `float` package was used, and the placement option was changed to `[H]` to ensure proper placement of text and figures in rendered document.

* **Broken links in PDF:** After rendering the document with `PyLaTeX`, I've noticed broken links displaying as **??** instead of numbers. This issue persisted even with manual rendering using `xelatex` via terminal. I associated this problem with auxiliary files being cleaned up by `PyLaTeX`. By rendering the document twice, first round without deleting auxiliary files and then with deletion, I successfully resolved the broken links and ensured correct link rendering in the final document.
