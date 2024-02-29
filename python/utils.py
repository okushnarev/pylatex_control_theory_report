import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from pylatex import Command, Document, NoEscape, Package, Math
from sympy import Eq, latex, sympify


# Plots
def init_matplotlib():
    matplotlib.rcParams.update({
        'pdf.fonttype': 42,
    })

    plt.rcParams.update({
        "font.family":           "Serif",
        "font.serif":            "Times New Roman",
        "font.size":             14,

        "axes.titlepad":         20,
        "axes.titlesize":        24,
        "axes.titleweight":      500,

        "axes.labelpad":         10,
        "axes.labelsize":        14,

        "xtick.major.pad":       5,
        "ytick.major.pad":       15,

        "xtick.labelsize":       14,
        "ytick.labelsize":       14,

        "figure.subplot.top":    0.8,
        "figure.subplot.bottom": 0.14,
        "figure.titlesize":      14,
    })


def my_plot(x, y, xlabel="", ylabel="", xlog=False, **kwargs):
    plt.plot(x, y, linewidth=3, color='k', **kwargs)
    ax = plt.gca()
    ax.set_axisbelow(True)
    ax.grid(color='k', alpha=0.5, linewidth=0.5)

    ax.set_xlabel(xlabel, loc='right', labelpad=-50)
    ax.set_ylabel(ylabel)

    if not xlog:
        # set the x-spine (see below for more info on `set_position`)
        ax.spines['left'].set_position('zero')

        # Move the left spine to x = 0 and y = 0, respectively.
        ax.spines["left"].set_position(("data", 0))

        # arrow tip on Y axis
        ax.plot(0, 1, "^k", transform=ax.get_xaxis_transform(), clip_on=False)

    # arrow tip on X axis
    ax.plot(1, 0, ">k", transform=ax.get_yaxis_transform(), clip_on=False)

    # Move the bottom spine to x = 0 and y = 0, respectively.
    ax.spines["bottom"].set_position(("data", 0))

    # Hide the top and right spines.
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    return ax


def plot_bode_vertical_lines(ax, exp_id, T_c_arr, wc, top=False):
    if top:
        text_pos = ax.get_ylim()[1] * 0.5
    else:
        text_pos = ax.get_ylim()[0] * 0.9

    for idx in range(len(T_c_arr[exp_id])):

        ax.axvline(x=1 / T_c_arr[exp_id][idx], color='k', linestyle='--', linewidth=1.5)
        ax.text(x=1 / T_c_arr[exp_id][idx] * 0.7, y=text_pos, s=f'$\\omega_{{c{idx + 1}}}$', rotation=90,
                usetex=True)

    ax.axvline(x=wc, color='k', linestyle=':', linewidth=2)
    ax.text(x=wc * 0.7, y=text_pos, s=r'$\omega_{cp}$', rotation=90, usetex=True)


def save_fig(filename, filepath='./figs'):
    plt.tight_layout()
    plt.savefig(f'{filepath}/{filename}.pdf')
    plt.close()


# Latex
def init_pylatex_doc():
    doc = Document('basic',
                   page_numbers=True,
                   documentclass='extarticle',
                   document_options=['14pt', 'a4paper', ],
                   geometry_options={'left':    '3cm',
                                     'right':   '1.5cm',
                                     'vmargin': '2cm'},
                   )

    # Add packages
    doc.packages.add(Package('float'))
    doc.packages.add(Package('indentfirst'))
    doc.packages.add(Package('titlesec'))
    doc.packages.add(Package('subcaption'))
    doc.packages.add(Package('babel', options='english, russian'))
    doc.packages.add(Package('pdfpages', 'final'))

    # Set font to Times
    doc.preamble.append(Command('babelfont', ['rm', 'Times']))

    # Custom Style for Section, SubSection, etc.
    for element in ['section', 'subsection', 'subsubsection', 'paragraph', 'subparagraph']:
        doc.preamble.append(NoEscape('\\titleformat*{\\' + element + r'}{\fontsize{14}{16}\bfseries}'))
        doc.preamble.append(NoEscape('\\titlespacing*{\\' + element + r'}{\parindent}{7pt}{7pt}'))

    # Set line parameters
    doc.preamble.append(Command('linespread', '1.5'))
    doc.preamble.append(Command('setlength', [Command('parindent'), '1.25cm']))

    # Set Russian Caption names
    doc.preamble.append(Command('addto\captionsrussian', Command('renewcommand',
                                                                 [Command('figurename'), 'Рисунок'])))

    # Set picture caption parameters
    doc.preamble.append(Command('DeclareCaptionFormat', ['custom', NoEscape('#1 #2 #3')]))
    doc.preamble.append(Command('DeclareCaptionLabelSeparator', ['custom', '–']))
    doc.preamble.append(Command('captionsetup', 'format=custom, labelsep=custom, justification=centering'))

    return doc


def to_latex(left, right):
    return latex(Eq(sympify(left), right))


def math_symbol(symbol):
    return Math(data=symbol, escape=False, inline=True).dumps()


# Data Operation
def check_if_one_element(data, d_range):
    try:
        len(data)
    except TypeError:
        data = [data for _ in range(len(d_range))]
    return data


def np_arr_to_str(arr):
    latex_expr = []
    for elem in arr:
        latex_expr.append(latex(sympify(np.array2string(elem, precision=2, suppress_small=True))))
    return ', '.join(latex_expr)
