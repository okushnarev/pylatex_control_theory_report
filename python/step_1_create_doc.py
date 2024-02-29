from pathlib import Path

import numpy as np
import pandas as pd
from pylatex import NoEscape, Section, Subsection, Figure, Math, Label, Ref, Marker, NewPage, Command
from sympy import latex, fraction, lambdify
from sympy.abc import s, omega
from sympy.physics.control.control_plots import pole_zero_numerical_data
from sympy.physics.control.lti import TransferFunction

from step_0_create_figures import TFs, VarSubs, T_c, inertial_tf_id
from utils import init_matplotlib, init_pylatex_doc, np_arr_to_str, math_symbol

# Path to export PDF
export_path = Path('../')

# Figures folder path related to export path
figs_folder = Path('./figs')

# Cover-page options
include_cover_page = False
cover_page_path = Path('./misc/cover_page.pdf')

if __name__ == '__main__':


    init_matplotlib()

    # ### Report Data

    TFs_latex = [
        'k',
        r'\frac{k}{T_1 s + 1}',
        r'\frac{k}{T_2^2 s^2 + T_1 s + 1}',
        r'\frac{k}{T_2^2 s^2 + T_1 s + 1}',
        r'\frac{k}{s}',
        r'\frac{k}{T_2^2 s^2 + T_1 s}',
        'ks',
        r'\frac{ks}{s + 1}'
    ]

    names = [
        'Безынерционное звено (пропорциональное)',
        'Инерционное звено 1-го порядка (апериодическое)',
        'Инерционное звено 2-го порядка (апериодическое)',
        'Инерционное звено 2-го порядка (колебательное)',
        'Идеальное интегрирующее звено',
        'Реальное интегрирующее звено',
        'Идеальное дифференцирующее звено',
        'Реальное дифференцирующее звено',
    ]

    commentaries = pd.read_csv('../misc/conclusions.csv', index_col=0)
    commentaries = commentaries.fillna('')

    # ### Begin doc
    # Using pylatex

    doc = init_pylatex_doc()

    if include_cover_page:
        doc.append(Command('includepdf', NoEscape(cover_page_path), NoEscape('pages=-')))
    else:
        doc.preamble.append(Command('setcounter', ['page', 2]))

    with doc.create(Section('Цель', numbering=False)):
        doc.append('Изучить модели и характеристики основных типовых динамических звеньев систем управления.')

    doc.append(Section('Ход работы', numbering=False))

    # ### Going through IDs

    for tf_id in range(8):

        with doc.create(Section(names[tf_id], numbering=True)):

            W_s_raw = TFs[tf_id]
            W_s = W_s_raw.subs({
                'k':  VarSubs['k'][tf_id],
                'T1': VarSubs['T1'][tf_id],
                'T2': VarSubs['T2'][tf_id]
            })
            W_s2 = W_s_raw.subs({
                'k':  VarSubs['k'][tf_id] * 2,
                'T1': VarSubs['T1'][tf_id],
                'T2': VarSubs['T2'][tf_id]
            })

            if tf_id == 1:
                W_sT = W_s_raw.subs({
                    'k':  VarSubs['k'][tf_id],
                    'T1': VarSubs['T1'][tf_id] * 3,
                    'T2': VarSubs['T2'][tf_id]
                })

            doc.append('Передаточная функция:')
            doc.append(Math(data=['W(s)', '=', TFs_latex[tf_id], '=', latex(W_s), '. '], escape=False))

            # ### Pole-Zero

            tf = TransferFunction(*fraction(W_s), s)
            z, p = pole_zero_numerical_data(tf)

            doc.append('Определим нули и полюса передаточной функции. ')
            if p.size == 0 and z.size == 0:
                doc.append('Нулей и полюсов нет.')

            else:

                if z.size > 0:
                    z_math_expr = Math(data=np_arr_to_str(z), inline=True, escape=False).dumps()
                    doc.append(NoEscape(f'Нули: {z_math_expr}. '))
                else:
                    doc.append('Нулей нет. ')

                if p.size > 0:
                    p_math_expr = Math(data=np_arr_to_str(p), inline=True, escape=False).dumps()
                    doc.append(NoEscape(f'Полюса: {p_math_expr}. '))
                else:
                    doc.append('Полюсов нет. ')

                mkr = Marker(prefix='polezero', name=str(tf_id))

                doc.append('Имеющиеся нули и полюса показаны на рисунке ')
                doc.append(Ref(mkr))
                doc.append('. ')

                doc.append(NoEscape(commentaries.loc['polezero'][str(tf_id)]))

                with doc.create(Figure(position='H')) as fig:
                    fig.add_image(f'{figs_folder}/pole_zero_{tf_id}.pdf')
                    fig.add_caption(f'График нулей и полюсов. {names[tf_id]}')
                    fig.append(Label(mkr))

            # ## Time Domain

            with doc.create(Subsection('Временные характеристики', numbering=False)):
                doc.append('Рассмотрим временные характеристики звена. ')

                # ### Step resopnse

                mkr = Marker(prefix='step', name=str(tf_id))

                doc.append('График переходного процесса показан на рисунке ')
                doc.append(Ref(mkr))
                doc.append('. ')

                doc.append(NoEscape(commentaries.loc['step'][str(tf_id)]))

                with doc.create(Figure(position='H')) as fig:
                    fig.add_image(f'{figs_folder}/step_{tf_id}.pdf')
                    fig.add_caption(f'Переходная характеристика. {names[tf_id]}')
                    fig.append(Label(mkr))

                # ### Impulse response

                doc.append('График импульсной характеристики показан на рисунке ')
                mkr = Marker(prefix='impulse', name=str(tf_id))
                doc.append(Ref(mkr))
                doc.append('. ')

                doc.append(NoEscape(commentaries.loc['impulse'][str(tf_id)]))

                with doc.create(Figure(position='H')) as fig:
                    fig.add_image(f'{figs_folder}/impulse_{tf_id}.pdf')
                    fig.add_caption(f'Импульсная характеристика. {names[tf_id]}')
                    fig.append(Label(mkr))

            # ## Frequency Domain

            with doc.create(Subsection('Частотные характеристики', numbering=False)):
                doc.append('Рассмотрим частотные характеристики звена. ')

                # ### Nyquist plot

                W_jw = W_s.subs(s, 1j * omega)
                W = lambdify(omega, W_jw, 'numpy')
                w_start = 0.01
                w_end = 200
                w_range = np.arange(w_start, w_end, 0.01)
                W = W(w_range)

                mkr = Marker(prefix='nyq', name=str(tf_id))

                doc.append('График АФЧХ для исследуемого звена показан на рисунке ')
                doc.append(Ref(mkr))
                doc.append('. ')

                doc.append(NoEscape(commentaries.loc['nyq'][str(tf_id)]))

                with doc.create(Figure(position='H')) as fig:
                    fig.add_image(f'{figs_folder}/nyq_{tf_id}.pdf')
                    # I know this is ^^^^^^ a dumb way to use pathlib Path in f-string instead of Path / string
                    # But here's a problem with .add_image() in pylatex
                    # It has additional fix_filename() function that ruins that nice way (Path / string) to pass a path

                    fig.add_caption(f'АФЧХ. {names[tf_id]}')
                    fig.append(Label(mkr))

                # ### Frequency resopnse

                A = np.abs(W)

                mkr = Marker(prefix='freq', name=str(tf_id))

                doc.append('График АЧХ и ФЧХ для исследуемого звена построены на рисунке ')
                doc.append(Ref(mkr))
                doc.append('. ')

                doc.append(NoEscape(commentaries.loc['freq'][str(tf_id)]))

                with doc.create(Figure(position='H')) as fig:
                    fig.add_image(f'{figs_folder}/freq_{tf_id}.pdf')
                    fig.add_caption(f'АЧХ и ФЧХ. {names[tf_id]}')
                    fig.append(Label(mkr))

                # ### Bode plot

                L = 20 * np.log(A)

                mkr = Marker(prefix='bode', name=str(tf_id))

                doc.append('Графики ЛАЧХ и ЛФЧХ для исследуемого звена построены на рисунке ')
                doc.append(Ref(mkr))
                doc.append('. ')

                if tf_id in inertial_tf_id:

                    if tf_id == 3:
                        wc = w_range[np.where(L <= L[0] - 3)[0][0]]
                    elif tf_id == 7:
                        wc = w_range[np.where(L <= np.max(L) - 3)[0][-1]]
                    else:
                        wc = w_range[np.where(L <= np.max(L) - 3)[0][0]]

                    doc.append(f'Данное звено является инерциальным, определим его характеристики. Частота среза ')
                    doc.append(NoEscape(math_symbol([r'\omega_{cp}', '=', wc]) + ' рад/с. Частоты сопряжения: '))

                    for idx in range(len(T_c[tf_id])):
                        w_cp = 1 / T_c[tf_id][idx]
                        doc.append(NoEscape(math_symbol([f'\omega_{{c{idx + 1}}}', '=', w_cp]) + ' рад/с. '))

                doc.append(NoEscape(commentaries.loc['bode'][str(tf_id)]))

                with doc.create(Figure(position='H')) as fig:
                    fig.add_image(f'{figs_folder}/bode_{tf_id}.pdf')
                    fig.add_caption(f'ЛАЧХ и ЛФЧХ. {names[tf_id]}')
                    fig.append(Label(mkr))

    doc.append(NewPage())
    with doc.create(Section('Заключение', numbering=False)):
        doc.append('В ходе работы были изучены основные типы динамических звеньев систем управления. '
                   'Экспериментальное исследование позволило получить представление о '
                   'поведении звеньев. '
                   'При исследовании каждого звена были подчеркнуты характерные ему особенности.')

    # ### Export doc

    # double generate_pdf for correct references
    doc.generate_pdf(export_path / 'report', clean_tex=False, clean=False, compiler='xelatex', silent=True)
    doc.generate_pdf(export_path / 'report', clean_tex=False, clean=True, compiler='xelatex', silent=True)
