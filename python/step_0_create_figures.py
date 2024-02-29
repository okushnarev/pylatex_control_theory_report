import numpy as np
from matplotlib import pyplot as plt
from sympy import symbols, fraction, lambdify, DiracDelta, re, im, inverse_laplace_transform
from sympy.abc import s, t, omega, k
from sympy.physics.control.control_plots import pole_zero_numerical_data
from sympy.physics.control.lti import TransferFunction

from utils import init_matplotlib, my_plot, save_fig, check_if_one_element, \
    plot_bode_vertical_lines

from pathlib import Path

# ### Report Data

T1, T2 = symbols('T1 T2')
TFs = [
    k,
    k / (T1 * s + 1),
    k / (T2 * s ** 2 + T1 * s + 1),
    k / (T2 * s ** 2 + T1 * s + 1),
    k / s,
    k / (T2 * s ** 2 + T1 * s),
    k * s,
    (k * s) / (T1 * s + 1)
]

VarSubs = {
    'k':  [10, 10, 10, 10, 10, 10, 10, 10],
    'T1': [0, 0.1, 0.1, 0.1, 1, 1, 0, 0.1],
    'T2': [0, 0, 1.6e-3, 0.04, 0, 0.1, 0, 0]
}

T_c = [[], [0.1], [0.08, 0.02], [0.2], [], [0.1], [], [0.1]]

inertial_tf_id = [1, 2, 3, 5, 7]

figs_folder = Path('../figs')

if __name__ == '__main__':
    
    figs_folder.mkdir(exist_ok=True, parents=True)

    init_matplotlib()

    # ### Begin doc
    # Using pylatex

    # ### Going through IDs

    for tf_id in range(8):

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

        # ### Pole-Zero

        tf = TransferFunction(*fraction(W_s), s)
        z, p = pole_zero_numerical_data(tf)

        if p.size != 0 or z.size != 0:

            plt.figure(figsize=(6.4, 3.5))
            ax = my_plot(0, 0, 'Re')
            ax.set_ylabel('Im', rotation=0, loc='top', labelpad=-70)
            handles, labels = ax.get_legend_handles_labels()

            if z.size > 0:
                scat_zeros = ax.scatter([re(x) for x in z], [im(x) for x in z], facecolors='none', edgecolors='k',
                                        marker='o',
                                        s=200, linewidths=2)
                handles.append(scat_zeros)
                labels.append('Нули')

            if p.size > 0:
                scat_poles = ax.scatter([re(x) for x in p], [im(x) for x in p], c='k', marker='X', s=200)
                handles.append(scat_poles)
                labels.append('Полюса')

            ax.legend(handles, labels, loc='lower center', framealpha=1)

            fig_name = f'pole_zero_{tf_id}'
            save_fig(fig_name, filepath=figs_folder)

        # ## Time Domain

        # ### Step resopnse

        h_t = inverse_laplace_transform(W_s / s, s, t)
        h_t = h_t.subs(DiracDelta(t), 0)
        h_t = lambdify(t, h_t, 'numpy')

        h_t2 = inverse_laplace_transform(W_s2 / s, s, t)
        h_t2 = h_t2.subs(DiracDelta(t), 0)
        h_t2 = lambdify(t, h_t2, 'numpy')

        if tf_id == 1:
            h_t3 = inverse_laplace_transform(W_sT / s, s, t)
            h_t3 = h_t3.subs(DiracDelta(t), 0)
            h_t3 = lambdify(t, h_t3, 'numpy')

        t_stop = 1
        t_st = 0.1
        while np.abs(h_t(t_stop + t_st) - h_t(t_stop - t_st)) > 0.001:
            t_stop += 1
            if t_stop >= 10:
                t_stop = 1
                break
        t_range = np.arange(0, t_stop, 0.01)
        plt.figure(figsize=(6.4, 4))
        ax = my_plot(0, 0, 'Время, с', 'Амплитуда')
        line1, = ax.plot(t_range, check_if_one_element(h_t(t_range), t_range), linewidth=3, color='k')
        line2, = ax.plot(t_range, check_if_one_element(h_t2(t_range), t_range), linestyle='--', linewidth=3,
                         color='k')
        if tf_id == 1:
            line3, = ax.plot(t_range, check_if_one_element(h_t3(t_range), t_range), linestyle=':', linewidth=3,
                             color='k')
            plt.legend([line1, line2, line3],
                       [
                           f'k = {VarSubs["k"][tf_id]}, T = {VarSubs["T1"][tf_id]}',
                           f'k = {VarSubs["k"][tf_id] * 2}, T = {VarSubs["T1"][tf_id]}',
                           f'k = {VarSubs["k"][tf_id]}, T = {VarSubs["T1"][tf_id] * 3:.1f}'
                       ],
                       loc='upper right',
                       framealpha=1)
        else:
            plt.legend([line1, line2],
                       [f'k = {VarSubs["k"][tf_id]}', f'k = {VarSubs["k"][tf_id] * 2}'],
                       loc='upper right',
                       framealpha=1)

        save_fig(f'step_{tf_id}', filepath=figs_folder)

        # ### Impulse response

        g_t = inverse_laplace_transform(W_s, s, t)
        g_t = g_t.subs(DiracDelta(t), 0)
        g_t = lambdify(t, g_t, 'numpy')

        plt.figure(figsize=(6.4, 4))
        my_plot(t_range, check_if_one_element(g_t(t_range), t_range), 'Время, с', 'Амплитуда')

        save_fig(f'impulse_{tf_id}', filepath=figs_folder)

        # ## Frequency Domain

        # ### Nyquist plot

        W_jw = W_s.subs(s, 1j * omega)
        W = lambdify(omega, W_jw, 'numpy')
        w_start = 0.01
        w_end = 200
        w_range = np.arange(w_start, w_end, 0.01)
        W = W(w_range)

        plt.figure(figsize=(6.4, 4))
        ax = my_plot(np.real(W), np.imag(W), 'Re', 'Im')
        if type(W) is int:
            ax.scatter(np.real(W), np.imag(W), color='k', linewidths=5)

        if tf_id == 4 or tf_id == 5:
            ax.set_ylim(-15, 0.01)
        elif tf_id == 6:
            ax.set_ylim(0.01, 15)

        save_fig(f'nyq_{tf_id}', filepath=figs_folder)

        # ### Frequency resopnse

        A = np.abs(W)
        phi = np.angle(W)

        plt.figure(figsize=(6.4, 4))
        plt.subplot(211)
        ax = my_plot(w_range, check_if_one_element(A, w_range), 'Частота, рад/с', 'Амплитуда')
        if tf_id == 4 or tf_id == 5:
            ax.set_ylim(0, 10)

        plt.subplot(212)
        ax = my_plot(w_range, check_if_one_element(np.rad2deg(phi), w_range), '', 'Фаза, град')
        ax.set_yticks([-180, -90, 0, 90])

        save_fig(f'freq_{tf_id}', filepath=figs_folder)

        # ### Bode plot

        L = 20 * np.log(A)

        fig = plt.figure(figsize=(6.4, 4))

        plt.subplot(211)
        ax = my_plot(w_range, check_if_one_element(L, w_range), '', 'Амплитуда, дБ', xlog=True)
        ax.set_xscale('log')
        ax.grid(True, which='minor', ls='--', linewidth=0.5, color='k', alpha=0.15)
        ax.set_xlim(w_start, w_end)

        if tf_id in inertial_tf_id:
            if tf_id == 3:
                wc = w_range[np.where(L <= L[0] - 3)[0][0]]
            elif tf_id == 7:
                wc = w_range[np.where(L <= np.max(L) - 3)[0][-1]]
            else:
                wc = w_range[np.where(L <= np.max(L) - 3)[0][0]]

            plot_bode_vertical_lines(ax, tf_id, T_c, wc)

        plt.subplot(212)
        ax = my_plot(w_range, check_if_one_element(np.rad2deg(phi), w_range), 'Частота, рад/с', 'Фаза, град',
                     xlog=True)
        ax.set_xscale('log')
        ax.grid(True, which='minor', ls='--', linewidth=0.5, color='k', alpha=0.15)
        ax.set_yticks(np.linspace(90, -180, 4))
        ax.set_xlim(w_start, w_end)

        if tf_id in inertial_tf_id:
            plot_bode_vertical_lines(ax, tf_id, T_c, wc, top=True)

        save_fig(f'bode_{tf_id}', filepath=figs_folder)
