import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def bfz_plot(b_at_det, freqs, det_pos):
    n_det, n_f, n_dim = np.shape(b_at_det)
    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    spec = gridspec.GridSpec(ncols=12, nrows=2, figure=fig)

    ax1 = fig.add_subplot(spec[0, 0:4])
    ax2 = fig.add_subplot(spec[0, 4:8])
    ax3 = fig.add_subplot(spec[0, 8:12])
    ax4 = fig.add_subplot(spec[1, 0:3])
    ax5 = fig.add_subplot(spec[1, 3:6])
    ax6 = fig.add_subplot(spec[1, 6:9])
    ax7 = fig.add_subplot(spec[1, 9:12])
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

    for ax in axes:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.set_xlabel(r'x[m]')
        ax.grid()
    ax1.set_ylabel(r'$|B_x|/|\mu_0 I_W/ 2\pi|$')
    ax2.set_ylabel(r'$|B_y|/|\mu_0 I_W/ 2\pi|$')
    ax3.set_ylabel(r'$|B_z|/|\mu_0 I_W/ 2\pi|$')
    ax4.set_ylabel(r'$|B_z/B_x|$')
    ax5.set_ylabel(r'$arg(B_x)/\pi$')
    ax6.set_ylabel(r'$arg(B_y)/\pi$')
    ax7.set_ylabel(r'$arg(B_z)/\pi$')

    for i in range(n_f):
        ax1.plot(det_pos[:, 0], abs(b_at_det[:, i, 0]))
        ax2.plot(det_pos[:, 0], abs(b_at_det[:, i, 1]))
        ax3.plot(det_pos[:, 0], abs(b_at_det[:, i, 2]), label='f=' + str(freqs[i]) + "Hz")

        ax4.plot(det_pos[:, 0], abs(b_at_det[:, i, 2]) / abs(b_at_det[:, i, 0]))
        ax5.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 0])), linestyle='-',
                 label='Bx, f=' + str(freqs[i]) + "Hz")
        ax6.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 1])), linestyle='--',
                 label='By, f=' + str(freqs[i]) + "Hz")
        ax7.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 2])), linestyle=':',
                 label='Bz, f=' + str(freqs[i]) + "Hz")
        # print(1 / np.pi * (np.angle(b_at_det[:, i, 0])))
        # exit()

        # ax5.plot(det_pos[:, 0], 1/np.pi*(np.angle(b_at_det[:, i, 2])-np.angle(b_at_det[:, i, 0])))
    ax5.set_ylim(-1, 1.1)
    ax6.set_ylim(-1, 1.1)
    ax7.set_ylim(-1, 1.1)
    # ax1.plot(det_pos[:, 1], 1/det_pos[:, 1])
    ax3.legend()
    # ax4.legend()
    # ax5.legend()
    # ax6.legend()
    # ax7.legend()


# def bfbz_plot(beta_vec, freqs, det_pos):
def bfbz_plot(b_at_det, beta_vec, freqs, det_pos):
    n_det, n_f, n_beta, n_dim = np.shape(b_at_det)
    print(np.shape(b_at_det))
    fig = plt.figure(constrained_layout=True, figsize=(8, 12))
    spec = gridspec.GridSpec(ncols=12, nrows=4, figure=fig)

    ax1 = fig.add_subplot(spec[0, 0:3])
    ax2 = fig.add_subplot(spec[0, 3:6])
    ax3 = fig.add_subplot(spec[0, 6:9])
    ax4 = fig.add_subplot(spec[0, 9:12])

    ax5 = fig.add_subplot(spec[1, 0:3])
    ax6 = fig.add_subplot(spec[1, 3:6])
    ax7 = fig.add_subplot(spec[1, 6:9])
    ax8 = fig.add_subplot(spec[1, 9:12])

    ax9 = fig.add_subplot(spec[2, 0:3])
    ax10 = fig.add_subplot(spec[2, 3:6])
    ax11 = fig.add_subplot(spec[2, 6:9])
    ax12 = fig.add_subplot(spec[2, 9:12])

    ax13 = fig.add_subplot(spec[3, 0:3])
    ax14 = fig.add_subplot(spec[3, 3:6])
    ax15 = fig.add_subplot(spec[3, 6:9])
    ax16 = fig.add_subplot(spec[3, 9:12])

    axes1 = [ax1, ax2, ax3, ax4]
    axes2 = [ax5, ax6, ax7, ax8]
    axes3 = [ax9, ax10, ax11, ax12]
    axes4 = [ax13, ax14, ax15, ax16]
    axes_list = [axes1, axes2, axes3, axes4]
    for axes in axes_list:
        for ax in axes:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                                useMathText=None)
            ax.grid()
    alphas = [r'$\pi/32$', r'$\pi/16$', r'$\pi/8$', r'$\pi/4$']
    for ax in axes1:
        ax_ind = axes1.index(ax)
        ax.title.set_text(r'$\alpha=$' + alphas[ax_ind])

    for ax in axes3:
        ax.set_ylim(-1, 1.1)
    for ax in axes4:
        ax.set_xlabel(r'x[m]')
        ax.set_ylim(-1, 1.1)

    ax1.set_ylabel(r'$|B_x|/|\mu_0 I_W/ 2\pi|$')
    # ax5.set_ylabel(r'$|B_y|/|\mu_0 I_W/ 2\pi|$')
    ax5.set_ylabel(r'$|B_z|/|\mu_0 I_W/ 2\pi|$')
    # ax9.set_ylabel(r'$|B_z/B_x|$')
    # ax5.set_ylabel(r'$arg(B_x)/\pi$')
    ax9.set_ylabel(r'$arg(B_x)/\pi$')
    ax13.set_ylabel(r'$arg(B_z)/\pi$')

    for i in range(n_f):
        ax1.plot(det_pos[:, 0], abs(b_at_det[:, i, 0, 0]))
        ax2.plot(det_pos[:, 0], abs(b_at_det[:, i, 1, 0]))
        ax3.plot(det_pos[:, 0], abs(b_at_det[:, i, 2, 0]))
        ax4.plot(det_pos[:, 0], abs(b_at_det[:, i, 3, 0]),
                 label='f=' + str(freqs[i]) + "Hz")

        ax5.plot(det_pos[:, 0], abs(b_at_det[:, i, 0, 2]))
        ax6.plot(det_pos[:, 0], abs(b_at_det[:, i, 1, 2]))
        ax7.plot(det_pos[:, 0], abs(b_at_det[:, i, 2, 2]))
        ax8.plot(det_pos[:, 0], abs(b_at_det[:, i, 3, 2]))

        ax9.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 0, 0])))  # , linestyle=':',
        # label='Bz, f=' + str(freqs[i]) + "Hz")
        ax10.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 1, 0])))
        ax11.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 2, 0])))
        ax12.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 3, 0])))

        ax13.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 0, 2])))
        ax14.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 1, 2])))
        ax15.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 2, 2])))
        ax16.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 3, 2])))

    #    ax2.plot(det_pos[:, 0], abs(b_at_det[:, i, 1]))
    #    ax3.plot(det_pos[:, 0], abs(b_at_det[:, i, 2]), label='f=' + str(freqs[i]) + "Hz")
    #
    #    ax4.plot(det_pos[:, 0], abs(b_at_det[:, i, 2]) / abs(b_at_det[:, i, 0]))
    #    ax5.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 0])), linestyle='-',
    #             label='Bx, f=' + str(freqs[i]) + "Hz")
    #    ax6.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 1])), linestyle='--',
    #             label='By, f=' + str(freqs[i]) + "Hz")
    #    ax7.plot(det_pos[:, 0], 1 / np.pi * (np.angle(b_at_det[:, i, 2])), linestyle=':',
    #             label='Bz, f=' + str(freqs[i]) + "Hz")
    #    # print(1 / np.pi * (np.angle(b_at_det[:, i, 0])))
    #    # exit()
    #
    #    # ax5.plot(det_pos[:, 0], 1/np.pi*(np.angle(b_at_det[:, i, 2])-np.angle(b_at_det[:, i, 0])))
    # ax5.set_ylim(-1, 1.1)
    # ax6.set_ylim(-1, 1.1)
    # ax7.set_ylim(-1, 1.1)
    ## ax1.plot(det_pos[:, 1], 1/det_pos[:, 1])
    # ax3.legend()
    # ax4.legend()
    # ax5.legend()
    # ax6.legend()
    ax4.legend()


def bfz_plot_measure(b_at_det, freqs, det_pos, row_ind):
    n_det, n_f, n_dim = np.shape(b_at_det)
    fig = plt.figure(constrained_layout=True, figsize=(8, 8))
    spec = gridspec.GridSpec(ncols=12, nrows=2, figure=fig)

    ax1 = fig.add_subplot(spec[0, 0:4])
    ax2 = fig.add_subplot(spec[0, 4:8])
    ax3 = fig.add_subplot(spec[0, 8:12])
    ax4 = fig.add_subplot(spec[1, 0:3])
    ax5 = fig.add_subplot(spec[1, 3:6])
    ax6 = fig.add_subplot(spec[1, 6:9])
    ax7 = fig.add_subplot(spec[1, 9:12])
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]

    for ax in axes:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.set_xlabel(r'x[m]')
        ax.grid()
    ax1.set_ylabel(r'$|B_x|/|\mu_0 I_W/ 2\pi|$')
    ax2.set_ylabel(r'$|B_y|/|\mu_0 I_W/ 2\pi|$')
    ax3.set_ylabel(r'$|B_z|/|\mu_0 I_W/ 2\pi|$')
    ax4.set_ylabel(r'$|B_z/B_x|$')
    ax5.set_ylabel(r'$arg(B_x)/\pi$')
    ax6.set_ylabel(r'$arg(B_y)/\pi$')
    ax7.set_ylabel(r'$arg(B_z)/\pi$')

    for i in range(n_f):
        ax1.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 0]))
        ax2.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 1]))
        ax3.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 2]), label='f=' + str(freqs[i]) + "Hz")

        ax4.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 2]) / abs(b_at_det[:, i, 0]))
        ax5.plot(det_pos[:, row_ind], 1 / np.pi * (np.angle(b_at_det[:, i, 0])), linestyle='-',
                 label='Bx, f=' + str(freqs[i]) + "Hz")
        ax6.plot(det_pos[:, row_ind], 1 / np.pi * (np.angle(b_at_det[:, i, 1])), linestyle='--',
                 label='By, f=' + str(freqs[i]) + "Hz")
        ax7.plot(det_pos[:, row_ind], 1 / np.pi * (np.angle(b_at_det[:, i, 2])), linestyle=':',
                 label='Bz, f=' + str(freqs[i]) + "Hz")
        # print(1 / np.pi * (np.angle(b_at_det[:, i, 0])))
        # exit()

        # ax5.plot(det_pos[:, 0], 1/np.pi*(np.angle(b_at_det[:, i, 2])-np.angle(b_at_det[:, i, 0])))
    ax5.set_ylim(-1, 1.1)
    ax6.set_ylim(-1, 1.1)
    ax7.set_ylim(-1, 1.1)
    # ax1.plot(det_pos[:, 1], 1/det_pos[:, 1])
    ax3.legend()
    # ax4.legend()
    # ax5.legend()
    # ax6.legend()
    # ax7.legend()


def bfz_plot_measure_rel_x(b_at_det_pack, freqs, det_pos, alpha_ind, exp_config):
    row_ind = 1
    mu0 = 4 * np.pi * 10 ** (-7)
    b_at_det, b_at_det_ref = b_at_det_pack
    n_det, n_f, n_dim = np.shape(b_at_det)
    fig = plt.figure(constrained_layout=True, figsize=(8, 10))
    alphas = np.linspace(0, 1, 13)
    fig.suptitle(r'$\alpha/\pi=$' + str(np.round(alphas[alpha_ind],decimals=3)), fontsize=16)

    spec = gridspec.GridSpec(ncols=12, nrows=10, figure=fig)

    ax1 = fig.add_subplot(spec[1:4, 0:4])
    ax2 = fig.add_subplot(spec[1:4, 4:8])
    ax3 = fig.add_subplot(spec[1:4, 8:12])

    ax1s = fig.add_subplot(spec[4:7, 0:4])
    ax2s = fig.add_subplot(spec[4:7, 4:8])
    ax3s = fig.add_subplot(spec[4:7, 8:12])

    ax4 = fig.add_subplot(spec[7:10, 0:6])
    ax5 = fig.add_subplot(spec[7:10, 6:12])

    axes_r1 = [ax1, ax2, ax3]
    axes_r2 = [ax1s, ax2s, ax3s]
    axes_r3 = [ax4, ax5]

    for ax in axes_r1:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                            useMathText=None)
        ax.grid()

    for ax in axes_r2:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                            useMathText=None)
        ax.set_xlabel(r'$y[m]$')
        ax.grid()

    for ax in axes_r3:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                            useMathText=None)
        ax.set_xlabel(r'$y[m]$')
        ax.grid()

    ax1.set_ylabel(r'$(|B_x|-|B_x^0|)/(\mu_0 I/(2\pi))$')
    ax2.set_ylabel(r'$(|B_y|-|B_y^0|)/(\mu_0 I/(2\pi))$')
    ax3.set_ylabel(r'$(|B_z|-|B_z^0|)/(\mu_0 I/(2\pi))$')
    ax1s.set_ylabel(r'$|B_x|/(\mu_0 I/(2\pi))$')
    ax2s.set_ylabel(r'$|B_y|/(\mu_0 I/(2\pi))$')
    ax3s.set_ylabel(r'$|B_z|/(\mu_0 I/(2\pi))$')

    ax4.set_ylabel(r'$|arg(B_x/B_y)|/\pi$')
    ax5.set_ylabel(r'$|arg(B_z/B_y)|/\pi$')

    for i in range(n_f):
        # ax1.plot(det_pos[:, row_ind],
        #         (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])) / abs(b_at_det_ref[:, i, 0]))
        # ax2.plot(det_pos[:, row_ind],
        #         (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])) / abs(b_at_det_ref[:, i, 1]))
        # ax3.plot(det_pos[:, row_ind],
        #         (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])) / abs(b_at_det_ref[:, i, 2]),
        #         label='f=' + str(freqs[i]) + "Hz")

        ax1.plot(det_pos[:, row_ind],
                 (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])) / (mu0 / (2 * np.pi)))
        ax2.plot(det_pos[:, row_ind],
                 (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])) / (mu0 / (2 * np.pi)))
        ax3.plot(det_pos[:, row_ind],
                 (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])) / (mu0 / (2 * np.pi)),
                 label='f=' + str(freqs[i]) + "Hz")

        ax1s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 0]) / (mu0 / (2 * np.pi)))
        ax2s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 1]) / (mu0 / (2 * np.pi)))
        ax3s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 2]) / (mu0 / (2 * np.pi)))
        #ax3s.plot(det_pos[:, row_ind], (abs(b_at_det_ref[:, i, 2]))/ (mu0 / (2 * np.pi))+100)


        ax4.plot(det_pos[:, row_ind], abs(1 / np.pi * (np.angle(b_at_det[:, i, 0] / b_at_det[:, i, 1]))), linestyle='-')
        ax5.plot(det_pos[:, row_ind], abs(1 / np.pi * (np.angle(b_at_det[:, i, 2] / b_at_det[:, i, 1]))), linestyle='-',
                 label='f=' + str(freqs[i]) + "Hz")

    ax4.set_ylim(-0.1, 1.1)
    ax5.set_ylim(-0.1, 1.1)

    ax5.legend(loc='lower right')
    exp_number = exp_config['experiment #']
    fig.savefig("./save_plots/eddy_y_" + str(exp_number) + ".pdf", bbox_inches='tight')


def bfz_plot_measure_active(b_at_det_pack, freqs, exp_config):
    s_ind = exp_config['sensor index']
    mu0 = 4 * np.pi * 10 ** (-7)
    b_at_det, z_vec = b_at_det_pack

    print(np.shape(b_at_det))

    n_z, n_det, n_f, n_dim = np.shape(b_at_det)
    det_pos = exp_config['sensor_pos'][exp_config['sensor index']]
    # lydis mon
    fig = plt.figure(constrained_layout=False, figsize=(14, 4))
    # fig = plt.figure(constrained_layout=True)
    plt.subplots_adjust(top=0.88, bottom=0.15, right=0.99, left=0.06,
                        hspace=0, wspace=0.2)
    spec = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)

    # fig = plt.figure(constrained_layout=False, figsize=(6, 2))
    ## fig = plt.figure(constrained_layout=True)
    # plt.subplots_adjust(top=0.88, bottom=0.15, right=0.98, left=0.06,
    #                    hspace=0, wspace=0.2)
    # spec = gridspec.GridSpec(ncols=3, nrows=1, figure=fig)

    ax1 = fig.add_subplot(spec[0, 0])
    ax2 = fig.add_subplot(spec[0, 1])
    ax3 = fig.add_subplot(spec[0, 2])
    axes_r1 = [ax1, ax2, ax3]

    for ax in axes_r1:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                            useMathText=None)
        ax.set_xlabel(r'$z_{obj}[m]$', fontsize=16)
        ax.grid()
    if exp_config['data type'] == 'AC':
        ax1.set_ylabel(r'$(|B_x|-|B_x^0|)/(\mu_0/(2\pi))$', fontsize=16)
        ax2.set_ylabel(r'$(|B_y|-|B_y^0|)/(\mu_0/(2\pi))$', fontsize=16)
        ax3.set_ylabel(r'$(|B_z|-|B_z^0|)/(\mu_0/(2\pi))$', fontsize=16)
    else:
        ax1.set_ylabel(r'$|B_x-B_x^0|/(\mu_0/(2\pi))$', fontsize=16)
        ax2.set_ylabel(r'$|B_y-B_y^0|/(\mu_0/(2\pi))$', fontsize=16)
        ax3.set_ylabel(r'$|B_z-B_z^0|/(\mu_0/(2\pi))$', fontsize=16)

    # ax1.text(0.95, 0.8, r'$x_{det}= $' + str(np.round(det_pos[0], decimals=2)) + '$m$, ' +
    #         '$y_{det}= $' + str(np.round(det_pos[1], decimals=2)) + '$m$, ' +
    #         r'$z_{det}= $' + str(np.round(det_pos[2], decimals=2)) + '$m$' ,
    #         verticalalignment='bottom', horizontalalignment='right',
    #         transform=ax1.transAxes, fontsize=15)

    # EXP 1 EXP 1EXP 1EXP 1EXP 1EXP 1EXP 1EXP 1
    title_text_p1 = ''
    if exp_config['experiment #'] == 1:

        title_text_p1 = r'Cube : $x_{obj}= $' + str(0.0) + '$m$, ' + '$y_{obj}= $' + str(0.16) + '$m$;    '


    # EXP 2 EXP 2 EXP 2 EXP 2 EXP 2 EXP 2 EXP 2 EXP 2 EXP 2
    elif exp_config['experiment #'] == 2:
        title_text_p1 = r'GPS : $x_{obj}= $' + str(-0.15) + '$m$, ' + '$y_{obj}=\pm$' + str(0.55) + '$m$;    '

    # EXP 3 EXP 3 EXP 3 EXP 3 EXP 3 EXP 3 EXP 3
    elif exp_config['experiment #'] == 3:
        title_text_p1 = r'Akku : $y_{obj}= $' + str(0.13) + '$m$, ' + '$z_{obj}= $' + str(0.134) + '$m$;    '

        for ax in axes_r1:
            ax.set_xlabel(r'$x_{obj}[m]$', fontsize=16)

    # EXP 4 EXP 4 EXP 4 EXP 4 EXP 4 EXP 4 EXP 4
    elif exp_config['experiment #'] == 4:
        title_text_p1 = r'Laptop : $y_{obj}= $' + str(0.13) + '$m$, ' + '$z_{obj}= $' + str(0.05) + '$m$;    '

        for ax in axes_r1:
            ax.set_xlabel(r'$x_{obj}[m]$', fontsize=16)

    # EXP 5 EXP 5 EXP 5 EXP 5
    elif exp_config['experiment #'] == 5:
        # title_text_p1 = r'GPS + Eddy: $x_{obj}= $' + str(-0.15) + '$m$, ' + '$y_{obj}=\pm$' + str(0.55) + '$m$;    '
        title_text_p1 = r'GPS + Eddy: $x_{obj}= $' + '?' + '$m$, ' + '$y_{obj}=\pm$' + '?' + '$m$;    '

    # EXP 6 EXP 6 EXP 6 EXP 6 EXP 6 EXP 6 EXP 6 EXP 6
    elif exp_config['experiment #'] == 6:
        title_text_p1 = r'Akku + Eddy : $y_{obj}= $' + str(0.0) + '$m$, ' + '$z_{obj}= $' + str(0.374) + '$m$;    '

        for ax in axes_r1:
            ax.set_xlabel(r'$x_{obj}[m]$', fontsize=16)

    else:
        print("ERROR: no experiment selected")

    fig.suptitle(title_text_p1 + r'$x_{det}= $' + str(np.round(det_pos[0], decimals=2)) + '$m$, ' +
                 '$y_{det}= $' + str(np.round(det_pos[1], decimals=2)) + '$m$, ' +
                 r'$z_{det}= $' + str(np.round(det_pos[2], decimals=2)) + '$m$ ' +
                 r', detector array: ' + str(int(exp_config['det_array'])), fontsize=14)

    for i in range(n_f):
        ax1.plot(z_vec, b_at_det[:, s_ind, i, 0] / mu0)
        ax2.plot(z_vec, b_at_det[:, s_ind, i, 1] / mu0)
        ax3.plot(z_vec, b_at_det[:, s_ind, i, 2] / mu0, label='f=' + str(freqs[i]) + "Hz")

    ax3.legend(loc='upper right')
    # exp_number = exp_config['experiment #']

    fig.savefig("./save_plots/Active_" + exp_config['data type'] + '_' + exp_config['Object'] + '_' + str(
        exp_config['experiment #']) + '_array_' + str(exp_config['det_array']) + '_sensor_' +
                str(exp_config['sensor index']) + ".pdf", bbox_inches='tight')


def bfz_plot_measure_rel_alpha(b_at_det_pack, freqs, det_pos, det_ind, exp_config):
    mu0 = 4 * np.pi * 10 ** (-7)
    b_at_det, b_at_det_ref = b_at_det_pack
    n_det, n_f, n_dim = np.shape(b_at_det)
    alphas = np.linspace(0, 1, n_det)
    fig = plt.figure(constrained_layout=True, figsize=(8, 10))
    obs_pos = det_pos[det_ind, 1]

    fig.suptitle(r'$y_{det}=$' + str(np.round(obs_pos, decimals=3)) + r'$m$', fontsize=16)
    spec = gridspec.GridSpec(ncols=12, nrows=10, figure=fig)

    # ax1 = fig.add_subplot(spec[0, 0:4])
    # ax2 = fig.add_subplot(spec[0, 4:8])
    # ax3 = fig.add_subplot(spec[0, 8:12])

    # ax1s = fig.add_subplot(spec[1, 0:4])
    # ax2s = fig.add_subplot(spec[1, 4:8])
    # ax3s = fig.add_subplot(spec[1, 8:12])

    # ax4 = fig.add_subplot(spec[2, 0:6])
    # ax5 = fig.add_subplot(spec[2, 6:12])

    ax1 = fig.add_subplot(spec[1:4, 0:4])
    ax2 = fig.add_subplot(spec[1:4, 4:8])
    ax3 = fig.add_subplot(spec[1:4, 8:12])

    ax1s = fig.add_subplot(spec[4:7, 0:4])
    ax2s = fig.add_subplot(spec[4:7, 4:8])
    ax3s = fig.add_subplot(spec[4:7, 8:12])

    ax4 = fig.add_subplot(spec[7:10, 0:6])
    ax5 = fig.add_subplot(spec[7:10, 6:12])

    axes_r1 = [ax1, ax2, ax3]
    axes_r2 = [ax1s, ax2s, ax3s]
    axes_r3 = [ax4, ax5]

    for ax in axes_r1:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                            useMathText=None)
        ax.grid()

    for ax in axes_r2:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                            useMathText=None)
        ax.set_xlabel(r'$\alpha/\pi$')
        ax.grid()

    for ax in axes_r3:
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                            useMathText=None)
        ax.set_xlabel(r'$\alpha/\pi$')
        ax.grid()

    ax1.set_ylabel(r'$(|B_x|-|B_x^0|)/(\mu_0 I/(2\pi))$')
    ax2.set_ylabel(r'$(|B_y|-|B_y^0|)/(\mu_0 I/(2\pi))$')
    ax3.set_ylabel(r'$(|B_z|-|B_z^0|)/(\mu_0 I/(2\pi))$')
    ax1s.set_ylabel(r'$|B_x|/(\mu_0 I/(2\pi))$')
    ax2s.set_ylabel(r'$|B_y|/(\mu_0 I/(2\pi))$')
    ax3s.set_ylabel(r'$|B_z|/(\mu_0 I/(2\pi))$')

    ax4.set_ylabel(r'$|arg(B_x/B_y)|/\pi$')
    ax5.set_ylabel(r'$|arg(B_z/B_y)|/\pi$')

    for i in range(n_f):
        ax1.plot(alphas,
                 (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])) / (mu0 / (2 * np.pi)))
        ax2.plot(alphas,
                 (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])) / (mu0 / (2 * np.pi)))
        ax3.plot(alphas,
                 (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])) / (mu0 / (2 * np.pi)),
                 label='f=' + str(freqs[i]) + "Hz")

        ax1s.plot(alphas, abs(b_at_det[:, i, 0]) / (mu0 / (2 * np.pi)))
        ax2s.plot(alphas, abs(b_at_det[:, i, 1]) / (mu0 / (2 * np.pi)))
        #ax2s.plot(alphas, abs(b_at_det_ref[:, i, 2]) / (mu0 / (2 * np.pi)))
        ax3s.plot(alphas, abs(b_at_det[:, i, 2]) / (mu0 / (2 * np.pi)))


        ax4.plot(alphas, abs(1 / np.pi * (np.angle(b_at_det[:, i, 0] / b_at_det[:, i, 1]))), linestyle='-')
        ax5.plot(alphas, abs(1 / np.pi * (np.angle(b_at_det[:, i, 2] / b_at_det[:, i, 1]))), linestyle='-',
                 label='f=' + str(freqs[i]) + "Hz")

    ax4.set_ylim(-0.1, 1.1)
    ax5.set_ylim(-0.1, 1.1)

    ax5.legend(loc='upper right')

    exp_number = exp_config['experiment #']
    fig.savefig("./save_plots/eddy_alpha_" + str(exp_number) + ".pdf", bbox_inches='tight')

    # plt.subplots_adjust(top=0.8)
    # plt.subplots_adjust(left=0.3)

    # spec.tight_layout(fig,rect=[0.5, 0, 1, 1], h_pad=0.5)

    def bfz_plot_measure_rel_x(b_at_det_pack, freqs, det_pos, row_ind):
        mu0 = 4 * np.pi * 10 ** (-7)
        b_at_det, b_at_det_ref = b_at_det_pack
        n_det, n_f, n_dim = np.shape(b_at_det)
        fig = plt.figure(constrained_layout=True, figsize=(8, 10))
        spec = gridspec.GridSpec(ncols=12, nrows=3, figure=fig)

        ax1 = fig.add_subplot(spec[0, 0:4])
        ax2 = fig.add_subplot(spec[0, 4:8])
        ax3 = fig.add_subplot(spec[0, 8:12])

        ax1s = fig.add_subplot(spec[1, 0:4])
        ax2s = fig.add_subplot(spec[1, 4:8])
        ax3s = fig.add_subplot(spec[1, 8:12])

        # ax4 = fig.add_subplot(spec[2, 0:3])
        ax4 = fig.add_subplot(spec[2, 0:6])
        ax5 = fig.add_subplot(spec[2, 6:12])
        # ax7 = fig.add_subplot(spec[2, 9:12])
        axes_r1 = [ax1, ax2, ax3]
        axes_r2 = [ax1s, ax2s, ax3s]
        axes_r3 = [ax4, ax5]

        for ax in axes_r1:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                                useMathText=None)
            ax.grid()

        for ax in axes_r2:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                                useMathText=None)
            ax.set_xlabel(r'y[m]')
            ax.grid()

        for ax in axes_r3:
            ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None,
                                useMathText=None)
            ax.set_xlabel(r'y[m]')
            ax.grid()

        ax1.set_ylabel(r'$(|B_x|-|B_x^0|)/(\mu_0 I/(2\pi))$')
        ax2.set_ylabel(r'$(|B_y|-|B_y^0|)/(\mu_0 I/(2\pi))$')
        ax3.set_ylabel(r'$(|B_z|-|B_z^0|)/(\mu_0 I/(2\pi))$')
        ax1s.set_ylabel(r'$|B_x|/(\mu_0/(2\pi))$')
        ax2s.set_ylabel(r'$|B_y|/(\mu_0/(2\pi))$')
        ax3s.set_ylabel(r'$|B_z|/(\mu_0/(2\pi))$')

        # ax4.set_ylabel(r'$|B_x/B_y|$')
        ax4.set_ylabel(r'$|arg(B_x/B_y)|/\pi$')
        ax5.set_ylabel(r'$|arg(B_z/B_y)|/\pi$')
        # ax7.set_ylabel(r'$arg(B_z)/\pi$')

        for i in range(n_f):
            # ax1.plot(det_pos[:, row_ind],
            #         (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])) / abs(b_at_det_ref[:, i, 0]))
            # ax2.plot(det_pos[:, row_ind],
            #         (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])) / abs(b_at_det_ref[:, i, 1]))
            # ax3.plot(det_pos[:, row_ind],
            #         (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])) / abs(b_at_det_ref[:, i, 2]),
            #         label='f=' + str(freqs[i]) + "Hz")

            ax1.plot(det_pos[:, row_ind],
                     (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])))
            ax2.plot(det_pos[:, row_ind],
                     (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])))
            ax3.plot(det_pos[:, row_ind],
                     (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])),
                     label='f=' + str(freqs[i]) + "Hz")

            ax1s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 0]) / (mu0 / (2 * np.pi)))
            ax2s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 1]) / (mu0 / (2 * np.pi)))
            ax3s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 2]) / (mu0 / (2 * np.pi)))
            # ax5.plot(det_pos[:, row_ind], 1 / np.pi * (np.angle(b_at_det[:, i, 0])/np.angle(b_at_det[:, i, 1])), linestyle='-',
            #         label='Bx, f=' + str(freqs[i]) + "Hz")
            ax4.plot(det_pos[:, row_ind], abs(1 / np.pi * (np.angle(b_at_det[:, i, 0] / b_at_det[:, i, 1]))),
                     linestyle='-')
            ax5.plot(det_pos[:, row_ind], abs(1 / np.pi * (np.angle(b_at_det[:, i, 2] / b_at_det[:, i, 1]))),
                     linestyle='-',
                     label='f=' + str(freqs[i]) + "Hz")
            # ax7.plot(det_pos[:, row_ind], 1 / np.pi * (np.angle(b_at_det[:, i, 2])), linestyle=':',
            #         label='Bz, f=' + str(freqs[i]) + "Hz")
            # print(1 / np.pi * (np.angle(b_at_det[:, i, 0])))
            # exit()

            # ax5.plot(det_pos[:, 0], 1/np.pi*(np.angle(b_at_det[:, i, 2])-np.angle(b_at_det[:, i, 0])))
        ax4.set_ylim(-0.1, 1.1)
        ax5.set_ylim(-0.1, 1.1)
        # ax7.set_ylim(-0.1, 1.1)
        # ax1.plot(det_pos[:, 1], 1/det_pos[:, 1])
        ax5.legend()
