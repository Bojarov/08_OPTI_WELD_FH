import numpy as np
import matplotlib.pyplot as plt
import obs_calc.i_fh__calc
import obs_calc.b_xyz
import obs_calc_analytical.bc_wire_dc

from matplotlib import rc, font_manager

sizeOfFont = 12
fontProperties = {'family': 'sans-serif', 'sans-serif': ['Helvetica'],
                  'weight': 'normal', 'size': sizeOfFont}
ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
                                         size=sizeOfFont, weight='normal', stretch='normal')
rc('text', usetex=True)
rc('font', **fontProperties)


def b_wire_comp_plot(xyz0_proj_re, xyz0_proj_im, params):
    """Very dirty, UNFINISHED function that plots analytic and numeric results
    from the Ferreira paper, still errors! need to finish!
    normalization of analytic result by total current through wire is missing
    only real part of FT is plotted, connection to FH result unclear
    factor of pi from FT not included/unclear
    """

    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params  # unpack params
    # zc_vec=np.array([0,1.5,1.9,2.0])
    # zc_vec=np.array([0,1.5,1.9,2.0,2.1,2.5,3.0])

    i_fh = i_fh_calc(params, xyz0_proj_re, xyz0_proj_im)

    zc_vec = np.array([2.1])
    rc_vec = np.linspace(1, 3, 10)
    rc_vec2 = np.linspace(0.01, 3, 20)
    bc_vec = np.zeros(len(rc_vec))
    bc_vec2 = np.zeros(len(rc_vec2))
    fig_b_wire_comp = plt.figure(figsize=(8, 5))
    ax_b_wire_comp = fig_b_wire_comp.add_subplot(1, 1, 1)

    # extract FH results
    for j in range(len(zc_vec)):
        if zc_vec[j] <= 2:
            for i in range(len(rc_vec)):
                bc_vec[i] = 0.1 * b_xyz(rc_vec[i] * ro_t, 0, zc_vec[j] * ro_t + flen / 2,
                                        xyz0_proj_re, xyz0_proj_im, params)

            ax_b_wire_comp.plot(rc_vec, bc_vec / (mu_0 * i_fh) * (2 * np.pi),
                                label=r'$\frac{z_C}{R}=$' + str(zc_vec[j]))

        if zc_vec[j] > 2:
            for i in range(len(rc_vec2)):
                bc_vec2[i] = ro_t * b_xyz(rc_vec2[i] * ro_t, 0, zc_vec[j] * ro_t + flen / 2,
                                          xyz0_proj_re, xyz0_proj_im, params)

            ax_b_wire_comp.plot(rc_vec2, bc_vec2 / (mu_0 * i_fh) * (2 * np.pi),
                                label=r'$\frac{z_C}{R}=$' + str(zc_vec[j]))
    # zc_vec3=np.array([0,1.5,1.9,1.99,2.1,2.5,3.0])
    zc_vec3 = np.array([2.1])
    rc_vec3 = np.linspace(0.01, 3, 50)
    bc_vec3 = np.zeros(len(rc_vec3))
    for i in zc_vec3:
        for j in range(len(rc_vec3)):
            bc_vec3[j] = bc_wire(rc_vec3[j], i, 1, flen / ro_t)[0]
        ax_b_wire_comp.plot(rc_vec3, bc_vec3, linestyle="--")

    ax_b_wire_comp.set_xlabel(r'$R_C/R$')
    ax_b_wire_comp.set_ylabel(r'$B_C/\frac{\mu_0 I_0}{2\pi}$')
    ax_b_wire_comp.legend(loc='upper left', fontsize=10)
    # ax31.set_xlim(-0.1, max(rc_vec)+0.1)
    # ax31.set_ylim(0, 1.1*max(bc_vec/(mu_0*i_fh)*(2*np.pi)))
