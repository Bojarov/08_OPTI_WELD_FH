import numpy as np
import matplotlib.pyplot as plt
import code.obs_calc_para_point as ocpp
import code.obs_calc_j as ocj
# from matplotlib import cm
# from matplotlib import rc, font_manager
import code.obs_calc_analytical as ocana
import code.obs_calc_b_field as ocb
import itertools
import matplotlib.gridspec as gridspec


def I_wire_plot(freq_vec, freq_vec_a, params, damage_params, filament_params, output_params,
                r_sub_vec, l_sub_vec, node_dens_vec, loop_list, sub_con_list):
    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

    i_fh_vec = np.zeros(len(freq_vec))
    i_b_vec = np.zeros(len(freq_vec_a))
    i_k_vec = np.zeros(len(freq_vec_a))

    for i in range(len(freq_vec)):
        print(i)

        params_l = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq_vec[i]]

        # calcs current densities at point in parameter + config space
        tube_lists = ocpp.para_point_calc(params_l, damage_params, filament_params, output_params, r_sub_vec, l_sub_vec,
                                          node_dens_vec, loop_list, sub_con_list, mode=0)

        # unpacking the lists
        tube_node_lists, tube_segment_lists, tube_pts_lists = tube_lists

        i_fh_vec[i] = ocj.i_fh_calc(r_sub_vec, l_sub_vec, node_dens_vec, tube_segment_lists,
                                    loop_list, sub_con_list, params_l, damage_params, filament_params)

        # print(i_fh_vec[0])
        # exit()
    for i in range(len(freq_vec_a)):
        params_l = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq_vec_a[i]]

        i_k_vec[i] = ocana.I_Knight_calc(params_l)
        i_b_vec[i] = ocana.I_Bessel_calc(params_l)

    fig = plt.figure(figsize=(3, 2))
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(freq_vec_a, i_k_vec, label=r'$I_{Knight}$')
    ax.plot(freq_vec_a, i_b_vec, label=r'$I_{Bessel}$', linestyle='--')
    ax.plot(freq_vec, i_fh_vec, label=r'$I_{FH}$', marker='x', linestyle='')
    ax.set_xscale('log')
    ax.legend(loc='upper right', fontsize=20)
    ax.set_xlabel(r'$f$[Hz]')
    ax.set_ylabel(r'$I$[A]')


# def B_wire_comp_plot(xyz0_proj_re, xyz0_proj_im, params):
def B_wire_comp_plot(params, damage_params, filament_params, r_sub_vec,
                     l_sub_vec, node_dens_vec, loop_list, sub_con_list, ana_pts=100,
                     ana_ac_pts=10, num_pts=20):
    """Compares analytic and FH B field of a massive finite wire
    """

    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

    fig = plt.figure(figsize=(8, 5))

    ax = fig.add_subplot(1, 1, 1)
    L = flen / ro_t

    # zc_vec=np.array([0.0,1.5,1.9,1.99,2.1,2.5,3.0])
    zc_vec = np.array([0.0])  # , 1.5, 2.1, 3.0])
    # zc_vec=np.array([0, L/2*0.99, L/2*1.01, 3/4*L])

    # zc_vec    = np.array([0.0])

    rc_vec_a = np.linspace(0.01, 3.0, ana_pts)
    rc_vec_a_ac = np.linspace(0.01, 3.0, ana_ac_pts)
    rc_vec_fh_pre = np.linspace(0.01, 3.0, num_pts)
    rc_vec_fh_pre[0] = 0.05
    bc_vec_a = np.zeros(len(rc_vec_a))
    bc_vec_a_ac = np.zeros(len(rc_vec_a_ac))

    rc_vec_a_ac = np.linspace(0.01, 3.0, ana_ac_pts)

    for i in zc_vec:
        for j in range(len(rc_vec_a)):
            # print("zC = ", i ,"RC_a = ", rc_vec_a[j])

            bc_vec_a[j] = ocana.BC_wire(rc_vec_a[j], i, 1, L)[0]
        ax.plot(rc_vec_a, bc_vec_a, label=r'$\frac{z_C}{R}=$' + str(i))

    # for i in zc_vec:
    #    for j in range(len(rc_vec_a_ac)):
    #        print("zC = ", i ,"RC_a = ", rc_vec_a_ac[j])

    #        bc_vec_a_ac[j] = ocana.BC_wire_AC( rc_vec_a_ac[j], i, 1, L, params)[0]
    #    ax.plot(rc_vec_a_ac, bc_vec_a_ac)

    # calculates s current densities and B_fields at points in parameter + config space

    for i in zc_vec:

        rc_vec_fh = []
        bc_vec_fh = []
        for j in range(len(rc_vec_fh_pre)):
            print("zC = ", i, "rc = ", rc_vec_fh_pre[j])

            x_detect = rc_vec_fh_pre[j] * ro_t
            y_detect = 0
            z_detect = i * ro_t + L / 2 * ro_t
            r_detect = np.sqrt(x_detect ** 2 + y_detect ** 2)

            if r_detect / ro_t > 1.0 or i > L / 2:
                print("Detector position is outside conductor - calculating B field numerically.")
                rc_vec_fh.append(r_detect / ro_t)

                output_params = [0, x_detect, y_detect, z_detect]

                tube_lists = ocpp.para_point_calc(params, damage_params, filament_params, output_params, r_sub_vec,
                                                  l_sub_vec, node_dens_vec, loop_list, sub_con_list, mode=1)

                # unpacking the lists
                tube_node_lists, tube_segment_lists, tube_pts_lists = tube_lists

                bc_vec_fh.append(ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                                   node_dens_vec, l_sub_vec, tube_segment_lists, loop_list,
                                                   sub_con_list))

            else:
                print("Error: Defined detector position (" + str(x_detect) + str(y_detect) + str(z_detect) +
                      ") is inside conductor. Skipping numerical evaluation.")

        # ax.plot(rc_vec_a, bc_vec_a , linestyle="-", label = r'$\frac{z_C}{R}=$'+ str(zc_vec[j]))
        ax.plot(rc_vec_fh, ro_t * np.array(bc_vec_fh), linestyle="--")
        ax.scatter(rc_vec_fh, ro_t * np.array(bc_vec_fh), marker="x", color="green")
        ax.set_xlabel(r'$R_C/R$')
        ax.set_ylabel(r'$B_C/\frac{\mu_0 I_0}{2\pi}$')
        ax.legend(loc='upper right', fontsize=15)
        # ax.set_ylim(0.0, 1.0)
        # print("gugugugugu")
        # print("gugugugugu")
        # print("gugugugugu")
        # print("gugugugugu")
        # print("gugugugugu")
        # print(bc_vec_fh)
    plt.tight_layout()


def b_at_det_iter(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                  params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    for i in itertools.product(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v):
        print(i)

        # iterate through all the parameter combinations

        x_detect = i[5]
        y_detect = i[6]
        z_detect = i[7]

        output_params = [0, x_detect, y_detect, z_detect]

        #####################################################################

        ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

        freq = i[8]

        params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

        #####################################################################
        rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

        phi_c_weld = i[0]
        d_weld = i[1]
        rd_i = i[2]
        rd_o = i[3]
        l_weld = i[4]

        damage_params = [rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage]

        if l_weld == flen:
            l_sub_vec = np.array([0.0, flen])
        else:
            l_sub_vec = np.array([0.0, 0.5 * (flen - l_weld), 0.5 * (flen + l_weld), flen])

        tube_lists = ocpp.para_point_calc(params, damage_params, filament_params, output_params, r_sub_vec,
                                          l_sub_vec, node_dens_vec, [], [], mode=1)

        # unpacking the lists
        tube_node_lists, tube_segment_lists, tube_pts_lists = tube_lists

        ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                          node_dens_vec, l_sub_vec, tube_segment_lists, [], [])


def b_at_det_plot_f(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                    params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    fig = plt.figure(figsize=(8, 5))

    ax = fig.add_subplot(1, 1, 1)

    for i in f_v:
        b_vec = []
        f_ind = f_v.index(i)

        for j in xd_v:

            # updating output params
            x_detect = j
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = i

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            # updating damage params
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = 0.0
            d_weld = 0.0
            rd_i = ri_t
            rd_o = ro_t
            l_weld = lw_v[0]

            # updating damage params

            damage_params = [rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage]

            if l_weld == flen:
                l_sub_vec = np.array([0.0, flen])
            else:
                l_sub_vec = np.array([0.0, 0.5 * (flen - l_weld), 0.5 * (flen + l_weld), flen])

            tube_lists = ocpp.para_point_calc(params, damage_params, filament_params, output_params, r_sub_vec,
                                              l_sub_vec, node_dens_vec, [], [], mode=1)
            # unpacking the lists
            tube_node_lists, tube_segment_lists, tube_pts_lists = tube_lists
            b_vec.append(ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                           node_dens_vec, l_sub_vec, tube_segment_lists, [], []))

        if f_ind == 0:
            b_vec_0 = b_vec

        break

    for i in f_v:
        b_vec = []
        f_ind = f_v.index(i)

        for j in xd_v:

            x_detect = j
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            #####################################################################

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = i

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            #####################################################################
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = phi_cv[0]
            d_weld = dw_v[0]
            rd_i = rdi_v[0]
            rd_o = rdo_v[0]
            l_weld = lw_v[0]

            damage_params = [rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage]

            if l_weld == flen:
                l_sub_vec = np.array([0.0, flen])
            else:
                l_sub_vec = np.array([0.0, 0.5 * (flen - l_weld), 0.5 * (flen + l_weld), flen])

            tube_lists = ocpp.para_point_calc(params, damage_params, filament_params, output_params, r_sub_vec,
                                              l_sub_vec, node_dens_vec, [], [], mode=1)

            # unpacking the lists
            tube_node_lists, tube_segment_lists, tube_pts_lists = tube_lists
            # b_at_d=[]
            b_vec.append(ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                           node_dens_vec, l_sub_vec, tube_segment_lists, [], []))

        # if f_ind == 0:
        #    b_vec_0 = b_vec

        ax.plot(xd_v, (np.array(b_vec) - np.array(b_vec_0)) / (np.array(b_vec_0)), label=r'$f=$' + str(i) + r'Hz')
        ax.set_xlabel(r'$R_C/R$')
        ax.set_ylabel(r'$(\tilde{B}_f-\tilde{B}_0)/\tilde{B}_0$', fontsize=15)
        ax.legend(loc='upper right', fontsize=15)

    plt.show()


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


def bfz_plot_measure_rel_x(b_at_det_pack, freqs, det_pos, alpha_ind):
    row_ind = 1
    mu0 = 4 * np.pi * 10 ** (-7)
    b_at_det, b_at_det_ref = b_at_det_pack
    n_det, n_f, n_dim = np.shape(b_at_det)
    fig = plt.figure(constrained_layout=True, figsize=(4, 10))
    alphas = np.linspace(0, 1, 13)
    fig.suptitle(r'$\alpha/\pi=$' + str(alphas[alpha_ind]), fontsize=16)

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

    ax1.set_ylabel(r'$(|B_x|-|B_x^0|)/|B_x^0|$')
    ax2.set_ylabel(r'$(|B_y|-|B_y^0|)/|B_y^0|$')
    ax3.set_ylabel(r'$(|B_z|-|B_z^0|)/|B_z^0|$')
    ax1s.set_ylabel(r'$|B_x|/(\mu_0/(2\pi))$')
    ax2s.set_ylabel(r'$|B_y|/(\mu_0/(2\pi))$')
    ax3s.set_ylabel(r'$|B_z|/(\mu_0/(2\pi))$')

    ax4.set_ylabel(r'$|arg(B_x/B_y)|/\pi$')
    ax5.set_ylabel(r'$|arg(B_z/B_y)|/\pi$')

    for i in range(n_f):
        ax1.plot(det_pos[:, row_ind],
                 (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])) / abs(b_at_det_ref[:, i, 0]))
        ax2.plot(det_pos[:, row_ind],
                 (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])) / abs(b_at_det_ref[:, i, 1]))
        ax3.plot(det_pos[:, row_ind],
                 (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])) / abs(b_at_det_ref[:, i, 2]),
                 label='f=' + str(freqs[i]) + "Hz")

        ax1s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 0]) / (mu0 / (2 * np.pi)))
        ax2s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 1]) / (mu0 / (2 * np.pi)))
        ax3s.plot(det_pos[:, row_ind], abs(b_at_det[:, i, 2]) / (mu0 / (2 * np.pi)))
        ax4.plot(det_pos[:, row_ind], abs(1 / np.pi * (np.angle(b_at_det[:, i, 0] / b_at_det[:, i, 1]))), linestyle='-')
        ax5.plot(det_pos[:, row_ind], abs(1 / np.pi * (np.angle(b_at_det[:, i, 2] / b_at_det[:, i, 1]))), linestyle='-',
                 label='f=' + str(freqs[i]) + "Hz")

    ax4.set_ylim(-0.1, 1.1)
    ax5.set_ylim(-0.1, 1.1)

    ax5.legend(loc='lower right')




def bfz_plot_measure_rel_alpha(b_at_det_pack, freqs, det_pos, det_ind):
    mu0 = 4 * np.pi * 10 ** (-7)
    b_at_det, b_at_det_ref = b_at_det_pack
    n_det, n_f, n_dim = np.shape(b_at_det)
    alphas = np.linspace(0, 1, n_det)
    fig = plt.figure(constrained_layout=True, figsize=(4, 10))
    obs_pos = det_pos[det_ind, 1]


    fig.suptitle(r'$y_{det}=$'+str(obs_pos)+r'$m$', fontsize=16)
    spec = gridspec.GridSpec(ncols=12, nrows=10, figure=fig)


    #ax1 = fig.add_subplot(spec[0, 0:4])
    #ax2 = fig.add_subplot(spec[0, 4:8])
    #ax3 = fig.add_subplot(spec[0, 8:12])

    #ax1s = fig.add_subplot(spec[1, 0:4])
    #ax2s = fig.add_subplot(spec[1, 4:8])
    #ax3s = fig.add_subplot(spec[1, 8:12])

    #ax4 = fig.add_subplot(spec[2, 0:6])
    #ax5 = fig.add_subplot(spec[2, 6:12])

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

    ax1.set_ylabel(r'$(|B_x|-|B_x^0|)/|B_x^0|$')
    ax2.set_ylabel(r'$(|B_y|-|B_y^0|)/|B_y^0|$')
    ax3.set_ylabel(r'$(|B_z|-|B_z^0|)/|B_z^0|$')
    ax1s.set_ylabel(r'$|B_x|/(\mu_0/(2\pi))$')
    ax2s.set_ylabel(r'$|B_y|/(\mu_0/(2\pi))$')
    ax3s.set_ylabel(r'$|B_z|/(\mu_0/(2\pi))$')

    ax4.set_ylabel(r'$|arg(B_x/B_y)|/\pi$')
    ax5.set_ylabel(r'$|arg(B_z/B_y)|/\pi$')

    for i in range(n_f):
        ax1.plot(alphas,
                 (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])) / abs(b_at_det_ref[:, i, 0]))
        ax2.plot(alphas,
                 (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])) / abs(b_at_det_ref[:, i, 1]))
        ax3.plot(alphas,
                 (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])) / abs(b_at_det_ref[:, i, 2]),
                 label='f=' + str(freqs[i]) + "Hz")

        ax1s.plot(alphas, abs(b_at_det[:, i, 0]) / (mu0 / (2 * np.pi)))
        ax2s.plot(alphas, abs(b_at_det[:, i, 1]) / (mu0 / (2 * np.pi)))
        ax3s.plot(alphas, abs(b_at_det[:, i, 2]) / (mu0 / (2 * np.pi)))
        ax4.plot(alphas, abs(1 / np.pi * (np.angle(b_at_det[:, i, 0] / b_at_det[:, i, 1]))), linestyle='-')
        ax5.plot(alphas, abs(1 / np.pi * (np.angle(b_at_det[:, i, 2] / b_at_det[:, i, 1]))), linestyle='-',
                 label='f=' + str(freqs[i]) + "Hz")

    ax4.set_ylim(-0.1, 1.1)
    ax5.set_ylim(-0.1, 1.1)

    ax5.legend(loc='lower right')
    #plt.subplots_adjust(top=0.8)
    #plt.subplots_adjust(left=0.3)

    #spec.tight_layout(fig,rect=[0.5, 0, 1, 1], h_pad=0.5)



    def bfz_plot_measure_rel_x(b_at_det_pack, freqs, det_pos, row_ind):
        mu0 = 4 * np.pi * 10 ** (-7)
        b_at_det, b_at_det_ref = b_at_det_pack
        n_det, n_f, n_dim = np.shape(b_at_det)
        fig = plt.figure(constrained_layout=True, figsize=(4, 8))
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

        ax1.set_ylabel(r'$(|B_x|-|B_x^0|)/|B_x^0|$')
        ax2.set_ylabel(r'$(|B_y|-|B_y^0|)/|B_y^0|$')
        ax3.set_ylabel(r'$(|B_z|-|B_z^0|)/|B_z^0|$')
        ax1s.set_ylabel(r'$|B_x|/(\mu_0/(2\pi))$')
        ax2s.set_ylabel(r'$|B_y|/(\mu_0/(2\pi))$')
        ax3s.set_ylabel(r'$|B_z|/(\mu_0/(2\pi))$')

        # ax4.set_ylabel(r'$|B_x/B_y|$')
        ax4.set_ylabel(r'$|arg(B_x/B_y)|/\pi$')
        ax5.set_ylabel(r'$|arg(B_z/B_y)|/\pi$')
        # ax7.set_ylabel(r'$arg(B_z)/\pi$')

        for i in range(n_f):
            ax1.plot(det_pos[:, row_ind],
                     (abs(b_at_det[:, i, 0]) - abs(b_at_det_ref[:, i, 0])) / abs(b_at_det_ref[:, i, 0]))
            ax2.plot(det_pos[:, row_ind],
                     (abs(b_at_det[:, i, 1]) - abs(b_at_det_ref[:, i, 1])) / abs(b_at_det_ref[:, i, 1]))
            ax3.plot(det_pos[:, row_ind],
                     (abs(b_at_det[:, i, 2]) - abs(b_at_det_ref[:, i, 2])) / abs(b_at_det_ref[:, i, 2]),
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



