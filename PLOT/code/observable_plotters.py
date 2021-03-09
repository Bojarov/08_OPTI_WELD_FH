import numpy as np
import matplotlib.pyplot as plt
import code.obs_calc_para_point as ocpp
import code.obs_calc_j as ocj
import code.obs_calc_analytical as ocana
import code.obs_calc_b_field as ocb
import itertools
from tempfile import TemporaryFile


# sizeOfFont = 6
# fontProperties = {'family':'sans-serif','sans-serif':['Helvetica'],
#    'weight' : 'normal', 'size' : sizeOfFont}
# ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
#    size=sizeOfFont, weight='normal', stretch='normal')
# rc('text', usetex=True)
# rc('font',**fontProperties)

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

    fig1 = plt.figure(figsize=(6, 4))

    ax1 = fig1.add_subplot(1, 1, 1)
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
        ax1.plot(rc_vec_a, bc_vec_a, label=r'$\frac{z_C}{R}=$' + str(i))

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
        ax1.plot(rc_vec_fh, ro_t * np.array(bc_vec_fh), linestyle="--")
        ax1.scatter(rc_vec_fh, ro_t * np.array(bc_vec_fh), marker="x", color="green")
        ax1.set_xlabel(r'$R_C/R$')
        ax1.set_ylabel(r'$B_C/\frac{\mu_0 I_0}{2\pi}$')
        ax1.legend(loc='upper right', fontsize=15)
        # ax.set_ylim(0.0, 1.0)
        # print(bc_vec_fh)
    plt.tight_layout()


# def b_at_detector_plot():
# print("")

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
        # print(l_weld)
        # exit()

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


def b_at_det_plot_f1(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    fig = plt.figure(figsize=(3, 2))

    ax = fig.add_subplot(1, 1, 1)

    b_vec = []

    for j in xd_v:

        # updating output params
        x_detect = j
        y_detect = yd_v[0]
        z_detect = zd_v[0]

        output_params = [0, x_detect, y_detect, z_detect]

        ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

        freq = 0.001

        params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

        # updating damage params
        rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

        phi_c_weld = 0.0
        d_weld = 0.0
        rd_i = ri_t
        rd_o = ro_t
        l_weld = flen  # lw_v[0]

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

    b_vec_0 = b_vec

    for i in dw_v:
        b_vec = []
        dw_ind = dw_v.index(i)

        for j in xd_v:

            x_detect = j
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            #####################################################################

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = f_v[0]

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            #####################################################################
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = phi_cv[0]
            d_weld = i  # dw_v[0]
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
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.plot(xd_v, (np.array(b_vec) - np.array(b_vec_0)) / (np.array(b_vec_0)), marker='o', markersize=2,
                label=r'$d_d=$' + str(i) + r'm')
        ax.set_xlabel(r'$x[m]$', fontsize=20)
        ax.set_ylabel(r'$(\tilde{B}_f-\tilde{B}_0)/\tilde{B}_0$', fontsize=20)
        ax.legend(loc='upper right', fontsize=10)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
        plt.title(
            r'$\varphi_d=$' + str(phi_c_weld / np.pi) + r'$\pi$, ' + r'$f=$' + str(freq) + r'$Hz$' + r', $l_d=$' + str(
                l_weld) + r'$m$')
        plt.grid(linewidth=0.5)
        # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
        # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()

    plt.show()


def b_at_det_plot_f2(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    bx_vecs_ideal = []
    for i in f_v:
        bx_vec_ideal = []
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
            bx_vec_ideal.append(ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                                  node_dens_vec, l_sub_vec, tube_segment_lists, [], [])[0])

        bx_vecs_ideal.append(bx_vec_ideal)

    bx_vecs = []
    for i in f_v:
        bx_vec = []
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
            bx_vec.append(ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                            node_dens_vec, l_sub_vec, tube_segment_lists, [], [])[0])
        bx_vecs.append(bx_vec)

    fig = plt.figure(figsize=(3, 2))

    ax = fig.add_subplot(1, 1, 1)

    for i in f_v:
        f_ind = f_v.index(i)

        plot_arr = (np.array(bx_vecs[f_ind]) - np.array(bx_vecs_ideal[f_ind])) / (np.array(bx_vecs_ideal[f_ind]))
        ax.plot(xd_v, plot_arr, marker='o', markersize=2,
                label=r'$f=$' + str(i) + r'Hz')

    # for i in f_v:
    #    f_ind = f_v.index(i)
    #    ax.plot(xd_v, np.array(bx_vecs_ideal[f_ind])+f_ind, marker='o', markersize=2, label=r'ideal $f=$' + str(i) + r'Hz')

    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)

    ax.set_xlabel(r'$x[m]$', fontsize=20)
    ax.set_ylabel(r'$(\tilde{B}^d_x-\tilde{B}^i_x)/\tilde{B}_x^i$', fontsize=20)
    ax.legend(loc='upper right', fontsize=10)
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
    # plt.title(r'$\varphi_d=$' + str(phi_c_weld / np.pi) + r'$\pi$, ' + r'$d_d=$' + str(
    #    d_weld) + r'$m$' + r', $l_d=$' + str(l_weld) + r'$m$')
    # plt.grid(linewidth=0.5)
    # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
    # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()
    plt.grid(True)

    with open('bvec.npy', 'wb') as f1:

        np.save(f1, np.array(bx_vecs))
    with open('bvec_i.npy', 'wb') as f2:

        np.save(f2, np.array(bx_vecs_ideal))

    plt.show()


def b_at_det_plot_f3(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    fig = plt.figure(figsize=(3, 2))

    ax = fig.add_subplot(1, 1, 1)

    # get the the ideal pipeline case at low low frequency
    b_vec = []
    for j in xd_v:

        # updating output params
        x_detect = j
        y_detect = yd_v[0]
        z_detect = zd_v[0]

        output_params = [0, x_detect, y_detect, z_detect]

        ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

        freq = 0.001

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

    b_vec_0 = b_vec

    for i in rdo_v:
        b_vec = []
        rdo_ind = rdo_v.index(i)

        for j in xd_v:

            x_detect = j
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            #####################################################################

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = f_v[0]

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            #####################################################################
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = phi_cv[0]
            d_weld = dw_v[0]
            rd_i = rdi_v[0]
            rd_o = i
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

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.plot(xd_v, (np.array(b_vec) - np.array(b_vec_0)) / (np.array(b_vec_0)), marker='o', markersize=2,
                label=r'$R_o-r_{d_o}/d=$' + "%.2f" % round((ro_t - i) / (ro_t - ri_t), 2))
        ax.set_xlabel(r'$x[m]$', fontsize=20)
        ax.set_ylabel(r'$(\tilde{B}_f-\tilde{B}_0)/\tilde{B}_0$', fontsize=20)
        ax.legend(loc='upper right', fontsize=10)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
        plt.title(
            r'$\varphi_d=$' + str(phi_c_weld / np.pi) + r'$\pi$, ' + r'$f=$' + str(freq) + r'$Hz$' + r', $d_d=$' + str(
                d_weld) + r'$m$' + r', $l_d=$' + str(l_weld) + r'$m$')
        plt.grid(linewidth=0.5)
        # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
        # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()

    plt.show()


def b_at_det_plot_f4(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    fig = plt.figure(figsize=(6.025, 4.225))

    ax = fig.add_subplot(1, 1, 1)

    # get the the ideal pipeline case at low low frequency
    b_vec = []
    for j in xd_v:

        # updating output params
        x_detect = j
        y_detect = yd_v[0]
        z_detect = zd_v[0]

        output_params = [0, x_detect, y_detect, z_detect]

        ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

        freq = 0.001

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

    b_vec_0 = b_vec

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

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.plot(xd_v, (np.array(b_vec) - np.array(b_vec_0)) / (np.array(b_vec_0)), marker='o', markersize=2,
                label=r'$f=$' + str(i) + r'Hz')
        ax.set_xlabel(r'$x[m]$', fontsize=20)
        ax.set_ylabel(r'$(\tilde{B}_f-\tilde{B}_0)/\tilde{B}_0$', fontsize=20)
        ax.legend(loc='upper right', fontsize=10)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
        # plt.title(r'$\varphi_d=$'+str(phi_c_weld/np.pi)+r'$\pi$, '+r'$d_d=$'+str(d_weld)+r'$m$'+r', $l_d=$'+str(
        # l_weld)+r'$m$')

        plt.title(r'$\varphi_d=$' + str(phi_c_weld / np.pi) + r'$\pi$, ' + r'$R_o-r_{d_o}/d=$' + "%.2f" % round(
            (ro_t - rd_o) / (ro_t - ri_t), 2) + r', $d_d=$' + str(d_weld) + r'$m$')
        plt.grid(linewidth=0.5)
        # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
        # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()

    plt.show()


def b_at_det_plot_f5(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    fig = plt.figure(figsize=(6.025, 4.225))

    ax = fig.add_subplot(1, 1, 1)

    # get the the ideal pipeline case at low low frequency
    b_vec = []
    for j in xd_v:

        # updating output params
        x_detect = j
        y_detect = yd_v[0]
        z_detect = zd_v[0]

        output_params = [0, x_detect, y_detect, z_detect]

        ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

        freq = 0.001

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

    b_vec_0 = b_vec

    for i in rdi_v:
        b_vec = []
        rdi_ind = rdi_v.index(i)

        for j in xd_v:

            x_detect = j
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            #####################################################################

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = f_v[0]

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            #####################################################################
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = phi_cv[0]
            d_weld = dw_v[0]
            rd_i = i
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

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.plot(xd_v, (np.array(b_vec) - np.array(b_vec_0)) / (np.array(b_vec_0)), marker='o', markersize=2,
                label=r'$(r_{d_i}-R_i)/d=$' + "%.2f" % round((i - ri_t) / (ro_t - ri_t), 2))
        ax.set_xlabel(r'$x[m]$', fontsize=20)
        ax.set_ylabel(r'$(\tilde{B}_f-\tilde{B}_0)/\tilde{B}_0$', fontsize=20)
        ax.legend(loc='upper right', fontsize=10)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
        # plt.title(r'$\varphi_d=$'+str(phi_c_weld/np.pi)+r'$\pi$, '+r'$d_d=$'+str(d_weld)+r'$m$'+r', $l_d=$'+str(
        # l_weld)+r'$m$')

        # plt.title(r'$\varphi_d=$'+str(phi_c_weld/np.pi)+r'$\pi$, '+r'$r_{d_i}-R_i/d=$' +"%.2f" % round((rd_i-Ri)/(
        # ro_t-ri_t),2)+r', $d_d=$'+str(d_weld)+r'$m$')
        plt.title(
            r'$\varphi_d=$' + str(phi_c_weld / np.pi) + r'$\pi$, ' + r'$f=$' + str(freq) + r'$Hz$' + r', $d_d=$' + str(
                d_weld) + r'$m$')

        plt.grid(linewidth=0.5)
        # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
        # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()

    plt.show()


def b_at_det_plot_f6(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    fig = plt.figure(figsize=(6.025, 4.225))

    ax = fig.add_subplot(1, 1, 1)

    # get the the ideal pipeline case at low low frequency
    b_vec = []
    for j in xd_v:

        # updating output params
        x_detect = j
        y_detect = yd_v[0]
        z_detect = zd_v[0]

        output_params = [0, x_detect, y_detect, z_detect]

        ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

        freq = 0.001

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

    b_vec_0 = b_vec

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

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.plot(xd_v, (np.array(b_vec) - np.array(b_vec_0)) / (np.array(b_vec_0)), marker='o', markersize=2,
                label=r'$f=$' + str(i) + r'Hz')
        ax.set_xlabel(r'$x[m]$', fontsize=20)
        ax.set_ylabel(r'$(\tilde{B}_f-\tilde{B}_0)/\tilde{B}_0$', fontsize=20)
        ax.legend(loc='upper right', fontsize=10)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
        # plt.title(r'$\varphi_d=$'+str(phi_c_weld/np.pi)+r'$\pi$, '+r'$d_d=$'+str(d_weld)+r'$m$'+r', $l_d=$'+str(
        # l_weld)+r'$m$')
        print(rd_i)
        plt.title(r'$\varphi_d=$' + str(phi_c_weld / np.pi) + r'$\pi$, ' + r'$r_{d_i}-R_i/d=$' + "%.2f" % round(
            (rd_i - ri_t) / (ro_t - ri_t), 2) + r', $d_d=$' + str(d_weld) + r'$m$')

        plt.grid(linewidth=0.5)
        # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
        # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()

    plt.show()


def b_at_det_plot_f7(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    fig = plt.figure(figsize=(6.025, 4.225))

    ax = fig.add_subplot(1, 1, 1)

    # get the the ideal pipeline case at low low frequency
    b_vec = []
    for j in xd_v:

        # updating output params
        x_detect = j
        y_detect = yd_v[0]
        z_detect = zd_v[0]

        output_params = [0, x_detect, y_detect, z_detect]

        ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

        freq = 0.001

        params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

        # updating damage params
        rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

        phi_c_weld = 0.0
        d_weld = 0.0
        rd_i = ri_t
        rd_o = ro_t
        l_weld = flen  # 9.99

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

    b_vec_0 = b_vec

    for i in lw_v:
        b_vec = []
        l_weld_ind = lw_v.index(i)

        for j in xd_v:

            x_detect = j
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            #####################################################################

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = f_v[0]

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            #####################################################################
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = phi_cv[0]
            d_weld = dw_v[0]
            rd_i = rdi_v[0]
            rd_o = rdo_v[0]
            l_weld = i
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

        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)
        ax.plot(xd_v, (np.array(b_vec) - np.array(b_vec_0)) / (np.array(b_vec_0)), marker='o', markersize=2,
                label=r'$l_d=$' + str(i) + r'$m$')
        ax.set_xlabel(r'$x[m]$', fontsize=20)
        ax.set_ylabel(r'$(\tilde{B}_f-\tilde{B}_0)/\tilde{B}_0$', fontsize=20)
        ax.legend(loc='upper right', fontsize=10)
        plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
        # plt.title(r'$\varphi_d=$'+str(phi_c_weld/np.pi)+r'$\pi$, '+r'$d_d=$'+str(d_weld)+r'$m$'+r', $l_d=$'+str(
        # l_weld)+r'$m$')

        plt.title(
            r'$\varphi_d=$' + str(phi_c_weld / np.pi) + r'$\pi$, ' + r'$f=$' + str(freq) + r'$Hz$' + r', $d_d=$' + str(
                d_weld) + r'$m$')

        plt.grid(linewidth=0.5)
        # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
        # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()

    plt.show()


def b_at_det_plot_f8(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):
    bx_vecs_ideal = []

    b_mat_ideal = np.zeros((len(xd_v), len(phi_cv), len(f_v)))

    for i in phi_cv:
        bx_vec_ideal = []

        for j in f_v:
            f_ind = f_v.index(j)

            # updating output params
            x_detect = xd_v[0]
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = j

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            # updating damage params
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = 0
            d_weld = 0.0
            rd_i = ri_t
            rd_o = ro_t
            l_weld = lw_v[0]
            sigma_damage = 0.001 * 1 / freq

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
            bx_vec_ideal.append(ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                           node_dens_vec, l_sub_vec, tube_segment_lists, [], [])[0])

        bx_vecs_ideal.append(bx_vec_ideal)

    bx_vecs = []
    for i in phi_cv:
        bx_vec = []

        for j in f_v:
            f_ind = f_v.index(j)

            x_detect = xd_v[0]
            y_detect = yd_v[0]
            z_detect = zd_v[0]

            output_params = [0, x_detect, y_detect, z_detect]

            #####################################################################

            ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

            freq = j

            params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

            #####################################################################
            rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

            phi_c_weld = i
            d_weld = dw_v[0]
            rd_i = rdi_v[0]
            rd_o = rdo_v[0]
            l_weld = lw_v[0]
            sigma_damage = 0.001 * 1 / freq

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
            bx_vec.append(ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                            node_dens_vec, l_sub_vec, tube_segment_lists, [], [])[0])
        bx_vecs.append(bx_vec)

    fig = plt.figure(figsize=(3, 2))

    ax = fig.add_subplot(1, 1, 1)

    for i in phi_cv:
        phi_ind = phi_cv.index(i)

        plot_arr = (np.array(bx_vecs[phi_ind]) - np.array(bx_vecs_ideal[phi_ind])) / (np.array(bx_vecs_ideal[phi_ind]))
        # print(np.shape(np.array(bx_vecs[phi_ind])))
        # print(np.shape(np.array(bx_vecs_ideal[phi_ind])))
        ax.plot(f_v, plot_arr, marker='o', markersize=2,
                label=r'$\varphi/\pi=$' + str(i / np.pi))

    # for i in f_v:
    #    f_ind = f_v.index(i)
    #    ax.plot(xd_v, np.array(bx_vecs_ideal[f_ind])+f_ind, marker='o', markersize=2, label=r'ideal $f=$' + str(i) + r'Hz')

    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useOffset=None, useLocale=None, useMathText=None)

    ax.set_xlabel(r'$f[Hz]$', fontsize=20)
    ax.set_ylabel(r'$(\tilde{B}^d_x-\tilde{B}^i_x)/\tilde{B}_x^i$', fontsize=20)
    ax.set_xscale('log', base=2)
    ax.legend(loc='upper right', fontsize=10)
    plt.subplots_adjust(left=0.15, bottom=0.15, right=0.95, top=None, wspace=None, hspace=None)
    plt.title(r'$d_d=$' + str(d_weld) + r'$m$' + r', $x_{det}=$' + str(xd_v[0]) + r'$m$')
    # plt.grid(linewidth=0.5)
    # fig.savefig("2_2_2.pdf", dpi=gcf().dpi)
    # plt.savefig('plot.pdf', dpi=fig.dpi, bbox_inches='tight')
    # plt.tight_layout()
    plt.grid(True)

    with open('bvec.npy', 'wb') as f1:

        np.save(f1, np.array(bx_vecs))
    with open('bvec_i.npy', 'wb') as f2:

        np.save(f2, np.array(bx_vecs_ideal))

    plt.show()


def b_at_det_plot_f9(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                     params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec):

    b_mat_ideal = np.zeros((len(phi_cv),len(xd_v), len(f_v)))

    b_mat = np.zeros((len(phi_cv),len(xd_v), len(f_v)))

    for k in xd_v:
        x_ind = xd_v.index(k)
        for i in phi_cv:
            phi_ind = phi_cv.index(i)
            for j in f_v:
                f_ind = f_v.index(j)

                # updating output params
                x_detect = k
                y_detect = yd_v[0]
                z_detect = zd_v[0]

                output_params = [0, x_detect, y_detect, z_detect]

                ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

                freq = j

                params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

                # updating damage params
                rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

                phi_c_weld = 0
                d_weld = 0.0
                rd_i = ri_t
                rd_o = ro_t
                l_weld = lw_v[0]
                sigma_damage = 0.001 * 1 / freq

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

                b_mat_ideal[phi_ind, x_ind,  f_ind] = ocb.b_at_detector(params, damage_params, filament_params, output_params,
                                                                 r_sub_vec, node_dens_vec, l_sub_vec,
                                                                 tube_segment_lists, [], [])[0]

    for k in xd_v:
        x_ind = xd_v.index(k)
        for i in phi_cv:
            phi_ind = phi_cv.index(i)

            for j in f_v:
                f_ind = f_v.index(j)

                x_detect = k
                y_detect = yd_v[0]
                z_detect = zd_v[0]

                output_params = [0, x_detect, y_detect, z_detect]

                #####################################################################

                ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

                freq = j

                params = [ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq]

                #####################################################################
                rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

                phi_c_weld = i
                d_weld = dw_v[0]
                rd_i = rdi_v[0]
                rd_o = rdo_v[0]
                l_weld = lw_v[0]
                sigma_damage = 0.001 * 1 / freq

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
                b_mat[ phi_ind, x_ind, f_ind] = \
                    ocb.b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec,
                                      node_dens_vec, l_sub_vec, tube_segment_lists, [], [])[0]

    with open('bvec.npy', 'wb') as f1:

        np.save(f1, np.array(b_mat))

    with open('bvec_i.npy', 'wb') as f2:

        np.save(f2, np.array(b_mat_ideal))
