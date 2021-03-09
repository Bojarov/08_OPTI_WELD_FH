##########################################################################
#
# Program writes input file for FASTHENRY (FH), runs simulation inside FH
# and collects output for a tube made of cuboids.
# Output is visualized within Python.
#
##########################################################################
import numpy as np
import matplotlib.pyplot as plt
# import code.obs_calc_j as ocj
import code.geometry_plotters_2D as gp2D
import code.observable_plotters as op

# Wire/Pipe specifications

Ro = 0.1  # outer radius
pipe_thickness = 0.0065
# pipe_thickness = 0.1
Ri = 0.0001  # inner radius
Ri = Ro - pipe_thickness
flen = 20  # length of wire/pipe
z0 = 0.0  # start of tube on z axis
node_dens = 80  # number of nodes on radius

# Physical FH parameters

units = "M"  # chose from km, m ,cm, mm, um, in , mils
sigma = 10 ** 6  # specify conductivity in 1/(Units*Ohms),
mu_0 = 4 * np.pi * 10 ** (-7)
mu_r = 10 ** 3  # 100.0
freq = 10 ** (0)  # #frequency of the AC
mu_r_max = 10 ** 3
max_freq = 10 ** 2

ndec = 1  # how many sub frequencies per decade should be calculated
nhinc = 1  # how many sub filaments in height of segment
nwinc = 1  # how many sub filaments in width of segment
sub_div_auto = 0  # FH internal auto subdivision of filaments off=0, on=1

min_skin_depth = np.sqrt(1 / (np.pi * sigma * mu_0 * mu_r_max * max_freq))
max_skin_depth = np.sqrt(1 / (np.pi * sigma * mu_0 * mu_r_max * 1))
skin_depth = np.sqrt(1 / (np.pi * sigma * mu_0 * mu_r * freq))

# print(skin_depth)
# min_res = Ro / (12 * node_dens)
max_res = Ro / node_dens

print("min_skin_depth")
print(min_skin_depth)
print("max_skin_depth")
print(max_skin_depth)
print("skin_depth")
print(skin_depth)
print("min_grid_width/min_skin_depth")
print(max_res / min_skin_depth)

# rc1=Ri+0.1*(Ro)                  #good for tubes
# rc2 = Ro-0.1*(Ro)

rc1 = Ri + min_skin_depth  # *(Ro-Ri)                  #good for tubes
rc2 = Ro - min_skin_depth  # 9/10*(Ro-Ri)
# rc1 = Ri + min_res

# rc1 = Ri+0.1*(Ro)                  #good for tubes
# rc2 = Ro-0.1*(Ro)

# damage parameters
rd_o = Ro
rd_i = Ri
d_weld = 0.0
l_weld = flen  # 0.1
phi_c_weld = np.pi / 2
# CAUTION in order to sim ferromagnets we use effective frequencies
# to mimics the materials
# freq_int = freq*mu_r
sigma_damage = 0.001 * 1 / freq  # therefore sigma in damage part has to be kept
# at same value when changing the frequency


r_sub_vec = np.array([Ri, rc1, rc2, Ro])  # vector of inner shell boundaries
# print(r_sub_vec)
# exit()

node_dens_vec = np.array(
    [int(node_dens * 1 / 2), int(node_dens * 1 / 2), int(node_dens)])  # vector of node densities in shells

# l_sub_vec = np.array([0.0, (flen-l_seg)/2, (flen+l_seg)/2, flen])              #vector of length wise internal
# subdivisions of the tube
# l_sub_vec = np.array([0.0, flen/4, flen*3/4, flen])

if l_weld == flen:
    l_sub_vec = np.array([0.0, flen])
else:
    l_sub_vec = np.array([0.0, 0.5 * (flen - l_weld), 0.5 * (flen + l_weld), flen])

# l_sub_vec = np.array([0.0, flen])

loop_list = []
sub_con_list = []
sub_con_node_list = []

###############################################################################
#
#   B-FIELD OUTPUT
#
###############################################################################
# Output specifications for the B-field plot

z_curr = z0  # flen*0.5
x_detect = 1.0  # width of the detector array
y_detect = 2.0  # height detector position above tube
z_detect = flen * 0.5  # detector position along tube axis
# Npts_detect = 1                                                                #amount of detector points


# packing the params
params = [Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq]
damage_params = [rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage]
filament_params = [nhinc, nwinc, ndec, units, sub_div_auto]
output_params = [z_curr, x_detect, y_detect, z_detect]

# mode = 0  # mode= 0 only current densities are calculated
# mode= 1 current densities and field at given detector point


freq_vec = np.logspace(0, 3, base=10, num=5)
freq_vec_a = np.logspace(0, 3, base=10, num=20)


def main():
    # phi_cv = [0.0]  # , np.pi/4, np.pi/2]
    # phi_cv = [np.pi/4]#, np.pi/2]
    # phi_cv = [np.pi / 2]
    # phi_cv = [0.0, np.pi]
    phi_cv = [0.0, np.pi/4, np.pi/2, 3/4*np.pi, np.pi]
    dw_v = [0.02]
    lw_v = [flen]
    xd_v = [-0.5, -0.16667, 0, 0.16667, 0.5]
    # xd_v = [-0.5]
    # xd_v = [-0.16667]
    # xd_v = [0.0]
    # xd_v = [0.16667]
    # xd_v = [0.5]

    yd_v = [2.0]

    zd_v = [flen / 2]

    #f_v = [2, 4, 8, 16, 32, 64, 128]
    f_v = [3, 6, 12, 24, 48, 96]

    # f_v = [10**0, 100]
    thick = Ro - Ri
    rdi_v = [Ri]  # ,Ri+0.25*thick, Ri+ 0.5*thick, Ri+0.75*thick, Ro]
    rdo_v = [Ro]

    # op.b_at_det_iter( phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
    #    params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec)

    op.b_at_det_plot_f9(phi_cv, dw_v, rdi_v, rdo_v, lw_v, xd_v, yd_v, zd_v, f_v,
                        params, damage_params, output_params, filament_params, node_dens_vec, r_sub_vec)

    # gp2D.wire_mesh_2D_plot_dyn(Ro, Ri, r_sub_vec, node_dens_vec, params,
    #                          filament_params, damage_params, l_sub_vec)

    # plt.show()
    # op.I_wire_plot(freq_vec, freq_vec_a, params, damage_params, filament_params,
    #               output_params, r_sub_vec, l_sub_vec, node_dens_vec, loop_list, sub_con_list)

    # visualize the cross section of the tube
    # vis_params = [ Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, 455]             #freq is dummy argument 455

    # op.B_wire_comp_plot(params, damage_params, filament_params, r_sub_vec,
    #                    l_sub_vec, node_dens_vec, loop_list, sub_con_list)

    plt.show()


if __name__ == '__main__':
    main()
