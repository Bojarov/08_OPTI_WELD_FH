import numpy as np
import matplotlib.pyplot as plt
import code.geometry_builders as gb
import code.geometry_plotters_3D as gp3D
import code.geometry_plotters_2D as gp2D
import code.FH_run_ZC as fhzc
import code.ZC_output_helpers as ohz
import code.obs_calc_ZC as ocz
import code.observable_plotters as op
import code.observable_plotters_ZC as opZC
import code.data_load_empit as dle

np.set_printoptions(linewidth=200)

# Physical FH parameters

units = "M"  # chose from km, m ,cm, mm, um, in , mils
sigma = 10 ** 6  # specify conductivity in 1/(Units*Ohms),
mu_0 = 4 * np.pi * 10 ** (-7)
mu_r = 10 ** 0  # 100.0
# freqs = [4, 1024]#, 8, 12, 24, 32, 64, 128, 256, 512, 1024]
freqs = [4, 16, 64, 256, 512, 1024]

phys_params = {"sigma": sigma, "mu0": mu_0, "mur": mu_r, "freqs": freqs}

ndec = 1  # how many sub frequencies per decade should be calculated
nhinc = 1  # how many sub filaments in height of segment
nwinc = 1  # how many sub filaments in width of segment
sub_div_auto = 0  # FH internal auto subdivision of filaments off=0, on=1
fil_params = [nhinc, nwinc, ndec]

# wire geometry
l_wire = 30  # length of wire/pipe
p1_wire = np.array([0.0, 0.0, 0.0])
p2_wire = np.array([0.0, 0.0, l_wire])
w_wire = 0.1
h_wire = 0.1

# Detector positions
n_det = 11  # number of detector loops
w_det = 1.0  # width of the detector array
det_pos = np.zeros((n_det, 3))

det_pos[:, 0] = np.linspace(-w_det / 2, w_det / 2, n_det)
# det_pos[:, 0] = np.linspace(0, w_det / 2, n_det)
det_pos[:, 1] = 1.0 * np.ones(n_det)
det_pos[:, 2] = l_wire / 2 * np.ones(n_det)

# det_pos[:, 0] = 0.0 * np.ones(n_det)
# det_pos[:, 1] = np.linspace(0.2, 1.0, n_det)
# det_pos[:, 2] = l_wire/2 * np.ones(n_det)


# Detector loop dimensions and direction of measured field component
i_xyz = 0  # direction index 0...x, 1...y, 2...z
w_l = 0.05  # width and height of loop
h_l = 0.05

# Filaments parameters of the detector loops
sigma_l = sigma * 0.00001
wf = 0.001
hf = 0.001
det_loop_fil_params = [wf, hf, nhinc, nwinc, sigma_l]

# plane 1 definition
p1 = tuple(np.array([1, 1 + 0.1, l_wire / 2 - 0.1]))
p2 = tuple(np.array([1, 1 + 0.1, l_wire / 2 + 0.1]))
p3 = tuple(np.array([1, 1 - 0.1, l_wire / 2 + 0.1]))
m_grid = 10
sigma_p = 10 ** 6
thick = 0.005

# 3D visualization
viso_point = [0, 0, l_wire / 2]
viso_dist = 1.5


def main():
    # dict of lists to gather defined objects
    geo_objects = { "det_loops": []}
    #geo_objects = {"wires": [], "det_loops": []}


    wire_fil_params = {"width": 0.1, "height": 0.1,
                    "width subs": 1, "height subs": 1, "conductivity":  10 ** 6}

    wire_build_params = {"start_point": p1_wire, 'end_point': p2_wire, "node count": 2,
                         "filament parameters": wire_fil_params}
    gb.wire_builder(wire_build_params, geo_objects)


    #gb.wire_builder(p1_wire, p2_wire, w_wire, h_wire, phys_params, fil_params, geo_objects["wires"], external=True)

    r_fil_params = {"width": 0.05, "height": 0.002,
                    "width subs": 1, "height subs": 1, "conductivity": 2.5 * 10 ** 6}

    r_build_params = {'center_pos': np.array([0.0, 0.0, l_wire / 2]).reshape(1, 3), "yaw": 0.0 * np.pi,
                      "pitch": 0.0 * np.pi,
                      "roll": 0.0 * np.pi, 'a_r': 0.5, 'b_r': 0.5, "node count": 5, "contact distance": 10 ** (-3),
                      "filament parameters": r_fil_params}

    plane_build_params = {'center_pos': np.array([0.0, 0.0, l_wire / 2]).reshape(1, 3),
                          "yaw": 0.0 * np.pi, "pitch": 0.0 * np.pi, "roll": 0.0 * np.pi,
                          'edge a': 0.5, 'edge b': 0.5, "loop count": 5, "contact distance": 10 ** (-3),
                          "filament parameters": r_fil_params}

    gb.loop_plane_builder(plane_build_params, geo_objects)

    # gb.rectangle_loop_builder(r_build_params, geo_objects)



    # adding the detector loops (for visualization only)
    gb.det_loop_builder(det_pos, 0, w_l, h_l, det_loop_fil_params, geo_objects["det_loops"])

    # adding a passive loop
    w_p = 0.05
    h_p = 0.05
    pass_pos = np.zeros((1, 3))
    pass_pos[0, 0] = 0.51
    pass_pos[0, 1] = 1.0
    pass_pos[0, 2] = l_wire / 2 + 0.0

    alpha_p = 0.0 * np.pi
    beta_p = 0.0 * np.pi
    gamma_p = 1.0 * np.pi
    pass_loop_fil_params = [w_p, h_p, nhinc, nwinc, sigma_p]

    # gb.loop_builder(pass_pos, alpha_p, beta_p, gamma_p, w_p, h_p, pass_loop_fil_params, geo_objects["passive_loops"])

    # adding circular loop
    circ_fil_params = {"width": 0.05, "height": 0.002,
                       "width subs": 1, "height subs": 1, "conductivity": 2.5 * 10 ** 6}

    circ_build_params = {"center_pos": pass_pos, "yaw": 0.0 * np.pi, "pitch": 0.25 * np.pi, "roll": 0.0 * np.pi,
                         "radius": 0.125, "node count": 41, "contact distance": 10 ** (-3),
                         "filament parameters": circ_fil_params}

    #gb.circular_loop_builder(circ_build_params, geo_objects)

    beta_vec = np.array([np.pi / 32, np.pi / 16, np.pi / 8, np.pi / 4])
    #b_at_det = ocz.b_det_f_beta_Z(beta_vec, phys_params, det_pos, w_l, h_l, det_loop_fil_params, geo_objects)
    #b_at_det = ocz.b_det_f_Z(phys_params, det_pos, w_l, h_l, det_loop_fil_params, geo_objects)
    # print(b_at_det)
    gp3D.ZC_viso(geo_objects, viso_point, viso_dist)
    #opZC.bfz_plot(b_at_det, freqs, det_pos)
    #op.bfbz_plot(b_at_det, beta_vec, freqs, det_pos)




    #dle.plot_empit_data()

    plt.show()


if __name__ == '__main__':
    main()
