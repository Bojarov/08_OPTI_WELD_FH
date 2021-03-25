import numpy as np
import matplotlib.pyplot as plt
import code.geometry_plotters_3D as gp3D
import code.geometry_plotters_2D as gp2D
import code.FH_run_ZC as fhzc
import code.ZC_output_helpers as ohz
import code.obs_calc_ZC as ocz
import code.observable_plotters as op
import code.data_load_empit as dle
import code.experiment_notebooks as en

np.set_printoptions(linewidth=200)

# TODO include effective frequency
# TODO check where mu_eff definition enters the code - guess: only in writer
# TODO rewrite the whole writer business - change all geo objects to dictionaries and standardize the objects properties
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
p2_wire = np.array([l_wire, 0.0, 0.0])
w_wire = 0.1
h_wire = 0.1

# Detector positions
n_det = 3  # number of detector loops
w_det = 1.3  # width of the detector array
det_pos = np.zeros((n_det, 9))

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
sigma_l = sigma * 0.00001  # TODO investigate how low conductivity of det loop can be & how result depends on it
wf = 0.001
hf = 0.001
det_loop_fil_params = [wf, hf, nhinc, nwinc, sigma_l]


def main():
    # dict of lists to gather defined objects
    geo_objects = {"wires": [], "det_loops": []}










    # 3D visualization


    viso_point = list((p2_wire-p1_wire)/2)

    viso_dist = 1.5

    #gp3D.ZC_viso(geo_objects, viso_point, viso_dist)

    dle.plot_empit_data()

    plt.show()


if __name__ == '__main__':
    main()
