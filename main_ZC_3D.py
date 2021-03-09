import numpy as np
import matplotlib.pyplot as plt
import code.geometry_builders as gb
import code.geometry_plotters_3D as gp3D
import code.FH_run_ZC as fhzc
import code.ZC_output_helpers as ozc
# TODO include effective frequency
# Physical FH parameters

units = "M"  # chose from km, m ,cm, mm, um, in , mils
sigma = 10 ** 0  # specify conductivity in 1/(Units*Ohms),
mu_0 = 4 * np.pi * 10 ** (-7)
mu_r = 10 ** 0  # 100.0
freq = 10 ** (0)  # #frequency of the AC
mu_r_max = 10 ** 3
max_freq = 10 ** 2

phys_params = [sigma, mu_0, mu_r, freq]

ndec = 1  # how many sub frequencies per decade should be calculated
nhinc = 1  # how many sub filaments in height of segment
nwinc = 1  # how many sub filaments in width of segment
sub_div_auto = 0  # FH internal auto subdivision of filaments off=0, on=1
fil_params = [nhinc, nwinc, ndec]

# wire geometry
l_wire = 20  # length of wire/pipe
p1_wire = np.array([0.0, 0.0, 0.0])
p2_wire = np.array([0.0, 0.0, l_wire])
w_wire = 0.1
h_wire = 0.1

# Detector positions
n_det = 1  # number of detector loops
w_det = 1.0  # width of the detector array
det_pos = np.zeros((n_det, 3))
det_pos[:, 0] = np.linspace(-w_det / 2, w_det / 2, n_det)
det_pos[:, 1] = 1.0 * np.ones(n_det)
det_pos[:, 2] = l_wire / 2 * np.ones(n_det)

# Detector loop dimensions and direction of measured field component
i_xyz = 0  # direction index 0...x, 1...y, 2...z
w_l = 0.01  # width and height of loop
h_l = 0.01

# Filaments parameters of the detector loops
sigma_l = sigma
wf = 0.001
hf = 0.001
det_loop_fil_params = [wf, hf, nhinc, nwinc, sigma_l]

# plane definition
p1 = tuple(np.array([1, 1 + 0.1, 10 - 0.1]))
p2 = tuple(np.array([1, 1 + 0.1, 10 + 0.1]))
p3 = tuple(np.array([1, 1 - 0.1, 10 + 0.1]))
m_grid = 10
thick = 0.005

# visualization
viso_point = [0, 1, 10]
viso_dist = 1.5

# dict of lists to gather defined objects
geo_objects = {"wires": [], "loops": [], "det_loops": [], "planes": []}


def main():

    gb.wire_builder(p1_wire, p2_wire, w_wire, h_wire, phys_params, fil_params, geo_objects["wires"], external=True)
    #gb.detector_builder(det_pos, i_xyz, w_l, h_l, det_loop_fil_params, geo_objects["loops"])
    gb.detector_builder(det_pos, i_xyz, w_l, h_l, det_loop_fil_params, geo_objects["det_loops"])
    #gb.plane_builder(p1, p2, p3, thick, m_grid, geo_objects["planes"])

    gp3D.ZC_viso(geo_objects, viso_point, viso_dist)

    fhzc.run_FH_ZC(phys_params, geo_objects, sub_div_auto=True)
    Z_mat = ozc.ZC_mat_extract(geo_objects)
    # TODO check Z mat
    print(Z_mat)


    plt.show()


if __name__ == '__main__':
    main()
