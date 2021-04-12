import numpy as np
import matplotlib.pyplot as plt
#import empit_base
#import empit_lab

#import code.observable_plotters as op

#import code.geometry_plotters_3D as gp3D
#import code.geometry_plotters_2D as gp2D
#import code.FH_run_ZC as fhzc
#import code.ZC_output_helpers as ohz
#import code.obs_calc_ZC as ocz
#import code.observable_plotter_ZC as opZC
import code.data_load_empit as dle
#import code.experiment_notebooks as en

#import code.geometry_helpers as gh

#import code.geometry_builders as gb
np.set_printoptions(linewidth=200)

# TODO include effective frequency
# TODO check where mu_eff definition enters the code - guess: only in writer
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
det_pos = np.zeros((n_det, 3))

det_pos[:, 0] = 0 * np.ones(n_det)
det_pos[:, 1] = np.linspace(-w_det / 2, w_det / 2, n_det)
det_pos[:, 2] = 1.0 * np.ones(n_det)

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
    #exp_config = {'experiment #': 8, 'rot angle index': 3, 'sensor index': 4,
    #              'det_array': 1, 'sensor_pos': None, 'data type': 'AC', 'Object': 'test'}


    #exp_config = {'experiment #': 1, 'rot angle index': 3, 'sensor index': 5,
    #              'det_array': 1, 'sensor_pos': None, 'data type': 'AC', 'Object': 'test'}  # choose experiment config and output

    # CONFIG FOR THE EXP 1 and long back det array
    #exp_config = {'experiment #': 1, 'rot angle index': 3, 'sensor index': 5, 'det_array': 1,
    #              'data type': 'DC', 'sensor_pos': None, 'Object': 'Cube'}
    # CONFIG FOR THE EXP 1 and short front det array
    #exp_config = {'experiment #': 1, 'rot angle index': 3, 'sensor index': 3, 'det_array': 2,
    #              'data type': 'DC', 'sensor_pos': None, 'Object': 'Cube'}

    # CONFIG FOR THE EXP 2 and long back array
    #exp_config = {'experiment #': 2, 'rot angle index': 3, 'sensor index': 0, 'data type': 'AC',
    #              'det_array': 1, 'sensor_pos': None, 'Object': 'GPS'}

    # CONFIG FOR THE EXP 3 and short front array
    #exp_config = {'experiment #': 3, 'rot angle index': 3, 'sensor index': 2,'data type': 'DC',
    #              'det_array': 3, 'sensor_pos': None, 'Object': 'Battery'}  # choose experiment config and output

    #CONFIG FOR THE EXP 3 and long back array
    #exp_config = {'experiment #': 3, 'rot angle index': 3, 'sensor index': 5,'data type': 'DC',
    #              'det_array': 1, 'sensor_pos': None, 'Object': 'Battery'}  # choose experiment config and output

    # CONFIG FOR THE EXP 4 and long back array
    #exp_config = {'experiment #': 4, 'rot angle index': 3, 'sensor index': 4, 'data type': 'AC',
    #              'det_array': 1, 'sensor_pos': None, 'Object': 'Laptop'}  # choose experiment config and output

    #exp_config = {'experiment #': 4, 'rot angle index': 3, 'sensor index': 2,
    #              'det_array': 2, 'sensor_pos': None, 'data type': 'DC', 'Object': 'Laptop'}  # choose experiment config and output


    # CONFIG FOR THE EXP 5 and long back array
    #exp_config = {'experiment #': 5, 'rot angle index': 3, 'sensor index': -1, 'data type': 'AC',
    #              'det_array': 1, 'sensor_pos': None, 'Object': 'GPS'}

    # CONFIG FOR THE EXP 6 and short arrays
    #exp_config = {'experiment #': 6, 'rot angle index': 3, 'sensor index': 2, 'data type': 'DC',
    #              'det_array': 3, 'sensor_pos': None, 'Object': 'Battery'}  # choose experiment config and output

    # CONFIG FOR THE EXP 6 and long array
    exp_config = {'experiment #': 6, 'rot angle index': 3, 'sensor index': 4,'data type': 'DC',
                  'det_array': 1, 'sensor_pos': None, 'Object': 'Battery'}  # choose experiment config and output


    #dle.load_exp_data_active(exp_config)  # load experimental data and plot it

    #dle.load_exp_data_active_DC(exp_config)


    dle.load_exp_data(exp_config)

    # simulation data

    # b_at_det = ocz.b_det_f_Z(phys_params, det_pos, w_l, h_l, det_loop_fil_params, geo_objects)
    plt.show()








if __name__ == '__main__':
    main()









#def main():

#TODO check if minimal simulations work (the FH simulations)
#TODO why is the old observable plotter not working?
#    exp_config = {'experiment #': 9, 'rot angle index': 3, 'sensor index': 4}  # choose experiment config and output
#
#    # build the wire
#    geo_objects = {"det_loops": []}
#
#    wire_fil_params = {"width": 0.1, "height": 0.1, "width subs": 1, "height subs": 1, "conductivity": 1.0 * 10 ** 6}
#
#    wire_build_params = {"start_point": np.array([-l_wire / 2, 0.0, -0.5]),
#                         'end_point': np.array([l_wire / 2, 0.0, -0.5]),
#                         "node count": 4,
#                         'filament parameters': wire_fil_params}
#
#    gb.wire_builder(wire_build_params, geo_objects)
#
#    gb.det_loop_builder(det_pos, 0, w_l, h_l, det_loop_fil_params, geo_objects["det_loops"])  # build the detectors
#
#    dle.load_exp_data(exp_config)  # load experimental data and plot it
#
#    en.exp_builder(geo_objects, exp_config)
#
#    # 3D visualization
#    wire_center = (wire_build_params["start_point"] + wire_build_params["end_point"]) / 2
#    viso_point = list(wire_center)
#    viso_dist = 1.5
#    gp3D.ZC_viso(geo_objects, viso_point, viso_dist)
#
#    # simulation data
#
#    # b_at_det = ocz.b_det_f_Z(phys_params, det_pos, w_l, h_l, det_loop_fil_params, geo_objects)
#    plt.show()