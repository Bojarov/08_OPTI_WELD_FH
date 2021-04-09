import empit_base as eb
import numpy as np
import code.observable_plotters_ZC as opZC
import code.fit_funcs as ff
from scipy.optimize import curve_fit


def load_exp_data(exp_config):
    i_alpha, i_y, i_exp = exp_config['rot angle index'], exp_config['sensor index'], exp_config['experiment #']

    data_path = '/home/buxbaum/EMPIT_DATA/Plane/eddy_currents_' + str(int(i_exp)) + '/'
    c_data_file_name = 'data_c.eyc'
    b_data_file_name = 'b_field_ac_vector_cplx.eyv'

    sa_c = eb.SensorArrayCalVector(filename=data_path + c_data_file_name)
    b_field_data = eb.DataVector(filename=data_path + b_data_file_name)
    freqs = 2 ** b_field_data.dim_values[1]
    det_pos_m = sa_c.sensor_coords_vector[0:9, :]  # only use sensor values of first sensor array
    b_field_array = b_field_data['vector'][:, 0:9, :, :]
    b_at_det_ref = (abs(b_field_array[0:3, :, :, :]).sum(axis=0) + abs(b_field_array[-3:, :, :, :]).sum(
        axis=0)) / 6  # average over reference

    b_at_det_angles = b_field_array[3:-3, :, :, :]

    b_at_det_ref = np.repeat(b_at_det_ref[np.newaxis, :, :, :], np.shape(b_at_det_angles)[0], axis=0)

    y_det = np.array(det_pos_m[:, 1])

    i_f = np.zeros(len(freqs))
    # print(np.shape(i_f))
    # mu0 = 4 * np.pi*10**(-7)
    for i in range(len(i_f)):
        by = b_at_det_ref[0, :, i, 1]
        p_opt, _ = curve_fit(ff.f_b_field, y_det, by)
        i_f[i] = p_opt[0]

    b_at_det_x = b_at_det_angles[i_alpha, :, :, :]
    b_at_det_x_ref = b_at_det_ref[i_alpha, :, :, :]

    b_at_det_pack_x = [b_at_det_x / i_f[np.newaxis, :, np.newaxis], b_at_det_x_ref / i_f[np.newaxis, :, np.newaxis]]

    opZC.bfz_plot_measure_rel_x(b_at_det_pack_x, freqs, det_pos_m, i_alpha, exp_config)

    b_at_det_alpha = b_at_det_angles[:, i_y, :, :]
    b_at_det_alpha_ref = b_at_det_ref[:, i_y, :, :]

    b_at_det_pack_alpha = [b_at_det_alpha / i_f[np.newaxis, :, np.newaxis],
                           b_at_det_alpha_ref / i_f[np.newaxis, :, np.newaxis]]
    opZC.bfz_plot_measure_rel_alpha(b_at_det_pack_alpha, freqs, det_pos_m, i_y, exp_config)


def load_exp_data_active(exp_config):
    i_alpha, i_y, i_exp, i_array = exp_config['rot angle index'], exp_config['sensor index'], \
                                   exp_config['experiment #'], exp_config['det_array']

    data_path = '/home/buxbaum/EMPIT_DATA/Active/' + str(int(i_exp)) + '/'
    c_data_file_name = 'data_c.eyc'
    b_data_file_name = 'b_field_ac_vector_cplx.eyv'

    sa_c = eb.SensorArrayCalVector(filename=data_path + c_data_file_name)
    b_field_data = eb.DataVector(filename=data_path + b_data_file_name)
    freqs = 2 ** b_field_data.dim_values[1]
    det_pos_m_9 = sa_c.sensor_coords_vector[0:9, :]  # use sensor values of first sensor array
    det_pos_m_5f = sa_c.sensor_coords_vector[14:, :]
    # print(det_pos_m_9)
    # print(det_pos_m_5f)
    # print(np.shape(b_field_data['vector']))
    # exit()

    # measurement number, sensor number, frequency, component

    if i_array == 1:
        b_field_array = b_field_data['vector'][:, 0:9, :, :]
        det_pos_array = sa_c.sensor_coords_vector[0:9, :]
    elif i_array == 2:
        b_field_array = b_field_data['vector'][:, 14:, :, :]
        det_pos_array = sa_c.sensor_coords_vector[14:, :]
    elif i_array == 3:
        b_field_array = b_field_data['vector'][:, 9:14, :, :]
        det_pos_array = sa_c.sensor_coords_vector[9:14, :]

    exp_config['sensor_pos'] = det_pos_array
    print(exp_config['sensor_pos'][exp_config['sensor index']])

    b_at_det_measurement = b_field_array[3:-3, :, :, :]
    # print(np.shape(b_at_det_measurement))
    # exit()
    print(b_at_det_measurement)

    b_at_det_ref = (abs(b_field_array[0:3, :, :, :]).sum(axis=0) + abs(b_field_array[-3:, :, :, :]).sum(
        axis=0)) / 6  # average over reference
    # b_at_det_ref = abs(b_field_array[1, :, :, :])
    b_at_det_ref = np.repeat(b_at_det_ref[np.newaxis, :, :, :], np.shape(b_at_det_measurement)[0], axis=0)

    if exp_config['experiment #'] == 1:
        z_obj = np.array([0.0, 0.12, 0.24, 0.47, 0.60, 0.75, 1.02, 1.4]) + 0.064

        b_at_det_pack = [abs(b_at_det_measurement[1:, :, :, :]) - b_at_det_ref[1:, :, :, :], z_obj]

    if exp_config['experiment #'] == 2:
        z_obj = np.array([-0.05, 0.01, 0.12, 0.335, 0.45, 0.63, 0.74, 1.02, 1.435]) + 0.05
        b_at_det_pack = [abs(b_at_det_measurement[2:, :, :, :]) - b_at_det_ref[2:, :, :, :], z_obj]

    if exp_config['experiment #'] == 3:
        x_obj = np.array([-0.06, -0.30, -0.50, -0.70, -1.0, -1.3])
        b_at_det_pack = [abs(b_at_det_measurement[1:, :, :, :]) - b_at_det_ref[1:, :, :, :], x_obj]

    if exp_config['experiment #'] == 4:
        x_obj = np.array([-0.06, -0.30, -0.50, -0.70, -1.0, -1.3])
        b_at_det_pack = [abs(b_at_det_measurement) - b_at_det_ref, x_obj]
        # print(len(x_obj))
        # print(np.shape(b_at_det_pack[0]))
        # exit()

    # bfz_plot_measure_active(b_at_det_pack, freqs, det_pos, alpha_ind, exp_config)
    opZC.bfz_plot_measure_active(b_at_det_pack, freqs, exp_config)


def load_exp_data_active_DC(exp_config):
    i_alpha, i_y, i_exp, i_array = exp_config['rot angle index'], exp_config['sensor index'], \
                                   exp_config['experiment #'], exp_config['det_array']

    data_path = '/home/buxbaum/EMPIT_DATA/Active/' + str(int(i_exp)) + '/'

    c_data_file_name = 'data_c.eyc'

    if exp_config['data type'] == 'AC':

        b_data_file_name = 'b_field_ac_vector_cplx.eyv'
        sa_c = eb.SensorArrayCalVector(filename=data_path + c_data_file_name)
        b_field_data = eb.DataVector(filename=data_path + b_data_file_name)
        b_field_array_full = b_field_data['vector']
        freqs = 2 ** b_field_data.dim_values[1]

    else:
        b_data_file_name = 'b_field_dc_vector.eyv'
        sa_c = eb.SensorArrayCalVector(filename=data_path + c_data_file_name)
        b_field_data = eb.DataVector(filename=data_path + b_data_file_name)
        b_field_array_full = b_field_data['sensor_signal']
        b_field_array_full = np.expand_dims(b_field_array_full, axis=2)
        freqs = np.array([0])

    # det_pos_m_9 = sa_c.sensor_coords_vector[0:9, :]  # use sensor values of first sensor array
    # det_pos_m_5f = sa_c.sensor_coords_vector[14:, :]
    # print(det_pos_m_9)
    # print(det_pos_m_5f)

    # measurement number, sensor number, frequency, component

    if i_array == 1:
        b_field_array = b_field_array_full[:, 0:9, :, :]
        det_pos_array = sa_c.sensor_coords_vector[0:9, :]
    elif i_array == 2:
        b_field_array = b_field_array_full[:, 14:, :, :]
        det_pos_array = sa_c.sensor_coords_vector[14:, :]
    elif i_array == 3:
        b_field_array = b_field_array_full[:, 9:14, :, :]
        det_pos_array = sa_c.sensor_coords_vector[9:14, :]
    else:
        print("ERROR: No valid detector array selected")

    exp_config['sensor_pos'] = det_pos_array

    b_at_det_measurement = b_field_array[3:-3, :, :, :]
    b_at_det_ref = (abs(b_field_array[0:3, :, :, :]).sum(axis=0) + abs(b_field_array[-3:, :, :, :]).sum(
        axis=0)) / 6  # average over reference

    b_at_det_ref = np.repeat(b_at_det_ref[np.newaxis, :, :, :], np.shape(b_at_det_measurement)[0], axis=0)

    if exp_config['data type'] == 'AC':
        b_at_det_object = abs(b_at_det_measurement) - b_at_det_ref

    else:
        # b_at_det_object = abs(b_at_det_measurement - b_at_det_ref)
        b_at_det_object = (b_at_det_measurement - b_at_det_ref)




    if exp_config['experiment #'] == 1:
        # CUBE
        z_obj = np.array([0.0, 0.12, 0.24, 0.47, 0.60, 0.75, 1.02, 1.4]) + 0.064

        b_at_det_pack = [b_at_det_object[1:, :, :, :], z_obj]

    if exp_config['experiment #'] == 2:
        # GPS
        z_obj = np.array([-0.05, 0.01, 0.12, 0.335, 0.45, 0.63, 0.74, 1.02, 1.435]) + 0.05
        b_at_det_pack = [b_at_det_object[2:, :, :, :], z_obj]

    if exp_config['experiment #'] == 3:
        # Battery
        x_obj = np.array([-0.06, -0.30, -0.50, -0.70, -1.0, -1.3])
        b_at_det_pack = [b_at_det_object[1:, :, :, :], x_obj]

    if exp_config['experiment #'] == 4:
        # Laptop
        x_obj = np.array([-0.06, -0.30, -0.50, -0.70, -1.0, -1.3])
        b_at_det_pack = [b_at_det_object, x_obj]


    if exp_config['experiment #'] == 5:
        # GPS + Eddycurrents
        z_obj = np.array([0.095, 0.18, 0.28, 0.41, 0.61, 0.81, 1.01, 1.3, 1.5]) + 0.05
        #b_at_det_pack = [b_at_det_object[2:, :, :, :], z_obj]
        b_at_det_pack = [b_at_det_object, z_obj]


    if exp_config['experiment #'] == 6:
        # Battery
        x_obj = np.array([-0.17, -0.32, -0.70, -1.0, -1.2, -1.4, -1.6])
        b_at_det_pack = [b_at_det_object, x_obj]

    opZC.bfz_plot_measure_active(b_at_det_pack, freqs, exp_config)
