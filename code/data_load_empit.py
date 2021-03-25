import empit_base as eb
import numpy as np
import code.observable_plotters as op
import code.fit_funcs as ff
from scipy.optimize import curve_fit

def plot_empit_data():
    data_path = '/home/buxbaum/EMPIT_DATA/Plane/eddy_currents_02/'
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

    i_f = np.zeros(len(y_det))
    # mu0 = 4 * np.pi*10**(-7)
    for i in range(len(i_f)):
        by = b_at_det_ref[0, :, i, 1]
        p_opt, _ = curve_fit(ff.f_b_field, y_det, by)
        i_f[i] = p_opt[0]
    #print(i_f)
    #exit()

    alpha_ind = 3

    b_at_det_x = b_at_det_angles[alpha_ind, :, :, :]
    b_at_det_x_ref = b_at_det_ref[alpha_ind, :, :, :]

    b_at_det_pack_x = [b_at_det_x/i_f[:, np.newaxis, np.newaxis], b_at_det_x_ref/i_f[:, np.newaxis, np.newaxis]]

    op.bfz_plot_measure_rel_x(b_at_det_pack_x, freqs, det_pos_m, alpha_ind)

    x_ind = 4

    b_at_det_alpha = b_at_det_angles[:, x_ind, :, :]
    b_at_det_alpha_ref = b_at_det_ref[:, x_ind, :, :]
    #TODO need to normalize the alpha plots too
    b_at_det_pack_alpha = [b_at_det_alpha, b_at_det_alpha_ref]
    op.bfz_plot_measure_rel_alpha(b_at_det_pack_alpha, freqs, det_pos_m, x_ind)

    # op.testplot()
