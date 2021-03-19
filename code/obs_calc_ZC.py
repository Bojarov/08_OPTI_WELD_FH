import code.geometry_builders as gb
import code.FH_run_ZC as fhzc
import numpy as np
import code.FH_run_ZC as fhz
import code.ZC_output_helpers as ohz
from numpy.linalg import inv


def b_det_f_Z(phys_params, det_pos, w_l, h_l, det_loop_fil_params, geo_objects):
    n_det, _ = np.shape(det_pos)
    freqs = phys_params["freqs"]
    mu0 = phys_params["mu0"]
    n_f = len(freqs)
    b_at_det_f = np.zeros((n_det, n_f, 3), dtype=complex)

    geo_objects_local = geo_objects.copy()
    wf_l = det_loop_fil_params[0]

    for i in range(n_det):
        print(i)
        det_pos_i = det_pos[i, :].reshape(1, 3)
        for j in range(n_f):
            freq = freqs[j]
            for k in range(3):
                i_xyz = k

                geo_objects_local["det_loops"] = []
                gb.det_loop_builder(det_pos_i, i_xyz, w_l, h_l, det_loop_fil_params, geo_objects_local["det_loops"])
                phys_params_run = [phys_params["sigma"], phys_params["mu0"], phys_params["mur"], freq]
                fhz.run_FH_ZC(phys_params_run, geo_objects_local, sub_div_auto=True)
                ZC_mat = ohz.ZC_mat_extract(geo_objects_local)  # impedance matrix
                Y = inv(ZC_mat)  # admittance matrix
                L12 = (ZC_mat[0, 1] / (2 * np.pi * freq)).imag  # mutual inductance between wire and detector loop
                Y11 = Y[0, 0]  # admittance of the wire
                Y_pass_loops = Y[2:, 0]
                L_pass_loops = (ZC_mat[1, 2:] / (2 * np.pi * freq)).imag
                b_pass_loops = np.sum(L_pass_loops * Y_pass_loops / Y11)
                b_at_det_f[i, j, k] = (L12 + b_pass_loops) / ((w_l - wf_l) ** 2)
                norm = mu0 / (2 * np.pi)

    return b_at_det_f / norm


def b_det_f_beta_Z(beta_vec, phys_params, det_pos, w_l, h_l, det_loop_fil_params, geo_objects):
    n_det, _ = np.shape(det_pos)
    freqs = phys_params["freqs"]
    mu0 = phys_params["mu0"]
    n_f = len(freqs)
    n_beta = len(beta_vec)
    b_at_det_f = np.zeros((n_det, n_f, n_beta, 3), dtype=complex)

    geo_objects_local = geo_objects.copy()
    wf_l = det_loop_fil_params[0]

    for i in range(n_det):
        print(i)
        det_pos_i = det_pos[i, :].reshape(1, 3)
        for j in range(n_f):
            freq = freqs[j]
            for k in range(len(beta_vec)):
                for l in range(3):
                    i_xyz = l
                    geo_objects_local["det_loops"] = []

                    circ = geo_objects_local["circ_pass_loops"][0].copy()
                    circ["build parameters"]["pitch"] = beta_vec[k]
                    geo_objects_local["circ_pass_loops"] = []
                    gb.circular_loop_builder(circ["build parameters"], geo_objects_local)

                    gb.det_loop_builder(det_pos_i, i_xyz, w_l, h_l, det_loop_fil_params, geo_objects_local["det_loops"])
                    phys_params_run = [phys_params["sigma"], phys_params["mu0"], phys_params["mur"], freq]
                    fhz.run_FH_ZC(phys_params_run, geo_objects_local, sub_div_auto=True)
                    ZC_mat = ohz.ZC_mat_extract(geo_objects_local)  # impedance matrix
                    Y = inv(ZC_mat)  # admittance matrix
                    L12 = (ZC_mat[0, 1] / (2 * np.pi * freq)).imag  # mutual inductance between wire and detector loop
                    Y11 = Y[0, 0]  # admittance of the wire
                    Y_pass_loops = Y[2:, 0]
                    L_pass_loops = (ZC_mat[1, 2:] / (2 * np.pi * freq)).imag
                    b_pass_loops = np.sum(L_pass_loops * Y_pass_loops / Y11)
                    b_at_det_f[i, j, k, l] = (L12 + b_pass_loops) / ((w_l - wf_l) ** 2)
                    norm = mu0 / (2 * np.pi)

    return b_at_det_f / norm
