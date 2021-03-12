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
                # TODO here we define a det loop only because we give the program det positions
                phys_params_run = [phys_params["sigma"], phys_params["mu0"], phys_params["mur"], freq]
                fhz.run_FH_ZC(phys_params_run, geo_objects_local, sub_div_auto=True)
                ZC_mat = ohz.ZC_mat_extract(geo_objects_local)
                tol = 10 ** (-16)
                #ZC_mat.real[abs(ZC_mat.real) < tol] = 0.0
                #ZC_mat.imag[abs(ZC_mat.imag) < tol] = 0.0
                Y = inv(ZC_mat)
                L12 = (ZC_mat[0, 1] / (2 * np.pi * freq)).imag          #mutual inductance wire, det loop
                L23 = (ZC_mat[1, 2] / (2 * np.pi * freq)).imag
                Y31 = Y[2, 0]
                #Y21 = Y[2, 1]
                Y11 = Y[0, 0]
                # print(L23)
                # B_pass_loops
                b_at_det_f[i, j, k] = (L12 + L23 *(Y31/Y11))/((w_l-wf_l)**2)
                #print(np.angle(Y[0, 0]) / np.pi, np.angle(Y[0, 1]) / np.pi, np.angle(Y[0, 1] / Y[0, 0]) / np.pi,
                #      (np.angle(Y[0, 1]) - np.angle(Y[0, 0])) / np.pi)
                # print(np.angle(Y[0, 0]), np.angle(Y[2, 0]), np.angle(Y[2, 0] / Y[0, 0] / (np.pi)))
                # exit()
                # print(Y[0, 0], Y[0, 2], np.angle(Y[2, 0] / Y[0, 0])/(np.pi))

                # TODO current in both loops should be same due to symmetry! check how this is done!
                # print('')
                # print(Y[1,0]/Y[0,0])
                # b_at_det_f[i, j, k] = L12/((w_l-wf_l/2)**2)           #only flux due to wire (should be the real field?)
                # b_at_det_f[i, j, k] = (L12+L23*(Y31/Y11))/((w_l-wf_l/2)**2)           #only flux due to wire (should be the real field?)
                # L22 = (ZC_mat[1, 1] / (2 * np.pi * freq)).imag
                # b_at_det_f[i, j, k] = (L12 + L22 * (ZC_mat[0, 1]/ZC_mat[1, 1]))/((w_l-wf_l)**2)

                norm = mu0 / (2 * np.pi)
                #print(b_at_det_f/norm)
    return b_at_det_f /norm
