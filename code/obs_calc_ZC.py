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
    wf_l=det_loop_fil_params[0]

    for i in range(n_det):
        print(i)
        det_pos_i = det_pos[i, :].reshape(1, 3)
        for j in range(n_f):
            freq = freqs[j]
            for k in range(3):
                i_xyz = k
                geo_objects_local["det_loops"] = []
                gb.detector_builder(det_pos_i, i_xyz, w_l, h_l, det_loop_fil_params, geo_objects_local["det_loops"])
                phys_params_run = [phys_params["sigma"], phys_params["mu0"], phys_params["mur"], freq]
                fhz.run_FH_ZC(phys_params_run, geo_objects_local, sub_div_auto=True)
                ZC_mat = ohz.ZC_mat_extract(geo_objects_local)
                #print(ZC_mat)
                #exit()
                #AZ = inv(ZC_mat)
                #b_at_det_f[i, j, k] = 4 * np.sqrt(2) / (w_l-wf_l) * AZ[0, 1] / AZ[1, 1]
                #b_at_det_f[i, j, k] = ZC_mat[0, 1].imag / (2 * np.pi * freq * w_l ** 2)
                L12 = (ZC_mat[0, 1] / (2 * np.pi * freq)).imag
                b_at_det_f[i, j, k] = L12/((w_l-wf_l)**2)           #only flux due to wire (should be the real field?)
                #L22 = (ZC_mat[1, 1] / (2 * np.pi * freq)).imag
                #b_at_det_f[i, j, k] = (L12 + L22 * (ZC_mat[0, 1]/ZC_mat[1, 1]))/((w_l-wf_l)**2)

                norm = mu0/(2*np.pi)
    return b_at_det_f/norm
