import numpy as np
import code.FH_output_helpers as FHout


# current from FH

def i_fh_calc_dyn(tube_segment_lists, l_sub_vec, xyz0_proj_re, xyz0_proj_im):
    """Calculates the total current through the wire from current densities and 
    pipeline geometry.
    """
    # extract the widths of the different filaments as array with same order as
    # the current density

    w_plaquete_vec = np.array(tube_segment_lists[0], dtype=object)[0::(len(l_sub_vec) - 1), 2]

    for segment_list in tube_segment_lists:
        seg_list_ind = tube_segment_lists.index(segment_list)

        if seg_list_ind > 0:
            w_plaquete_vec = np.append(w_plaquete_vec,
                                       np.array(segment_list, dtype=object)[0::(len(l_sub_vec) - 1), 2])

    i_fh_re = np.dot(xyz0_proj_re[:, 5], w_plaquete_vec ** 2)
    i_fh_im = np.dot(xyz0_proj_im[:, 5], w_plaquete_vec ** 2)
    i_fh = np.sqrt(i_fh_re ** 2 + i_fh_im ** 2)
    return i_fh


def i_fh_calc(r_sub_vec, l_sub_vec, node_dens_vec, tube_segment_lists,
              loop_list, sub_con_list, params, damage_params, filament_params):
    """Scans data base for given parameters extracts corresponding current
    densities and initiates I_FH calculation.
    """

    j_dens_ind = FHout.j_dens_storage(r_sub_vec, l_sub_vec, node_dens_vec, loop_list, sub_con_list, params,
                                      damage_params, filament_params, message=False)

    currents_path = 'FH_output_files/current_densities/current_dens.npy'

    with open(currents_path, 'rb') as currents:
        currents_list = list(np.load(currents, allow_pickle=True))

    j_dens_pack = currents_list[j_dens_ind]
    j_dens_cuts = j_dens_pack[5]
    xyz0_proj_re = j_dens_cuts[2]  # current in the xy plane at z=0
    xyz0_proj_im = j_dens_cuts[5]

    i_fh = i_fh_calc_dyn(tube_segment_lists, l_sub_vec, xyz0_proj_re, xyz0_proj_im)

    return i_fh
