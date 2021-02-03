# import numpy as np
import time
import code.obs_calc_b_field as ocb
import code.FH_run_dyn as FHrun
import code.geometry_builders as gb


def para_point_calc(params, damage_params, filament_params, output_params,
                    r_sub_vec, l_sub_vec, node_dens_vec, loop_list, sub_con_list, mode):
    # hierarchy of geometry objects
    # tube->shell->segments->nodes
    tube_node_lists = []  # lists to carry the node and segment lists for each shell
    tube_segment_lists = []
    tube_pts_lists = []

    start_time = time.time()

    # unpacking the params
    Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq = params
    rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params
    nhinc, nwinc, ndec, units, sub_div_auto = filament_params
    z0, x_detect, y_detect, z_detect = output_params

    # build the geometry in python from the input params
    gb.tube_builder(Ro, Ri, r_sub_vec, l_sub_vec, node_dens_vec,
                    params, damage_params, filament_params, tube_node_lists,
                    tube_segment_lists, tube_pts_lists)

    FHrun.run_FH_dyn(params, damage_params, filament_params, output_params,
                     r_sub_vec, node_dens_vec, l_sub_vec,
                     tube_segment_lists, tube_pts_lists,
                     loop_list, sub_con_list)

    if mode == 1:
        b_at_det = \
            ocb.b_at_detector(params, damage_params, filament_params, output_params,
                              r_sub_vec, node_dens_vec, l_sub_vec, tube_segment_lists, loop_list,
                              sub_con_list)

    # USER INFO
    no_of_fils = len(tube_segment_lists[0]) + len(tube_segment_lists[1]) \
                 + len(tube_segment_lists[2]) + 4 * len(loop_list) \
                 + len(sub_con_list)

    print("Number of basic segments: " + str(no_of_fils))
    print("Of these there are: " + str(4 * len(loop_list)) + " loop filaments")
    print("Of these there are: " + str(len(sub_con_list)) + " sub connection filaments")
    print("--- %s seconds ---" % (time.time() - start_time))
    tube_lists = [tube_node_lists, tube_segment_lists, tube_pts_lists]  # packing

    return tube_lists
