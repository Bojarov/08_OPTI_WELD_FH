import numpy as np
# import os.path
import code.obs_calc_j as ocj
import code.obs_calc_a_field as oca
import code.FH_output_helpers as FHout


# import code.obs_calc_analytical as ocana

# B field calculation

def b_point_calc(a_cross_re_z, a_cross_im_z, params):
    """Determines the tangential B-field in a point (x,y,z) 
       (field in direction of the unit vector e_phi)
       from the vector potential at that point,
       must be given for not only the point but for a cross formed plaquete
       around the point x,y,z

       CAUTION: ONLY MAGNITUDE OF B PHASOR IS CALCULATED,
       STILL UNCLEAR/NOT CLARIFIED HOW TO DEAL WITH THE PHASE

       BIT DIRTY DUE CAUSE REUSING OLD CODE BUT WORKS :)
    """
    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params  # unpack
    skin_depth = np.sqrt(2 / (2 * np.pi * freq * sigma * mu_0))

    # it is important to define the same h as for the A field calculation here!!!
    # need to check if the discritisation is fine !!!!!
    h = min(0.01 * ro_t / node_dens, 0.1 * skin_depth)

    # create 3x3 plaquete in which A field is calculated for each spatial dimension
    axyz0_re = np.zeros((3, 3, 3))
    axyz0_im = np.zeros((3, 3, 3))

    # for now fill only z :)
    # z-cross
    axyz0_re[2, 1, 2] = a_cross_re_z[0]
    axyz0_re[0, 1, 2] = a_cross_re_z[1]
    axyz0_re[1, 1, 2] = a_cross_re_z[2]
    axyz0_re[1, 2, 2] = a_cross_re_z[3]
    axyz0_re[1, 0, 2] = a_cross_re_z[4]

    axyz0_im[2, 1, 2] = a_cross_im_z[0]
    axyz0_im[0, 1, 2] = a_cross_im_z[1]
    axyz0_im[1, 1, 2] = a_cross_im_z[2]
    axyz0_im[1, 2, 2] = a_cross_im_z[3]
    axyz0_im[1, 0, 2] = a_cross_im_z[4]

    # B_z can be calculated from the cross of 5 points which are stored
    # in A-dpts_cross_xyz_re/im in the following order:
    # y               ______
    # A              |     |
    # |              |  4  |
    # |         _____|     |_____
    # |        |                 |
    # |        |  2     3     1  |
    # |        |_____      ______|
    # |              |     |
    # |              |  5  |
    # |              |_____|
    # |
    # ----------------------------------> x

    az_xy_re = axyz0_re[:, :, 2]  # real part of A_z in the z=z0 plane
    az_xy_im = axyz0_im[:, :, 2]  # im   part of A_z in the z=z0 plane

    az_tot_xy = np.sqrt(az_xy_re ** 2 + az_xy_im ** 2)

    # calculate the tangential B field in center of plaquete
    # dx_az_re = (az_xy_re[2, 1] - az_xy_re[0, 1]) / (2 * h)
    # dy_az_re = (az_xy_re[1, 2] - az_xy_re[1, 0]) / (2 * h)
    # dx_az_im = (az_xy_im[2, 1] - az_xy_im[0, 1]) / (2 * h)
    # dy_az_im = (az_xy_im[1, 2] - az_xy_im[1, 0]) / (2 * h)
    dx_az_tot = \
        (az_tot_xy[2, 1] - az_tot_xy[0, 1]) / (2 * h)
    dy_az_tot = \
        (az_tot_xy[1, 2] - az_tot_xy[1, 0]) / (2 * h)

    b_mag = np.sqrt(dx_az_tot ** 2 + dy_az_tot ** 2)
    # for now only B field magnitude is of interest
    # needs to change when symmetry is lower then cylindrical
    return b_mag


def b_at_detector(params, damage_params, filament_params, output_params, r_sub_vec, node_dens_vec, l_sub_vec,
                  tube_segment_lists, loop_list, sub_con_list):
    """Calculates B field at detector point above ground
    """

    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params  # unpack params

    # l_weld, phi_c_weld, sigma_damage = damage_params

    # z, x_detect, y_detect, z_detect = output_params

    # extract the current densities from the data base
    j_dens_ind = FHout.j_dens_storage(r_sub_vec, l_sub_vec, node_dens_vec, loop_list, sub_con_list, params,
                                      damage_params, filament_params, message=False)

    currents_path = 'FH_output_files/current_densities/current_dens.npy'
    with open(currents_path, 'rb') as currents:
        currents_list = list(np.load(currents, allow_pickle=True))

    j_dens_pack = currents_list[j_dens_ind]
    j_dens_cuts = j_dens_pack[5]
    j_tube_re = j_dens_cuts[1]  # current densities at tube segment nodes
    j_tube_im = j_dens_cuts[4]

    xyz0_proj_re = j_dens_cuts[2]  # current in the xy plane at z=0
    xyz0_proj_im = j_dens_cuts[5]

    i_fh = ocj.i_fh_calc_dyn(tube_segment_lists, l_sub_vec, xyz0_proj_re, xyz0_proj_im)
    #print("i_fh ", i_fh)
    a_dpts_cross_xyz_re_z, a_dpts_cross_xyz_im_z = oca.A_int_dyn(
        j_tube_re, j_tube_im, r_sub_vec, l_sub_vec, node_dens_vec,
        params, output_params, tube_segment_lists, loop_list)

    normalization = mu_0 * i_fh / (2 * np.pi)  # normalize by current

    a_cross_re_z = a_dpts_cross_xyz_re_z[0, :]
    a_cross_im_z = a_dpts_cross_xyz_im_z[0, :]

    b_at_det = \
        b_point_calc(a_cross_re_z, a_cross_im_z, params) / normalization

    return b_at_det
