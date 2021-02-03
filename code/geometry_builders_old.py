import numpy as np
import code.geometry_helpers as gh


# from itertools import product

# dyn specific code

def node_points_filaments_dyn(ro, ri, params, damage_params, filament_params,
                              segment_list, node_list, r_sub_vec, node_dens_vec, l_sub_vec, shell_ind):
    """generates the node points and filaments for a shell
    """

    # extract physical parameters and 
    Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq = params
    nhinc, nwinc, ndec, units, sub_div_auto = filament_params

    # extract the parameters to model damaged parts
    rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

    del_phi_weld = d_weld / (Ro)
    phi_s = (phi_c_weld - 0.5 * del_phi_weld) % (2 * np.pi)
    phi_e = (phi_c_weld + 0.5 * del_phi_weld) % (2 * np.pi)

    # if we are in one of the surface rings
    if ro == Ro or ri == Ri:
        # CAUTION: THIS CODE IS NOT GENERAL, for now only one inner ring works ok
        if ri == Ri:
            node_dens = node_dens_vec[0]
            # if Ri==0:
            #    node_dens=node_dens_vec[1]
        if ro == Ro:
            node_dens = node_dens_vec[-1]

        pts = list(gh.points_in_ring_dyn(Ro, ro, ri, node_dens))

        # define the cuboids representing the filaments
        for i in pts:
            fil_ind = pts.index(i)
            for k in range(len(l_sub_vec) - 1):
                sub_ind = k
                Ns = i + (l_sub_vec[k],)
                x = i[0]
                y = i[1]

                Ne = i + (l_sub_vec[k + 1],)
                w = Ro / node_dens
                h = Ro / node_dens

                phi_i = gh.phi_def(x, y)

                name = [shell_ind, fil_ind, sub_ind]

                if phi_e > phi_i > phi_s:
                    gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc,
                              nwinc, sigma_damage, name)  # defines the cuboids
                    # and appends them to list
                elif phi_s > phi_e > phi_i and phi_i >= 0:
                    gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc, nwinc, sigma_damage, name)
                elif phi_s > phi_e and phi_i > phi_s and phi_i < 2 * np.pi:
                    gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc, nwinc, sigma_damage, name)
                else:
                    gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc, nwinc, sigma, name)
    # if we are in the inner ring
    else:

        pts_lists = gh.ring_cover_classify(Ro, r_sub_vec, node_dens_vec)

        all_pts = []
        for i in range(len(pts_lists)):
            all_pts.extend(pts_lists[i])
        fil_ind = 0
        for pts in pts_lists:

            j = pts_lists.index(pts)
            w = Ro / node_dens_vec[j]
            h = w

            # define the cuboids representing the filaments
            for i in pts:

                for k in range(len(l_sub_vec) - 1):
                    sub_ind = k
                    Ns = i + (l_sub_vec[k],)

                    x = i[0]
                    y = i[1]
                    Ne = i + (l_sub_vec[k + 1],)

                    sub_ind = k

                    name = [shell_ind, fil_ind, sub_ind]

                    phi_i = gh.phi_def(x, y)

                    if phi_e > phi_i > phi_s:
                        gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc,
                                  nwinc, sigma_damage, name)  # defines the cuboids
                        # and appends them to list
                    elif phi_s > phi_e > phi_i >= 0:
                        gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc, nwinc, sigma_damage, name)
                    elif phi_e < phi_s < phi_i < 2 * np.pi:
                        gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc, nwinc, sigma_damage, name)
                    else:
                        gh.cuboid(segment_list, node_list, Ns, Ne, w, h, nhinc, nwinc, sigma, name)
                fil_ind = fil_ind + 1

        pts = all_pts

    return pts


def tube_builder(Ro, Ri, r_sub_vec, l_sub_vec, node_dens_vec, params, damage_params, filament_params, tube_node_lists,
                 tube_segment_lists, tube_pts_lists):
    # tube_node_lists, tube_segment_lists, tube_pts_lists = tube_lists    #unpacking

    Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq = params
    rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params
    nhinc, nwinc, ndec, units, sub_div_auto = filament_params

    for i in range(len(r_sub_vec) - 1):
        shell_segment_list = []
        shell_node_list = []
        # build the segments and nodes from the pts in the shell,
        # later get rid of the nodes list if possible as it is a doubling
        ro = r_sub_vec[i + 1]
        ri = r_sub_vec[i]
        node_dens = node_dens_vec[i]

        shell_ind = i

        shell_pts_list = node_points_filaments_dyn(ro, ri, params, damage_params, filament_params, shell_segment_list,
                                                   shell_node_list, r_sub_vec, node_dens_vec, l_sub_vec, shell_ind)

        tube_node_lists.append(shell_node_list)
        tube_segment_lists.append(shell_segment_list)
        tube_pts_lists.append(shell_pts_list)


def sub_connection(s1, s2, phi_con, l_sub_vec, r_sub_vec, sub_con_list, sub_con_node_list, tube_lists, sub_fil_params):
    """
    - s1, s2 are the indices of the requested sub connections in l_sub_vec
    and thus determine where along the tube axis the nodes are ("z coordinate")
    - phi_con is the polar angle of connection of s1
    s2 is located at phi_con + pi
    this function finds the nodes of the tube which are closest to the desired 
    sub connection nodes and creates a cuboid between them
    """
    """PROBLEMS: shell index of the nodes for the subcon is wrong if inner
        shell is empty
    """

    ws, hs, nhinc_s, nwinc_s, sigma_s = sub_fil_params

    xs_pre = r_sub_vec[0] * np.cos(phi_con)  # requested start node
    ys_pre = r_sub_vec[0] * np.sin(phi_con)
    zs = l_sub_vec[s1]

    xe_pre = r_sub_vec[0] * np.cos(phi_con + np.pi)  # requested end node
    ye_pre = r_sub_vec[0] * np.sin(phi_con + np.pi)
    ze = l_sub_vec[s2]
    tube_node_lists, tube_segment_lists, tube_pts_lists = tube_lists

    connection_candidates_e = []
    connection_candidates_s = []

    for tube_pts_list in tube_pts_lists:  # determine the closest
        # nodes of the tube
        shell_ind = tube_pts_lists.index(tube_pts_list)

        if len(tube_pts_list) > 0:
            tube_pts_array = np.array(tube_pts_list)
            tube_x_pts = tube_pts_array[:, 0]
            tube_y_pts = tube_pts_array[:, 1]

            dist_shell_s = gh.dist_2D(xs_pre, ys_pre, tube_x_pts, tube_y_pts)
            dist_shell_e = gh.dist_2D(xe_pre, ye_pre, tube_x_pts, tube_y_pts)

            min_dist_ind_s = np.argmin(dist_shell_s)
            min_dist_ind_e = np.argmin(dist_shell_e)

            connection_candidates_s. \
                append([dist_shell_s[min_dist_ind_s],
                        tube_pts_array[min_dist_ind_s, :], min_dist_ind_s])

            connection_candidates_e. \
                append([dist_shell_e[min_dist_ind_e],
                        tube_pts_array[min_dist_ind_e, :], min_dist_ind_e])

            if len(connection_candidates_s) == 2:
                break

    con_s = (np.array(connection_candidates_s, dtype=object))
    node_ind_s = np.argmin(con_s[:, 0])
    node_xy_s = con_s[node_ind_s, 1]
    fil_ind_s = con_s[node_ind_s, 2]

    con_e = (np.array(connection_candidates_e, dtype=object))
    node_ind_e = np.argmin(con_e[:, 0])
    node_xy_e = con_e[node_ind_e, 1]
    fil_ind_e = con_e[node_ind_e, 2]

    Ns = (node_xy_s[0], node_xy_s[1], zs)
    Ne = (node_xy_e[0], node_xy_e[1], ze)

    shell_ind_s = node_ind_s
    shell_ind_e = node_ind_e

    name = [shell_ind_s, fil_ind_s, s1, shell_ind_e, fil_ind_e, s2]  # indices for the node names

    gh.cuboid(sub_con_list, sub_con_node_list, Ns, Ne, ws, hs, nhinc_s, nwinc_s, sigma_s,
              name)  # create the subcon segment
