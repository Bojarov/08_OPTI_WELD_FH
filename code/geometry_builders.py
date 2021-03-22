import numpy as np
import code.geometry_helpers as gh
from collections import deque


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

    del_phi_weld = d_weld / Ro
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
                r_node = np.sqrt(x ** 2 + y ** 2)

                Ne = i + (l_sub_vec[k + 1],)
                w = Ro / node_dens
                h = Ro / node_dens

                phi_i = gh.phi_def(x, y)

                name = [shell_ind, fil_ind, sub_ind]
                if l_weld == flen:
                    gh.damage_sort_line(r_node, rd_o, rd_i, phi_e, phi_i, phi_s, segment_list, node_list, Ns,
                                        Ne, w, h, nhinc, nwinc, sigma, sigma_damage, name)
                else:
                    gh.damage_sort_patch(r_node, rd_o, rd_i, sub_ind, phi_e, phi_i, phi_s, segment_list, node_list, Ns,
                                         Ne, w, h, nhinc, nwinc, sigma, sigma_damage, name)
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

                    r_node = np.sqrt(x ** 2 + y ** 2)

                    name = [shell_ind, fil_ind, sub_ind]

                    phi_i = gh.phi_def(x, y)

                    if l_weld == flen:
                        gh.damage_sort_line(r_node, rd_o, rd_i, phi_e, phi_i, phi_s, segment_list, node_list, Ns,
                                            Ne, w, h, nhinc, nwinc, sigma, sigma_damage, name)
                    else:
                        gh.damage_sort_patch(r_node, rd_o, rd_i, sub_ind, phi_e, phi_i, phi_s, segment_list, node_list,
                                             Ns,
                                             Ne, w, h, nhinc, nwinc, sigma, sigma_damage, name)
                fil_ind = fil_ind + 1

        pts = all_pts

    return pts


def tube_builder(Ro, Ri, r_sub_vec, l_sub_vec, node_dens_vec, params,
                 damage_params, filament_params, tube_node_lists,
                 tube_segment_lists, tube_pts_lists):
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


def sub_connection(s1, s2, phi_con, l_sub_vec, r_sub_vec, sub_con_list,
                   sub_con_node_list, tube_lists, sub_fil_params):
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
            # f_v = [0.001, 1, 4, 16, 64, 100]

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


# ZC code

def node_builder(my_func, build_params):
    alpha, beta, gamma = build_params['yaw'], build_params['pitch'], build_params['roll']
    center_pos = build_params['center_pos']
    R_mat = gh.r_ypr(alpha, beta, gamma)
    nodes_pos = my_func(build_params)
    nodes_pos = np.einsum('ij,kj->ik', R_mat, nodes_pos) + np.transpose(center_pos)
    nodes_pos = np.transpose(nodes_pos)

    return nodes_pos


def det_loop_builder(det_pos, i_xyz, w_l, h_l, loop_fil_params, loop_list):
    # if i_xyz == 0:
    alpha, beta, gamma = np.array([0, np.pi / 2, 0])
    if i_xyz == 1:
        alpha, beta, gamma = np.array([0, 0, np.pi / 2])
    elif i_xyz == 2:
        alpha, beta, gamma = np.zeros(3)

    wf, hf, nhinc_f, nwinc_f, sigma_l = loop_fil_params
    n_det, _ = np.shape(det_pos)

    for i in range(n_det):
        p = det_pos[i, :]
        loop_list.append([p, w_l, h_l, alpha, beta, gamma, sigma_l, wf, hf, nhinc_f, nwinc_f])


def loop_builder(loop_pos, alpha, beta, gamma, w_l, h_l, loop_fil_params, loop_list):
    wf, hf, nhinc_f, nwinc_f, sigma_l = loop_fil_params
    n_det, _ = np.shape(loop_pos)

    for i in range(n_det):
        p = loop_pos[i, :]
        loop_list.append([p, w_l, h_l, alpha, beta, gamma, sigma_l, wf, hf, nhinc_f, nwinc_f])


def circular_loop_builder(build_params, geo_objects):
    if not ('circ_pass_loops' in geo_objects):
        geo_objects['circ_pass_loops'] = []

    alpha, beta, gamma = build_params["yaw"], build_params["pitch"], build_params["roll"]
    r_l, loop_pos, n_nodes = build_params["radius"], build_params['center_pos'], build_params["node count"]
    c_d = build_params["contact distance"]
    loop_fil_params = build_params["filament parameters"]

    R_mat = gh.r_ypr(alpha, beta, gamma)
    nodes_pos = np.zeros((n_nodes, 3))
    # TODO change design so that nodes_pos has a method that builds these from the geometry params to reuse it for others

    nodes_pos[:, 0] = r_l * np.cos(np.linspace(0, 2 * np.pi - c_d, n_nodes))
    nodes_pos[:, 1] = r_l * np.sin(np.linspace(0, 2 * np.pi - c_d, n_nodes))
    nodes_pos = np.einsum('ij,kj->ik', R_mat, nodes_pos) + np.transpose(loop_pos)
    nodes_pos = np.transpose(nodes_pos)
    nodes_pos_shifted = np.vstack((nodes_pos[1:n_nodes, :], nodes_pos[0, :]))

    n_loop = len(geo_objects["circ_pass_loops"])

    node_names = np.chararray((n_nodes,), unicode=True) + "N_CL_" + str(n_loop) + "_" + list(
        map(str, list(np.arange(n_nodes))))

    segments = [nodes_pos, nodes_pos_shifted]

    segment_names = np.chararray((n_nodes - 1,), unicode=True) + "E_CL_" + list(
        map(str, list(np.arange(n_nodes - 1)))) + " " + node_names[0:n_nodes - 1] + " " + node_names[1:n_nodes]

    seg_dir = (nodes_pos_shifted - nodes_pos)
    seg_centers = nodes_pos[:-1] + 0.5 * seg_dir[:-1]
    seg_w_vec = seg_centers - loop_pos

    seg_params = [seg_centers, seg_w_vec, loop_fil_params]

    gate_list = [[node_names[0], node_names[-1]]]
    # TODO think if category "passive loops" or only "loops" is enough and the rest can be done by keyword "type"
    circ_loop = {"name": 'CL_' + str(n_loop), "type": "circular loop",
                 "super object": None, "nodes": nodes_pos, "node_names": list(node_names), "segments": segments,
                 "segment_names": list(segment_names), "seg_params": seg_params, "external": gate_list,
                 "build parameters": build_params}

    geo_objects["circ_pass_loops"].append(circ_loop)


def rectangle_loop_builder(r_build_params, geo_objects):
    if not ('r_pass_loops' in geo_objects):
        geo_objects['r_pass_loops'] = []

    n_nodes = 5
    a_r, b_r, center_pos = r_build_params['a_r'], r_build_params['b_r'], r_build_params['center_pos']
    fil_params = r_build_params["filament parameters"]
    n_loop = len(geo_objects['r_pass_loops'])

    node_names = np.chararray((5,), unicode=True) + "N_RL_" + str(n_loop) + "_" + list(
        map(str, list(np.arange(5))))
    nodes_pos = node_builder(gh.f_rlb, r_build_params)

    nodes_pos_shifted = np.vstack((nodes_pos[1:n_nodes, :], nodes_pos[0, :]))

    segments = [nodes_pos, nodes_pos_shifted]
    segment_names = np.chararray((n_nodes - 1,), unicode=True) + "E_RL_" + list(
        map(str, list(np.arange(n_nodes - 1)))) + " " + node_names[0:n_nodes - 1] + " " + node_names[1:n_nodes]

    seg_dir = (nodes_pos_shifted - nodes_pos)
    seg_centers = nodes_pos[:-1] + 0.5 * seg_dir[:-1]
    seg_w_vec = seg_centers - center_pos

    seg_params = [seg_centers, seg_w_vec, fil_params]

    gate_list = [[node_names[0], node_names[-1]]]

    rect_loop = {"name": 'RL_' + str(n_loop), "type": "rectangle loop",
                 "super object": None, "nodes": nodes_pos, "node_names": list(node_names), "segments": segments,
                 "segment_names": list(segment_names), "seg_params": seg_params, "external": gate_list,
                 "build parameters": r_build_params}

    geo_objects["r_pass_loops"].append(rect_loop)


def wire_builder(p1, p2, w_wire, h_wire, phys_params, fil_params, wire_list, external=False):
    # TODO change wire def to the style of the circular loop

    sigma = phys_params["sigma"]
    nhinc, nwinc, _ = fil_params
    name = {"Node_1": None, "Node_2": None, "Segment": None, "external": external}
    wire_list.append([p1, p2, w_wire, h_wire, nhinc, nwinc, sigma, name])


def plane_builder_angle(p, alpha_p, beta_p, gamma_p, a, b, thick, sigma_p, m_grid, plane_list):
    """defines a plane made of segments,
    p is array with the center of weight coordinates
    alpha_p, beta_p, gamma_p are yaw pitch roll angles of the plane normal
    initial orientation of the normal, when all angles are zero, is parallel to positive z axis
    so that edge with length a is parallel to x axis
    m_grid is the amount of filaments along an edge
    we use the same amount along the two edges
    thick is the thickness of the segments and thus the plane thickness
    """
    p1 = np.array([-a / 2, b / 2, 0])
    p2 = np.array([-a / 2, -b / 2, 0])
    p3 = np.array([a / 2, -b / 2, 0])

    r_mat = gh.r_ypr(alpha_p, beta_p, gamma_p)  # rotation matrix
    q1 = p + np.matmul(r_mat, p1)  # rotate the initial points
    q2 = p + np.matmul(r_mat, p2)
    q3 = p + np.matmul(r_mat, p3)

    if m_grid >= 2:
        plane_list.append([q1, q2, q3, thick, m_grid, sigma_p])
    else:
        print("Invalid plane grid, increase number of segments to m_grid>=2.")


def loop_plane_builder(build_params, geo_objects):
    """defines a plane made of loops,
    p is array with the center of weight coordinates
    alpha_p, beta_p, gamma_p are yaw pitch roll angles of the plane normal
    initial orientation of the normal, when all angles are zero, is parallel to positive z axis
    so that edge with length a is parallel to x axis
    n_loops is the amount of loops the plane is made of, each loop has same filament parameters
    thick is the thickness of the segments and thus the plane thickness
    """

    if not ('planes' in geo_objects):
        geo_objects['planes'] = []



    plane_pos = build_params["center_pos"]
    alpha, beta, gamma = build_params["yaw"], build_params["pitch"], build_params["roll"]
    cd = build_params["contact distance"]
    a, b = build_params["edge a"], build_params["edge b"]
    n_loops = build_params["loop count"]

    if a >= b:
        w_f = 0.5 * b / n_loops
        a_loops = a - 2 * w_f * np.arange(n_loops)
        b_loops = b - 2 * w_f * np.arange(n_loops)

    else:
        w_f = 0.5 * a / n_loops
        a_loops = a - w_f * np.arange(n_loops)
        b_loops = b - w_f * np.arange(n_loops)

    for i in range(n_loops):
        r_fil_params = {"width": w_f, "height": 0.002,
                        "width subs": 1, "height subs": 1, "conductivity": 2.5 * 10 ** 6}
        r_build_params = {'center_pos': plane_pos, "yaw": alpha, "pitch": beta, "roll": gamma,
                          'a_r': a_loops[i], 'b_r': b_loops[i],
                          'node count': 5, 'contact distance': cd, "filament parameters": r_fil_params}

        rectangle_loop_builder(r_build_params, geo_objects)

    plane = {'build_params': build_params}

    geo_objects['planes'].append(plane)



    # TODO plane fil params depend on the definition of the plane thickness and
    #plane_fil_params = plane_build_params['filament parameters']

    # 1. build loops to fill plane from plane params
    # plane rotations are done in the loops!

    # 2. save plane characteristics to a dictionary

    # TODO 3. think if necessary to label the impedance entries in Z_mat so
    #   that field contributions from different objects can be separated!

    #plane = {"build parameters": plane_build_params}




