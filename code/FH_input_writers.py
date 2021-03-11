# writer functions
import code.geometry_helpers as gh
import numpy as np
from textwrap import wrap


def write_header_input(units, sigma, freq, mu_r, ndec, out_inp):
    """puts the physical parameters in the header of the fasthenry input file
    """
    out_inp.writelines("*Input to calculate current distribution in a wire\n")
    out_inp.writelines("*Conductivity" + "\n")
    # sigma_str = '%.2E' % sigma
    freq_int = freq * mu_r
    freq_int_str = '%.2E' % freq_int  # this makes FH use the effective frequency

    # out_inp.writelines(".Default "+ "sigma="+sigma_str+"\n")
    out_inp.writelines(".Units " + units + "\n")
    out_inp.writelines("*Frequency range of interest" + "\n")
    out_inp.writelines("*Note that this is an effective frequency to include FMs" + "\n")
    out_inp.writelines(".freq " + "fmin=" + freq_int_str + " fmax=" + freq_int_str + " ndec="
                       + str(ndec) + "\n")
    out_inp.close()


def write_header_ZC(phys_params_run, input_filename, ndec=1, units="M"):
    sigma, mu_0, mu_r, freq = phys_params_run
    """puts the physical parameters in the header of the fasthenry input file
    """
    out_inp = open(input_filename, "a")
    out_inp.writelines("*Input to calculate current distribution in a wire\n")
    out_inp.writelines("*Conductivity" + "\n")
    # sigma_str = '%.2E' % sigma
    freq_int = freq * mu_r
    freq_int_str = '%.2E' % freq_int  # this makes FH use the effective frequency

    # out_inp.writelines(".Default "+ "sigma="+sigma_str+"\n")
    out_inp.writelines(".Units " + units + "\n")
    out_inp.writelines("*Frequency range of interest" + "\n")
    out_inp.writelines("*Note that this is an effective frequency to include FMs" + "\n")
    out_inp.writelines(".freq " + "fmin=" + freq_int_str + " fmax=" + freq_int_str + " ndec="
                       + str(ndec) + "\n")
    out_inp.close()


def write_node_seg_wire(geo_objects, input_filename):
    """
    writes all defined cuboids into the input file
    """
    # TODO finish the wire writer
    # TODO define wire with multiple subdivs...requires input format as series of node points

    wire_list = geo_objects['wires']

    out_inp = open(input_filename, "a")

    # write all the nodes into the input file

    out_inp.writelines("\n" + "*The nodes of the wires\n")

    for wire in wire_list:  # loop through the shells
        wire_ind = wire_list.index(wire)

        out_inp.writelines("N_" + str(wire_ind) + "_" + str(1)
                           + " x=" + str(wire[0][0])
                           + " y=" + str(wire[0][1])
                           + " z=" + str(wire[0][2]) + "\n")

        # add the end node of each filament to input file

        out_inp.writelines("N_" + str(wire_ind) + "_" + str(2)
                           + " x=" + str(wire[1][0])
                           + " y=" + str(wire[1][1])
                           + " z=" + str(wire[1][2]) + "\n")
        wire[7]["Node_1"] = "N_" + str(wire_ind) + "_" + str(1)
        wire[7]["Node_2"] = "N_" + str(wire_ind) + "_" + str(2)

    # write all the segments between the nodes into input file

    out_inp.writelines("\n" + "*The segments of the cuboids\n")

    for wire in wire_list:
        wire_ind = wire_list.index(wire)

        out_inp.writelines(
            "E_" + str(wire_ind)
            + " N_" + str(wire_ind) + "_" + str(1)
            + " N_" + str(wire_ind) + "_" + str(2)
            + " w=" + str(wire[2]) + " h=" + str(wire[3]) + " sigma= " + str(wire[6])
            + " nhinc=" + str(wire[4]) + " nwinc=" + str(wire[5]) + "\n")
        wire[7]["Segment"] = "E_" + str(wire_ind)

        if wire[7]["external"]:
            node_in = wire[7]["Node_1"]
            node_out = wire[7]["Node_2"]
            external = True
    if external:
        out_inp.writelines("\n" + "*Define in and out" + "\n")  # define in and out nodes
        out_inp.writelines("\n" + ".External " + node_in + " " + node_out)


def write_node_seg_input_dyn(tube_segment_lists, input_filename, l_sub_vec, out_inp):
    """
    writes all defined cuboids into the input file
    """
    out_inp = open(input_filename, "a")

    # write all the nodes into the input file

    out_inp.writelines("\n" + "*The nodes of the cuboids\n")

    for segment_list in tube_segment_lists:  # loop throug the shells
        seg_list_ind = tube_segment_lists.index(segment_list)

        for segment in segment_list:  # loop through segments in a shell
            seg_ind = segment_list.index(segment)

            seg_in_fil_ind = int(seg_ind % (len(l_sub_vec) - 1))  # index of a segment in a filament
            fil_ind = int((seg_ind / (len(l_sub_vec) - 1)))  # index of a filament

            out_inp.writelines("N_" + str(seg_list_ind)
                               + "_" + str(fil_ind) + "_" + str(seg_in_fil_ind)
                               + " x=" + str(segment[0][0])
                               + " y=" + str(segment[0][1])
                               + " z=" + str(segment[0][2]) + "\n")

            # add the end node of each filament to input file
            if seg_in_fil_ind % (len(l_sub_vec) - 1) == (len(l_sub_vec) - 2):
                out_inp.writelines("N_" + str(seg_list_ind)
                                   + "_" + str(fil_ind) + "_" + str(seg_in_fil_ind + 1)
                                   + " x=" + str(segment[1][0])
                                   + " y=" + str(segment[1][1])
                                   + " z=" + str(segment[1][2]) + "\n")

    # write all the segments between the nodes into input file

    out_inp.writelines("\n" + "*The segments of the cuboids\n")

    for segment_list in tube_segment_lists:
        seg_list_ind = tube_segment_lists.index(segment_list)

        for segment in segment_list:
            seg_ind = segment_list.index(segment)
            seg_in_fil_ind = int(seg_ind % (len(l_sub_vec) - 1))  # index of a segment in a filament
            fil_ind = int((seg_ind / (len(l_sub_vec) - 1)))  # index of a filament

            out_inp.writelines(
                "E_" + str(seg_list_ind) + "_" + str(fil_ind) + "_" + str(seg_in_fil_ind)
                + " N_" + str(seg_list_ind) + "_" + str(fil_ind) + "_" + str(seg_in_fil_ind)
                + " N_" + str(seg_list_ind) + "_" + str(fil_ind) + "_" + str(seg_in_fil_ind + 1)
                + " w=" + str(segment[2]) + " h=" + str(segment[3]) + " sigma= " + str(segment[6])
                + " nhinc=" + str(segment[4]) + " nwinc=" + str(segment[5]) + "\n")


def write_node_ports_input_dyn(tube_segment_lists, input_filename, l_sub_vec, out_inp):
    out_inp = open(input_filename, "a")

    out_inp.writelines("\n" + "*Equalize nodes at sub divisions" + "\n")

    n_sub = len(l_sub_vec)  # number of filament subdivisions

    equi_strings = [''] * n_sub  # list to collect name strings of equal
    # subdivisions

    sub_string = np.chararray((n_sub, 1), itemsize=n_sub, unicode=True)
    sub_string[:, 0] = list(map(str, list(np.arange(n_sub))))  # char array with
    # subdivision index

    for segment_list in tube_segment_lists:

        seg_list_ind = tube_segment_lists.index(segment_list)  # shell index

        n_fil_seg = int(len(segment_list) / (len(l_sub_vec) - 1))  # number of filaments in the shell
        fil_string = np.chararray((1, n_fil_seg),
                                  itemsize=n_fil_seg, unicode=True)

        fil_string[0, :] = list(map(str, list(np.arange(n_fil_seg))))  # char array with filament index

        seg_str_mat = ("N_" + str(seg_list_ind) + "_"
                       + fil_string + "_" + sub_string)  # string matrix of all
        # nodes in this shell

        for i in range(n_sub):  # gather the name strings
            # of all nodes at same subdiv

            str_vec = seg_str_mat[i, :]
            str_full = ' '.join(str_vec)

            equi_strings[i] = equi_strings[i] + ' ' + str_full
            str_wrap = wrap(str_full, 60)

            str_full = '\n+ '.join(str_wrap)

    for i in range(n_sub):  # write equilisation of nodes
        # at same subdiv into file
        equi_str_wrap = wrap(equi_strings[i], 60)
        equi_str_full = '\n+ '.join(equi_str_wrap)[1:]
        equi_strings[i] = equi_str_full

        out_inp.writelines("\n" + "*Equipotential nodes at filament sub division "
                           + str(i + 1) + " out of " + str(n_sub) + " sub divisions" + "\n")

        out_inp.writelines(".Equiv " + equi_strings[i])

    out_inp.writelines("\n")

    out_inp.writelines("\n" + "*Define tube in and out" + "\n")  # define tube in and out nodes
    out_inp.writelines("\n" + ".External N_0_0_0 N_0_0_" + str(n_sub - 1))
    out_inp.close()


def write_loop_input(geo_objects, input_filename):
    loop_list = geo_objects['loops']
    out_inp = open(input_filename, "a")
    if len(loop_list) > 0:
        for loop_ind in range(len(loop_list)):
            # unpacking loop parameters
            p, wl, hl, alpha, beta, gamma, sigma_l, wf, hf, nhinc_f, nwinc_f = \
                loop_list[loop_ind]

            out_inp.writelines("\n" + "*The nodes and segments of loop # " + str(loop_ind) + " \n")

            corners = list(gh.loop_corners(loop_list[loop_ind]))
            for i in range(len(corners)):
                corner = corners[i]
                # write corner node coordinates into file
                out_inp.writelines("N_L_" + str(loop_ind) + "_" + str(i)
                                   + " x=" + str(corner[0])
                                   + " y=" + str(corner[1])
                                   + " z=" + str(corner[2]) + "\n")

            for i in range(len(corners)):
                out_inp.writelines("E_L_" + str(loop_ind) + "_" + str(i)
                                   + " N_L_" + str(loop_ind) + "_" + str(i) + " N_L_" + str(loop_ind) + "_" + str(
                    (i + 1) % len(corners))
                                   + " w=" + str(wf) + " h=" + str(hf) + " sigma= " + str(sigma_l)
                                   + " nhinc=" + str(nhinc_f) + " nwinc=" + str(nwinc_f) + "\n")

    out_inp.close()


def write_det_loop_input(geo_objects, input_filename):
    loop_list = geo_objects['det_loops']

    out_inp = open(input_filename, "a")
    if len(loop_list) > 0:
        for loop_ind in range(len(loop_list)):
            # unpacking loop parameters
            p, wl, hl, alpha, beta, gamma, sigmal, wf, hf, nhinc_f, nwinc_f = \
                loop_list[loop_ind]

            out_inp.writelines("\n" + "*The nodes and segments of detector loops # " + str(loop_ind) + " \n")

            corners = list(gh.det_loop_corners(loop_list[loop_ind]))
            for i in range(len(corners)):
                corner = corners[i]
                # write corner node coords into file
                out_inp.writelines("N_DL_" + str(loop_ind) + "_" + str(i)
                                   + " x=" + str(corner[0])
                                   + " y=" + str(corner[1])
                                   + " z=" + str(corner[2]) + "\n")

            for i in range(len(corners) - 1):
                out_inp.writelines("E_DL_" + str(loop_ind) + "_" + str(i)
                                   + " N_DL_" + str(loop_ind) + "_" + str(i) + " N_DL_" + str(loop_ind) + "_" + str(
                    (i + 1) % len(corners))
                                   + " w=" + str(wf) + " h=" + str(hf) + " sigma= " + str(sigmal)
                                   + " nhinc=" + str(nhinc_f) + " nwinc=" + str(nwinc_f) + "\n")

            out_inp.writelines("\n" + "*Define in and out" + "\n")  # define in and out nodes
            out_inp.writelines("\n" + ".External " + "N_DL_" + str(loop_ind)
                               + "_" + "0" + " " + "N_DL_" + str(loop_ind) + "_" + str(len(corners) - 1))
            out_inp.writelines("\n")
    out_inp.close()


def write_subcon_seg_input(sub_con_list, input_filename):
    out_inp = open(input_filename, "a")

    for sub_con in sub_con_list:
        sub_con_ind = sub_con_list.index(sub_con)

        shell_ind_s, fil_ind_s, s1, shell_ind_e, fil_ind_e, s2 = sub_con[7]

        out_inp.writelines("\n")
        out_inp.writelines("\n" + "*Sub connections between tube filaments" + "\n")  # define tube in and out nodes
        out_inp.writelines(
            "E_SUB_" + str(sub_con_ind)
            + " N_" + str(shell_ind_s) + "_" + str(fil_ind_s) + "_" + str(s1)
            + " N_" + str(shell_ind_e) + "_" + str(fil_ind_e) + "_" + str(s2)
            + " w=" + str(sub_con[2]) + " h=" + str(sub_con[3]) + " sigma= " + str(sub_con[6])
            + " nhinc=" + str(sub_con[4]) + " nwinc=" + str(sub_con[5]) + "\n")

    out_inp.close()


def write_plane_input(geo_objects, input_filename):
    plane_list = geo_objects['planes']
    """
    writes all defined planes into the input file
    """

    out_inp = open(input_filename, "a")
    out_inp.writelines("\n" + "*The planes defined...\n")
    for i in plane_list:
        plane_index = plane_list.index(i)
        sigma_plane = i[5]

        out_inp.writelines("g_" + str(plane_index)  # write the corners of the plane
                           + " x1=" + str(i[0][0])
                           + " y1=" + str(i[0][1])
                           + " z1=" + str(i[0][2])
                           + " x2=" + str(i[1][0])
                           + " y2=" + str(i[1][1])
                           + " z2=" + str(i[1][2])
                           + " x3=" + str(i[2][0])
                           + " y3=" + str(i[2][1])
                           + " z3=" + str(i[2][2])
                           + " thick=" + str(i[3]) + " sigma=" + str(sigma_plane) + " file = NONE"
                           + "\n"
                           + "+ contact initial_grid (" + str(int(i[4])) + "," + str(int(i[4])) + ")")


def write_end_input(input_filename):
    out_inp = open(input_filename, "a")

    out_inp.writelines("\n")
    out_inp.writelines("\n" + ".End")

    out_inp.close()
