# writer functions
import code.geometry_helpers as gh


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
    wire_list = geo_objects['wires']

    out_inp = open(input_filename, "a")

    # write all the nodes into the input file

    for wire in wire_list:
        wire_ind = wire_list.index(wire)
        out_inp.writelines("\n" + "*The nodes of the wire # " + str(wire_ind) + "\n")

        nodes_pos, node_names = wire["nodes"], wire["node_names"]

        for node_name in node_names:
            node_ind = node_names.index(node_name)

            out_inp.writelines(node_name
                               + " x=" + str(nodes_pos[node_ind, 0])
                               + " y=" + str(nodes_pos[node_ind, 1])
                               + " z=" + str(nodes_pos[node_ind, 2]) + "\n")

        segments, segment_names = wire["segments"], wire["segment_names"]
        seg_centers, seg_w_vec, loop_fil_params = wire["seg_params"]

        wire_fil_params = wire["build parameters"]["filament parameters"]
        wf = wire_fil_params["width"]
        hf = wire_fil_params["height"]
        nhinc = wire_fil_params["height subs"]
        nwinc = wire_fil_params["width subs"]
        sigma = wire_fil_params["conductivity"]

        for segment_name in segment_names:
            seg_ind = segment_names.index(segment_name)
            wx = seg_w_vec[seg_ind, 0]
            wy = seg_w_vec[seg_ind, 1]
            wz = seg_w_vec[seg_ind, 2]
            out_inp.writelines(segment_name
                               + " w=" + str(wf) + " h=" + str(hf) + " sigma= " + str(sigma)
                               #                               + " wx=" + str(wx) + " wy=" + str(wy) + " wz=" + str(wz)
                               + " nhinc=" + str(nhinc) + " nwinc=" + str(nwinc) + "\n")

        out_inp.writelines("\n" + "*Define in and out" + "\n")  # define in and out nodes
        gate_list = wire["external"]
        for gates in gate_list:
            out_inp.writelines("\n" + ".External " + gates[0] + " " + gates[1])
            out_inp.writelines("\n")

    # write all the segments between the nodes into input file

    out_inp.writelines("\n" + "*The segments of the cuboids\n")


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


def write_pass_loop_input(geo_objects, input_filename):
    loop_list = geo_objects['pass_loops']

    out_inp = open(input_filename, "a")
    if len(loop_list) > 0:
        for loop_ind in range(len(loop_list)):
            loop = geo_objects['pass_loops'][loop_ind]

            nodes_pos, node_names = loop["nodes"], loop["node_names"]

            out_inp.writelines("\n" + "*The nodes and segments of passive loops # " + str(loop_ind) + " \n")

            for node_name in node_names:
                node_ind = node_names.index(node_name)

                out_inp.writelines(node_name
                                   + " x=" + str(nodes_pos[node_ind, 0])
                                   + " y=" + str(nodes_pos[node_ind, 1])
                                   + " z=" + str(nodes_pos[node_ind, 2]) + "\n")

            segments, segment_names = loop["segments"], loop["segment_names"]
            seg_centers, seg_w_vec, loop_fil_params = loop["seg_params"]

            circ_fil_params = loop["build parameters"]["filament parameters"]
            wf = circ_fil_params["width"]
            hf = circ_fil_params["height"]
            nhinc = circ_fil_params["height subs"]
            nwinc = circ_fil_params["width subs"]
            sigma = circ_fil_params["conductivity"]

            for segment_name in segment_names:
                seg_ind = segment_names.index(segment_name)
                wx = seg_w_vec[seg_ind, 0]
                wy = seg_w_vec[seg_ind, 1]
                wz = seg_w_vec[seg_ind, 2]
                out_inp.writelines(segment_name
                                   + " w=" + str(wf) + " h=" + str(hf) + " sigma= " + str(sigma)
                                   + " wx=" + str(wx) + " wy=" + str(wy) + " wz=" + str(wz)
                                   + " nhinc=" + str(nhinc) + " nwinc=" + str(nwinc) + "\n")
            out_inp.writelines("\n" + "*Define in and out" + "\n")  # define in and out nodes
            gate_list = loop["external"]
            for gates in gate_list:
                out_inp.writelines("\n" + ".External " + gates[0] + " " + gates[1])
                out_inp.writelines("\n")
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


