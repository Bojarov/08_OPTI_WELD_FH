import os
import os.path
import re

import code.obs_calc_a_field as oca
import numpy as np


def tube_geometry_compare(r_sub_vec, l_sub_vec, node_dens_vec, params, damage_params, filament_params, j_dens_pack):
    """compares the tube geometry to a already existant geometry that is saved in a data pack
    with the corresponding current density (j_dens_pack)
    """
    if ((j_dens_pack[0] == r_sub_vec).all() and
            (j_dens_pack[1] == node_dens_vec).all() and
            (j_dens_pack[2] == params) and
            (j_dens_pack[3] == np.array(damage_params)).all() and
            (j_dens_pack[4] == np.array(filament_params)).all() and
            len(j_dens_pack[8]) == len(l_sub_vec)):
        if (j_dens_pack[8] == l_sub_vec).all():
            same = True
        else:
            same = False
    else:
        same = False
    return same


def loop_compare(loop_list, j_dens_pack):
    """Compares the loops of a geometry saved in j_dens_pack. The order of the
    loops matters. If order of loops in list is different to saved order the
    new geometry will be declared different and a new j_dens_pack appended to
    the currents_list. 
    """
    if len(loop_list) == len(j_dens_pack[6]):
        if len(loop_list) == 0:
            same = True
        else:
            for loop_ind in range(len(loop_list)):
                loop = loop_list[loop_ind]
                if (j_dens_pack[6][loop_ind][0] == loop[0]).all() and \
                        (j_dens_pack[6][loop_ind][1:len(loop)] == loop[1:len(loop)]):
                    same = True
                else:
                    same = False
                    break

    else:
        same = False

    return same


def sub_con_compare(sub_con_list, j_dens_pack):
    """compares the subconnections in the geometry,
    order of subcons is important! if order in list idffers geometries are
    recognized as different - might need to fix this!
    """
    same = False

    if len(sub_con_list) == len(j_dens_pack[7]):
        if len(sub_con_list) == 0:
            same = True
        else:
            for sub in sub_con_list:
                sub_ind = sub_con_list.index(sub)
                if (np.array(j_dens_pack[7][sub_ind][0]) == np.array(sub[0])).all() and \
                        (np.array(j_dens_pack[7][sub_ind][1]) == np.array(sub[1])).all() and \
                        (j_dens_pack[7][sub_ind][2:len(sub) - 1] == sub[2:len(sub) - 1]):
                    same = True
                else:
                    same = False
                    break

    else:
        same = False

    return same


def j_dens_storage(r_sub_vec, l_sub_vec, node_dens_vec, loop_list, sub_con_list,
                   params, damage_params, filament_params, message=True):
    """checks if current density has been produced for this set of params
    """

    if not os.path.isdir('FH_output_files/current_densities'):
        os.makedirs('FH_output_files/current_densities')
    currents_path: str = 'FH_output_files/current_densities/current_dens.npy'

    # make an empty currents storage list if currents folder is empty
    if not os.path.exists(currents_path):
        with open(currents_path, 'wb') as currents:
            np.save(currents, [])

    with open(currents_path, 'rb') as currents:
        currents_list = list(np.load(currents, allow_pickle=True))

    if message:
        print("")
        print("Checking current density database")
        print("")

    # THIS IS A TOO SIMPLE WAY OF COMPARISON BECAUSE THE ORDER OF THE LOOPS IS IMPORTANT
    # in the LOOP LIST even when the resulting for inerchanged order stays the same
    # as does the current

    j_dens_ind_final = -1  # just a control flag, -1 means no currents found

    if len(currents_list) > 0:
        for j_dens_pack in currents_list:
            j_dens_ind = currents_list.index(j_dens_pack)

            # compare the bare pipe geometry
            if (tube_geometry_compare(r_sub_vec, l_sub_vec, node_dens_vec,
                                      params, damage_params, filament_params, j_dens_pack)):
                # print("bare tube geometry agrees")
                same_geo = True
            else:

                same_geo = False

            # compare loops
            if loop_compare(loop_list, j_dens_pack):
                same_loops = True
            else:
                same_loops = False

            # compare sub cons
            if sub_con_compare(sub_con_list, j_dens_pack):
                same_subs = True
            else:
                same_subs = False

            if same_geo and same_loops and same_subs:
                j_dens_ind_final = j_dens_ind
                # print("Found matching geometry")
                break

    return j_dens_ind_final


def j_dens_catogorize(j_dens_full, tube_segment_list, loop_list, sub_con_list):
    """seperates the total current density from the FH output into 
       current density of the tube and current density of the defined loops
       and defined sub connections
    """
    Ns_t = 0
    for segment_list in tube_segment_list:
        Ns_t = Ns_t + len(segment_list)

    Ns_l = 4 * len(loop_list)
    Ns_s = len(sub_con_list)

    j_dens_tube = j_dens_full[0: Ns_t, :]
    j_dens_loops = j_dens_full[Ns_t: Ns_t + Ns_l, :]
    j_dens_subs = j_dens_full[Ns_t + Ns_l: Ns_t + Ns_l + Ns_s, :]

    return j_dens_tube, j_dens_loops, j_dens_subs


def j_dens_extract(r_sub_vec, l_sub_vec, node_dens_vec, params, damage_params,
                   filament_params, output_params, tube_segment_lists, loop_list, sub_con_list):
    """Extracts the current density for the pipe segments from FH output
       CAUTION: USING ONLY CURRENT DENSITIES AT ONE END F THE WIRE
       THIS IS ONLY OK AS LONG AS THERE ARE NO RADIAL OR TANGENTIAL
       CURRENTS POSSIBLE in the tube
    """

    if not os.path.isdir('FH_output_files/current_densities'):
        os.makedirs('FH_output_files/current_densities')

    if os.path.exists('./Jreal1_0.mat') and \
            os.path.exists('./Jimag1_0.mat'):
        print("Extracting current density at desired points")

        # prepare the real part of the current density
        with open('./Jreal1_0.mat') as fp_real:
            # 1. iterate over file line-by-line
            # 2. split line by spaces into list (of number strings)
            # 3. convert number substrings to float values
            # 4. convert map object to list
            data_re = [list(map(float, re.split(' +', line))) for line in fp_real]

        x_re = []
        y_re = []
        z_re = []
        u_re = []
        v_re = []
        w_re = []

        for i in range(len(data_re)):
            # print(len(data_re))
            x_re.append(data_re[i][0])
            y_re.append(data_re[i][1])
            z_re.append(data_re[i][2])
            u_re.append(data_re[i][3])
            v_re.append(data_re[i][4])
            w_re.append(data_re[i][5])

        # full real part data
        xyz_re = np.array(list(zip(x_re, y_re, z_re, u_re, v_re, w_re)))

        j_tube_re, j_loops_re, j_subs_re = j_dens_catogorize(xyz_re, tube_segment_lists, loop_list, sub_con_list)

        # extract the j real at z=0
        j_tube_z_proj_re = np.array([i for i in list(j_tube_re) if i[2] == 0])

        # prepare the imag part of the current density
        with open('./Jimag1_0.mat') as fp_imag:
            # 1. iterate over file line-by-line
            # 2. split line by spaces into list (of number strings)
            # 3. convert number substrings to float values
            # 4. convert map object to list
            data_im = [list(map(float, re.split(' +', line))) for line in fp_imag]

        x_im = []
        y_im = []
        z_im = []
        u_im = []
        v_im = []
        w_im = []

        for i in range(len(data_im)):
            x_im.append(data_im[i][0])
            y_im.append(data_im[i][1])
            z_im.append(data_im[i][2])
            u_im.append(data_im[i][3])
            v_im.append(data_im[i][4])
            w_im.append(data_im[i][5])

        # full im part data
        xyz_im = np.array(list(zip(x_im, y_im, z_im, u_im, v_im, w_im)))

        j_tube_im, j_loops_im, j_subs_im = j_dens_catogorize(xyz_im, tube_segment_lists,
                                                             loop_list, sub_con_list)

        # extract the j real at z=0
        j_tube_z_proj_im = np.array([i for i in list(j_tube_im) if i[2] == 0])

        j_dens_cuts = [j_loops_re, j_tube_re, j_tube_z_proj_re,
                       j_loops_im, j_tube_im, j_tube_z_proj_im,
                       j_subs_re, j_subs_im]

        j_pack = [r_sub_vec, node_dens_vec,
                  params, damage_params, filament_params, j_dens_cuts, loop_list, sub_con_list, l_sub_vec]
    else:
        print("Could not find current density *.mat files, make sure to run FH")

    return j_pack


def A_ker_storage(r_sub_vec, l_sub_vec, node_dens_vec,
                  tube_segment_lists, loop_list, params, output_params):
    """Manages the storage of already simulated pipeline-detector 
       configurations/geometries:
       - checks if configuration already exists
       - for now we simply save a geometry pack numpy array in a geometry list
       - a loop goes through all generated 
         geometry runs and compares the relevant params
       - if geometry does not exists it is calculated and appended to the list
        return: index of current geometry in the geometry list
    """

    # z0, w_detect, h_detect, z_detect = output_params

    # create the geometry folder if not existent
    if not os.path.isdir('FH_output_files/geometry_files'):
        os.makedirs('FH_output_files/geometry_files')

    geometries_path = 'FH_output_files/geometry_files/geometries.npy'

    # make an empty geometries list if geometry folder is empty
    if not os.path.exists(geometries_path):
        with open(geometries_path, 'wb') as geometries:
            np.save(geometries, [])

    with open(geometries_path, 'rb') as geometries:
        geometries_list = list(np.load(geometries, allow_pickle=True))

    # check if geometry already exists
    print(" ")
    print("Scanning geometry files to check if requested conductor-detector setup already exists:")
    print(" ")

    geo_ind = -1

    if len(geometries_list) > 0:
        # go through all geometries and search for current one

        for geometry in geometries_list:

            same_geo = False
            if ((geometry[0] == r_sub_vec).all() and
                    (geometry[1] == node_dens_vec).all() and
                    (np.array(geometry[2]) == np.array(output_params)).all() and
                    (len(geometry[7]) == len(l_sub_vec))):

                if (geometry[7] == l_sub_vec).all():

                    same_geo = True
                else:
                    same_geo = False

            # if these match compare the defined loops in the geometry
            same_loop = False
            if len(loop_list) == len(geometry[6]):
                if len(loop_list) == 0:
                    same_loop = True
                else:
                    for loop in loop_list:
                        loop_ind = loop_list.index(loop)

                        if (geometry[6][loop_ind][0] == loop[0]).all() and \
                                (geometry[6][loop_ind][1:len(loop)] == loop[1:len(loop)]):

                            same_loop = True
                        else:
                            same_loop = False
                            break

            if same_geo and same_loop:
                geo_ind = geometries_list.index(geometry)
                break

    if geo_ind == -1:
        print("------>New setup - calculating integrals for field at detectors")
        A_ker_pre_pack = oca.A_ker_pre_pack_calc(r_sub_vec, node_dens_vec,
                                                 loop_list, tube_segment_lists, params, output_params, l_sub_vec)

        geometries_list.append(A_ker_pre_pack)
        geo_ind = len(geometries_list) - 1

        with open(geometries_path, 'wb') as geometries:
            np.save(geometries, np.array(geometries_list, dtype=object))

    else:
        print("Found conductor-detector setup in database - using stored field integrals")
        print("")

    return geo_ind
