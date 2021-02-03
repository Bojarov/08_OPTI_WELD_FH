import code.FH_input_writers as FHin
import code.geometry_builders as gb
import code.FH_output_helpers as FHout
import os.path
import os
import subprocess
import numpy as np


def run_FH_dyn(params, damage_params, filament_params, output_params,
               r_sub_vec, node_dens_vec, l_sub_vec, tube_segment_lists, tube_pts_lists,
               loop_list, sub_con_list):
    """ takes list with all the defined segments of the geometry
        writes FH input files and runs FH to obtain the current
        density in the conductor as a matlab file *.mat

        stores current densities in output folder


        all other FH output like impedance matrices are deleted at the end of
        the function
    """

    Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq = params
    rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params
    nhinc, nwinc, ndec, units, sub_div_auto = filament_params
    z0, x_detect, y_detect, z_detect = output_params

    # write the input file for FASTHENRY
    FH_input_filename = "input_Ro_" + str(Ro) + "_Ri_" + str(Ri) + "_nd_" + str(node_dens) + \
                        "_l_" + str(flen) + "_sigm_" + str(sigma) + "_mur_" + str(mu_r) + "_f_" + \
                        str(freq) + "_phi_c_" + str(phi_c_weld / np.pi) + \
                        "_lweld_" + str(l_weld) + ".inp"

    if not os.path.isdir('FH_input_files'):
        os.mkdir('FH_input_files')

    # check if current density already exists for given set of parameters
    j_dens_ind = FHout.j_dens_storage(r_sub_vec, l_sub_vec, node_dens_vec, loop_list,
                                      sub_con_list, params, damage_params, filament_params, message=True)

    if j_dens_ind == -1:
        print("No current densities found for requested set of parameters - running FH and saving output")
        print("")

        out_inp = open(FH_input_filename, "a")
        out_inp.truncate(0)  # delete old file
        FHin.write_header_input(units, sigma, freq, mu_r, ndec, out_inp)  # write physical parameters
        # into header of input file

        FHin.write_node_seg_input_dyn(tube_segment_lists, FH_input_filename,
                                      l_sub_vec, out_inp)  # write the segments

        FHin.write_node_ports_input_dyn(tube_segment_lists, FH_input_filename,
                                        l_sub_vec, out_inp)  # write the ports and equalize
        # tube subdivisions

        FHin.write_loop_node_seg_input(loop_list, FH_input_filename)  # write loop segments

        FHin.write_subcon_seg_input(sub_con_list, FH_input_filename)  # write tube sub connections

        FHin.write_end_input(FH_input_filename)

        out_inp.close()

        os.rename('./' + FH_input_filename, './FH_input_files/' + FH_input_filename)

        # Run FASTHENRY
        if sub_div_auto == 1:
            p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
                                  './FH_input_files/' + FH_input_filename, '-d', 'grids'],
                                 stdout=subprocess.PIPE)

            output_fh = p.communicate()
            print(output_fh)

        elif sub_div_auto == 0:

            p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
                                  './FH_input_files/' + FH_input_filename, '-d', 'grids', '-a', 'off'],
                                 stdout=subprocess.PIPE)

            output_fh = p.communicate()
            print(output_fh)

        else:
            "wrong auto refinement chosen in FH_run, only 0-off or 1-on"
        print("j_dens_extract is called")
        j_pack = FHout.j_dens_extract(r_sub_vec, l_sub_vec, node_dens_vec, params,
                                      damage_params, filament_params, output_params,
                                      tube_segment_lists, loop_list, sub_con_list)

        # save the calculated current density
        currents_path = 'FH_output_files/current_densities/current_dens.npy'

        with open(currents_path, 'rb') as currents:
            currents_list = list(np.load(currents, allow_pickle=True))
        currents_list.append(j_pack)

        with open(currents_path, 'wb') as currents:
            np.save(currents, np.array(currents_list, dtype=object))

    else:
        print("Found FH current density for requested set of parameters")
        print("")

    if os.path.exists('b.mat'):
        os.remove('b.mat')
    if os.path.exists('Zc.mat'):
        os.remove('Zc.mat')
    if os.path.exists('Jmag1_0.mat'):
        os.remove('Jmag1_0.mat')
