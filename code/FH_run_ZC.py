import code.FH_input_writers as FHin
import os.path
import os
import subprocess
import textwrap
import numpy as np

FH_input_filename = "input.inp"


def run_FH_ZC(phys_params, geo_objects, sub_div_auto):
    if not os.path.isdir('ZC_input_files'):
        os.mkdir('ZC_input_files')


    FHin.write_header_ZC(phys_params, FH_input_filename)
    FHin.write_node_seg_wire(geo_objects, FH_input_filename)
    FHin.write_det_loop_input(geo_objects, FH_input_filename)
    #FHin.write_loop_input(geo_objects, FH_input_filename)
    FHin.write_plane_input(geo_objects, FH_input_filename)


    FHin.write_end_input(FH_input_filename)
    os.rename('./' + FH_input_filename, './ZC_input_files/' + FH_input_filename)

    # Run FASTHENRY
    if sub_div_auto:
        p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
                              './ZC_input_files/' + FH_input_filename, '-d', 'grids'],
                             stdout=subprocess.PIPE)
        #p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
        #                      './ZC_input_files/' + FH_input_filename],
        #                     stdout=subprocess.PIPE)

        output_fh = p.communicate()
        print(output_fh)

    else:

        p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
                              './ZC_input_files/' + FH_input_filename, '-a', 'off'],
                             stdout=subprocess.PIPE)

        output_fh = p.communicate()
        print(output_fh)

