import code.FH_input_writers_ZC as FHin
import os.path
import os
import subprocess
import textwrap
import numpy as np

FH_input_filename = "input.inp"


def run_FH_ZC(phys_params_run, geo_objects, sub_div_auto, message=False):
    if not os.path.isdir('ZC_input_files'):
        os.mkdir('ZC_input_files')

    FHin.write_header_ZC(phys_params_run, FH_input_filename)
    FHin.write_node_seg_wire(geo_objects, FH_input_filename)
    FHin.write_det_loop_input(geo_objects, FH_input_filename)
    if 'pass_loops' in geo_objects:
        FHin.write_pass_loop_input(geo_objects, FH_input_filename)

    FHin.write_end_input(FH_input_filename)
    os.rename('./' + FH_input_filename, './ZC_input_files/' + FH_input_filename)
    # Run FASTHENRY
    if sub_div_auto:
        p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
                              './ZC_input_files/' + FH_input_filename, '-d', 'grids'],
                             stdout=subprocess.PIPE)
        # p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
        #                      './ZC_input_files/' + FH_input_filename],
        #                     stdout=subprocess.PIPE)

        output_fh = p.communicate()
        if message:
            print(output_fh)

    else:

        p = subprocess.Popen(['/usr/local/xictools/bin/fasthenry',
                              './ZC_input_files/' + FH_input_filename, '-a', 'off'],
                             stdout=subprocess.PIPE)

        output_fh = p.communicate()
        if message:
            print(output_fh)


    dir_files = os.listdir()



    for file in dir_files:
        if file.startswith("J"):
            os.remove(file)

