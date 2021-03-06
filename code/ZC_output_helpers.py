import os.path
import re
import numpy as np


def ZC_mat_extract(geo_objects):
    #TODO works for now, but I would like to get the number of gates from the keys of the dicts
    n_gates = len(geo_objects["det_loops"]) + len(geo_objects["pass_loops"]) + len(geo_objects["wires"])

    if os.path.exists('./Zc.mat'):
        Z_mat = np.zeros((n_gates, n_gates), dtype=complex)
        with open('./Zc.mat') as ZC_file:
            for i in range(n_gates + 1):
                next(ZC_file)  # skip header
            rows = [np.array(re.split(' +', line)) for line in ZC_file]

            for i in range(len(rows)):
                row = rows[i]
                full_row = np.array(row[1:2 * n_gates+1].astype(np.complex))
                Z_mat[i, :] = full_row[::2] + full_row[1::2]
    return Z_mat
