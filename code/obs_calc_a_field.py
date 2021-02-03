import numpy as np
import scipy
import code.FH_output_helpers as FHout
from scipy.integrate import tplquad
from code.c_code import KerA


# A field calculations

def A_ker_pre_point_calc(x1, y1, z1, tube_segment_lists, params):
    """Calculates the contribution of every filament to the Vector potential
    divided by the current density of the particular filament and stores it
    
    the matrix of contributions is (5, 3, n_seg) large and
    can be used to efficiently recalculate the B field in the point 
    x1 ,y1, z1 in the case when the physical params change but the geometry of
    the conductor detector setup remains the same

    B_z can be calculated from the cross of 5 points which are stored
    in the following order:
    y               ______
    A              |     |
    |              |  4  |
    |         _____|     |_____
    |        |                 |
    |        |  2     3     1  |
    |        |_____      ______|
    |              |     |
    |              |  5  |
    |              |_____|
    |              
    ----------------------------------> x

    """
    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params  # unpack

    skin_depth = np.sqrt(2 / (2 * np.pi * freq * sigma * mu_0))

    integrand = scipy.LowLevelCallable.from_cython(KerA, 'integrand')  # load the C
    # implemented
    # kernel

    h_del = min(0.01 * ro_t / node_dens, 0.1 * skin_depth)

    cross_xy = np.array([[x1 + h_del, y1],
                         [x1 - h_del, y1],
                         [x1, y1],
                         [x1, y1 + h_del],
                         [x1, y1 - h_del]])

    # calculate the total number of filaments
    n_seg = 0
    for segment_list in tube_segment_lists:
        n_seg = n_seg + len(segment_list)

    # create cross (5) for each spatial dimension (3) in which A field from each
    # of the (N_fil) filaments are calculated
    a_ker_mat = np.zeros((5, 3, n_seg))
    #print("a_ker_mat does not need to have 3 spatial dimensions!")
    # i define it so that additional spatial dimensions can be included, for npw only z

    #print("CAUTION in obs_calc, when defined shell has no segments -> error in loop through segment_list")
    #print("CAUTION in obs_calc, overthink the h definition")
    for i in range(len(cross_xy)):
        a_ker_vec = []
        point = cross_xy[i]

        for segment_list in tube_segment_lists:
            for segment in segment_list:
                x1 = point[0]  # a point of the cross at the detector
                y1 = point[1]

                xn = segment[0][0]  # node position at one face of
                yn = segment[0][1]  # of a filament
                zn = segment[0][2]

                # xe = segment[1][0]  # opposite face
                # ye = segment[1][1]
                ze = segment[1][2]

                w = segment[2]
                h = segment[3]

                a_ker_vec.append(np.array(tplquad(integrand, -w / 2, w / 2,
                                                  lambda x2: -h / 2, lambda x2: h / 2,
                                                  lambda x2, y2: zn,
                                                  lambda x2, y2: ze, args=(x1 - xn, y1 - yn, z1 - zn)))[0])

        a_ker_mat[i, 2, :] = np.array(a_ker_vec[:])
    return a_ker_mat


def A_ker_pre_pack_calc(r_sub_vec, node_dens_vec, loop_list, tube_segment_lists, params, output_params, l_sub_vec):
    """Creates the geometry pack which contains the characteristics
        of the conductor-detector setup and the purely geometrical A field
        contributions from each filament for all detector positions,
        the later must be multiplied by the current density of the filament to
        obtain the A field
    """

    z, x1, y1, z1 = output_params

    a_ker_mat_pack = []

    detector_positions = [[x1, y1, z1]]
    # calculate the cross of A field contributions for each detector position
    a_ker_mat = A_ker_pre_point_calc(x1, y1, z1, tube_segment_lists, params)

    a_ker_mat_pack.append(a_ker_mat)

    a_ker_pre_pack = np.array([r_sub_vec, node_dens_vec, output_params, detector_positions, tube_segment_lists,
                               a_ker_mat_pack, loop_list, l_sub_vec], dtype=object)

    return a_ker_pre_pack


def A_int_dyn(j_tube_re, j_tube_im, r_sub_vec, l_sub_vec, node_dens_vec, params, output_params, tube_segment_lists,
              loop_list):
    """Calculates the z component of the vector potential A in the points
    defined in the output specifications
    
    uses precalculated geometry if setup already exists
    """
    #print("CAUTION: only z directed currents possible in obs_calc.A_int_dyn !")
    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params  # unpack
    z0, w_detect, h_detect, _ = output_params

    # search for geometry in geometries database and load it
    # will extend database by new setup if not found
    # get index of current geometry in database

    geo_ind = FHout.A_ker_storage(r_sub_vec, l_sub_vec, node_dens_vec, tube_segment_lists, loop_list, params,
                                  output_params)

    # open database
    geometries_path = 'FH_output_files/geometry_files/geometries.npy'

    with open(geometries_path, 'rb') as geometries:
        geometries_list = list(np.load(geometries, allow_pickle=True))

    geometry_pack = geometries_list[geo_ind]

    a_ker_mat_pack = np.array(geometry_pack[5], dtype=object)
    # multiplying the filament integrals at all spatial points
    # with individual current density in z-direction
    # need to include other directions
    # PUT IN FULL TUBE CURRENT LISTS HERE!

    a_mat_pack_re_z = a_ker_mat_pack[:, :, 2, :] * j_tube_re[:, 5]
    a_mat_pack_im_z = a_ker_mat_pack[:, :, 2, :] * j_tube_im[:, 5]

    # matrices with Vector potential at the (dpts) detector points
    # for each dpt there is a cross (5)
    # for each spatial dimension (1)
    # so matrices have shape (N_pts_detect, 5)
    a_dpts_cross_xyz_re_z = mu_0 * mu_r / (4 * np.pi) * np.sum(a_mat_pack_re_z, axis=-1)
    a_dpts_cross_xyz_im_z = mu_0 * mu_r / (4 * np.pi) * np.sum(a_mat_pack_im_z, axis=-1)

    return a_dpts_cross_xyz_re_z, a_dpts_cross_xyz_im_z
