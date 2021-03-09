from itertools import product
import numpy as np



def r_ypr(alpha, beta, gamma):
    """Rotation matrix - yaw/pitch/roll angles
    """
    r_mat = np.zeros((3, 3))
    r_mat[0, 0] = np.cos(alpha) * np.cos(beta)
    r_mat[0, 1] = np.cos(alpha) * np.sin(beta) * np.sin(gamma) - np.sin(alpha) * np.cos(gamma)
    r_mat[0, 2] = np.cos(alpha) * np.sin(beta) * np.cos(gamma) + np.sin(alpha) * np.sin(gamma)
    r_mat[1, 0] = np.sin(alpha) * np.cos(beta)
    r_mat[1, 1] = np.sin(alpha) * np.sin(beta) * np.sin(gamma) + np.cos(alpha) * np.cos(gamma)
    r_mat[1, 2] = np.sin(alpha) * np.sin(beta) * np.cos(gamma) - np.cos(alpha) * np.sin(gamma)
    r_mat[2, 0] = - np.sin(beta)
    r_mat[2, 1] = np.cos(beta) * np.sin(gamma)
    r_mat[2, 2] = np.cos(beta) * np.cos(gamma)
    return r_mat


def phi_def(x, y):
    """returns the arg of a point in a plane with the branch cut
       at x=0 towards x=infinity
       essentially its a converter from cartesian to polar coordinates
       can made better but at least it works correctly 
    """
    if x != 0 and x > 0 and y >= 0:
        phi_i = np.arctan(y / x)
    if x != 0 and x < 0 <= y:
        phi_i = np.pi + np.arctan(y / x)
    if x != 0 and x < 0 and y < 0:
        phi_i = np.pi + np.arctan(y / x)
    if x != 0 and x > 0 > y:
        phi_i = 2 * np.pi + np.arctan(y / x)
    if x == 0 and y >= 0:
        phi_i = np.pi / 2
    if x == 0 and y < 0:
        phi_i = 3 * np.pi / 2

    return phi_i


def dist_2d(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

def norm(a):
    return np.sqrt(np.sum(a ** 2))

def loc_sub_con_coord(sub_con_seg, detect_point):
    """Takes in and out nodes of a sub connection and calculates cylinder
    coordinates zC and rc for the detector position
    the form is required by the analytical solution
    in BC_wire_DC in module obs_calc_analytical

    returns zC, rc and the tangential vector in the detector position
    the tangential vector is defined so that it points along the magnetic field 
    vector for a current flow from p1 to p2
    """
    # print(sub_con_seg)
    p1 = np.array(sub_con_seg[0])
    p2 = np.array(sub_con_seg[1])
    p3 = detect_point  # np.array

    # L_sub_con = np.linalg.norm(p2-p1)
    l_wire = np.linalg.norm(p2 - p1)  # length of the wire
    OL = 0.5 * (p2 - p1)  # center of the wire,
    # local coordinate origin
    b_vec = (p2 - p1) / np.linalg.norm(p2 - p1)  # vector in fil direction
    a_vec = p1
    lamb = np.dot(b_vec, (p3 - a_vec))
    S = p1 + lamb * b_vec  # lot on the wire from detector

    rc = np.linalg.norm(p3 - S)
    zc = np.linalg.norm(S - OL)
    tang_vec = np.cross(p2 - p1, p3 - S) / np.linalg.norm(np.cross(p2 - p1, p3 - S))

    return rc, zc, l_wire, tang_vec


# filaments

def cuboid(segment_list, node_list, ns, ne, w, h, nhinc, nwinc, sigma, name):
    """defines a cuboid conducting wire between nodes Ns and Ne (3-tuples) 
    w, h are the width and height of the wire
    Ns and Ne sit on the area focus of the w*h surfaces
    nhinc, nwinc define the amount of sub filaments the wire is divided
    during simulation
    """
    # list to store relevant params
    # for naming the cuboid in FH
    segment_list.append([ns, ne, w, h, nhinc, nwinc, sigma, name])
    node_list.append(ns)
    node_list.append(ne)


def damage_sort_line(r_node, rd_o, rd_i, phi_e, phi_i, phi_s, segment_list,
                node_list, ns, ne, w, h, nhinc, nwinc, sigma, sigma_damage, name):

    if phi_e > phi_i > phi_s and rd_i < r_node < rd_o :
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc,
                  nwinc, sigma_damage, name)  # defines the cuboids and appends them to list
    elif phi_s > phi_e > phi_i >= 0 and rd_i < r_node < rd_o :
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc, nwinc, sigma_damage, name)
    elif phi_e < phi_s < phi_i < 2 * np.pi and rd_i < r_node < rd_o :
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc, nwinc, sigma_damage, name)
    else:
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc, nwinc, sigma, name)

def damage_sort_patch(r_node, rd_o, rd_i, sub_ind, phi_e, phi_i, phi_s, segment_list,
                node_list, ns, ne, w, h, nhinc, nwinc, sigma, sigma_damage, name):

    if phi_e > phi_i > phi_s and rd_i < r_node < rd_o and sub_ind == 1:
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc,
                  nwinc, sigma_damage, name)  # defines the cuboids and appends them to list
    elif phi_s > phi_e > phi_i >= 0 and rd_i < r_node < rd_o and sub_ind == 1:
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc, nwinc, sigma_damage, name)
    elif phi_e < phi_s < phi_i < 2 * np.pi and rd_i < r_node < rd_o and sub_ind == 1:
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc, nwinc, sigma_damage, name)
    else:
        cuboid(segment_list, node_list, ns, ne, w, h, nhinc, nwinc, sigma, name)

# loops

def loop(p, wl, hl, alpha, beta, gamma, loop_fil_params, loop_list):
    """Creates a square loop with center at point p = (x, y, z)
       initial position is in the x-y plane
       side that is hl wide is parallel to x axis
       side that is wl wide is parallel to y axis
       can be rotated by yaw, pitch and roll angles around point p

    """
    wf, hf, nhinc_f, nwinc_f, sigma_l = loop_fil_params

    loop_list.append([p, wl, hl, alpha, beta, gamma, sigma_l, wf, hf, nhinc_f, nwinc_f])


def loop_corners(loop):
    """returns the corners of a loop, these will be nodes in FH
    """
    p, wl, hl, alpha, beta, gamma, sigma_l, wf, hf, nhinc_f, nwinc_f = loop
    p1 = np.array([hl / 2, -wl / 2, 0])
    p2 = np.array([hl / 2, wl / 2, 0])
    p3 = np.array([-hl / 2, wl / 2, 0])
    p4 = np.array([-hl / 2, -wl / 2, 0])

    r_mat = r_ypr(alpha, beta, gamma)  # rotation matrix
    q1 = p + np.matmul(r_mat, p1)  # rotate the initial points
    q2 = p + np.matmul(r_mat, p2)
    q3 = p + np.matmul(r_mat, p3)
    q4 = p + np.matmul(r_mat, p4)
    return q1, q2, q3, q4

def det_loop_corners(loop):
    """returns the node points of a detector loop
    """
    p, wl, hl, alpha, beta, gamma, sigma_l, wf, hf, nhinc_f, nwinc_f = loop
    p1 = np.array([hl / 2, -wl / 2, 0])
    p2 = np.array([hl / 2, wl / 2, 0])
    p3 = np.array([-hl / 2, wl / 2, 0])
    p4 = np.array([-hl / 2, -wl / 2, 0])
    p5 = np.array([0.999 * hl / 2,   -wl / 2, 0])
    r_mat = r_ypr(alpha, beta, gamma)  # rotation matrix
    q1 = p + np.matmul(r_mat, p1)  # rotate the initial points
    q2 = p + np.matmul(r_mat, p2)
    q3 = p + np.matmul(r_mat, p3)
    q4 = p + np.matmul(r_mat, p4)
    q5 = p + np.matmul(r_mat, p5)
    return q1, q2, q3, q4, q5

# dyn_mesh specific code

def points_in_ring_dyn(ro_t, ro_shell, ri_shell, pt_dens):
    """determines points on a grid with spacing Ro/pt_dens
       that lie inside a ring of width ro_shell-ri_shell
       @param ro_t: 
       @param ri_shell: 
       @param ro_shell: 
       @type pt_dens: int
    """
    wh = 0.5 * ro_t / pt_dens

    if ri_shell > 0.0:
        for x, y in product(np.linspace(wh, ro_t - wh, pt_dens, endpoint=True),
                            np.linspace(wh, ro_t - wh, pt_dens, endpoint=True)):

            if (x ** 2 + y ** 2 < ro_shell ** 2) and (x ** 2 + y ** 2 > ri_shell ** 2):
                yield from {(x, -y), (x, y), (-x, y), (-x, -y)}

    elif ri_shell == 0.0:
        for x, y in product(np.linspace(0.0, ro_t, pt_dens, endpoint=True),
                            repeat=2):

            if x ** 2 + y ** 2 < ro_shell ** 2:
                yield from {(x, -y), (x, y), (-x, y), (-x, -y)}


def points_ring_cover(ro_t, ro_shell, ri_shell, pt_dens):
    """determines points on a grid with spacing ro_t/pt_dens
       patches with centers on these points and width height of wh
       cover a ring of width ro_shell-ri_shell
    """
    wh = 0.5 * ro_t / pt_dens

    if ri_shell > 0.0:
        for x, y in product(np.linspace(wh, ro_t - wh, pt_dens, endpoint=True),
                            np.linspace(wh, ro_t - wh, pt_dens, endpoint=True)):

            if ((x - wh) ** 2 + (y - wh) ** 2 <= ro_shell ** 2) and ((x + wh) ** 2 + (y + wh) ** 2 > ri_shell ** 2):
                yield from {(x, y)}


def patch_refine(cover_patches, r_sub_vec, node_dens_vec, ro_t, edge_type):
    """sub divides the patches covering the inner or the outer edge of the ring
    """
    # outer_node_dens=node_dens_vec[1+1]
    # node_dens=node_dens_vec[1]
    # inner_node_dens=node_dens_vec[0]

    ri_shell = r_sub_vec[1]
    ro_shell = r_sub_vec[2]
    # can get rid of the edge type later!
    # can be defined by the ri, ro index!???
    node_dens = node_dens_vec[1]
    if edge_type == 'i':
        node_dens_sub = node_dens_vec[0]
    elif edge_type == 'o':
        node_dens_sub = node_dens_vec[2]

    refined_cover_patches = []
    inner_patches_from_cover = []
    for patch in cover_patches:
        x, y = patch
        whc = 0.5 * ro_t / node_dens
        whs = 0.5 * ro_t / node_dens_sub
        ax = np.linspace(x - whc + whs, x + whc - whs, int(whc / whs))
        ay = np.linspace(y - whc + whs, y + whc - whs, int(whc / whs))

        sub_patches = list(product(ax, ay))

        sub_patches_removed = 0
        for sub_patch in reversed(sub_patches):
            xs, ys = sub_patch
            if edge_type == 'i':
                if xs ** 2 + ys ** 2 <= ri_shell ** 2 or xs ** 2 + ys ** 2 >= ro_shell ** 2:
                    sub_patches_removed = sub_patches_removed + 1
                    sub_patches.remove(sub_patch)
            if edge_type == 'o':
                if xs ** 2 + ys ** 2 > ro_shell ** 2 or xs ** 2 + ys ** 2 < ri_shell ** 2:
                    sub_patches_removed = sub_patches_removed + 1
                    sub_patches.remove(sub_patch)

        # when all centers of the sub patches are inside the ring we keep the
        # patch to avoid unnecessary refinement
        if sub_patches_removed == 0:
            sub_patches = []
            inner_patches_from_cover.append(patch)

        refined_cover_patches.extend(sub_patches)
    return refined_cover_patches, inner_patches_from_cover


def patches_sector_to_ring(patches):
    """takes patches in the 1st quadrant and mirrors
       them to all other quadrants
    """
    full_ring_patches = []
    if len(patches) != 0:
        refine_sec1_array = np.array(patches)
        refine_sec2 = list(map(tuple, refine_sec1_array * [-1, 1]))
        refine_sec3 = list(map(tuple, refine_sec1_array * [-1, -1]))
        refine_sec4 = list(map(tuple, refine_sec1_array * [1, -1]))

        full_ring_patches.extend(patches)
        full_ring_patches.extend(refine_sec2)
        full_ring_patches.extend(refine_sec3)
        full_ring_patches.extend(refine_sec4)
    return full_ring_patches


def ring_cover_classify(ro_t, r_sub_vec, node_dens_vec):
    """
        CAUTION FOR NOW WORKS ONLY FOR ONE MIDDLE RING!
        CAUTION FOR NOW WORKS ONLY FOR ONE MIDDLE RING!
       takes the patches that cover a ring in the first quadrant and separates
       them into 3 categories, the ones:
                 - covering the inner radius
                 - lying fully inside ring
                 - covering the outer radius
       then it refines the ones on the radii to match the node density of the 
       surrounding rings to make a transition between potential different mesh sizes
       on the different rings
    """
    ri_shell = r_sub_vec[1]
    ro_shell = r_sub_vec[2]
    inner_edge_cover = []
    inner_patches = []
    outer_edge_cover = []

    node_dens = node_dens_vec[1]
    patches = list(points_ring_cover(ro_t, ro_shell, ri_shell, node_dens))

    wh = 0.5 * ro_t / node_dens
    for patch in patches:
        x, y = patch
        if ((x + wh) ** 2 + (y + wh) ** 2 <= ro_shell ** 2) and ((x - wh) ** 2 + (y - wh) ** 2 > ri_shell ** 2):
            inner_patches.append(patch)
        elif ((x - wh) ** 2 + (y - wh) ** 2 <= ro_shell ** 2) and ((x + wh) ** 2 + (y + wh) ** 2 > ro_shell ** 2):
            outer_edge_cover.append(patch)
        elif ((x - wh) ** 2 + (y - wh) ** 2 <= ri_shell ** 2) and ((x + wh) ** 2 + (y + wh) ** 2 > ri_shell ** 2):
            inner_edge_cover.append(patch)

    outer_node_dens = node_dens_vec[1 + 1]
    node_dens = node_dens_vec[1]
    inner_node_dens = node_dens_vec[0]

    refined_inner_cover_patches, inner_patches_from_inner_cover = \
        patch_refine(inner_edge_cover, r_sub_vec, node_dens_vec, ro_t, 'i')

    refined_outer_cover_patches, inner_patches_from_outer_cover = \
        patch_refine(outer_edge_cover, r_sub_vec, node_dens_vec, ro_t, 'o')

    inner_patches.extend(inner_patches_from_inner_cover)
    inner_patches.extend(inner_patches_from_outer_cover)

    refined_inner_cover_patches = patches_sector_to_ring(refined_inner_cover_patches)
    refined_outer_cover_patches = patches_sector_to_ring(refined_outer_cover_patches)
    inner_patches = patches_sector_to_ring(inner_patches)

    return refined_inner_cover_patches, inner_patches, refined_outer_cover_patches
