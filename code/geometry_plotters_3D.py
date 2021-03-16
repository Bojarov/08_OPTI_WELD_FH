import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import code.geometry_helpers as gh
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def data_for_cylinder_along_z(center_x, center_y, radius, height_z):
    z = np.linspace(0, height_z, 50)
    theta = np.linspace(0, 2 * np.pi, 50)
    theta_grid, z_grid = np.meshgrid(theta, z)
    x_grid = radius * np.cos(theta_grid) + center_x
    y_grid = radius * np.sin(theta_grid) + center_y
    return x_grid, y_grid, z_grid


def plane_verts(plane):
    """Plots the planes with contacts
    """
    P1 = plane[0]
    P2 = plane[1]
    P3 = plane[2]
    P4 = tuple(np.array(P3) - (np.array(P2) - np.array(P1)))
    thick = plane[3]

    P1P2 = np.array(P1) - np.array(P2)
    P2P3 = np.array(P2) - np.array(P3)
    normal = np.cross(P1P2, P2P3) / gh.norm(np.cross(P1P2, P2P3))

    P1v_t = tuple(P1 + 0.5 * thick * normal)
    P2v_t = tuple(P2 + 0.5 * thick * normal)
    P3v_t = tuple(P3 + 0.5 * thick * normal)
    P4v_t = tuple(np.array(P3) - (np.array(P2) - np.array(P1)) + 0.5 * thick * normal)

    P1v_c = P1
    P2v_c = P2
    P3v_c = P3
    P4v_c = tuple(np.array(P3) - (np.array(P2) - np.array(P1)))

    P1v_b = tuple(P1 - 0.5 * thick * normal)
    P2v_b = tuple(P2 - 0.5 * thick * normal)
    P3v_b = tuple(P3 - 0.5 * thick * normal)
    P4v_b = tuple(np.array(P3) - (np.array(P2) - np.array(P1)) - 0.5 * thick * normal)

    verts_t = [[P1v_t, P2v_t, P3v_t, P4v_t]]
    verts_c = [[P1v_c, P2v_c, P3v_c, P4v_c]]
    verts_b = [[P1v_b, P2v_b, P3v_b, P4v_b]]

    return verts_t, verts_c, verts_b


def ZC_viso(geo_objects, view_pos, dist):
    x, y, z = view_pos

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('$x$', fontsize=20, rotation=0)
    ax.set_ylabel('$y$', fontsize=20)
    ax.set_zlabel('$z$', fontsize=20, rotation=0)

    ax.set_xlim3d(x - dist, x + dist)
    ax.set_ylim3d(y - dist, y + dist)
    ax.set_zlim3d(z - dist, z + dist)

    passive_loop_list = geo_objects['passive_loops']
    det_loop_list = geo_objects['det_loops']
    wire_list = geo_objects['wires']
    plane_list = geo_objects['planes']
    circular_loops = geo_objects['circ_pass_loops']

    for loop in passive_loop_list:
        # plot the node points
        Q_list = list(gh.det_loop_corners(loop))
        ax.scatter(*zip(*Q_list), color='blue', s=2)
        # connect the nodes to represent the loop
        ax.plot(*zip(*Q_list), color='blue')
        ax.scatter(Q_list[-1][0], Q_list[-1][1], Q_list[-1][2], color='black', marker='x', s=2)

    for loop in det_loop_list:
        # plot the node points
        Q_list = list(gh.det_loop_corners(loop))
        ax.scatter(*zip(*Q_list), color='orange', s=2)
        # connect the nodes to represent the loop
        ax.plot(*zip(*Q_list), color='orange')
        ax.scatter(Q_list[-1][0], Q_list[-1][1], Q_list[-1][2], color='black', marker='x', s=5)

    for wire in wire_list:
        p1, p2, w_wire, h_wire, nhinc, nwinc, sigma, name = wire
        Q_list = list([p1, p2])
        ax.scatter(*zip(*Q_list), color='g')
        ax.plot(*zip(*Q_list), color='g')

    for plane in plane_list:
        verts_t, verts_c, verts_b = plane_verts(plane)

        plane_poly_t = Poly3DCollection(verts_t)
        plane_poly_t.set_alpha(0.1)
        plane_poly_t.set_edgecolor("black")
        plane_poly_c = Poly3DCollection(verts_c)
        plane_poly_c.set_alpha(0.3)
        plane_poly_c.set_color('red')
        plane_poly_c.set_edgecolor("black")
        plane_poly_b = Poly3DCollection(verts_b)
        plane_poly_b.set_alpha(0.1)
        plane_poly_b.set_edgecolor("black")
        ax.add_collection3d(plane_poly_t)
        ax.add_collection3d(plane_poly_c)
        ax.add_collection3d(plane_poly_b)

    for circ_loop in circular_loops:
        nodes = circ_loop["nodes"]
        ax.scatter(nodes[1:-1, 0], nodes[1:-1, 1], nodes[1:-1, 2], color='g')
        ax.plot(nodes[:, 0], nodes[:, 1], nodes[:, 2], color='g')
        ax.scatter(nodes[0, 0], nodes[0, 1], nodes[0, 2], color='r')
        ax.scatter(nodes[-1, 0], nodes[-1, 1], nodes[-1, 2], color='r')

        # seg_centers = circ_loop["seg_params"][0]
        # seg_w_vec = circ_loop["seg_params"][1]
        # ax.scatter(seg_centers[:, 0], seg_centers[:, 1], seg_centers[:, 2], color='red')
        # ax.quiver(seg_centers[:, 0], seg_centers[:, 1], seg_centers[:, 2],
        #          seg_w_vec[:, 0], seg_w_vec[:, 1], seg_w_vec[:, 2])

