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

    passive_loop_list = geo_objects['pass_loops']
    det_loop_list = geo_objects['det_loops']
    wire_list = geo_objects['wires']
    #geo_objects['circ_pass_loops']

    #r_loops = geo_objects['r_pass_loops']

    #for loop in passive_loop_list:
    #    # plot the node points
    #    Q_list = list(gh.det_loop_corners(loop))
    #    ax.scatter(*zip(*Q_list), color='blue', s=2)
    #    # connect the nodes to represent the loop
    #    ax.plot(*zip(*Q_list), color='blue')
    #    ax.scatter(Q_list[-1][0], Q_list[-1][1], Q_list[-1][2], color='black', marker='x', s=2)

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

    if 'planes' in geo_objects:
        plane_list = geo_objects['planes']

        for plane in plane_list:
            corners = gh.corners(plane['build_params'])
            verts = gh.plane_verts(corners)
            plane_poly = Poly3DCollection(verts)
            plane_poly.set_alpha(0.2)
            plane_poly.set_color('blue')
            plane_poly.set_edgecolor("black")
            ax.add_collection3d(plane_poly)

    if 'pass_loops' in geo_objects:
        pass_loops = geo_objects['pass_loops']
        for pass_loop in pass_loops:
            nodes = pass_loop["nodes"]
            ax.scatter(nodes[1:-1, 0], nodes[1:-1, 1], nodes[1:-1, 2], color='g')
            ax.plot(nodes[:, 0], nodes[:, 1], nodes[:, 2], color='g')
            ax.scatter(nodes[0, 0], nodes[0, 1], nodes[0, 2], color='r')
            ax.scatter(nodes[-1, 0], nodes[-1, 1], nodes[-1, 2], color='r')
    # TODO add segment viso if requested


