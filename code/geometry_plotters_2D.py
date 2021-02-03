import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import code.geometry_builders as gb
from matplotlib import rc, font_manager

sizeOfFont = 12
fontProperties = {'family': 'sans-serif', 'sans-serif': ['Helvetica'],
                  'weight': 'normal', 'size': sizeOfFont}
ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
                                         size=sizeOfFont, weight='normal', stretch='normal')
rc('text', usetex=True)
rc('font', **fontProperties)


def wire_mesh_2D_plot(params, damage_params, segment_list, pts):
    """
    make a 2D plot of the mesh in the cross section of the wire/tube
    for better understanding
    input: list of segments forming the tube which is build in 
           geometry_builders
    """

    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params  # unpack

    rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

    fig_mesh2d, ax_mesh2D = plt.subplots()

    ax_mesh2D.set_aspect('equal')
    plt.xlim(-(ro_t + 0.1 * ro_t), (ro_t + 0.1 * ro_t))
    plt.ylim(-(ro_t + 0.1 * ro_t), (ro_t + 0.1 * ro_t))

    # plot outer and inner surface of the pipe

    circle_o2 = plt.Circle((0, 0), ro_t, color='b', fill=False)
    circle_i2 = plt.Circle((0, 0), ri_t, color='b', fill=False)
    ax_mesh2D.add_patch(circle_o2)
    ax_mesh2D.add_patch(circle_i2)

    for i in pts:
        j = pts.index(i)
        sigma_segment = segment_list[j][6]
        if sigma_segment == sigma_damage:
            color = 'red'
        else:
            color = 'blue'

        rect = patches.Rectangle(
            tuple(np.subtract(i, (0.5 * ro_t / (node_dens - 1), 0.5 * ro_t / (node_dens - 1)))),
            ro_t / (node_dens - 1), ro_t / (node_dens - 1), linewidth=1.5 * (ro_t - ri_t) / node_dens,
            edgecolor='black', facecolor=color, alpha=0.5
        )

        ax_mesh2D.add_patch(rect)

    plt.scatter(*zip(*pts), marker='o', s=10, color='green', zorder=2)
    ax_mesh2D.title.set_text(r'Conductor cross section and filaments')
    ax_mesh2D.set_xlabel(r'$x[m]$')
    ax_mesh2D.set_ylabel('$y[m]$')

    plt.tight_layout()

    plt.grid()

    plt.show()


def wire_mesh_2D_plot_dyn(ro_t, ri_t, r_sub_vec, node_dens_vec, params,
                          filament_params, damage_params, l_sub_vec):
    """
    make a 2D plot of the mesh in the cross section of the wire
    for better understanding
    """

    # hierarchy of geometry objects
    # tube->shell->segments->nodes
    tube_node_lists = []  # lists to carry the node and segment lists for each shell
    tube_segment_lists = []
    tube_pts_lists = []

    # build the geometry in python from the input params
    gb.tube_builder(ro_t, ri_t, r_sub_vec, l_sub_vec, node_dens_vec,
                    params, damage_params, filament_params, tube_node_lists,
                    tube_segment_lists, tube_pts_lists)

    ro_t, ri_t, flen, node_dens, sigma, mu_0, mu_r, freq = params

    rd_i, rd_o, d_weld, l_weld, phi_c_weld, sigma_damage = damage_params

    n_segments = len(tube_segment_lists[0]) + len(tube_segment_lists[1]) + len(tube_segment_lists[2])

    fig2, ax2 = plt.subplots()

    ax2.set_aspect('equal')
    plt.xlim(-(ro_t + 0.1 * ro_t), (ro_t + 0.1 * ro_t))
    plt.ylim(-(ro_t + 0.1 * ro_t), (ro_t + 0.1 * ro_t))

    shell_colors = ['blue', 'green']

    # plot outer and inner surface of the pipe
    for i in r_sub_vec:
        circle = plt.Circle((0, 0), i, color='b', fill=False)
        ax2.add_patch(circle)

    for pts in tube_pts_lists:
        i = tube_pts_lists.index(pts)

        segment_list = tube_segment_lists[i]
        node_dens = node_dens_vec[i]

        for pt in pts:
            j = pts.index(pt)
            #
            # if j==0:
            # print(len(l_sub_vec))
            # print(len(pts))
            # print(len(segment_list)/3)
            #    print(len(segment_list))
            #    print(segment_list[117])
            # print(j*(len(l_sub_vec)-1))

            fil_ind = (j * (len(l_sub_vec) - 1))

            # sigma_segment = segment_list[j][6]
            # w = segment_list[j][2]
            sigma_segment = segment_list[fil_ind][6]
            w = segment_list[fil_ind][2]
            if sigma_segment == sigma_damage and sigma_damage != sigma:
                color = 'red'
            else:
                color = shell_colors[i % 2]

            rect = patches.Rectangle(
                tuple(np.subtract(pt, (0.5 * w, 0.5 * w))),
                w, w, linewidth=1.5 * (ro_t - ri_t) / node_dens,
                edgecolor='black', facecolor=color, alpha=0.5
            )

            ax2.add_patch(rect)

        if len(pts) > 0:
            plt.scatter(*zip(*pts), marker='o', s=10, color=shell_colors[i % 2], zorder=2)

    ax2.title.set_text(r'Conductor cross section and filaments')
    ax2.set_xlabel(r'$x[m]$')
    ax2.set_ylabel('$y[m]$')

    plt.tight_layout()

    plt.grid()

    print("The geometry is made of " + str(n_segments) + " segments.")
