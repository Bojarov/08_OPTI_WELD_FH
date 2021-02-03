import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import numpy.cos as cos
#import numpy.sin as sin
def R_ypr(alpha, beta, gamma):
    """Rotation matrix - yaw/pitch/roll angles
    """
    R_mat=np.zeros((3,3))
    R_mat[0,0] =   np.cos(alpha)*np.cos(beta)
    R_mat[0,1] =   np.cos(alpha)*np.sin(beta)*np.sin(gamma) - np.sin(alpha)*np.cos(gamma)
    R_mat[0,2] =   np.cos(alpha)*np.sin(beta)*np.cos(gamma) + np.sin(alpha)*np.sin(gamma)
    R_mat[1,0] =   np.sin(alpha)*np.cos(beta)
    R_mat[1,1] =   np.sin(alpha)*np.sin(beta)*np.sin(gamma) + np.cos(alpha)*np.cos(gamma)
    R_mat[1,2] =   np.sin(alpha)*np.sin(beta)*np.cos(gamma) - np.cos(alpha)*np.sin(gamma)
    R_mat[2,0] = - np.sin(beta)
    R_mat[2,1] =   np.cos(beta)*np.sin(gamma)
    R_mat[2,2] =   np.cos(beta)*np.cos(gamma)
    return R_mat

def data_for_cylinder_along_z(center_x, center_y, radius, z_1, z_2):
    z = np.linspace(z_1, z_2, 50)
    theta = np.linspace(0, 2*np.pi, 50)
    theta_grid, z_grid=np.meshgrid(theta, z)
    x_grid = radius*np.cos(theta_grid) + center_x
    y_grid = radius*np.sin(theta_grid) + center_y
    return x_grid,y_grid,z_grid

def loop(P, w, h, alpha, beta, gamma, loop_fil_params, loop_list):
    """Creates a square loop with center at point P = (x, y, z)
       initial position is in the x-y plane
       side that is h wide is parallel to x axis
       side that is w wide is parallel to y axis
       can be rotated by yaw, pitch and roll angles around point P

    """
    wf, hf, nhinc_f, nwinc_f, sigmal = loop_fil_params

    loop_list.append([P, w, h, alpha, beta, gamma, wf, hf, nhinc_f, nwinc_f])
    
def loop_builder(loop):
    """returns the corners of a loop, these will be nodes in FH
    """
    P, w, h, alpha, beta, gamma, wf, hf, nhinc_f, nwinc_f = loop
    P1 = np.array([ h/2, -w/2 , 0])
    P2 = np.array([ h/2,  w/2 , 0])
    P3 = np.array([-h/2,  w/2 , 0])
    P4 = np.array([-h/2, -w/2 , 0])


    R_mat = R_ypr(alpha, beta, gamma)       #rotation matrix
    Q1 = P + np.matmul(R_mat, P1)                #rotate the initial points
    Q2 = P + np.matmul(R_mat, P2)
    Q3 = P + np.matmul(R_mat, P3)
    Q4 = P + np.matmul(R_mat, P4)

    print(Q1)
    print(Q2)
    print(Q3)
    print(Q4)
    return Q1, Q2, Q3, Q4

#cylinder representing the tube:
x1   =  0.0
y1   =  0.0
flen = 1.0
Ro   =  0.1
thick = 0.01
Ri=Ro-thick

Pl1 = np.array([0,0,0.05])
Pl2 = np.array([0.06,0,0.05])
wl = 0.1*0.5
hl = 0.1*0.5


#loop parameters
alpha = np.pi*0.0
beta  = np.pi*0.0
gamma = np.pi*0.0


#parameters to define the filament properties of a loop
wf = 0.01
hf = 0.01
nhinc_f = 1     #number of FH internal sub divisions
nwinc_f = 1
sigmal = 100


loop_fil_params = [wf, hf, nhinc_f, nwinc_f, sigmal]

loop_list = []


#loop(Pl1, wl, hl, alpha, beta, gamma, loop_fil_params, loop_list)
#loop(Pl2, wl, hl, alpha, beta, gamma, loop_fil_params, loop_list)

#loop1 = loop_list[0]
#loop2 = loop_list[1]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

Xo1,Yo1,Zo1 = data_for_cylinder_along_z(x1, y1, Ro, 0, flen*0.4)
ax.plot_surface(Xo1, Yo1, Zo1, alpha=0.25, color='green')
Xo2,Yo2,Zo2 = data_for_cylinder_along_z(x1, y1, Ro, flen*0.4, flen*0.6)
ax.plot_surface(Xo2, Yo2, Zo2, alpha=0.25, color ='red')
Xo3,Yo3,Zo3 = data_for_cylinder_along_z(x1, y1, Ro, flen*0.6, flen)
ax.plot_surface(Xo3, Yo3, Zo3, alpha=0.25, color='green')

Xi1,Yi1,Zi1 = data_for_cylinder_along_z(x1, y1, Ri, 0, flen*0.4)
ax.plot_surface(Xi1, Yi1, Zi1, alpha=0.25, color='green')
Xi2,Yi2,Zi2 = data_for_cylinder_along_z(x1, y1, Ri, flen*0.4, flen*0.6)
ax.plot_surface(Xi2, Yi2, Zi2, alpha=0.25, color='red')
Xi3,Yi3,Zi3 = data_for_cylinder_along_z(x1, y1, Ri, flen*0.6, flen)
ax.plot_surface(Xi3, Yi3, Zi3, alpha=0.25, color = 'green')


ax.set_xlabel('$x$', fontsize=20, rotation=0)
ax.set_ylabel('$y$', fontsize=20)
ax.set_zlabel('$z$', fontsize=20, rotation=0)
#ax.set_zlim3d(-0.0, 0.1)                    # viewrange for z-axis should be [-4,4] 
#ax.set_ylim3d(-0.1, 0.1)                    # viewrange for y-axis should be [-2,2] 
ax.set_xlim3d(-0.5, 0.5)                    # viewrange for x-axis should be [-2,2] 
ax.set_ylim3d(-0.5, 0.5)                    # viewrange for x-axis should be [-2,2] 
#ax.set_aspect('auto')
#plot the node points
#Q_list1 = list(loop_builder(loop1))
#Q_list2 = list(loop_builder(loop2))
#ax.scatter(*zip(*Q_list1))
#ax.scatter(*zip(*Q_list2))
#connect the nodes to represent the loop
#Q_list1.append(Q_list1[0])
#Q_list2.append(Q_list2[0])

#ax.plot(*zip(*Q_list1), color='b')
#ax.plot(*zip(*Q_list2), color='orange')

plt.show()
