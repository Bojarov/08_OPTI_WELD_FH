import numpy as np
#import obs_calc.I_FH_calc
import obs_calc.I_Knight
#import sys, os
import os.path
import matplotlib.pyplot as plt


"""
THIS MODULE IS DAMAGED!
THIS MODULE IS DAMAGED!
THIS MODULE IS DAMAGED!
THIS MODULE IS DAMAGED!
THIS MODULE IS DAMAGED!
THIS MODULE IS DAMAGED!
"""






def A_xyz0_calc(x_win, y_win, z0, B_pts, xyz0_proj_re, xyz0_proj_im, params):
    """A field in a x_win*y_win window above the wire where in the z0 plane
       the field is calculated on a grid of B_pts**2 points
    """
    Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq = params   #unpack

    x_val=np.linspace( -x_win/2,  x_win/2, B_pts)
    y_val=np.linspace( -y_win/2,  y_win/2, B_pts)
    hx = abs(x_val[1]-x_val[0])
    hy = abs(y_val[1]-y_val[0])

    w=Ro/(node_dens - 1)
    h=Ro/(node_dens - 1)


    #check if plot data has been already produced
    filenameAzre='Axy_re'+str(len(x_val))+'_'+str(len(y_val))+'_z='+str(z0)+\
                str(Ro)+"_Ri_"+str(Ri)+"_nd_"+str(node_dens)+\
               "_l_"+str(flen)+"_sigm"+str(sigma)+"_mur_"+str(mu_r)+"_f_"+\
               str(freq)+'.npy' 

    filenameAzim='Axy_im'+str(len(x_val))+'_'+str(len(y_val))+'_z='+str(z0)+\
                str(Ro)+"_Ri_"+str(Ri)+"_nd_"+str(node_dens)+\
               "_l_"+str(flen)+"_sigm"+str(sigma)+"_mur_"+str(mu_r)+"_f_"+\
               str(freq)+'.npy' 

    filenameAzre_err='Axy_re_err'+str(len(x_val))+'_'+str(len(y_val))+'_z='+str(z0)+\
                str(Ro)+"_Ri_"+str(Ri)+"_nd_"+str(node_dens)+\
               "_l_"+str(flen)+"_sigm"+str(sigma)+"_mur_"+str(mu_r)+"_f_"+\
               str(freq)+'.npy' 

    filenameAzim_err='Axy_im_err_'+str(len(x_val))+'_'+str(len(y_val))+'_z='+str(z0)+\
                str(Ro)+"_Ri_"+str(Ri)+"_nd_"+str(node_dens)+\
               "_l_"+str(flen)+"_sigm"+str(sigma)+"_mur_"+str(mu_r)+"_f_"+\
               str(freq)+'.npy' 


    if os.path.exists(filenameAzre):
        print ("using existing Az_xy files")
    
        with open(filenameAzre, 'rb') as f_Are:
    
            Axyz0_re = np.load(f_Are)
        with open(filenameAzim, 'rb') as f_Aim:
    
            Axyz0_im = np.load(f_Aim)
        with open(filenameAzre_err, 'rb') as f_Are_err:
    
            Axyz0_re_err = np.load(f_Are_err)
        with open(filenameAzim_err, 'rb') as f_Aim_err:
    
            Axyz0_im_err = np.load(f_Aim_err)
    else:
        print ("Az_xy files do not yet exist")
        Axyz0_re = np.zeros((B_pts, B_pts, 3))
        Axyz0_im = np.zeros((B_pts, B_pts, 3))
        Axyz0_re_err = np.zeros((B_pts, B_pts, 3))
        Axyz0_im_err = np.zeros((B_pts, B_pts, 3))



        for i in range(B_pts):
            for j in range(B_pts):
                print(i,j)
                #if(x_val[i]**2+y_val[j]**2>Ro**2):
                if(np.sqrt(np.round(x_val[i]**2,14)+np.round(y_val[j]**2,14))>=Ro):

                #if(abs(x_val[i])>w/2.0 or abs(y_val[j])>h/2.0):


                    Axyz0_re[i,j,:], Axyz0_im[i,j,:],\
                    Axyz0_re_err[i,j,:], Axyz0_im_err[i,j,:] =\
                        A_int(x_val[i], y_val[j], z0,\
                        xyz0_proj_re,xyz0_proj_im, params)


        #Q:!!!why is z0 irrelevant here?????
        #print("z0 problem here in B_xyz0_calc!")

        with open(filenameAzre, 'wb') as Arefile:
            np.save(Arefile, Axyz0_re)
        with open(filenameAzim, 'wb') as Aimfile:
            np.save(Aimfile, Axyz0_im)
    
        with open(filenameAzre_err, 'wb') as Arefile_err:
            np.save(Arefile_err, Axyz0_re_err)
        with open(filenameAzim_err, 'wb') as Aimfile_err:
            np.save(Aimfile_err, Axyz0_im_err)

    return Axyz0_re, Axyz0_im, Axyz0_re_err, Axyz0_im_err

def B_calc_plane(B_pts, x_win, y_win, z0, xyz0_proj_re, xyz0_proj_im, params):
    """Calculates B = curl A in a the xy-plane at z=z0
    """
    Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq = params   #unpack params
    Axyz0_re, Axyz0_im, Axyz0_re_err, Axyz0_im_err =\
    A_xyz0_calc(x_win, y_win, z0, B_pts, xyz0_proj_re, xyz0_proj_im, params)
    
    x_val=np.linspace( -x_win/2 , x_win/2, B_pts)
    y_val=np.linspace( -y_win/2 , y_win/2, B_pts)
    hx = abs(x_val[1]-x_val[0])
    hy = abs(y_val[1]-y_val[0])
    
    w=Ro/(node_dens - 1)
    h=Ro/(node_dens - 1)
    
    dxAz_re = np.zeros((np.shape(Axyz0_re)[0]-2,np.shape(Axyz0_re)[0]-2))
    dyAz_re = np.zeros((np.shape(Axyz0_re)[0]-2,np.shape(Axyz0_re)[0]-2))
    
    dxAz_im = np.zeros((np.shape(Axyz0_im)[0]-2,np.shape(Axyz0_im)[0]-2))
    dyAz_im = np.zeros((np.shape(Axyz0_im)[0]-2,np.shape(Axyz0_im)[0]-2))
    
    dxAz_tot = np.zeros((np.shape(Axyz0_re)[0]-2,np.shape(Axyz0_re)[0]-2))
    dyAz_tot = np.zeros((np.shape(Axyz0_re)[0]-2,np.shape(Axyz0_re)[0]-2))
    
    Az_xy_re=Axyz0_re[:,:,2]    # real part of A_z in the z=z0 plane
    Az_xy_im=Axyz0_im[:,:,2]    # im   part of A_z in the z=z0 plane

    Az_tot_xy=np.sqrt(Az_xy_re**2+Az_xy_im**2)

    Az_xy_re_err=Axyz0_re_err[:,:,2]    # error real part of A_z in the z=z0 plane
    Az_xy_im_err=Axyz0_im_err[:,:,2]    # error im   part of A_z in the z=z0 plane


    x_val_new=x_val[1:B_pts-1]
    y_val_new=y_val[1:B_pts-1]


    for i in np.arange(len(x_val)-2):
        for j in np.arange(len(y_val)-2):
            #if(abs(x_val_new[i])>w/2.0+1.0*hx  or abs(y_val_new[j])>h/2.0+1.0*hy ):
            if((abs(x_val_new[i]))**2+(abs(y_val_new[j]))**2 >(Ro+np.sqrt((hx/2)**2+(hy/2)**2))**2):

                dxAz_re[i,j]=\
                (Az_xy_re[ i + 2 , j + 1]   - Az_xy_re[ i     , j + 1 ])/(2*hx)
                dyAz_re[i,j]=\
                (Az_xy_re[ i + 1 , j + 2 ]  - Az_xy_re[ i + 1 , j     ])/(2*hy)
    
                dxAz_im[i,j]=\
                (Az_xy_im[ i + 2 , j + 1]   - Az_xy_im[ i     , j + 1 ])/(2*hx)
                dyAz_im[i,j]=\
                (Az_xy_im[ i + 1 , j + 2 ]  - Az_xy_im[ i + 1 , j     ])/(2*hy)

                dxAz_tot[i,j]=\
                (Az_tot_xy[ i + 2 , j + 1]   - Az_tot_xy[ i     , j + 1 ])/(2*hx)
                dyAz_tot[i,j]=\
                (Az_tot_xy[ i + 1 , j + 2 ]  - Az_tot_xy[ i + 1 , j     ])/(2*hy)


    return Az_tot_xy, dxAz_tot, dyAz_tot, Az_xy_re, Az_xy_im, dxAz_re, dxAz_im, dyAz_re, dyAz_im,\
            Az_xy_re_err, Az_xy_im_err

def B_plot_plane(B_pts, x_win, y_win, z0, pts, xyz0_proj_re, xyz0_proj_im, params):
    """Plots absolute value of the magntiude of A and B field around the wire
    in the x-y plane at z=z0
    """
    Ro, Ri, flen, node_dens, sigma, mu_0, mu_r, freq = params   #unpack params
    Az_tot_xy, dxAz_tot, dyAz_tot, Az_xy_re, Az_xy_im, dxAz_re, dxAz_im, dyAz_re, dyAz_im, \
    Az_xy_re_err, Az_xy_im_err =\
            B_calc_plane(B_pts, x_win, y_win, z0, xyz0_proj_re, xyz0_proj_im, params)


    x_val=np.linspace( -x_win/2 , x_win/2, B_pts)
    y_val=np.linspace( -y_win/2 , y_win/2, B_pts)
    hx=x_win/(B_pts-1)
    hy=y_win/(B_pts-1)

    fig = plt.figure(figsize=(12,3))

    ax = fig.add_subplot(131)

#    q1=np.min(Az_xy_re[np.nonzero(Az_xy_re)])
#    q2=np.max(Az_xy_re[np.nonzero(Az_xy_re)])
#    AZ_XY_RE=ax.pcolormesh(Az_xy_re,vmin=q1,vmax=q2)
#    ax.set_aspect('equal')
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes("right", size="5%", pad=0.05)
#    fig.colorbar(AZ_XY_RE, ax=ax, cax=cax)
#    ax.title.set_text(r'${A}_z(x,y,z=$'+str(z0)+r'$)$')
#    ax.set_xticks(np.linspace(0,B_pts,5))
#    ax.set_xticklabels(-x_win/2+np.linspace(0,B_pts,5)*x_win/B_pts)
#    ax.set_yticks(np.linspace(0,B_pts,5))
#    ax.set_yticklabels(-y_win/2+np.linspace(0,B_pts,5)*y_win/B_pts)
#    ax.set_xlabel('x[m]')
#    ax.set_ylabel('y[m]')
#
#
#    ax = fig.add_subplot(232)
#
#    q5=np.min(Az_xy_re_err[np.nonzero(Az_xy_re_err)])
#    q6=np.max(Az_xy_re_err[np.nonzero(Az_xy_re_err)])
#    AZ_XY_RE_ERR=ax.pcolormesh(Az_xy_re_err,vmin=q5,vmax=q6)
#    ax.set_aspect('equal')
#    divider = make_axes_locatable(ax)
#    cax = divider.append_axes("right", size="5%", pad=0.05)
#    fig.colorbar(AZ_XY_RE_ERR, ax=ax, cax=cax)
#    ax.title.set_text(r'$\delta{A}_z(x,y,z=$'+str(z0)+r'$)$')
#    ax.set_xticks(np.linspace(0,B_pts,5))
#    ax.set_xticklabels(-x_win/2+np.linspace(0,B_pts,5)*x_win/B_pts)
#    ax.set_yticks(np.linspace(0,B_pts,5))
#    ax.set_yticklabels(-y_win/2+np.linspace(0,B_pts,5)*y_win/B_pts)
#    ax.set_xlabel('x[m]')
#    ax.set_ylabel('y[m]')
#
#    ax = fig.add_subplot(233)
#
#
    q3=np.min(Az_tot_xy[np.nonzero(Az_tot_xy)])
    q4=np.max(Az_tot_xy[np.nonzero(Az_tot_xy)])
    AZ_XY_TOT=ax.pcolormesh(Az_tot_xy,vmin=q3,vmax=q4)
    ax.set_aspect('equal')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(AZ_XY_TOT, ax=ax, cax=cax)
    ax.title.set_text(r'${A}_z(x,y,z=$'+str(z0)+r'$)$')
    ax.set_xticks(np.linspace(0,B_pts,5))
    ax.set_xticklabels(-x_win/2+np.linspace(0,B_pts,5)*x_win/B_pts)
    ax.set_yticks(np.linspace(0,B_pts,5))
    ax.set_yticklabels(-y_win/2+np.linspace(0,B_pts,5)*y_win/B_pts)
    ax.set_xlabel('x[m]')
    ax.set_ylabel('y[m]')

    ax = fig.add_subplot(132)
    
    #B_mag_tot=np.sqrt(dxAz_tot**2+dyAz_tot**2)
    B_mag_tot=np.sqrt(dxAz_tot**2+dyAz_tot**2)
    ZZZ=ax.pcolormesh(B_mag_tot)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(ZZZ, ax=ax,cax=cax)
    arrow_dist=4
    
    #print("arrow_dist")
    #print(arrow_dist)
    for i in np.arange(len(x_val)-2):
        for j in np.arange(len(y_val)-2):
            if (i%int(arrow_dist)==0 and j %(int(arrow_dist))==0):
                ax.quiver( i+1-0.5 , j+1-0.5, dyAz_tot[i,j], -dxAz_tot[i,j] )


    ax.set_aspect('equal')
    ax.title.set_text(r'$\mathbf{B}(x,y,z=$'+str(z0)+r'$)$')

    ax.set_xticks(np.linspace(0,B_pts-2,5))
    ax.set_xticklabels(np.round(-x_win/2+hx+np.linspace(0,B_pts-2,5)*(x_win-2*hx)/(B_pts-2),2))
    ax.set_yticks(np.linspace(0,B_pts-2,5))
    ax.set_yticklabels(np.round(-y_win/2+hy+np.linspace(0,B_pts-2,5)*(y_win-2*hy)/(B_pts-2),2))
    #ax.ticklabel_format(axis='both', style='sci',scilimits=(0,0))
    ax.set_xlabel('x[m]')
    ax.set_ylabel('y[m]')


    ax = fig.add_subplot(133)
    
    B_diag=np.diagonal(B_mag_tot)[0:int(B_pts/2)]
    B_diag=B_diag[::-1]


    x_val_diag=x_val[int(B_pts/2):len(x_val)-1]
    y_val_diag=y_val[int(B_pts/2):len(y_val)-1]
    xy_diag=np.sqrt(x_val_diag**2+y_val_diag**2)

    #with HiddenPrints():
    I_tot_analytic=I_Knight_calc(params) 
    #I_FH = I_FH_calc(params, xyz0_proj_re, xyz0_proj_im)

    xy_diag_theo=np.linspace(Ro,max(xy_diag))

    B_theo_diag=mu_0*I_tot_analytic/(2*np.pi*xy_diag_theo)

    ax.plot(xy_diag,B_diag,label="num. FH for filament")
    ax.plot(xy_diag_theo,B_theo_diag,'--',label="theory for wire")
    ax.ticklabel_format(axis='y', style='sci',scilimits=(0,0))
    ax.title.set_text(r'${B}(x=y,z=$'+str(z0)+r'$)$')
    ax.set_xlabel('r[m]')
    ax.set_ylabel('B[T]')
    ax.set_aspect(aspect=1.0/(max(B_theo_diag)/max(xy_diag)))
    ax.legend(loc='upper right')

    #ax = fig.add_subplot(236)

    #ax.plot(np.array(xy_diag), B_diag/(mu_0*I_fil)*(2*np.pi))
    #ax.plot(np.linspace(0.1,max(xy_diag),50),1/np.linspace(0.1,max(xy_diag),50))
    plt.tight_layout()

