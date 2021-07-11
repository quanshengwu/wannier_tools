#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Changming Yue, yuechangming8@gmail.com


import numpy as np
import sys, os
from myplane import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mp3d

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def get_brillouin_zone_3d(cell):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """

    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (3, 3)

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        # WHY 13 ????
        # The Voronoi ridges/facets are perpendicular to the lines drawn between the
        # input points. The 14th input point is [0, 0, 0].
        if(pid[0] == 13 or pid[1] == 13):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges, bz_facets

# run it as " python 3D_bz.py Nodes.dat_symm  "
if __name__ == "__main__":
    cell = np.zeros((3,3),dtype=np.float64)
# change the cell to your own 
    cell[0]=[  -2.126699999999998,      2.126699999999998,      7.32050000000000] # a_vec
    cell[1]=[   2.126699999999998,     -2.126699999999998,      7.32050000000000] # b_vec
    cell[2]=[   2.126699999999998,      2.126699999999999,     -7.32050000000000] # c_vec

    if os.path.exists("POSCAR"):
        with open("POSCAR","r") as rd:
           rd.readline()
           rd.readline()
           for i in range(3):
               a1,a2,a3=rd.readline().strip().split()
               cell[i,0]=np.float64(a1) 
               cell[i,1]=np.float64(a2) 
               cell[i,2]=np.float64(a3) 

#   icell is reciprocal lattice 
    icell = 2.0*np.pi*np.linalg.inv(cell).T                
# 
#  -0.000000    1.477215    0.429150 # b1_vec
#   1.477215    0.000000    0.429150 # b2_vec
#   1.477215    1.477215    0.000000 # b3_vec
 
# length of b1 b2 b3
    b1, b2, b3 = np.linalg.norm(icell, axis=1)   

# vortices, edges and faces of brillouin zone
    v, e, f = get_brillouin_zone_3d(icell)

    fig = plt.figure( figsize=(6, 6), dpi=300)
    ax = plt.subplot(111, projection='3d',proj_type = 'ortho')

    for xx in e:
        ax.plot(xx[:, 0], xx[:, 1], xx[:, 2], color='k', lw=1.0)

#------- set kxky plane (kz=0) by finding all the vortices 
    pts=[]
    lens=[]
    phis=[]

    for xx in e:
        for i in range(len(xx)-1):  
            if abs(xx[i,2])<1.0E-6: # vortices
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*1.0+xx[i+1,0]*0.0, xx[i,1]*1.0+xx[i+1,1]*0.0, xx[i,2]*1.0+xx[i+1,2]*0.0)
		t=(xx[i,0]*1.0+xx[i+1,0]*0.0, xx[i,1]*1.0+xx[i+1,1]*0.0, xx[i,2]*1.0+xx[i+1,2]*0.0)
		c=xx[i,0]*1.0+xx[i+1,0]*0.0+(xx[i,1]*1.0+xx[i+1,1]*0.0)*1.0j
            elif abs(xx[i+1,2])<1.0E-6: # vortices
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*0.0+xx[i+1,0]*1.0, xx[i,1]*0.0+xx[i+1,1]*1.0, xx[i,2]*0.0+xx[i+1,2]*1.0)
  		t=(xx[i,0]*0.0+xx[i+1,0]*1.0, xx[i,1]*0.0+xx[i+1,1]*1.0, xx[i,2]*0.0+xx[i+1,2]*1.0)
		c=xx[i,0]*0.0+xx[i+1,0]*1.0 + (xx[i,1]*0.0+xx[i+1,1]*1.0)*1.0j
            elif abs(xx[i,2]*0.5+xx[i+1,2]*0.5)<1.0E-6: # centers 
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*0.5+xx[i+1,0]*0.5, xx[i,1]*0.5+xx[i+1,1]*0.5, xx[i,2]*0.5+xx[i+1,2]*0.5)
		t=(xx[i,0]*0.5+xx[i+1,0]*0.5, xx[i,1]*0.5+xx[i+1,1]*0.5, xx[i,2]*0.5+xx[i+1,2]*0.5)
		c=(xx[i,0]*0.5+xx[i+1,0]*0.5)+(xx[i,1]*0.5+xx[i+1,1]*0.5)*1.0j
	    else:
	        continue

            if s not in lens:
                pts.append(t)
                ang=np.angle(c)
                phis.append(ang)
                lens.append(s)

    phis_= np.sort(phis)
    kxky_plane=[]
    for i in phis_:
        idx=phis.index(i) 
        kxky_plane.append(pts[idx])	
    kxky_plane.append(kxky_plane[0])	

#------- set kykz plane (kx=0) by finding all the vortices 
    pts=[]
    lens=[]
    phis=[]

    for xx in e:
        for i in range(len(xx)-1):  
            if abs(xx[i,0])<1.0E-6: # vortices
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*1.0+xx[i+1,0]*0.0, xx[i,1]*1.0+xx[i+1,1]*0.0, xx[i,2]*1.0+xx[i+1,2]*0.0)
		t=(xx[i,0]*1.0+xx[i+1,0]*0.0, xx[i,1]*1.0+xx[i+1,1]*0.0, xx[i,2]*1.0+xx[i+1,2]*0.0)
		c=xx[i,1]*1.0+xx[i+1,1]*0.0+(xx[i,2]*1.0+xx[i+1,2]*0.0)*1.0j
            elif abs(xx[i+1,0])<1.0E-6: # vortices
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*0.0+xx[i+1,0]*1.0, xx[i,1]*0.0+xx[i+1,1]*1.0, xx[i,2]*0.0+xx[i+1,2]*1.0)
  		t=(xx[i,0]*0.0+xx[i+1,0]*1.0, xx[i,1]*0.0+xx[i+1,1]*1.0, xx[i,2]*0.0+xx[i+1,2]*1.0)
		c=xx[i,1]*0.0+xx[i+1,1]*1.0 + (xx[i,2]*0.0+xx[i+1,2]*1.0)*1.0j
            elif abs(xx[i,0]*0.5+xx[i+1,0]*0.5)<1.0E-6: # centers 
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*0.5+xx[i+1,0]*0.5, xx[i,1]*0.5+xx[i+1,1]*0.5, xx[i,2]*0.5+xx[i+1,2]*0.5)
		t=(xx[i,0]*0.5+xx[i+1,0]*0.5, xx[i,1]*0.5+xx[i+1,1]*0.5, xx[i,2]*0.5+xx[i+1,2]*0.5)
		c=(xx[i,1]*0.5+xx[i+1,1]*0.5)+(xx[i,2]*0.5+xx[i+1,2]*0.5)*1.0j
	    else:
	        continue

            if s not in lens:
                pts.append(t)
                ang=np.angle(c)
                phis.append(ang)
                lens.append(s)

    phis_= np.sort(phis)
    kykz_plane=[]
    for i in phis_:
        idx=phis.index(i) 
        kykz_plane.append(pts[idx])	
    kykz_plane.append(kykz_plane[0])	

#------- set kzkx plane (ky=0) by finding all the vortices 
    pts=[]
    lens=[]
    phis=[]

    for xx in e:
        for i in range(len(xx)-1):  
            if abs(xx[i,1])<1.0E-6: # vortices
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*1.0+xx[i+1,0]*0.0, xx[i,1]*1.0+xx[i+1,1]*0.0, xx[i,2]*1.0+xx[i+1,2]*0.0)
		t=(xx[i,0]*1.0+xx[i+1,0]*0.0, xx[i,1]*1.0+xx[i+1,1]*0.0, xx[i,2]*1.0+xx[i+1,2]*0.0)
		c=xx[i,2]*1.0+xx[i+1,2]*0.0+(xx[i,0]*1.0+xx[i+1,0]*0.0)*1.0j
            elif abs(xx[i+1,1])<1.0E-6: # vortices
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*0.0+xx[i+1,0]*1.0, xx[i,1]*0.0+xx[i+1,1]*1.0, xx[i,2]*0.0+xx[i+1,2]*1.0)
  		t=(xx[i,0]*0.0+xx[i+1,0]*1.0, xx[i,1]*0.0+xx[i+1,1]*1.0, xx[i,2]*0.0+xx[i+1,2]*1.0)
		c=xx[i,2]*0.0+xx[i+1,2]*1.0 + (xx[i,0]*0.0+xx[i+1,0]*1.0)*1.0j
            elif abs(xx[i,1]*0.5+xx[i+1,1]*0.5)<1.0E-6: # centers 
                s="{:16.8f}{:16.8f}{:16.8f}".format(xx[i,0]*0.5+xx[i+1,0]*0.5, xx[i,1]*0.5+xx[i+1,1]*0.5, xx[i,2]*0.5+xx[i+1,2]*0.5)
		t=(xx[i,0]*0.5+xx[i+1,0]*0.5, xx[i,1]*0.5+xx[i+1,1]*0.5, xx[i,2]*0.5+xx[i+1,2]*0.5)
		c=(xx[i,2]*0.5+xx[i+1,2]*0.5)+(xx[i,0]*0.5+xx[i+1,0]*0.5)*1.0j
	    else:
	        continue

            if s not in lens:
                pts.append(t)
                ang=np.angle(c)
                phis.append(ang)
                lens.append(s)

    phis_= np.sort(phis)
    kzkx_plane=[]
    for i in phis_:
        idx=phis.index(i) 
        kzkx_plane.append(pts[idx])	
    kzkx_plane.append(kzkx_plane[0])	

    face1 = mp3d.art3d.Poly3DCollection([kxky_plane], alpha=0.35, linewidth=1)
    face2 = mp3d.art3d.Poly3DCollection([kykz_plane], alpha=0.35, linewidth=1)
    face3 = mp3d.art3d.Poly3DCollection([kzkx_plane], alpha=0.35, linewidth=1)
    face1.set_facecolor('tab:gray') 
    face2.set_facecolor('tab:gray') 
    face3.set_facecolor('tab:gray') 
    ax.add_collection3d(face1)
    ax.add_collection3d(face2)
    ax.add_collection3d(face3)

# read in Nodes Points and plot it
    fl=sys.argv[1] # the Nodes.dat file name
    nl = sum(1 for line in open(fl))
    print "number of Nodes are ", nl-2 

    num_node= nl-2
    xp=np.zeros((num_node),dtype=np.float64)
    yp=np.zeros((num_node),dtype=np.float64)
    zp=np.zeros((num_node),dtype=np.float64)
    with open(fl,"r") as rd:
       rd.readline()
       rd.readline()
       for i in range(num_node):
           kx,ky,kz,gp,en,k1,k2,k3=rd.readline().strip().split()
           xp[i]=np.float64(kx)
           yp[i]=np.float64(ky)
           zp[i]=np.float64(kz)

# Modify here if you want plot different Node points with different color
    for i in range(num_node):
       if (abs(zp[i])) < 0.02: 
          ax.scatter3D(xp[i], yp[i], zp[i], s=1, color = "red", alpha=1.0) # s is for size, alpha for transparency
       elif (abs(xp[i])) < 0.02: 
          ax.scatter3D(xp[i], yp[i], zp[i], s=1, color = "blue", alpha=1.0) # s is for size, alpha for transparency
       elif (abs(yp[i])) < 0.02: 
          ax.scatter3D(xp[i], yp[i], zp[i], s=1, color = "green", alpha=1.0) # s is for size, alpha for transparency

    #ax.scatter3D(xp[0:num_node-40], yp[0:num_node-40], zp[0:num_node-40], s=1, color = "black", alpha=1.0) # s is for size, alpha for transparency
    #ax.scatter3D(xp[num_node-40:num_node-20], yp[num_node-40:num_node-20], zp[num_node-40:num_node-20], s=1, color = "red", alpha=1.0) # s is for size, alpha for transparency
    #ax.scatter3D(xp[num_node-20:num_node], yp[num_node-20:num_node], zp[num_node-20:num_node], s=1, color = "blue", alpha=1.0) # s is for size, alpha for transparency

    ax.set_xlim(-b1/2, b1/2)
    ax.set_ylim(-b2/2, b2/2)
    ax.set_zlim(-b3/2, b3/2)

    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.zaxis.set_ticklabels([])
    plt.axis('off')

# set equal xyz axis
    ax.set_aspect('equal')
    set_axes_equal(ax)

    plt.show()


