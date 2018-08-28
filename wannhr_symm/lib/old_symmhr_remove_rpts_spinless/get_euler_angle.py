#!/usr/bin/env python
import sys
import numpy as np

def rmat2euler(rmat):
    """
    Given normalized rotation matrix, return the Euler angle
    rmat should be right hand 
    Parameters
    ----------
    rmat:   3*3 matrix

    Returns
    ---------
    alpha:  Euler angle, in radian, [0,2*pi]
    beta:   Euler angle, in radian, [0,  pi]
    gamma:  Euler angle, in radian, [0,2*pi]
   
    """

# check orthogonality
    if abs(np.dot(np.transpose(rmat),rmat)-np.eye(3,dtype=np.float64)).sum() > 1.0E-6:
        print "rmat is not orthogonal, please check it!"
        sys.exit(0)

# now rmat is orthogonal. replace y axis with z-axis cross product x-axis to 
# make a right-hand coordinate. */
    rrmt = rmat * 1.0
    rrmt[:,1] = np.cross(rrmt[:,2],rrmt[:,0])
    xx = rrmt[0,0]; yx = rrmt[0,1]; zx = rrmt[0,2];
    xy = rrmt[1,0]; yy = rrmt[1,1]; zy = rrmt[1,2];
    xz = rrmt[2,0]; yz = rrmt[2,1]; zz = rrmt[2,2];

    if np.abs(np.abs(zz)-1.0) < 1.0E-6: # new z-axis is along old z-axis
        if zz > 0: # along the positive direction
            beta = 0.0    
        else:             # has a 180-degree rotation
            beta = np.pi 

        gamma = 0.0 # /* in this case, alpha and gamma are the same, rotating around  the same z-axis.*/

        # cos(alpha)cos(beta)    -sin(alpha)     0
        # sin(alpha)cos(beta)     cos(alpha)     0
        #          0                   0     cos(beta)
        alpha = np.arcsin(xy/zz)
        if xx/zz < 0.0:
             alpha = np.pi - alpha
        if alpha<0.0:
             alpha = 2.0 * np.pi + alpha

    else:  # general case
        beta = np.arccos(zz) # beta is always in [0, PI] 
        alpha= np.arcsin(zy/np.sin(beta)) # zy = sin(alpha)sin(beta) asin() gives vale between [-PI/2, PI/2] 
        if  zx < 0.0: # zx =sin(beta)*cos(alpha) and sin(beta)>=0.0, this means cos(alpha)<0.0
            alpha = np.pi - alpha
        if  alpha < 0.0:
            alpha = 2.0 * np.pi + alpha

        # determine gamma now
        if np.abs(np.abs(-xz/np.sin(beta))-1.0) < 1.0E-5:
            if -xz/np.sin(beta) < 0.0:
                gamma = np.pi  
            else: 
                gamma = 0.0 
        else:
            gamma = np.arccos(-xz/np.sin(beta)) # xz=-cos(gamma)*sin(beta). acos() gieve a value between [0, PI]

        # /* we need sin(gamma) to finally determin gamma. xx=cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma) */
        tmp = -(xx-(-xz/np.sin(beta))*zz*np.cos(alpha))/np.sin(alpha) # /* tmp is sin(gamma) */ 
        if  tmp < 0.0:
            gamma = 2.0 * np.pi-gamma;
    return alpha, beta, gamma
