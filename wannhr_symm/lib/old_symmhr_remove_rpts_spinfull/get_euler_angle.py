#!/usr/bin/env python
import sys
import numpy as np
from numpy import cos,sin

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
        print zz, beta/np.pi

# gamma can be zero in this case since when beta=0 or pi, dmat only rely on alpha+gamma 
# or alpha-gamma, respectively
# actuallly we are calculating (alpha + gamm) when zz=1 or alpha-gamma
        gamma = 0.0 # /* in this case, alpha and gamma are the same, rotating around  the same z-axis.*/

        # cos(alpha)cos(beta)    -sin(alpha) cos(alpha)*cos(beta)
        # sin(alpha)cos(beta)     cos(alpha)        0
        #          0                   0        cos(beta)

        tmp = xy/zz
        if np.abs(np.abs(tmp)-1.0)<1.0E-6: tmp = tmp/np.abs(tmp)
        alpha = np.arcsin(tmp)

# move to -pi to pi
        if (-yx) >= 0.0 and yy >= 0.0: 
            alpha = alpha  
        elif (-yx) >= 0.0 and yy < 0.0: 
            alpha = np.pi - alpha  
        elif (-yx) < 0.0 and yy >= 0.0: 
            alpha = alpha
        elif (-yx) < 0.0 and yy < 0.0: 
            alpha = -np.pi - alpha 

        if np.abs(zx/zz)<1.0E-6:  # alpha is actually pi/2 or -pi/2
            if zz>=0.0: # alpha+gamma is 
               if (alpha - np.pi/2)< - np.pi:
                   gamma = alpha + np.pi/2
                   alpha = np.pi/2
               elif (alpha + np.pi/2)> np.pi:
                   gamma = alpha - np.pi/2
                   alpha = np.pi/2
               else: 
                   gamma = alpha - np.pi/2
                   alpha = np.pi/2

            else: # alpha-gamma is 
               if (-np.pi/2-alpha)< - np.pi:
                   gamma = np.pi/2 - alpha
                   alpha = np.pi/2  
               elif (np.pi/2-alpha)> np.pi:
                   gamma =-np.pi/2 - alpha
                   alpha =-np.pi/2  
               else: 
                   gamma =-np.pi/2 - alpha
                   alpha =-np.pi/2  

    else:  # general case
        beta = np.arccos(zz) # beta is always in [0, PI] 

        # determine gamma first since sin(beta)/=0
        csg = -xz/np.sin(beta) 
        sng =  yz/np.sin(beta)

        if np.abs(np.abs(sng)-1.0)<1.0E-6: sng = sng/np.abs(sng) # incase of sna > 1.0
        gamma = np.arcsin(sng)
 
# move to -pi to pi
        if sng >= 0.0 and csg >= 0.0: 
            gamma = gamma  
        elif sng >= 0.0 and csg < 0.0: 
            gamma = np.pi - gamma  
        elif sng < 0.0 and csg >= 0.0: 
            gamma = gamma
        elif sng < 0.0 and csg < 0.0: 
            gamma = -np.pi - gamma 

# determine alpha finally
        sna = zy/np.sin(beta) # sin(alpha)

        if np.abs(np.abs(sna)-1.0)<1.0E-6: sna = sna/np.abs(sna) # incase of sna > 1.0
        alpha= np.arcsin(sna) # zy = sin(alpha)sin(beta) asin() gives vale between [-PI/2, PI/2] 
 
# two case: zz=0 and zz/=0
        if np.abs(zz)>1.0E-6:
            csa = zx/zz           # cos(alpha)
   
        else: # zz=0
            # two case:  sin(r)=0 and /=0
            if abs(yz)>1.0E-6: 
                csa = xy / np.sin(gamma)#/=0
            else: 
                csa = yy / np.cos(gamma)#/=0

# move to -pi to pi
        if sna >= 0.0 and csa >= 0.0: 
            alpha = alpha  
        elif sna >= 0.0 and csa < 0.0: 
            alpha = np.pi - alpha  
        elif sna < 0.0 and csa >= 0.0: 
            alpha = alpha
        elif sna < 0.0 and csa < 0.0: 
            alpha = -np.pi - alpha 

# check whether this alpha beta gamma gives the same rmat
    rmtn = np.zeros((3,3),dtype = np.float64)
    rmtn[0,0] = cos(alpha)*cos(beta)*cos(gamma)-sin(alpha)*sin(gamma)
    rmtn[1,0] = sin(alpha)*cos(beta)*cos(gamma)+cos(alpha)*sin(gamma)
    rmtn[2,0] =-cos(gamma)*sin(beta)

    rmtn[0,1] =-sin(gamma)*cos(alpha)*cos(beta)-sin(alpha)*cos(gamma)
    rmtn[1,1] =-sin(gamma)*sin(alpha)*cos(beta)+cos(alpha)*cos(gamma)
    rmtn[2,1] = sin(beta)*sin(gamma)
    
    rmtn[0,2] = cos(alpha)*cos(beta)
    rmtn[1,2] = sin(alpha)*sin(beta)
    rmtn[2,2] = cos(beta)
 
    for i in range(3):
        for j in range(3):
            if abs(rmtn[i,j])<1.0E-9: rmtn[i,j] = 0.0
    print "rot_glb_new\n", rmtn
    print "rot_glb minus rot_glb_new\n", rmat-rmtn 
    sm = sum(sum(abs(rmat-rmtn)))
    if sm < 1.0E-9: sm = 0.0
    print "sum abs minus", sm
     
    return alpha, beta, gamma

rmat = np.zeros((3,3),dtype = np.float64)
rmat[0,0] = -2.0/3
rmat[1,0] =  1.0/3
rmat[2,0] = -2.0/3

rmat[0,1] =  1.0/3
rmat[1,1] = -2.0/3
rmat[2,1] = -2.0/3

rmat[0,2] = -2.0/3
rmat[1,2] = -2.0/3
rmat[2,2] =  1.0/3


pi = np.pi
a,b,c=rmat2euler(rmat)
print a/pi,b/pi,c/pi





#0 local axis
#[[-0.6666666667  0.3333333333 -0.6666666667]
# [ 0.3333333333 -0.6666666667 -0.6666666667]
# [-0.6666666667 -0.6666666667  0.3333333333]]
#1 local axis
#[[ 0.6666666667 -0.3333333333  0.6666666667]
# [-0.3333333333  0.6666666667  0.6666666667]
# [-0.6666666667 -0.6666666667  0.3333333333]]
#2 local axis
#[[ 0.6666666667  0.6666666667 -0.3333333333]
# [-0.6666666667  0.3333333333 -0.6666666667]
# [-0.3333333333  0.6666666667  0.6666666667]]
#3 local axis
#[[-0.6666666667 -0.6666666667  0.3333333333]
# [ 0.6666666667 -0.3333333333  0.6666666667]
# [-0.3333333333  0.6666666667  0.6666666667]]
#
