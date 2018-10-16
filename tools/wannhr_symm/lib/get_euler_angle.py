#!/usr/bin/env python
import sys
import numpy as np
import numpy.linalg as LA

def get_angle(sn,co):
    '''
    get angle in 0-2pi given sin(angle) and cos(angle)
    '''
    if np.abs(np.abs(sn)-1.0)<1.0E-3: sn = sn/np.abs(sn)
    angle = np.arcsin(sn)

    # move to 0-2pi
    if   sn >= 0.0 and co >= 0.0: 
        angle = angle  
    elif sn >= 0.0 and co  < 0.0: 
        angle = np.pi - angle  
    elif sn  < 0.0 and co >= 0.0: 
        angle = angle + 2*np.pi
    elif sn  < 0.0 and co  < 0.0: 
        angle = np.pi - angle 
    else:
        print "Error in get euler angle. Please contact the author! "
        sys.exit(0)
    if np.abs(np.abs(sn))<1.0E-3: 
       if co>0.0: angle = 0.0
       if co<0.0: angle = np.pi
    return angle

def rmat2euler(rmat):
    """
    Given normalized rotation matrix, return the Euler angle
    rmat should be right hand 
    For point goups, elements of rmat are sort of special. This function 
    does not apply to general cases. Please force precision before using it.
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

    if np.abs(np.abs(zz)-1.0) < 1.0E-3: # new z-axis is along old z-axis
        if zz > 0: # along the positive direction
            beta = 0.0    
        else:             # has a 180-degree rotation
            beta = np.pi 
        gamma = 0.0 # /* in this case, alpha and gamma can be combined since rotating around  the same z-axis.*/

        # cos(alpha)cos(beta)    -sin(alpha) cos(alpha)
        # sin(alpha)cos(beta)     cos(alpha)     0
        #          0                   0     cos(beta)
        sna = xy/zz
	coa = yy
	alpha = get_angle(sna,coa)

    else:  # general case
        beta = np.arccos(zz) # beta is always in (0, PI)

        # determine gamma first since sin(beta)/=0
        sng =  yz/np.sin(beta)
        cog = -xz/np.sin(beta) 
        gamma = get_angle(sng,cog)

# determine alpha finally
        sna = zy/np.sin(beta) # sin(alpha)

# two case: cos(beta) = zz != 0
        if np.abs(zz)>1.0E-3:
            csa = zx/zz 
   
        else: # cos(beta) = zz = 0, sin(beta) != 0
            if abs(yz)>1.0E-3: # sin(gamma)!=0
                csa = xy / np.sin(gamma)
            else: # sin(gamma)=0, cos(gamma)!=0
                csa = yy / np.cos(gamma)#/=0

        alpha = get_angle(sna,csa)

    return alpha, beta, gamma


def spin_reps(prep):
    """
    Copied from G Winker's symmetry_soc.py.
    Calculates the spin rotation matrices. The formulas to determine the rotation axes and angles
    are taken from `here <http://scipp.ucsc.edu/~haber/ph116A/rotation_11.pdf>`_.

    :param prep:   List that contains 3d rotation matrices.
    :type prep:    list(array)
    """
    # general representation of the D1/2 rotation about the axis (l,m,n) around the
    # angle phi
    D12 = lambda l,m,n,phi: np.array([[np.cos(phi/2.) - 1j*n*np.sin(phi/2.), (-1j*l -m)*np.sin(phi/2.)],
                               [(-1j*l + m)*np.sin(phi/2.), np.cos(phi/2.) + 1j*n*np.sin(phi/2.)]])

    #print "prep\n", prep
    n = np.zeros(3)
    tr = np.trace(prep)
    det = np.round(np.linalg.det(prep),5)
    if  det == 1.: #rotations
        theta = np.arccos(0.5*(tr-1.))
        if theta != 0:
            n[0] = prep[2,1]-prep[1,2]
            n[1] = prep[0,2]-prep[2,0]
            n[2] = prep[1,0]-prep[0,1]               
            if np.round(np.linalg.norm(n),5) == 0.: # theta = pi, that is C2 rotations
                e,v = LA.eig(prep)
                n = v[:,list(np.round(e,10)).index(1.)]
                spin=np.round(D12(n[0],n[1],n[2],np.pi),15)
            else:
                n /= np.linalg.norm(n)
                spin=np.round(D12(n[0],n[1],n[2],theta),15)
        else: # case of unitiy
            spin=D12(0,0,0,0)
    elif det == -1.: #improper rotations and reflections
        theta = np.arccos(0.5*(tr+1.)) 
        if np.round(theta,5) != np.round(np.pi,5):                 
            n[0] = prep[2,1]-prep[1,2]
            n[1] = prep[0,2]-prep[2,0]
            n[2] = prep[1,0]-prep[0,1]                
            if np.round(np.linalg.norm(n),5)== 0.: # theta = 0 (reflection)
                e,v = LA.eig(prep)
                n = v[:,list(np.round(e,10)).index(-1.)] # normal vector is eigenvector to eigenvalue -1
                spin=np.round(D12(n[0],n[1],n[2],np.pi),15) #spin is a pseudovector!
            else:
                n /= np.linalg.norm(n)
                # rotation followed by reflection:
                spin=np.round(np.dot(D12(n[0],n[1],n[2],np.pi),D12(n[0],n[1],n[2],theta)),15)
        else: # case of inversion (does not do anything to spin)
            spin=D12(0,0,0,0)
    return np.array(spin)


