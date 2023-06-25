#!/usr/bin/env python

import numpy as np

def euler_to_rmat(alpha, beta, gamma):
    """
    Given Euler angle: alpha, beta, gamma, generate the rotation matrix R

    Parameters
    ----------
    alpha:  Euler angle, in radian, [0,2*pi]
    beta:   Euler angle, in radian, [0,  pi]
    gamma:  Euler angle, in radian, [0,2*pi]
   
    Returns
    ----------
    rmat:   The rotation matrix
    """

    rmat=np.zeros((3,3), dtype=np.float64)
    rmat[0,0] =  np.cos(alpha) * np.cos(beta) * np.cos(gamma) - np.sin(alpha) * np.sin(gamma)
    rmat[0,1] = -np.sin(gamma) * np.cos(alpha) * np.cos(beta) - np.sin(alpha) * np.cos(gamma)
    rmat[0,2] =  np.cos(alpha) * np.sin(beta)
    rmat[1,0] =  np.sin(alpha) * np.cos(beta) * np.cos(gamma) + np.cos(alpha) * np.sin(gamma)
    rmat[1,1] = -np.sin(gamma) * np.sin(alpha) * np.cos(beta) + np.cos(alpha) * np.cos(gamma)
    rmat[1,2] =  np.sin(alpha) * np.sin(beta)
    rmat[2,0] = -np.cos(gamma) * np.sin(beta)
    rmat[2,1] =  np.sin(beta) * np.sin(gamma)
    rmat[2,2] =  np.cos(beta)
    return rmat

def rmat_to_euler(rmat):
    """
    Given rotation matrix, return the Euler angle

    Parameters
    ----------
    rmat:   3*3 matrix

    Returns
    ---------
    alpha:  Euler angle, in radian, [0,2*pi]
    beta:   Euler angle, in radian, [0,  pi]
    gamma:  Euler angle, in radian, [0,2*pi]
   
    """

    if np.abs(rmat[2,2]) < 1.0:
        beta = np.arccos(rmat[2,2]) 

        cos_gamma = -rmat[2,0] / np.sin(beta)
        sin_gamma =  rmat[2,1] / np.sin(beta)
        gamma = where_is_angle(sin_gamma, cos_gamma)
        
        cos_alpha = rmat[0,2] / np.sin(beta)
        sin_alpha = rmat[1,2] / np.sin(beta)
        alpha = where_is_angle(sin_alpha, cos_alpha)
    else:
        if rmat[2,2] > 0: # cos(beta) = 1, beta = 0, sin(beta/2)=0.0
            beta = 0.0    
#       cos(alpha+gamma)     -sin(alpha+gamma)     0
#      -sin(alpha+gamma)      cos(alpha+gamma)     0
#                0                   0             1
            gamma = 0.0
            alpha = np.arccos(rmat[1,1])
            if   -rmat[0,1] < 0.0:
                 alpha = -1.0*alpha

        else:             # cos(beta) =-1, beta =pi, sin(beta/2)=1.0
            beta = np.pi 
#      -cos(alpha-gamma)     -sin(alpha-gamma)     0
#      -sin(alpha-gamma)      cos(alpha-gamma)     0
#                0                   0            -1
            gamma = 0.0
            alpha = np.arccos(rmat[1,1]) # 0~pi pi/2 if rmat[0,1] 
            if -rmat[0,1] < 0.0: alpha = -1.0*alpha
    return alpha, beta, gamma

def where_is_angle(sina, cosa):
    """
    Given sin and cos of an angle, return the angle range from [0,2*pi]

    Parameters
    ----------
    sina:   sin(alpha)
    cosa:   cos(alpha)

    Returns
    ----------
    alpha:  range from [0,2*pi]
    """
    if cosa > 1.0:
        cosa =1.0
    elif cosa < -1.0:
        cosa = -1.0
    alpha = np.arccos(cosa)
    if sina < 0.0:
        alpha = 2.0 * np.pi - alpha
    return alpha

def dmat_spinor(alpha,beta,gamma):
    """
    Given three Euler angle alpha, beta, gamma, return the transformation
    matrix for 1/2 spinor
   
    Parameters
    ----------
    alpha:   Euler angle [0,2*pi]
    beta :   Euler angle [0,  pi]
    gamma:   Euler angle [0,2*pi]
    
    Returns
    ----------
    dmat:    np.ndarray((2,2), dtype=np.complex128)
    """
    dmat = np.zeros((2,2), dtype=np.complex128)
    dmat[0,0] =  np.exp(-(alpha+gamma)/2.0 * 1j) * np.cos(beta/2.0) 
    dmat[0,1] = -np.exp(-(alpha-gamma)/2.0 * 1j) * np.sin(beta/2.0) 
    dmat[1,0] =  np.exp( (alpha-gamma)/2.0 * 1j) * np.sin(beta/2.0) 
    dmat[1,1] =  np.exp( (alpha+gamma)/2.0 * 1j) * np.cos(beta/2.0) 
    return dmat

def zx_to_rmat(z,x):
    """
    Given z vector and x vector, return y vector satisfy right-hand coordinate,
    and normalize them if needed

    Parameters
    ----------
    z:    np.ndarray((3,1), dtype=np.float64)
    x:    np.ndarray((3,1), dtype=np.float64)
    
    Returns
    ---------
    rmat:   [x,y,z] np.ndarray((3,3), dtype=np.float64)
    """
    z = np.array(z, dtype=np.float64) 
    x = np.array(x, dtype=np.float64) 
    xx = x / np.sqrt(np.dot(x,x))
    zz = z / np.sqrt(np.dot(z,z))
    yy = np.cross(zz,xx)
    print yy

    rmat = np.zeros((3,3), dtype=np.float64)
    rmat[:,0] = xx
    rmat[:,1] = yy
    rmat[:,2] = zz

    return rmat
