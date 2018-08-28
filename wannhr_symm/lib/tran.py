#!/usr/bin/env python

import numpy as np

def tran_op(op, tmat):

    """
    transform quantum operator from representation A to 
    another representation B
    
    Args:
        op:   the matrix form of operator in representation A
        tmat: the unitary transform matrix
    """

    return np.dot(np.dot(np.conj(np.transpose(tmat)), op), tmat)

def tmat_c2r(case, ispin=False):

    """
    the transform matrix from complex shperical harmonics to 
    real spherical harmonics
    
    Args:
        case: label for different systems 
        ispin: whether to include spin or not 
    """
    sqrt2 = np.sqrt(2.0)
    ci    = np.complex128(0.0+1.0j)
    cone  = np.complex128(1.0+0.0j)

    if case.strip() == 's':
        nband = 1
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        t_c2r[0,0] =  cone

    elif case.strip() == 'p':
        nband = 3
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        # px=1/sqrt(2)( |1,-1> - |1,1> )
        t_c2r[0,0] =  cone/sqrt2
        t_c2r[2,0] = -cone/sqrt2
        # py=i/sqrt(2)( |1,-1> + |1,1> )
        t_c2r[0,1] =  ci/sqrt2
        t_c2r[2,1] =  ci/sqrt2
        # pz=|1,0>
        t_c2r[1,2] =  cone

    elif case.strip() == 'pwien': # in wien by default px py pz
        nband = 3
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        # px=1/sqrt(2)( |1,-1> - |1,1> )
        t_c2r[0,0] =  cone/sqrt2
        t_c2r[2,0] = -cone/sqrt2
        # py=i/sqrt(2)( |1,-1> + |1,1> )
        t_c2r[0,1] =  ci/sqrt2
        t_c2r[2,1] =  ci/sqrt2
        # pz=|1,0>
        t_c2r[1,2] =  cone

    elif case.strip() == 'pwann': # p in wannier basis vasp
        nband = 3
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        t_r2w = np.zeros((nband, nband), dtype=np.complex128)
        t_c2w = np.zeros((nband, nband), dtype=np.complex128)
        # px=1/sqrt(2)( |1,-1> - |1,1> )
        t_c2r[0,0] =  cone/sqrt2
        t_c2r[2,0] = -cone/sqrt2
        # py=i/sqrt(2)( |1,-1> + |1,1> )
        t_c2r[0,1] =  ci/sqrt2
        t_c2r[2,1] =  ci/sqrt2
        # pz=|1,0>
        t_c2r[1,2] =  cone
    
        # pz = (px,py,pz) (0,0,1)^T 
        t_r2w[2,0] = cone
        # px = (px,py,pz) (1,0,0)^T 
        t_r2w[0,1] = cone
        # py = (px,py,pz) (0,1,0)^T 
        t_r2w[1,2] = cone
        t_c2w = np.dot(t_c2r,t_r2w)
        t_c2r = t_c2w

    elif case.strip() == 't2g':
        nband = 3
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        # dzx --> py=i/sqrt(2)( |1,-1> + |1,1> )
        t_c2r[0,0] =  ci/sqrt2
        t_c2r[2,0] =  ci/sqrt2
        # dzy --> px=1/sqrt(2)( |1,-1> - |1,1> )
        t_c2r[0,1] =  cone/sqrt2
        t_c2r[2,1] = -cone/sqrt2
        # dxy --> pz=|1,0>
        t_c2r[1,2] =  cone
    elif case.strip() == 'd':
        nband = 5
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        # dz2=|2,0>
        t_c2r[2,0] =  cone
        # dzx=1/sqrt(2)( |2,-1> - |2,1> )
        t_c2r[1,1] =  cone/sqrt2 
        t_c2r[3,1] = -cone/sqrt2
        # dzy=i/sqrt(2)( |2,-1> + |2,1> )
        t_c2r[1,2] =  ci/sqrt2
        t_c2r[3,2] =  ci/sqrt2
        # dx2-y2=1/sqrt(2)( |2,-2> + |2,2> )
        t_c2r[0,3] =  cone/sqrt2
        t_c2r[4,3] =  cone/sqrt2
        # dxy=i/sqrt(2)( |2,-2> - |2,2> )
        t_c2r[0,4] =  ci/sqrt2
        t_c2r[4,4] = -ci/sqrt2

    elif case.strip() == 'dwien': # by default wien2k: dxy dzx dyz dx2-y2 dz2
        nband = 5
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        # dz2=|2,4>
        t_c2r[2,4] =  cone
        # dzx=1/sqrt(2)( |2,-1> - |2,1> )
        t_c2r[1,1] =  cone/sqrt2 
        t_c2r[3,1] = -cone/sqrt2
        # dzy=i/sqrt(2)( |2,-1> + |2,1> )
        t_c2r[1,2] =  ci/sqrt2
        t_c2r[3,2] =  ci/sqrt2
        # dx2-y2=1/sqrt(2)( |2,-2> + |2,2> )
        t_c2r[0,3] =  cone/sqrt2
        t_c2r[4,3] =  cone/sqrt2
        # dxy=i/sqrt(2)( |2,-2> - |2,2> )
        t_c2r[0,0] =  ci/sqrt2
        t_c2r[4,0] = -ci/sqrt2

    elif case.strip() == 'f':
        nband = 7
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        # fz3 = |3,0>
        t_c2r[3, 0] =  cone 
        # fxz2 = 1/sqrt(2)( |3,-1> - |3,1> )
        t_c2r[2, 1] =  cone/sqrt2
        t_c2r[4, 1] = -cone/sqrt2
        # fyz2 = i/sqrt(2)( |3,-1> + |3,1> )
        t_c2r[2, 2] =  ci/sqrt2
        t_c2r[4, 2] =  ci/sqrt2
        # fz(x2-y2) = 1/sqrt(2)( |3,-2> + |3,2> )
        t_c2r[1, 3] =  cone/sqrt2
        t_c2r[5, 3] =  cone/sqrt2
        # fxyz = i/sqrt(2)( |3,-2> - |3,2> )
        t_c2r[1, 4] =  ci/sqrt2
        t_c2r[5, 4] = -ci/sqrt2
        # fx(x2-3y2) = 1/sqrt(2) ( |3,-3> - |3,3> )
        t_c2r[0, 5] =  cone/sqrt2
        t_c2r[6, 5] = -cone/sqrt2
        # fy(3x2-y2) = i/sqrt(2) ( |3,-3> + |3,3> )
        t_c2r[0, 6] =  ci/sqrt2
        t_c2r[6, 6] =  ci/sqrt2

    elif case.strip() == 'fwien': # fxz2 fyz2 fz3 fx(x2-3y2) fy(3x2-y2) fz(x2-y2) fxyz
        nband = 7
        t_c2r = np.zeros((nband, nband), dtype=np.complex128)
        # fz3 = |3,0>
        t_c2r[3, 2] =  cone 
        # fxz2 = 1/sqrt(2)( |3,-1> - |3,1> )
        t_c2r[2, 0] =  cone/sqrt2
        t_c2r[4, 0] = -cone/sqrt2
        # fyz2 = i/sqrt(2)( |3,-1> + |3,1> )
        t_c2r[2, 1] =  ci/sqrt2
        t_c2r[4, 1] =  ci/sqrt2
        # fz(x2-y2) = 1/sqrt(2)( |3,-2> + |3,2> )
        t_c2r[1, 5] =  cone/sqrt2
        t_c2r[5, 5] =  cone/sqrt2
        # fxyz = i/sqrt(2)( |3,-2> - |3,2> )
        t_c2r[1, 6] =  ci/sqrt2
        t_c2r[5, 6] = -ci/sqrt2
        # fx(x2-3y2) = 1/sqrt(2) ( |3,-3> - |3,3> )
        t_c2r[0, 3] =  cone/sqrt2
        t_c2r[6, 3] = -cone/sqrt2
        # fy(3x2-y2) = i/sqrt(2) ( |3,-3> + |3,3> )
        t_c2r[0, 4] =  ci/sqrt2
        t_c2r[6, 4] =  ci/sqrt2

    else:
        print "don't support t_c2r for this case: ", case
        return

    if ispin:
        norbs=2*nband
        t_c2r_spin = np.zeros((norbs,norbs), dtype=np.complex128)
        t_c2r_spin[0:norbs:2,0:norbs:2] = t_c2r      
        t_c2r_spin[1:norbs:2,1:norbs:2] = t_c2r      
        return t_c2r_spin
    else:
        return t_c2r

def tmat_r2c(case, ispin=False):
    
    """
    the transform matrix from real spherical harmonics to
    complex shperical harmonics
    
    Args:
        case: label for different systems
        ispin: whether to include spin or not
    """

    return np.conj(np.transpose(tmat_c2r(case, ispin)))

def tmat_r2cub(ispin=False):
    
    """
    the transform matrix from real spherical harmonics to the cubic
    spherical harmonics, just for f system

    Args:
        ispin: whether to include spin or not
    """

    a = np.sqrt(10.0) / 4.0 + 0.0j
    b = np.sqrt(6.0) / 4.0  + 0.0j
    c = 1.0 + 0.0j

    nband = 7
    t_r2cub = np.zeros((nband,nband), dtype=np.complex128)
    # fx3 = -sqrt(6)/4 fxz2 + sqrt(10)/4 fx(x2-3y2) 
    t_r2cub[1, 0] = -b
    t_r2cub[5, 0] =  a
    # fy3 = -sqrt(6)/4 fyz2 - sqrt(10)/4 fy(3x2-y2) 
    t_r2cub[2, 1] = -b
    t_r2cub[6, 1] = -a
    # fz3 = fz3
    t_r2cub[0, 2] =  c
    # fx(y2-z2) = -sqrt(10)/4 fxz2 - sqrt(6)/4 fx(x2-3y2)
    t_r2cub[1, 3] = -a
    t_r2cub[5, 3] = -b
    # fy(z2-x2) = sqrt(10)/4 fyz2 - sqrt(6)/4 fy(3x2-y2)
    t_r2cub[2, 4] =  a
    t_r2cub[6, 4] = -b
    # fz(x2-y2) = fz(x2-y2)
    t_r2cub[3, 5] =  c
    # fxyz = fxyz
    t_r2cub[4, 6] =  c

    if ispin:
        norbs = 2 * nband
        t_r2cub_spin = np.zeros((norbs, norbs), dtype=np.complex128)
        t_r2cub_spin[0:norbs:2,0:norbs:2] = t_r2cub
        t_r2cub_spin[1:norbs:2,1:norbs:2] = t_r2cub
        return t_r2cub_spin
    else:
        return t_r2cub

def tmat_cub2r(ispin=False):

    """
    the transform matrix from cubic spherical harmonics to
    real spherical harmonics, just for f system

    Args:
        ispin: whether to include spin or not
    """

    return np.conj( np.transpose( tmat_r2cub(ispin) ) )

def tmat_c2j(l):
    """
    the transform matrix from complex shperical harmonics to 
    j2,jz basis
    
    Args:
        case: label for different systems 
    """

    if l == 1:
        t_c2j = np.zeros((6, 6), dtype=np.complex128)
        t_c2j[0,0] = -np.sqrt(2.0/3.0) 
        t_c2j[3,0] =  np.sqrt(1.0/3.0) 
        t_c2j[2,1] = -np.sqrt(1.0/3.0) 
        t_c2j[5,1] =  np.sqrt(2.0/3.0) 
        t_c2j[1,2] =  1.0
        t_c2j[0,3] =  np.sqrt(1.0/3.0)
        t_c2j[3,3] =  np.sqrt(2.0/3.0)
        t_c2j[2,4] =  np.sqrt(2.0/3.0)
        t_c2j[5,4] =  np.sqrt(1.0/3.0)
        t_c2j[4,5] =  1.0
        return t_c2j
    elif l == 2:
        t_c2j = np.zeros((10, 10), dtype=np.complex128)
        t_c2j[0,0] =  -np.sqrt(4.0/5.0)
        t_c2j[3,0] =   np.sqrt(1.0/5.0)
        t_c2j[2,1] =  -np.sqrt(3.0/5.0)
        t_c2j[5,1] =   np.sqrt(2.0/5.0)
        t_c2j[4,2] =  -np.sqrt(2.0/5.0)
        t_c2j[7,2] =   np.sqrt(3.0/5.0)
        t_c2j[6,3] =  -np.sqrt(1.0/5.0)
        t_c2j[9,3] =   np.sqrt(4.0/5.0)
        t_c2j[1,4] =   1.0
        t_c2j[0,5] =   np.sqrt(1.0/5.0)
        t_c2j[3,5] =   np.sqrt(4.0/5.0)
        t_c2j[2,6] =   np.sqrt(2.0/5.0)
        t_c2j[5,6] =   np.sqrt(3.0/5.0)
        t_c2j[4,7] =   np.sqrt(3.0/5.0)
        t_c2j[7,7] =   np.sqrt(2.0/5.0)
        t_c2j[6,8] =   np.sqrt(4.0/5.0)
        t_c2j[9,8] =   np.sqrt(1.0/5.0)
        t_c2j[8,9] =   1.0
        return t_c2j
    elif l == 3:
        t_c2j = np.zeros((14,14), dtype=np.complex128)
        t_c2j[0,0] = -np.sqrt(6.0/7.0)
        t_c2j[3,0] =  np.sqrt(1.0/7.0)
        t_c2j[2,1] = -np.sqrt(5.0/7.0)
        t_c2j[5,1] =  np.sqrt(2.0/7.0)
        t_c2j[4,2] = -np.sqrt(4.0/7.0)
        t_c2j[7,2] =  np.sqrt(3.0/7.0)
        t_c2j[6,3] = -np.sqrt(3.0/7.0)
        t_c2j[9,3] =  np.sqrt(4.0/7.0)
        t_c2j[8,4] = -np.sqrt(2.0/7.0)
        t_c2j[11,4] =  np.sqrt(5.0/7.0)
        t_c2j[10,5] = -np.sqrt(1.0/7.0)
        t_c2j[13,5] =  np.sqrt(6.0/7.0)
        t_c2j[1,6] =  1.0
        t_c2j[0,7] =  np.sqrt(1.0/7.0)
        t_c2j[3,7] =  np.sqrt(6.0/7.0)
        t_c2j[2,8] =  np.sqrt(2.0/7.0)
        t_c2j[5,8] =  np.sqrt(5.0/7.0)
        t_c2j[4,9] =  np.sqrt(3.0/7.0)
        t_c2j[7,9] =  np.sqrt(4.0/7.0)
        t_c2j[6,10] =  np.sqrt(4.0/7.0)
        t_c2j[9,10] =  np.sqrt(3.0/7.0)
        t_c2j[8,11] =  np.sqrt(5.0/7.0)
        t_c2j[11,11] =  np.sqrt(2.0/7.0)
        t_c2j[10,12] =  np.sqrt(6.0/7.0)
        t_c2j[13,12] =  np.sqrt(1.0/7.0)
        t_c2j[12,13] =  1.0
        return t_c2j
    else:
        print "NOT Implemented !!!"


def fourier_hr2hk(norbs, nkpt, kvec, nrpt, rvec, deg_rpt, hr):

    """
    Fourier transform from R-space to K-space

    Args:
        norbs:   number of orbitals
        nkpt:    number of K-points
        kvec:    fractional coordinate for K-points
        nrpt:    number of R-points
        rvec:    fractional coordinate for R-points
        deg_rpt: the degenerate for each R-point
        hr:      Hamiltonian in R-space

    Return:
        hk:      Hamiltonian in K-space 
    """

    hk = np.zeros((nkpt, norbs, norbs), dtype=np.complex128)
    for i in range(nkpt):
        #print "kvec", i, kvec[i,:]
        for j in range(nrpt):
            coef  = 2*np.pi*np.dot(kvec[i,:], np.float64(rvec[j,:]))
            ratio = (np.cos(coef) + np.sin(coef) * 1j) / np.float64(deg_rpt[j])
            hk[i,:,:] = hk[i,:,:] + ratio * hr[j,:,:]
    return hk

def fourier_hr2h1k(norbs, kvec, nrpt, rvec, deg_rpt, hr):

    """
    Fourier transform from R-space to K-space

    Args:
        norbs:   number of orbitals
        nkpt:    number of K-points
        kvec:    fractional coordinate for K-points
        nrpt:    number of R-points
        rvec:    fractional coordinate for R-points
        deg_rpt: the degenerate for each R-point
        hr:      Hamiltonian in R-space

    Return:
        hk:      Hamiltonian in K-space 
    """

    hk = np.zeros((norbs, norbs), dtype=np.complex128)
    for i in range(nrpt):
        coef =  2*np.pi*np.dot(kvec, np.float64(rvec[i,:]))
        ratio = (np.cos(coef) + np.sin(coef) * 1.0j) / np.float64(deg_rpt[i])
        hk[:,:] = hk[:,:] + ratio * hr[i,:,:]
    return hk

def myfourier_hr2hk(norbs, nkpt, kvec, nrpt, rvec, deg_rpt, hr):

    """
    Fourier transform from R-space to K-space

    Args:
        norbs:   number of orbitals
        nkpt:    number of K-points
        kvec:    fractional coordinate for K-points
        nrpt:    number of R-points
        rvec:    fractional coordinate for R-points
        deg_rpt: the degenerate for each R-point
        hr:      Hamiltonian in R-space

    Return:
        hk:      Hamiltonian in K-space 
    """
    print "ALERT: Gauge changed: not exp(-ikR) but exp(ikR)*exp(ik*tau_mu)"

    hk = np.zeros((nkpt, norbs, norbs), dtype=np.complex128)
    for i in range(nkpt):
        for j in range(nrpt):
            coef = 2*np.pi*np.dot(kvec[i,:], rvec[j,:])
            ratio = (np.cos(coef) + np.sin(coef) * 1j) / float(deg_rpt[j])
            hk[i,:,:] = hk[i,:,:] + ratio * hr[j,:,:]
    return hk

def fourier_hr2hk_gauge(norbs, nkpt, kvec, nrpt, rvec, deg_rpt, hr):

    """
    Fourier transform from R-space to K-space

    Args:
        norbs:   number of orbitals
        nkpt:    number of K-points
        kvec:    fractional coordinate for K-points
        nrpt:    number of R-points
        rvec:    fractional coordinate for R-points
        deg_rpt: the degenerate for each R-point
        hr:      Hamiltonian in R-space

    Return:
        hk:      Hamiltonian in K-space 
    """

    hk = np.zeros((nkpt, norbs, norbs), dtype=np.complex128)
    for i in range(nkpt):
        print "kvec", i, kvec[i,:]
        for j in range(nrpt):
            coef = 2*np.pi*np.dot(kvec[i,:], rvec[j,:])
            ratio = (np.cos(coef) + np.sin(coef) * 1j) / float(deg_rpt[j])
            hk[i,:,:] = hk[i,:,:] + ratio * hr[j,:,:]
    return hk


