#!/usr/bin/env python
'''
get rotation matrix for a given rotation matrix in basis of s, p, d, f and t2g
orbitals.
'''

import sympy, sys
from numpy.linalg import inv, det
import numpy as np
from sympy import *

x = sympy.Symbol('x')
y = sympy.Symbol('y')
z = sympy.Symbol('z')

def ss(x,y,z):
    return 1+0*(x+y+z)

def pz(x,y,z):
    return z

def px(x,y,z):
    return x

def py(x,y,z):
    return y

def dz2(x,y,z):
    return (2*z*z-x*x-y*y)/(2*sympy.sqrt(3))

def dxz(x,y,z):
    return x*z

def dyz(x,y,z):
    return y*z

def dx2_y2(x,y,z):
    return (x*x-y*y)/2

def dxy(x,y,z):
    return x*y

def fz3(x,y,z):
    return z*(2*z*z-3*x*x-3*y*y)/(2*sympy.sqrt(15))

def fxz2(x,y,z):
    return x*(4*z*z-x*x-y*y)/(2*sympy.sqrt(10))

def fyz2(x,y,z):
    return y*(4*z*z-x*x-y*y)/(2*sympy.sqrt(10))

def fzx2_zy2(x,y,z):
    return z*(x*x-y*y)/2

def fxyz(x,y,z):
    return x*y*z

def fx3_3xy2(x,y,z):
    return x*(x*x-3*y*y)/(2*sympy.sqrt(6))

def f3yx2_y3(x,y,z):
    return y*(3*x*x-y*y)/(2*sympy.sqrt(6))

ft2g=[dxz,dyz,dxy]
fs=[ss]
fp=[pz,px,py]
fd=[dz2,dxz,dyz,dx2_y2,dxy]
ff=[fz3,fxz2,fyz2,fzx2_zy2,fxyz,fx3_3xy2,f3yx2_y3]
forb=[fs,fp,fd,ff,ft2g]
#      0  1  2  3  4

def get_orb_map_s(xp,yp,zp,ndim,orbi):
    if ndim!=1 or orbi!=0: 
       print "Error: local rot matrix for s orbital is of dim 1 !"
       sys.exit(0)

    rmat=np.zeros((ndim,ndim),dtype=np.float64)

    for j in range(ndim):
        e=(forb[orbi][j](xp,yp,zp)).expand()
        rmat[j,0] = e.subs(x,0).subs(y,0).subs(z,0).evalf()
    return rmat

def get_orb_map_p(xp,yp,zp,ndim,orbi):
    if ndim!=3 or orbi!=1: 
       print "Error: local rot matrix for p orbital is of dim 3 !"
       sys.exit(0)

    fst=[x,y,z]
    #    0 1 2

    rmat=np.zeros((ndim,ndim),dtype=np.float64)
    for j in range(ndim):
        e=(forb[orbi][j](xp,yp,zp)).expand()

        # " get coefficients of pz "
        t=e
        for st in fst:
            if st != z: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c0=t
        rmat[0,j]=c0.evalf()
        
        # " get coefficients of px "
        t=e
        for st in fst:
            if st != x: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c1=t
        rmat[1,j]=c1.evalf()
        
        # " get coefficients of py "
        t=e
        for st in fst:
            if st != y: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c2=t
        rmat[2,j]=c2.evalf()
        
        #print pretty([c0,c1,c2])
    return rmat

def get_orb_map_d(xp,yp,zp,ndim,orbi):
    if ndim!=5 or orbi!=2: 
       print "Error: local rot matrix for d orbital is of dim 5 !"
       sys.exit(0)

    x2 =x*x; xy=x*y; xz=x*z; yz=y*z; y2=y*y; z2=z*z;
    fst=[x2,xy,xz,y2,yz,z2]
    #    0  1  2  3  4  5

    rmat=np.zeros((ndim,ndim),dtype=np.float64)
    for j in range(ndim):
        e=(forb[orbi][j](xp,yp,zp)).expand()

        # " get coefficients of dz2 "
        t=e
        for st in fst:
            if st != z2: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c0=t*(sympy.sqrt(3))
        rmat[0,j]=c0.evalf()
        # " get coefficients of dxz "
        t=e
        for st in fst:
            if st != xz: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c1=t
        rmat[1,j]=c1.evalf()
        
        # " get coefficients of dyz "
        t=e
        for st in fst:
            if st != yz: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c2=t
        rmat[2,j]=c2.evalf()

        # " get coefficients of dx2-y2 "
        t=e
        for st in fst:
            if st != x2: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c3=2*t+c0/(sympy.sqrt(3))
        rmat[3,j]=c3.evalf()
        
        # " get coefficients of dxy "
        t=e
        for st in fst:
            if st != xy: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c4=t
        rmat[4,j]=c4.evalf()
        
        #print pretty([c0,c1,c2,c3,c4])
    return rmat

def get_orb_map_f(xp,yp,zp,ndim,orbi):
    if ndim!=7 or orbi!=3: 
       print "Error: local rot matrix for f orbital is of dim 7 !"
       sys.exit(0)

    x3 =x*x*x; xy2=x*y*y; xz2=x*z*z;
    yx2=y*x*x; y3 =y*y*y; yz2=y*z*z;
    zx2=z*x*x; zy2=z*y*y; z3 =z*z*z; 
    xyz=x*y*z;
    fst=[x3,xy2,xz2,yx2,y3,yz2,zx2,zy2,z3,xyz]
    #    0   1   2   3  4   5  6   7   8  9

    rmat=np.zeros((ndim,ndim),dtype=np.float64)
    for j in range(ndim):
        e=(forb[orbi][j](xp,yp,zp)).expand()

        # " get coefficients of fz3 "
        t=e
        for st in fst:
            if st != z3: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c0=t*sympy.sqrt(15)
        rmat[0,j]=c0.evalf()
        
        # " get coefficients of fxz2 "
        t=e
        for st in fst:
            if st != xz2: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c1=t*sympy.sqrt(10)/2
        rmat[1,j]=c1.evalf()
        
        # " get coefficients of fyz2 "
        t=e
        for st in fst:
            if st != yz2: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c2=t*sympy.sqrt(10)/2
        rmat[2,j]=c2.evalf()

        # " get coefficients of fz(x2-y2) "
        t=e
        for st in fst:
            if st != zx2 :
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c3=2*t+3*c0/(sympy.sqrt(15))
        rmat[3,j]=c3.evalf()
        
        # " get coefficients of fxyz "
        t=e
        for st in fst:
            if st != xyz: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c4=t
        rmat[4,j]=c4.evalf()
        
        # " get coefficients of fx(x2-3y2) "
        t=e
        for st in fst:
            if st != x3: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c5=(2*t+c1/(sympy.sqrt(10)))*sympy.sqrt(6) 
        rmat[5,j]=c5.evalf()
        
        # " get coefficients of fy(3x2-y2) "
        t=e
        for st in fst:
            if st != y3: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c6=(-2*t-c2/(sympy.sqrt(10)))*sympy.sqrt(6) 
        rmat[6,j]=c6.evalf()
        
        #print pretty([c0,c1,c2,c3,c4,c5,c6])
    return rmat

def get_orb_map_t2g(xp,yp,zp,ndim,orbi):
    if ndim!=3 or orbi!=4: 
       print "Error: local rot matrix for t2g orbital is of dim 3 !"
       sys.exit(0)

    xy=x*y; xz=x*z; yz=y*z; xx=x*x; yy=y*y; zz=z*z
    fst=[xz,yz,xy,xx,yy,zz]
    #    0  1  2

    rmat=np.zeros((ndim,ndim),dtype=np.float64)
    for j in range(ndim):
        e=(forb[orbi][j](xp,yp,zp)).expand()

        # " get coefficients of dxz "
        print "e", e
        t=e
        for st in fst:
            if st != xz: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c0=t
        #print "c0",c0
        rmat[0,j]=c0.evalf()
        
        # " get coefficients of dyz "
        t=e
        for st in fst:
            if st != yz: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c1=t
        rmat[1,j]=c1.evalf()
        
        # " get coefficients of dxy "
        t=e
        for st in fst:
            if st != xy: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c2=t
        rmat[2,j]=c2.evalf()
        
        #print pretty([c0,c1,c2])
    return rmat

def get_any_rot_orb_twostep(case,rot):
    """
    
    Args:
        case: label for different kinds of orbitals
        rot: rotation matrix in local coordinates defined as Rloc(wB,xA)=wB*e_wB*Rg*e_xA
    """

    invrot = np.zeros((3,3),dtype=np.float64)
    invrot = inv(np.float64(rot))
       
    # map of variables in function f(x,y,z)->f(xp,yp,zp)
    xp=np.dot(invrot[0,:],np.transpose([x,y,z]))
    yp=np.dot(invrot[1,:],np.transpose([x,y,z]))
    zp=np.dot(invrot[2,:],np.transpose([x,y,z]))

    if case.strip() == 's': 
        ndim = 1
        orbi = 0
        rmat = get_orb_map_s(xp,yp,zp,ndim,orbi)
    elif case.strip() == 'p': 
        ndim = 3
        orbi = 1
        rmat = get_orb_map_p(xp,yp,zp,ndim,orbi)
    elif case.strip() == 'd': 
        ndim = 5
        orbi = 2
        rmat = get_orb_map_d(xp,yp,zp,ndim,orbi)
    elif case.strip() == 'f': 
        ndim = 7
        orbi = 3
        rmat = get_orb_map_f(xp,yp,zp,ndim,orbi)
    elif case.strip() == 't2g':
        ndim = 3
        orbi = 4
        rmat = get_orb_map_t2g(xp,yp,zp,ndim,orbi)
    else:
        print "don't support orbitals for this case: ", case
        return
# make infinite small element zero
    for i in range(ndim):
        for j in range(ndim):
            if abs(rmat[j,i])<1.0E-6:
                rmat[j,i]=0.0
    return rmat
