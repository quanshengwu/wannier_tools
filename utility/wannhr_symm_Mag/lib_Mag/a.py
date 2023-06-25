def get_orb_map_d(xp,yp,zp,ndim,orbi):
    if ndim!=1 or orbi!=6: 
       print "Error: local rot matrix for dz2 orbital is of dim 1 !"
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
        c1=t.evalf()
        
        # " get coefficients of dyz "
        t=e
        for st in fst:
            if st != yz: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c2=t.evalf()

        # " get coefficients of dx2-y2 "
        t=e
        for st in fst:
            if st != x2: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c3=2*t+c0/(sympy.sqrt(3))
        c3=c3.evalf()
        
        # " get coefficients of dxy "
        t=e
        for st in fst:
            if st != xy: 
               t = t.subs(st,0)
            else:
               t = t.subs(st,1)
        c4=t.evalf()
        if abs(c1)+abs(c2)+abs(c3)+abs(c4)>1.0E-3: 
            print "ERROR ! dz2 only is not complete sinze it is rotated to other d orbitals! STOP"     
            sys.exit(0)
        
    return rmat


