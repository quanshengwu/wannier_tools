def get_angle(sn,co)
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


