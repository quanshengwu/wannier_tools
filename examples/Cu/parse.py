import re

def parse_fortran_namelist(namelist):
    # Remove comments
    namelist = re.sub(r'!.*', '', namelist)
    # Parse parameters
    pattern = r"(\w+)\s*=\s*([^\s]+)"
    matches = re.findall(pattern, namelist)
    parameters = {key: float(value) if '.' in value else int(value) for key, value in matches}
    return parameters

namelist = """
&PARAMETERS
OmegaNum = 1        ! omega number
OmegaMin =  0.0     ! energy interval
OmegaMax =  0.0     ! energy interval E_i= OmegaMin+ (OmegaMax-OmegaMin)/(OmegaNum-1)*(i-1)
Nk1 =41            ! Kmesh(1) for KCUBE_BULK
Nk2 =41            ! Kmesh(2) for KCUBE_BULK
Nk3 =41            ! Kmesh(3) for KCUBE_BULK
BTauNum= 100        ! Number of B*tau we calculate
BTauMax = 40.0      ! The maximum B*tau, starting from Btau=0.
Tmin = 30           ! Temperature in Kelvin
Tmax = 330          ! Temperature in Kelvin
NumT = 11           ! number temperature we calculate. T_i=Tmin+(Tmax-Tmin)*(i-1)/(NumT-1)
Nslice_BTau_Max = 20000 ! increase this number if negative magnetoresistance occurs, default =5000
/
"""

parameters = parse_fortran_namelist(namelist)
print(parameters)

