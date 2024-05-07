# Author: Hanqi Pi,
# Email : hqpi1999@gmail.com or hqpi@iphy.ac.cn
# Date  : 2023/06/25
# version : 0.2
# ( in testing )

'''
0. Required Packages: PYTHON3, numpy, matplotlib, scipy, re, os
   Required files : wt.in, sigma_bands_mu_*.dat
   Copy this file to your work folder where contains wt.in and sigma files 
   Usage : python post_sigma_OHE.py

1. Function Descriptions

1.1 Function readwtin()

    Description: Reads namelist PARAMETERS from wt.in
    Returns: Tmin, Tmax, NumT, OmegaMin, OmegaMax, OmegaNum, BTauMax, BTauNum

1.2 Function readsigma_writerho(band_list, tau_list)

    Description: Reads sigma/tau tensor data from sigma_bands_mu_***eV.dat and calculates the total resistivity.
    The resistivity is stored in rhotau_total_mu_***eV.dat and rho_all.npy.
    Parameters:
    - band_list: A list specifying the bands to be included, in accordance with the SELECTEDBANDS in wt.in
    - tau_list: A list specifying the ratio of relaxation time for the corresponding bands.
    Note: The values in tau_list cannot be zero. To exclude the contribution of certain bands,
    set the corresponding value in tau_list to a very small value (e.g., 0.00001).

1.3 Function plot_cbar(component, cmap)

    Description: Plots rhotau(Btau).
    Parameters:
    - component: A list specifying the tensor component to be plotted, such as ['xx']

'''

# from cProfile import label
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from scipy.interpolate import RegularGridInterpolator,CubicSpline
import matplotlib
import re
params={'font.family':'serif',
        'font.serif':'Times New Roman',
        'font.style':'normal',# or 'italic
        'font.weight':'normal', #or 'blod'
        'font.size':'15',#or large,small
        }
rcParams.update(params)
#set the font of mathtext same with the others
mpl.rcParams['mathtext.default'] = 'regular'

# read parameters from wt.in
def readwtin():
    with open('wt.in', 'r') as infile:
        file_content = infile.read()
        pattern = re.compile(r'&PARAMETERS\n(.*?)\n/', re.DOTALL)
        parameters_section = pattern.findall(file_content)[0]
        parameter_lines = parameters_section.split('\n')
        parameters = {}
        for line in parameter_lines:
            line = line.split('!')[0].strip()  # Ignore comments
            if '=' in line:
                key, value = line.split('=')
                key = key.upper().strip()  # Convert keys to uppercase
                value = value.strip()
                if '.' in value or 'e' in value.lower():  # The value is a float number
                    parameters[key] = float(value)
                elif ',' in value:  # The value is a tuple of integers
                    parameters[key] = tuple(int(i) for i in value.split(','))
                else:  # The value is an integer
                    parameters[key] = int(value)

    Tmin = parameters['TMIN']
    Tmax = parameters['TMAX']
    NumT = parameters['NUMT']
    OmegaMin = parameters['OMEGAMIN']
    OmegaMax = parameters['OMEGAMAX']
    OmegaNum = parameters['OMEGANUM']
    BTauMax = parameters['BTAUMAX']
    try:
        BTauNum = parameters['BTAUNUM']
    except:
        BTauNum = parameters['NBTAU']

    print('Tmin = ',Tmin, '\nTmax = ',Tmax, '\nNumT = ',NumT, '\nOmegaMin = ',OmegaMin, '\nOmegaMax = ',
            OmegaMax, '\nOmegaNum = ',OmegaNum, '\nBTauMax = ',BTauMax, '\nBTauNum = ',BTauNum)

    return Tmin, Tmax, NumT, OmegaMin, OmegaMax, OmegaNum, BTauMax, BTauNum

# read the data of the magnetic resistivity
def readsigma_writerho(band_list, tau_list):
    
    assert len(band_list) == len(tau_list), 'the length of band_list and tau_list are not equal'
    # number of bands
    Nbands = len(band_list)
    # obtain the max number in tau_list
    tau_max = max(tau_list)

    # read parameters from wt.in
    Tmin, Tmax, NumT, OmegaMin, OmegaMax, OmegaNum, BTauMax, BTauNum = readwtin()

    # generate the list of temperatures and chemical potentials
    Tlist = [Tmin+i*(Tmax-Tmin)/(NumT-1) for i in range(int(NumT))]
    if OmegaNum == 1:
        Omegalist = [OmegaMin]
    else:
        Omegalist = [OmegaMin+i*(OmegaMax-OmegaMin)/(OmegaNum-1) for i in range(int(OmegaNum))]
    # generate the list of btau which is the product of magnetic field and relaxation time
    btaulist= [ i*(BTauMax)/(BTauNum-1) for i in range(int(BTauNum))]
    new_btaumax = BTauMax/tau_max
    new_btaulist = [ i*(new_btaumax)/(BTauNum-1) for i in range(int(BTauNum))]
    
    # check if sigma_bands.npy exists
    if os.path.exists('sigma_bands.npy'):
        sigma_bands = np.load('sigma_bands.npy', allow_pickle=True).item()['sigma_bands']
    else:
        # read band-resolved data from the subdirectories
        # interpolate the data and store them in array 'sigma_bands'
        sigma_bands = np.zeros((OmegaNum, NumT, BTauNum, Nbands, 9))
        # read the conductivities at different potentials
        for imu, mu in enumerate(Omegalist):
            mu = mu + 0.000001
            with open(f'sigma_bands_mu_{mu:.2f}eV.dat', 'r') as rdata:
                print(f'read data at mu = {mu:.2f} eV')
                rdata = rdata.readlines()[4:]
                for ib, band in enumerate(band_list):
                    bandindex = int(float(rdata[ib*(2+NumT*(3+BTauNum+1)+1)].split()[-1]))
                    assert band == bandindex, 'band index is wrong'
                    band_data = rdata[ib*(2+NumT*(3+BTauNum+1)+1)+2 : (ib+1)*((2+NumT*(3+BTauNum+1)+1))-1]
                    print(f'read data at band {band:d}')
                    for iT, T in enumerate(Tlist):
                        # store the data at temperature Tlist[iT] in rdata_local
                        rdata_local = band_data[iT*(3+BTauNum+1)+3 : (iT+1)*(3+BTauNum+1)-1]
                        
                        # store the conductivity in sigma_bands
                        sigma_ori = np.zeros((9, BTauNum))
                        for ibtau, lines in enumerate(rdata_local):
                            for icomp in range(9):
                                # rescale the conductivity by the relaxation time
                                sigma_ori[icomp, ibtau] = float(lines.split()[icomp+1])*tau_list[ib]
                        
                        # rescale the btau by the relaxation time
                        # interpolate the 'sigma_ori' to obtain the rescaled conductivity 
                        # store them in array 'sigma_bands'
                        btau_band = np.array(btaulist)/tau_list[ib]
                        for icomp in range(9):
                            interpolator = CubicSpline(btau_band, sigma_ori[icomp])
                            sigma_bands[imu, iT, :, ib, icomp] = interpolator(new_btaulist)
            
            # # for debug, output the modified sigma*tau for every bands
            # sbandfile = open(f'sigma_bands_mu_{mu:.3f}eV.datnew', 'w')
            # for ib, band in enumerate(band_list):
            #     sbandfile.write(f'# band = {band:d}\n')
            #     for iT, T in enumerate(Tlist):
            #         sbandfile.write(f'# T = {T}K\n')
            #         sbandfile.write('#{:>15s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}\n'.format('BTau','xx','xy','xz','yx','yy','yz','zx','zy','zz'))
            #         for ibtau, btau in enumerate(new_btaulist):
            #             row = [btau] + [sigma_bands[imu, iT, ibtau, ib, i] for i in range(9)]
            #             sbandfile.write('{:>15.6f}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}\n'.format(*row))

        # store Tlist, Omegalist, btaulist and sigma_bands in sigma_bands.npy in the form of dictionary
        np.save('sigma_bands.npy',{'Tlist':Tlist, 'Omegalist':Omegalist, 'btaulist':new_btaulist, 'sigma_bands':sigma_bands}, allow_pickle=True)
        
    
    # write sigma and rhotau file
    rho_all = np.zeros((OmegaNum, NumT, BTauNum, 9))
    
    for imu, mu in enumerate(Omegalist):
        # sum the conductivity and resistivity of all bands and write them in 
        # sigma_total_mu_*.dat and rhotau_total_mu_*.dat

        sigmafile = open(f'sigma_total_mu_{mu:.3f}eV.dat', 'w')
        rhofile = open(f'rhotau_total_mu_{mu:.3f}eV.dat', 'w')
        for iT, T in enumerate(Tlist):
            sigmafile.write(f'# T = {T}K\n')
            sigmafile.write('#{:>15s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}\n'.format('BTau','xx','xy','xz','yx','yy','yz','zx','zy','zz'))
            rhofile.write(f'# T = {T}K\n')
            rhofile.write('#{:>15s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}\n'.format('BTau','xx','xy','xz','yx','yy','yz','zx','zy','zz'))
            for ibtau, btau in enumerate(new_btaulist):
                sigma  = sum(sigma_bands[imu, iT, ibtau])
                row = [btau] + [sigma[i] for i in range(9)]
                sigmafile.write('{:>15.6f}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}\n'.format(*row))
                if abs(np.linalg.det(sigma.reshape(3,3)))>1e-9:
                    rho = np.linalg.inv(sigma.reshape(3,3))
                    rho_all[imu, iT, ibtau] = rho.reshape(1,9)[0] 
                    row = [btau] + [rho.reshape(1,9)[0][i] for i in range(9)]
                    rhofile.write('{:>15.6f}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}\n'.format(*row))
                else:
                    rhofile.write('# error: sigma is zero since no k points contribute to the calculations of MR\n')
                    rho_all[imu, iT, ibtau] = np.NaN
            sigmafile.write('\n')
            rhofile.write('\n')
    # store Tlist, Omegalist, btaulist and rho_all in rho_all.npy in the form of dictionary
    np.save('rho_all.npy', {'Tlist':Tlist, 'Omegalist':Omegalist, 'btaulist':new_btaulist, 'rho_all':rho_all}, allow_pickle=True)

def plot_cbar(component=['xx'], cmap='rainbow'):
    
    Ncomp = len(component)    
    # read Tlist, Omegalist, btaulist and rho_all stored as dictionary from rho_all.npy 
    # Note that the true chemical potential is mu = Omega + 0.04 eV
    Omegalist = np.load('rho_all.npy', allow_pickle=True).item()['Omegalist']
    Tlist = np.load('rho_all.npy', allow_pickle=True).item()['Tlist']
    btaulist = np.load('rho_all.npy', allow_pickle=True).item()['btaulist']
    rho_all = np.load('rho_all.npy', allow_pickle=True).item()['rho_all']
    
    OmegaNum = len(Omegalist)
    NumT = len(Tlist)
    BTauNum = len(btaulist)

    rho = np.zeros((OmegaNum, Ncomp, NumT, BTauNum))
    for imu, mu in enumerate(Omegalist):
        for iT, T in enumerate(Tlist):
            icomp = 0
            if 'xx' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 0] 
                icomp += 1
            if 'xy' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 1] 
                icomp += 1
            if 'xz' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 2] 
                icomp += 1
            if 'yx' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 3] 
                icomp += 1
            if 'yy' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 4] 
                icomp += 1
            if 'yz' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 5] 
                icomp += 1
            if 'zx' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 6] 
                icomp += 1
            if 'zy' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 7] 
                icomp += 1
            if 'zz' in component:
                rho[imu, icomp, iT, :] = rho_all[imu, iT, :, 8] 
                icomp += 1
            assert icomp == Ncomp 
    
        # The subplots will have the same width
        # and there will be a horizontal space of 0.4 between the subplots
        fig, axs = plt.subplots(figsize=(12,5), nrows=1, ncols=Ncomp,gridspec_kw={'width_ratios': [1 for i in range(Ncomp)], 'wspace': 0.4})
        print('plot data at mu = {:.3f} eV'.format((mu+0.00001)))

        cmap = plt.get_cmap(cmap,len(Tlist))
        for icomp, comp in enumerate(component):
            for iT, T in enumerate(Tlist):
                if T:
                    if Ncomp >1:
                        if comp == 'xx' or comp == 'yy' or comp =='zz':
                            # axs[icomp].plot(btaulist, rho[imu, icomp, iT, :], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                            axs[icomp].plot(btaulist, (rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0])/rho[imu, icomp, iT, 0]*100, linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                        else:
                            axs[icomp].plot(btaulist, rho[imu, icomp, iT, :], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                            # axs[icomp].plot(btaulist, rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                    else:
                        if comp == 'xx' or comp == 'yy' or comp =='zz':
                            # axs.plot(btaulist, rho[imu, icomp, iT, :], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                            axs.plot(btaulist, (rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0])/rho[imu, icomp, iT, 0]*100, linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                        else:
                            axs.plot(btaulist, rho[imu, icomp, iT, :], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                            # axs.plot(btaulist, rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))

            if Ncomp >1:
                 axs[icomp].set_xlabel(r'$B\tau\ (T\cdot ps)$')
                 if comp == 'xx' or comp == 'yy' or comp =='zz':
                    #  axs[icomp].set_ylabel(r'$ \rho_{%s}\tau\ (\Omega\cdot m\cdot s) $'%(comp))
                     axs[icomp].set_ylabel(r'$ MR_{%s}\  (\%%)$'%(comp))
                 else:
                     axs[icomp].set_ylabel(r'$ \rho_{%s}\tau\ \  (\Omega\cdot m\cdot s) $'%(comp))
                     #  axs[icomp].set_ylabel(r'$ (\rho_{%s}-\rho_0)\ (\Omega\cdot m) $'%(comp))
            else:
                 axs.set_xlabel(r'$B\tau\ (T\cdot ps)$')
                 if comp == 'xx' or comp == 'yy' or comp =='zz':
                     axs.set_ylabel(r'$ MR_{%s} $'%(comp))
                 else:
                     axs.set_ylabel(r'$ (\rho_{%s}-\rho_0)\tau\ \ (\Omega\cdot m\cdot s) $'%(comp))
            
            
        # Add title
        fig.suptitle(r'$\mu=%.0f$ meV'%((mu+0.00001)*1000), ha='center', y=0.95)
   
        # Adjust the subplot margins
        plt.subplots_adjust(top=0.8)
        plt.subplots_adjust(bottom=0.15)

        norm = mpl.colors.Normalize(vmin=Tlist[0],vmax=Tlist[-1])
        sm  = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        clb = plt.colorbar(sm, ticks=np.linspace(Tlist[0], Tlist[-1],  5))
        clb.ax.set_title('T(K)')
        plt.savefig('%.0fmeV.pdf'%(mu*1000), bbox_inches='tight')
        #plt.show()
        plt.close()

                

if __name__ == '__main__':

    # this setting is for Cu example, modify the settings on your own systems
    
    # band list set in SELECTEDBANDS in wt.in. e.g. band_list=[6,7]
    band_list=[6]

    # tau list for each band, unit is ps. e.g. tau_list= [1.000, 2.000], tau_list=[0.0001, 1.000]
    # tau could be very small like 0.0001 to ignore the contribution from one band, but can not be zero.
    tau_list= [1]  

    # note: band_list and tau_list should have the same length

    component=['xx','yx','zz']

    # if you want to remove the existed 'sigma.npy' and 'rho.npy', please uncomment the following two lines
    # os.system('rm sigma_bands.npy')
    # os.system('rm rho_all.npy')

    # start to calculate the sigma and rho with provided band_list and tau_list
    os.system('rm *.npy')
    readsigma_writerho(band_list, tau_list)

    # plot MR or Hall 

   
    # cmap could be 'viridis', 'plasma', 'inferno', 'magma', 'cividis', 
    # 'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
    # 'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
    # 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'
    # flag', 'prism', 'ocean', 'gist_earth', 'terrain',
    # 'gist_stern', 'gnuplot', 'gnuplot2', 'CMRmap',
    # 'cubehelix', 'brg', 'gist_rainbow', 'rainbow', 'jet',
    # 'turbo', 'nipy_spectral', 'gist_ncar'
    # 'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
    # 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'
    # full list of cmap is shown in https://matplotlib.org/stable/tutorials/colors/colormaps.html
    cmap='rainbow'
    plot_cbar(component=component, cmap=cmap)
