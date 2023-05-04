# Author: Hanqi Pi,
# Email : hqpi1999@gmail.com or hqpi@iphy.ac.cn
# Date  : 2023/05/04
# version : 0.1
# ( in testing )

'''
1.function readsigma_writerho() 
  This function read sigma/tau tensor data from every single band, adds them together and calculate the total resistivity. 
  You can modify the script if you need electrons in different bands to have different relaxation time.
  This function reads paramenters from `wt.in` and sigma tensor from `sigma_bands_mu_***eV.dat`.
  The resistivity is stored in `rhotau_total_mu_***eV.dat` and `rho_all.npy`

2.function plot_cbar(component)
  This function is designed to plot rhotau(Btau).
  You should provide a list called `component` that includes the tensor component that you want to plot, such as ['xx']


'''

from cProfile import label
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
params={'font.family':'serif',
        'font.serif':'Times New Roman',
        'font.style':'normal',# or 'italic
        'font.weight':'normal', #or 'blod'
        'font.size':'15',#or large,small
        }
rcParams.update(params)
#set the font of mathtext same with the others
mpl.rcParams['mathtext.default'] = 'regular'

def readsigma_writerho():

    # read parameters from wt.in
    with open('wt.in', 'r') as rdata:
        rdata = rdata.readlines()
        for cnt, lines in enumerate(rdata):
            # read lines containing 'Tmin' or raise the error


            if 'Tmin' in lines:
                Tmin = float(lines.split()[2])
                print('Tmin = ',Tmin)
            if 'Tmax' in lines:
                Tmax = float(lines.split()[2])
                print('Tmax = ',Tmax)
            if 'NumT' in lines:
                NumT = int(float(lines.split()[2]))
                print('NumT = ',NumT)
            if 'BTauNum' in lines:
                BTauNum = int(float(lines.split()[2]))
                print('BTauNum = ',BTauNum)
            if 'BTauMax' in lines:
                BTauMax = float(lines.split()[2])
                print('BTauMax = ',BTauMax)
            if 'OmegaNum' in lines:
                OmegaNum = int(float(lines.split()[2]))
                print('OmegaNum = ',OmegaNum)
            if 'OmegaMin' in lines:
                OmegaMin = float(lines.split()[2])
                print('OmegaMin = ',OmegaMin)
            if 'OmegaMax' in lines:
                OmegaMax = float(lines.split()[2])
                print('OmegaMax = ',OmegaMax)
            if 'SELECTEDBANDS' in lines:
                Nbands = int(rdata[cnt+1].split()[0])
                if '!' in rdata[cnt+2]:
                    band_list = [int(i) for i in rdata[cnt+2].split('!')[0].split()]
                else:
                    band_list = [int(i) for i in rdata[cnt+2].split()]
                print('SELECTEDBANDS = ',band_list)

    Tlist = [Tmin+i*(Tmax-Tmin)/(NumT-1) for i in range(int(NumT))]
    if OmegaNum == 1:
        Omegalist = [OmegaMin]
    else:
        Omegalist = [OmegaMin+i*(OmegaMax-OmegaMin)/(OmegaNum-1) for i in range(int(OmegaNum))]
    btaulist= [ i*(BTauMax)/(BTauNum-1) for i in range(BTauNum)]
    
    # read data from sigma_total_mu_*.dat and skip all the comments
    sigma_bands = np.zeros((OmegaNum, Nbands, NumT, BTauNum, 9))
    for imu, mu in enumerate(Omegalist):
        with open(f'sigma_bands_mu_{mu:.2f}eV.dat', 'r') as rdata:
            print(f'#####read data of mu = {mu:.2f} eV#####')
            print('')
            rdata = rdata.readlines()[4:]
            # read data from different bands
            for ib, band in enumerate(band_list):
                data_band = rdata[ib*(2+NumT*(2+BTauNum+1)+1) : (ib+1)*(2+NumT*(2+BTauNum+1)+1)]
                # read data from different temperatures
                for iT, T in enumerate(Tlist):
                    data_T = data_band[2+iT*(2+BTauNum+1)+2 : 2+(iT+1)*(2+BTauNum+1)-1]
                    for ibtau, lines in enumerate(data_T):
                        for icomp in range(9):
                            sigma_bands[imu, ib, iT, ibtau, icomp] = float(lines.split()[icomp+2])

    # write  rhotau file
    rho_all = np.zeros((OmegaNum, NumT, BTauNum, 9))
    
    for imu, mu in enumerate(Omegalist):
        rhofile = open(f'rhotau_total_mu_{mu:.3f}eV.dat', 'w')
        for iT, T in enumerate(Tlist):
            rhofile.write(f'# T = {T}K\n')
            rhofile.write('#{:>15s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}{:>16s}\n'.format('BTau','xx','xy','xz','yx','yy','yz','zx','zy','zz'))
            for ibtau, btau in enumerate(btaulist):
                sigma  = sum(sigma_bands[imu, :, iT, ibtau])
                row = [btau] + [sigma[i] for i in range(9)]
                if abs(np.linalg.det(sigma.reshape(3,3)))>1e-9:
                    rho = np.linalg.inv(sigma.reshape(3,3))
                    rho_all[imu, iT, ibtau] = rho.reshape(1,9)[0] 
                    row = [btau] + [rho.reshape(1,9)[0][i] for i in range(9)]
                    rhofile.write('{:>15.6f}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}{:>16.6e}\n'.format(*row))
                else:
                    rhofile.write('# error: sigma is zero since no k points contribute to the calculations of MR\n')
                    rho_all[imu, iT, ibtau] = np.NaN
    
    np.save('rho_all.npy', {'Tlist':Tlist, 'Omegalist':Omegalist, 'btaulist':btaulist, 'rho_all':rho_all}, allow_pickle=True)

def plot_cbar(component=['xx']):
    
    Ncomp = len(component)    
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

        cmap = plt.get_cmap('jet',len(Tlist))
        for icomp, comp in enumerate(component):
            for iT, T in enumerate(Tlist):
                if Ncomp >1:
                    if comp == 'xx' or comp == 'yy' or comp =='zz':
                        axs[icomp].plot(btaulist, (rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0])/rho[imu, icomp, iT, 0]*100, linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                    else:
                        axs[icomp].plot(btaulist, rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                else:
                    if comp == 'xx' or comp == 'yy' or comp =='zz':
                        axs.plot(btaulist, (rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0])/rho[imu, icomp, iT, 0]*100, linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                    else:
                        axs.plot(btaulist, rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))

            if Ncomp >1:
                 axs[icomp].set_xlabel(r'$B\tau\ (T\cdot ps)$')
                 if comp == 'xx' or comp == 'yy' or comp =='zz':
                     axs[icomp].set_ylabel(r'$ MR_{%s} $'%(comp))
                 else:
                     axs[icomp].set_ylabel(r'$ (\rho_{%s}-\rho_0)\tau\ (\Omega\cdot m\cdot s) $'%(comp))
            else:
                 axs.set_xlabel(r'$B\tau\ (T\cdot ps)$')
                 if comp == 'xx' or comp == 'yy' or comp =='zz':
                     axs.set_ylabel(r'$ MR_{%s} $'%(comp))
                 else:
                     axs.set_ylabel(r'$ (\rho_{%s}-\rho_0)\tau\ (\Omega\cdot m\cdot s) $'%(comp))
            
            

        # Add title
        fig.suptitle(r'$\mu=%.0f$ meV'%(mu*1000), ha='center', y=0.95)
   
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

   
if __name__ == '__main__':
    # read the data from sigma_bands file, you can modify the relaxation time of bands and get the resistivity
    # stored in 'rhotau_total_mu_{mu:.3f}eV.dat' and 'rho_all.npy'
    readsigma_writerho()
    
    # plot the component of rhotau tensor
    component=['xx','xz','zz']
    plot_cbar(component=component)
