# pre
# install python3 numpy, matplotlib, scipy, cProfile
# pip install numpy scipy matplotlib cProfile
# 
# usage
# modify 
# python  
# python resol
from cProfile import label
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from scipy.interpolate import RegularGridInterpolator,CubicSpline
import matplotlib

# set the font of the figure 
params={'font.family':'serif',
        'font.serif':'Times New Roman',
        'font.style':'normal',# or 'italic
        'font.weight':'normal', #or 'blod'
        'font.size':'15',#or large,small
        }
rcParams.update(params)
#set the font of mathtext same with the others
matplotlib.rcParams['mathtext.default'] = 'regular'

# generate the directory to calculate the magnetic resistivity
def gen_band_dir(band_list):
    home = os.getcwd()
    for band in band_list:
        print(band)
        os.system('mkdir %s'%band)
        os.system('cp ./wt.in ./wt.slurm %s'%band)
        os.system('ln -s ./wannier90_hr.dat_nsymm8 %s'%band)
        os.chdir(band)
        os.system("sed -i '128s/.*/%s/g' wt.in"%(band.split('band')[1]))
        os.system('sbatch wt.slurm')
        os.chdir(home)

# read the data of the magnetic resistivity
def readsigma_writerho(band_list, tau_list):
    
    assert len(band_list) == len(tau_list), 'the length of band_list and tau_list are not equal'

    # number of bands
    Nbands = len(band_list)

    # obtain the max number in tau_list
    tau_max = max(tau_list)

    # read parameters from wt.in
    with open('wt.in', 'r') as rdata:
        rdata = rdata.readlines()
        for cnt, lines in enumerate(rdata):
            lines = lines.upper()

            # read the parameters in wt.in
            if 'TMIN' in lines:

                # split by '=' and '!'
                Tmin = float(lines.split('=')[1].split('!')[0])
                print('Tmin = ',Tmin)
            if 'TMAX' in lines:
                Tmax = float(lines.split('=')[1].split('!')[0])
                print('Tmax = ',Tmax)
            if 'NUMT' in lines:
                NumT = int(float(lines.split('=')[1].split('!')[0]))
                print('NumT = ',NumT)
            if 'BTAUNUM' in lines or 'NBTAU' in lines :
                BTauNum = int(float(lines.split('=')[1].split('!')[0]))
                print('BTauNum = ',BTauNum)
            if 'BTAUMAX' in lines:
                BTauMax = float(lines.split('=')[1].split('!')[0])
                print('BTauMax = ',BTauMax)
            if 'OMEGANUM' in lines:
                OmegaNum = int(float(lines.split('=')[1].split('!')[0]))
                print('OmegaNum = ',OmegaNum)
            if 'OMEGAMIN' in lines:
                OmegaMin = float(lines.split('=')[1].split('!')[0])
                print('OmegaMin = ',OmegaMin)
            if 'OMEGAMAX' in lines:
                OmegaMax = float(lines.split('=')[1].split('!')[0])
                print('OmegaMax = ',OmegaMax)

    # generate the list of temperatures and chemical potentials
    Tlist = [Tmin+i*(Tmax-Tmin)/(NumT-1) for i in range(int(NumT))]
    if OmegaNum == 1:
        Omegalist = [OmegaMin]
    else:
        Omegalist = [OmegaMin+i*(OmegaMax-OmegaMin)/(OmegaNum-1) for i in range(int(OmegaNum))]
    print('Omegalist:', Omegalist)
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
            mu = (mu+ 0.000001)
            with open(f'sigma_bands_mu_{mu:.2f}eV.dat', 'r') as rdata:
                print(f'read data at mu = {mu:.2f} eV')
                # skip the first 5 lines
                rdata = rdata.readlines()[4:]
                for ib, band in enumerate(band_list):
                    bandindex = int(float(rdata[ib*(3+NumT*(3+BTauNum+1)+1)].split()[-1]))
                    assert band == bandindex, 'band index is wrong'
                    band_data = rdata[ib*(3+NumT*(3+BTauNum+1)+1)+2 : (ib+1)*((3+NumT*(3+BTauNum+1)+1))-1]
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

def plot_cbar(component=['xx']):
    
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

    print('Omegalist', Omegalist)
    print('here debug', OmegaNum, NumT, BTauNum)
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
        print('plot data at mu = {:.3f} eV'.format(mu))

        cmap = plt.get_cmap('jet',len(Tlist))
        for icomp, comp in enumerate(component):
            for iT, T in enumerate(Tlist):
                if T:
                    if Ncomp >1:
                        if comp == 'xx' or comp == 'yy' or comp =='zz':

                            # uncomment this to plot resistivity 
                            axs[icomp].plot(btaulist, rho[imu, icomp, iT, :], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                            # plot MR
                            #axs[icomp].plot(btaulist, (rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0])/rho[imu, icomp, iT, 0]*100, linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                        else:

                            # uncomment this to plot resistivity 
                            axs[icomp].plot(btaulist, rho[imu, icomp, iT, :], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                            # plot rho_xy-rho_xy(0)
                            # axs[icomp].plot(btaulist, rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                    else:
                        if comp == 'xx' or comp == 'yy' or comp =='zz':
                            axs.plot(btaulist, (rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0])/rho[imu, icomp, iT, 0]*100, linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))
                        else:
                            axs.plot(btaulist, rho[imu, icomp, iT, :]-rho[imu, icomp, iT, 0], linewidth=2.5, color=cmap(iT), label='T=%.2f K'%(T))

            if Ncomp >1:
                 axs[icomp].set_xlabel(r'$B\tau\ (T\cdot ps)$')
                 if comp == 'xx' or comp == 'yy' or comp =='zz':
                     axs[icomp].set_ylabel(r'$ \rho_{%s}\tau\ (\Omega\cdot m\cdot s) $'%(comp))
                     #axs[icomp].set_ylabel(r'$ MR_{%s} $'%(comp))
                 else:
                     axs[icomp].set_ylabel(r'$ \rho_{%s}\tau\ (\Omega\cdot m\cdot s) $'%(comp))
                     #axs[icomp].set_ylabel(r'$ (\rho_{%s}-\rho_0)\ (\Omega\cdot m) $'%(comp))
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
        plt.savefig('%.0fmeV.pdf'%((mu+0.000001)*1000), bbox_inches='tight')
        #plt.show()
        plt.close()

def interpolate(component=['xx']):
    # component is a list of components to be plotted
    Ncomp = len(component)

    # read Tlist, Omegalist, btaulist, rho_all from rho_all.npy
    # rho_all is a 4D array with shape (len(Omegalist), len(Tlist), len(btaulist), 9)
    Omegalist = np.load('rho_all.npy', allow_pickle=True).item()['Omegalist']
    Tlist_interpolated = np.load('rho_all.npy', allow_pickle=True).item()['Tlist']
    btaulist = np.load('rho_all.npy', allow_pickle=True).item()['btaulist']
    rho_all = np.load('rho_all.npy', allow_pickle=True).item()['rho_all']
    NBtau = len(btaulist)
    # print(NumT, NOmega, NBtau)

    # read concentration, Tlist, mu from fixdoping_mu.npy
    # mulist is a 2D array with shape (len(concenlist), len(Tlist))
    concenlist = np.load('fixdoping_mu.npy', allow_pickle=True).item()['concentration']
    Tlist = np.load('fixdoping_mu.npy', allow_pickle=True).item()['Tlist']
    mulist = np.load('fixdoping_mu.npy', allow_pickle=True).item()['mu']
    # choose T in Tlist which is in the range of Tlist_interpolated
    Tlist = [T for T in Tlist if (T > min(Tlist_interpolated) and T < max(Tlist_interpolated))]
    Nconcen, NumT = len(concenlist), len(Tlist)


    for comp in component:
        # make sure the component is valid or print error message
        if comp not in ['xx', 'xy', 'xz', 'yx', 'yy', 'yz', 'zx', 'zy', 'zz']:
            print('Invalid component: %s'%(comp))
            return
        
        # generate the array for the component and use 'comp' as the prefix of the array name
        if comp == 'xx':
            rho_comp = rho_all[:, :, :, 0]
        elif comp == 'xy':
            rho_comp = rho_all[:, :, :, 1]
        elif comp == 'xz':
            rho_comp = rho_all[:, :, :, 2]
        elif comp == 'yx':
            rho_comp = rho_all[:, :, :, 3]
        elif comp == 'yy':
            rho_comp = rho_all[:, :, :, 4]
        elif comp == 'yz':
            rho_comp = rho_all[:, :, :, 5]
        elif comp == 'zx':
            rho_comp = rho_all[:, :, :, 6]
        elif comp == 'zy':
            rho_comp = rho_all[:, :, :, 7]
        elif comp == 'zz':
            rho_comp = rho_all[:, :, :, 8]
            
        # # replace the NaN in rho_comp with interpolated values
        # if np.isnan(rho_comp).any():
        #     print('NaN detected in {} component'.format(comp))
        #     nan_mask = np.isnan(rho_comp)
        #     # nan_idx = np.where(nan_mask)
        #     # for idx in zip(nan_idx[0],nan_idx[1],nan_idx[2]):
        #         # print ('Omega, T, btau = ', '{:.3f}'.format(Omegalist[idx[0]]), '{:.3f}'.format(Tlist[idx[1]]), '{:.3f}'.format(btaulist[idx[2]]))                 
            
        #     rho_comp[nan_mask] = np.interp(np.flatnonzero(nan_mask), np.flatnonzero(~nan_mask), rho_comp[~nan_mask])
        #     # print(rho_comp[nan_mask])

        
        ################ define the interpolation function#############
        # print(np.shape(rho_comp))
        interpolator = RegularGridInterpolator((Omegalist, Tlist_interpolated,  btaulist), rho_comp)
        ################ interpolate the resistivity at every concentration and temperature
        ################ store the data in file 'rho{}_fixdoping.dat'.format{comp}
        with open('rho{}_fixdoping.dat'.format(comp), 'w') as f:
            print('write data of rho{}'.format(comp))
            for icon in range(Nconcen):
                concentration = concenlist[icon]
                # write a comment line to 
                print('write rho{} at n = {:.3e}'.format(comp, concentration))
                f.write('# n = {:.3e}\n'.format(concentration))
                for it in range(NumT):
                    T = Tlist[it]
                    mu = mulist[icon, it]
                    f.write('# T = {:.3f} K, mu = {:.3f} eV\n'.format(T,mu))
                    f.write('#{:>15s}{:>16s}\n'.format('BTau',comp))
                    for btau in btaulist:
                        rho_interp= interpolator((mu, T, btau))
                        f.write('{:>15.6f}{:>16.6e}\n'.format(btau, rho_interp))
    

def plot_interpolate(component=['xx']):
    # component is a list of components to be plotted
    Ncomp = len(component)

    # read Tlist, Omegalist, btaulist, rho_all from rho_all.npy
    # rho_all is a 4D array with shape (len(Omegalist), len(Tlist), len(btaulist), 9)
    Omegalist = np.load('rho_all.npy', allow_pickle=True).item()['Omegalist']
    Tlist_interpolated = np.load('rho_all.npy', allow_pickle=True).item()['Tlist']
    btaulist = np.load('rho_all.npy', allow_pickle=True).item()['btaulist']
    rho_all = np.load('rho_all.npy', allow_pickle=True).item()['rho_all']
    NBtau = len(btaulist)
    # print(NumT, NOmega, NBtau)

    # read concentration, Tlist, mu from fixdoping_mu.npy
    # mulist is a 2D array with shape (len(concenlist), len(Tlist))
    concenlist = np.load('fixdoping_mu.npy', allow_pickle=True).item()['concentration']
    Tlist = np.load('fixdoping_mu.npy', allow_pickle=True).item()['Tlist']
    mulist = np.load('fixdoping_mu.npy', allow_pickle=True).item()['mu']
    # choose T in Tlist which is in the range of Tlist_interpolated
    Tlist = [T for T in Tlist if (T > min(Tlist_interpolated) and T < max(Tlist_interpolated))]
    Nconcen, NumT = len(concenlist), len(Tlist)

    
    ################# plot the resistivity at every concentration and temperature
    ################# store the data in file 'rho{}_fixdoping.pdf'.format{comp}
    cmap = plt.get_cmap('jet',len(Tlist))
    for icon in range(Nconcen):
        print('plot data at n = {:.3e}'.format(concenlist[icon]))
        concentration = concenlist[icon]
        # plt.figure()
        fig, axs = plt.subplots(figsize=(13,4), nrows=1, ncols=Ncomp+1 ,
                                gridspec_kw={'width_ratios': [1 for i in range(Ncomp+1)], 'wspace': 0.5})
        fig.suptitle(r'$n = {:.1e}$'.format(concentration)+' cm$^{-3}$')
        for icomp, comp in enumerate(component):
            muvsT = []
            for it in range(NumT):
                T = Tlist[it]
                muvsT.append(mulist[icon, it])
                # print(1+icon*(NumT*(NBtau+2)+1)+it*(NBtau+2))
                rho = np.loadtxt('rho{}_fixdoping.dat'.format(comp), comments='#', skiprows=1+icon*(NumT*(NBtau+2)+1)+it*(NBtau+2), max_rows=len(btaulist))
                ############### the first subplot is the resistivity vs BTau
                # axs[icomp].plot(rho[:,0], rho[:,1], color=cmap(it),label='T = {:.3f} K'.format(T))
                ############### the first subplot is the MR vs BTau
                if 0:# comp == 'xx' or comp =='yy' or comp =='zz':              
                    axs[icomp].plot(np.array(rho[:,0])*T, (rho[:,1]-rho[0,1])/rho[0,1], color=cmap(it),label='T = {:.3f} K'.format(T))
                else:
                    axs[icomp].plot(np.array(rho[:,0])*T, (rho[:,1]-rho[0,1]), color=cmap(it),label='T = {:.3f} K'.format(T))
            
            if comp == 'xx' or comp =='yy' or comp =='zz': 
                # the xlabel of the first subplot is BTau, the xlabel of the second subplot is T(K)
                axs[icomp].set_xlabel('BTau')
                # axs[icomp].set_ylabel('MR{}'.format(comp))
                axs[icomp].set_ylabel('MR{}'.format(comp))
                # set xlim to [0,1]
                axs[icomp].set_xlim([0,20])
                axs[icomp].set_ylim([0,1e-17])
            else:
                # the xlabel of the first subplot is BTau, the xlabel of the second subplot is T(K)
                axs[icomp].set_xlabel('BTau')
                # axs[icomp].set_ylabel('MR{}'.format(comp))
                axs[icomp].set_ylabel(r'$\rho_{}$'.format(comp))
                axs[icomp].set_ylim([0,1e-17])
                axs[icomp].set_xlim([0,10])


        # add colorbar to the first subplot
        norm = mpl.colors.Normalize(vmin=Tlist[0],vmax=Tlist[-1])
        sm  = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        clb = plt.colorbar(sm, ticks=np.linspace(Tlist[0], Tlist[-1],  5))
        clb.ax.set_title('T(K)')

        # the second subplot is the chemical potential vs T
        axs[Ncomp].plot(Tlist, muvsT, color='k')
        axs[Ncomp].set_ylabel(r'$\mu$ (eV)')
        axs[Ncomp].set_xlabel('T(K)')
        
        plt.savefig('rho_{:.3e}.pdf'.format( concentration), bbox_inches='tight')
        # plt.savefig('MR_{:.3e}.pdf'.format( concentration), bbox_inches='tight')
        #plt.show()
        plt.close()
        

def read_BoltzTrap():
    
    # open file 'case.trace_fixdoping' and read the chemical potential at each temperature at each concentration
    # the file contains the data of ten concentrations and every concentration are temperature dependent
    # the temperature varies from 5k to 300k with step 5k
    Tlist = np.arange(5, 305, 5)
    concentration = np.array([-10e18, -7e18, -5E18, -3e18, -1e18, 1e18, 3e18, 5e18, 7e18, 10e18])
    NumT = len(Tlist)
    Ncon = len(concentration)
    mu = np.zeros((Ncon, NumT))
    

    with open('case.trace_fixdoping', 'r') as f:
        # the first line is the comment line, skip it
        f = f.readlines()[1:]

        # generate the loop to read data of every concentration
        for icon in range(Ncon):

            # the temperature-dependent icon-th concentration is stored in the [icon*(NumT+1):icon*(NumT+1)+NumT] rows
            data = f[icon*(NumT+1):icon*(NumT+1)+NumT]

            # generate the loop to read data of every temperature
            # the first column is the temperature, the tenth column is the chemical potential
            # Note that the chemical potential is in unit of Rydberg, so we need to multiply it by 13.6057 to convert it to eV
            for it in range(NumT):
                mu[icon, it] = float(data[it].split()[9])*13.6057

    # store the concentration, Tlist, and mu as a dictionary in the file 'fixdoping_mu.npy'
    np.save('fixdoping_mu.npy', {'concentration':concentration, 'Tlist':Tlist, 'mu':mu}, allow_pickle=True)
                

# main program 
if __name__ == '__main__':
    # define the band list
    band_list=[6]
    tau_list= [1]

    # band_list=['band50']
    component=['xx','yx','yy']
    # component=['xx']
    # component=['xz']

    # remove all npy files in the folder
    os.system('rm -f *.npy')

    readsigma_writerho(band_list, tau_list)
    plot_cbar(component=component)
    # read_BoltzTrap()
    # interpolate(component=component)
    # plot_interpolate(component=component)
