import numpy as np
import matplotlib.pyplot as plt
from math import *
import cmath 
import functools
from numpy import *
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.cm import ScalarMappable


def orbital_HC(ind,colors, labels):

    #with open('E:\\python_code\\transport\\mulit-hall\\LOAM\\haldane\\m0.1_tnnn0.1_orbital_texture.txt', 'r') as phase: 
    #with open('E:\\python_code\\transport\\mulit-hall\\LOAM\\G_G_2311.06447\\GAGAOAM.txt', 'r') as phase: 
    #with open('E:\\python_code\\transport\\mulit-hall\\LOAM\\G_G_2311.06447\\Ga_Ga_0.1_sigma_OHE.txt', 'r') as phase: 
    #with open('E:\\python_code\\transport\\mulit-hall\\LOAM\\kagome\\orbital_cond_lambda0.00_eta0.0005.txt', 'r') as phase:
    #with open('E:\\python_code\\transport\\mulit-hall\\LOAM\\square\\t1_tin1_cond.txt', 'r') as phase: 
    
    with open('Intra_Orbital_Hall_Conductivity.dat', 'r') as phase: 
    #with open('E:\\DFT_calculation\\wannier_tools\\KIn\\sigma_ahc_eta10.00meV.txt', 'r') as phase: 
    
        orbital_value = phase.readlines()
    phase.close()

    for ii in range(len(orbital_value)):
        orbital_value[ii] = orbital_value[ii].split()
    #print(orbital_value)
    orbital_value = orbital_value[4:]
    #for ii in range(len(orbital_value)):
    #    print(len(orbital_value[ii]))
    #plt.figure(figsize=(8,5))

    #print(orbital_value.shape)
    texture = []
    omega_list = []
    print(len(orbital_value[0]))
    for ii in range(len(orbital_value)):
        
        if(abs(float(orbital_value[ii][ind]))>1000000):
            continue
        omega_list.append(float(orbital_value[ii][0]))
        texture.append(-float(orbital_value[ii][ind])*2*pi)
        #sigma_ohe_perE.append(float(orbital_value[ii][2]))
    para = [ii for ii in range(len(texture))]
    
    plt.plot(omega_list, texture,color = colors ,linewidth = 5,label=labels)




    


def main():
    fig = plt.figure(figsize=(8,5))
    sub = fig.add_subplot(111)
    plt.tick_params(size=8,width=1.5,labelsize = 15)
    sub.spines['bottom'].set_linewidth(1.6)
    sub.spines['left'].set_linewidth(1.6)
    sub.spines['right'].set_linewidth(1.6)
    sub.spines['top'].set_linewidth(1.6)

    sub.spines['bottom'].set_color('black')
    sub.spines['top'].set_color('black')
    sub.spines['right'].set_color('black')
    sub.spines['left'].set_color('black')

    color_list = ['#1a9c1a', '#1a9c4c', '#1a9c7f', '#1a7f9c', '#1a4c9c', '#1a1a9c', '#4c1a9c', '#7f1a9c', '#9c1a7f', '#9c1a4c', '#9c1a1a', '#9c4c1a', '#9c7f1a', '#7f9c1a', '#4c9c1a', '#1a9c1a', '#1a9c4c', '#1a9c7f', '#1a7f9c', '#1a4c9c', '#1a1a9c', '#4c1a9c', '#7f1a9c', '#9c1a7f', '#9c1a4c', '#9c1a1a', '#9c4c1a']

    label_list = ['\sigma_Lx_xx', '\sigma_Lx_xy', '\sigma_Lx_xz', '\sigma_Lx_yx', '\sigma_Lx_yy', '\sigma_Lx_yz', '\sigma_Lx_zx', '\sigma_Lx_zy', '\sigma_Lx_zz', 
                  '\sigma_Ly_xx', '\sigma_Ly_xy', '\sigma_Ly_xz', '\sigma_Ly_yx', '\sigma_Ly_yy', '\sigma_Ly_yz', '\sigma_Ly_zx', '\sigma_Ly_zy', '\sigma_Ly_zz',
                  '\sigma_Lz_xx', '\sigma_Lz_xy', '\sigma_Lz_xz', '\sigma_Lz_yx', '\sigma_Lz_yy', '\sigma_Lz_yz', '\sigma_Lz_zx', '\sigma_Lz_zy', '\sigma_Lz_zz']

    for ii in range(18):
        orbital_HC(1+2*ii, color_list[1*ii],label_list[0+1*ii])
    #orbital_HCd()
    #orbital_HCp()

    #plt.plot([-10,15],[0,0],'orange',linewidth=5,label='ssssssssssss')
    plt.legend(prop = {'size':18 })

    ### CrS2
    plt.xlim([-2.4,1.6])
    plt.ylim([-3,3])
    ### MoS2
    #plt.xlim([-2.2,1.9])
    #plt.xlim([-4,4])
    #plt.ylim([-3,2])

    ### WS2
    #plt.xlim([-6,6])
    #plt.ylim([-4,2])
    plt.legend()
    #plt.show()
    plt.savefig('OHC_Lx_and_Ly.pdf')
    
if __name__ == '__main__':
    main()
