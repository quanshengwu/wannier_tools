import numpy as np
import matplotlib.pyplot as plt
from math import *
from numpy import *
import os

current_path = os.getcwd()

path = file_path = os.path.join(current_path,  'wt.in')
with open(path, 'r') as file:
    orbital_value = file.readlines()
    file.close()
for ii in range(len(orbital_value)):
    orbital_value[ii] = orbital_value[ii].split()
#print(orbital_value)

# SOC == 0 : don't consider SOC
# SOC == 1 : consider SOC
SOC = 0
Num_wann = 0
atom_number = 0
atom_name_array = []
atom_orbital_number_array = []
orbital_dict = {"s":0, "pz":0, "px":0, "py":0, "dz2":0, "dxz":0, "dyz":0, "dx2-y2":0, "dxy":0}

for ii in range(len(orbital_value)):
    if(len(orbital_value[ii]) > 0):
        if(orbital_value[ii][0] == 'SOC'):
            temp_string = orbital_value[ii]
            for jj in range(len(temp_string)):
                if(temp_string[jj]=='0' or temp_string[jj]=='1'):
                    SOC = int(temp_string[jj])
                

        if(orbital_value[ii][0] == 'PROJECTORS'):
            
            atom_number_array = orbital_value[ii+1]
            for jj in range(len(atom_number_array)):
                string = atom_number_array[jj]
                # if it is a number: num_flag == 1, else: num_flag == 0 
                num_flag = 1
                for ss in range(len(string)):
                    if( ord(string[ss]) < 48 or ord(string[ss]) > 57 ):
                        num_flag = 0
                        break
                if(num_flag == 1):
                    atom_number = atom_number + 1
                    Num_wann = Num_wann + int(string)
                    atom_orbital_number_array.append(int(string))
            
            for aa in range(atom_number):
                string = orbital_value[ii+2+aa]
                #print(string)
                atom_name_array.append(string[0])

print("atom_number:", atom_number)
print("atom_name_array:", atom_name_array)
print("atom_orbital_number_array:", atom_orbital_number_array)



if(SOC == 1):
    Num_wann = Num_wann * 2

print("Num_wann:", Num_wann)
print("SOC:", SOC)



if(SOC == 1):
    Num_orb = Num_wann // 2
elif(SOC == 0):
    Num_orb = Num_wann
else:
    print("wrong!!!!!!!!!!!!!!!!!!!")


Lx = np.zeros((Num_wann, Num_wann),dtype=complex)
Ly = np.zeros((Num_wann, Num_wann),dtype=complex)
Lz = np.zeros((Num_wann, Num_wann),dtype=complex)


for ii in range(len(orbital_value)):
    if(len(orbital_value[ii]) > 0):
        if(orbital_value[ii][0] == 'PROJECTORS'):
            for aa in range(atom_number):
                orbital_index = {}

                string = orbital_value[ii+2+aa]
                #print(string)
                for ee in range(len(string)):
                    if string[ee] in orbital_dict:
                        #print(ee, ' ', string[ee])
                        orbital_index[string[ee]] = int(ee - 1 + sum(atom_orbital_number_array[:aa]))
                
                print(orbital_index)

                if( 'px' in orbital_index and 'py' in orbital_index and 'pz' in orbital_index ):
                    print('p orbital')
                    Lz[orbital_index['px'], orbital_index['py']] = 1j
                    Lz[orbital_index['py'], orbital_index['px']] =-1j

                    Lx[orbital_index['py'], orbital_index['pz']] = 1j
                    Lx[orbital_index['pz'], orbital_index['py']] =-1j

                    Ly[orbital_index['pz'], orbital_index['px']] = 1j
                    Ly[orbital_index['px'], orbital_index['pz']] =-1j


                if( 'dxy' in orbital_index and 'dyz' in orbital_index and 'dxz' in orbital_index and 'dz2' in orbital_index and 'dx2-y2' in orbital_index ):
                    print('d orbital')

                    Lx[orbital_index['dz2'], orbital_index['dyz']] = -1j*np.sqrt(3)
                    Lx[orbital_index['dxz'], orbital_index['dxy']] = -1j
                    Lx[orbital_index['dyz'], orbital_index['dx2-y2']] = 1j

                    Lx[orbital_index['dyz'], orbital_index['dz2']] = 1j*np.sqrt(3)
                    Lx[orbital_index['dxy'], orbital_index['dxz']] = 1j
                    Lx[orbital_index['dx2-y2'], orbital_index['dyz']] = -1j

                    Ly[orbital_index['dz2'], orbital_index['dxz']] = 1j*np.sqrt(3)
                    Ly[orbital_index['dxz'], orbital_index['dx2-y2']] = 1j
                    Ly[orbital_index['dyz'], orbital_index['dxy']] = 1j

                    Ly[orbital_index['dxz'], orbital_index['dz2']] =  -1j*np.sqrt(3)
                    Ly[orbital_index['dxy'], orbital_index['dyz']] = -1j
                    Ly[orbital_index['dx2-y2'], orbital_index['dxz']] = -1j

                    Lz[orbital_index['dxz'], orbital_index['dyz']] = 1j
                    Lz[orbital_index['dxy'], orbital_index['dx2-y2']] = -1j*2

                    Lz[orbital_index['dyz'], orbital_index['dxz']] = -1j
                    Lz[orbital_index['dx2-y2'], orbital_index['dxy']] = 1j*2
                
                






if(SOC == 1):
    Lx[Num_orb:Num_orb*2,Num_orb:Num_orb*2] = Lx[0:Num_orb,0:Num_orb]
    Ly[Num_orb:Num_orb*2,Num_orb:Num_orb*2] = Ly[0:Num_orb,0:Num_orb]
    Lz[Num_orb:Num_orb*2,Num_orb:Num_orb*2] = Lz[0:Num_orb,0:Num_orb]


print("Lx:\n",Lx)
print("Ly:\n",Ly)
print("Lz:\n",Lz)

path = file_path = os.path.join(current_path,  'oam.dat')

with open(path, 'w') as file:
    file.write(str(Num_wann)+'\n')
    file.write('Lx\n')
    for ii in range(Num_wann):
        for jj in range(Num_wann):
            if(abs(Lx[ii,jj]) > 0.0001):
                file.write( str(ii+1) + '\t' +  str(jj+1) + '\t' + "{:.12f}".format(Lx[ii,jj].real) + '\t' + "{:.12f}".format(Lx[ii,jj].imag) + '\n' )
    
    file.write('Ly\n')
    for ii in range(Num_wann):
        for jj in range(Num_wann):
            if(abs(Ly[ii,jj]) > 0.0001):
                file.write( str(ii+1) + '\t' +  str(jj+1) + '\t' + "{:.12f}".format(Ly[ii,jj].real) + '\t' + "{:.12f}".format(Ly[ii,jj].imag) + '\n' )


    file.write('Lz\n')
    for ii in range(Num_wann):
        for jj in range(Num_wann):
            if(abs(Lz[ii,jj]) > 0.0001):
                file.write( str(ii+1) + '\t' +  str(jj+1) + '\t' + "{:.12f}".format(Lz[ii,jj].real) + '\t' + "{:.12f}".format(Lz[ii,jj].imag) + '\n' )
    




file.close()



                    





                



