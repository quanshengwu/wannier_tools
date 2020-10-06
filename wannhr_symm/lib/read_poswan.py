#!/usr/bin/env python
'''
Define primitive lattice vectors, atoms' fractional coordinates and local 
axises for each atom
'''
from numpy import array,dot,cross,zeros
import numpy as np
from numpy.linalg import norm

class poswan():
    """ 
    define Bravis cell, wannier projections

    """
    def __init__(self, Amat, natom_pos, atoms_frac, atoms_cart, atoms_symbol, ispinor, natom_wan, type_orbs, num_orbs, atom_num, axis, nband, norbs):
        self.Amat=Amat
        self.natom_pos = natom_pos
        self.atoms_frac = atoms_frac
        self.atoms_cart = atoms_cart
        self.atoms_symbol = atoms_symbol
        self.ispinor = ispinor
        self.natom_wan = natom_wan
        self.type_orbs = type_orbs
        self.num_orbs = num_orbs
        self.atom_num = atom_num
        self.axis = axis
        self.nband = nband
        self.norbs = norbs

    @staticmethod
    def get_poswan(): 
        """
        read poscar.in adjusted from POSCAR
        output Amat, atoms_frac, atoms_cart, natom_wan_pos, atoms_symbol
        """
        # define types and size
        Amat=np.zeros((3,3),dtype=np.float64)
        fname = 'poscar.in'
        try: 
           with open(fname, 'r') as f:
               f.readline() # skip comment
               f.readline() # skip comment
               x,y,z = f.readline().strip().split() 
               Amat[:,0] = np.array([np.float64(x),np.float64(y),np.float64(z)],dtype=np.float64)
               x,y,z = f.readline().strip().split() 
               Amat[:,1] = np.array([np.float64(x),np.float64(y),np.float64(z)],dtype=np.float64)
               x,y,z = f.readline().strip().split() 
               Amat[:,2] = np.array([np.float64(x),np.float64(y),np.float64(z)],dtype=np.float64)
               f.readline() # skip comment
               f.readline() # skip comment
               nt = f.readline().strip().split() 
        
               #number of atoms in crystal cell
               natom_pos = np.int32(nt)[0]
               atoms_frac=np.zeros((3,natom_pos),dtype=np.float64)
               axis = np.zeros((3,3,natom_pos),dtype=np.float64)
               atoms_frac=np.zeros((3,natom_pos),dtype=np.float64)
               atoms_cart=np.zeros((3,natom_pos),dtype=np.float64)
        
               atoms_symbol = []
               for iatom in range(natom_pos):
                   x,y,z,s = f.readline().strip().split() 
                   atoms_symbol.append(s)
                   atoms_frac[:,iatom] = np.array([np.float64(x),np.float64(y),np.float64(z)],dtype=np.float64)
                   atoms_cart[:,iatom] = np.dot(Amat,atoms_frac[:,iatom])
        
               # for the following case, all {fij} 1>fij>=0 
               #begin projections
               #Sn:l=0;l=1
               #end projections
               #ALERT: If in your wannier90.win/projections blocks you set negative fractional coordinates
               #you should be careful. maybe I should define a new wannatom_frac.
	       print "ATTENTION!!! make poscar.in consistent with your wannier90.win and wannier90.wout. The poscar should gives right wannier center"
	       # uncoment the following lines if you do not set negetive fracs. move line # 65 belowing line #80
               #for iatom in range(natom_pos): 
               #     for j in range(3): 
               #         if atoms_frac[j,iatom] < 0.0:
               #            atoms_frac[j,iatom] = atoms_frac[j,iatom] + 1.0
               #         elif atoms_frac[j,iatom] > 0.9999: 
               #            atoms_frac[j,iatom] = atoms_frac[j,iatom] - 1.0
               #    print iatom+1, atoms_frac[:,iatom].T

        except IOError:
            print "File:" + "\"" + fname + "\"" +  " doesn't exist!"
        
        fname = 'wann.in'
        try: 
            # read in wann.in
            with open(fname,'r') as f:
               orb_dim = {'s':1, 'pz':1, 'p':3, 'd':5, 'f':7, 't2g':3}
               type_orbs = []
               num_orbs = []
               atom_num=[]
               nband = 0
               for i in range(23):
                   f.readline() # skip comment
               # numbe of atoms in wannier90.win being projected
               nt = f.readline().strip().split() 
               natom_wan = np.int32(nt)[0]
               # read in isinor
               tmp = f.readline().strip().split() 
               ispinor = tmp[0].upper()[0]=='T'
 
               for iatom in range(natom_wan):
                   t1 = f.readline().strip().split() 
                   atom_num.append(np.int32(t1[1]))
                   nt = np.int32(t1[1])-1
                   ot = []
                   no = []
                   for i in range(2,len(t1)):
                       ot.append(t1[i])
                       no.append(orb_dim[t1[i]])
                       nband = nband + orb_dim[t1[i]]
                   type_orbs.append(ot)
                   num_orbs.append(no)
           
            # spin full dim
            if ispinor:
               norbs = nband * 2
            else: 
               norbs = nband * 1

            # check consistency
            if len(num_orbs) != natom_wan: 
               print "number of items in num_orbs should be same to natom_wan in wannier projection"
               sys.exit(0)
            
            if len(type_orbs) != natom_wan: 
               print "number of items in type_orbs should be same to natom_wan in wannier projection"
               sys.exit(0)
            
            if len(atom_num) != natom_wan: 
               print "number of items in atom_num should be same to natom_wan in wannier projection"
               sys.exit(0)
            
            # number of orbitals for each atom 
            nors = np.zeros((natom_wan),dtype=np.int32)
            
            # index offset of each atom
            offs = np.zeros((natom_wan),dtype=np.int32)
            
            # get offset and num-orbs, and print out infos
            for jatom in range(natom_wan):
               nors[jatom] = sum(num_orbs[jatom][:])
               offs[jatom] = sum(nors[0:jatom])

        except IOError:
            print "File:" + "\"" + fname + "\"" +  " doesn't exist!"

        # read in local axis for atoms involved in wannier90.win/projection
        fname = 'locaxis.in'
        try: 
            # read in local axis
            with open('locaxis.in','r') as f:
                # local axices for atoms in wannier projection, only set z and x axis, and
                # y axis will be got by vector cross.
                for i in range(natom_wan):
                    axis[:,2,i] = [0.0,0.0,1.0]  # z 
                    axis[:,0,i] = [1.0,0.0,0.0]  # x
                for i in range(19):
                    f.readline() # skip comment
                # numbe of atoms in wannier90.win being projected
                nl = np.int32(f.readline().strip().split())
                if (nl!=0 and nl ==natom_wan): 
                   for i in range(nl):
                      nt,zx,zy,zz,xx,xy,xz = f.readline().strip().split() 
                      iw = np.int32(nt)-1
                      axis[:,2,iw] = np.array([zx,zy,zz],dtype=np.float64)
                      axis[:,0,iw] = np.array([xx,xy,xz],dtype=np.float64)
                for i in range(natom_wan):
                    axis[:,1,i] = cross(axis[:,2,i], axis[:,0,i]) # y
                    for j in range(3):
                         axis[:,j,i] = axis[:,j,i] / norm(axis[:,j,i]) # normalization, NECESSARY
                    print i, "local axis\n", axis[:,:,i]
            return poswan(Amat, natom_pos, atoms_frac, atoms_cart, atoms_symbol, ispinor, natom_wan, type_orbs, num_orbs, atom_num, axis, nband, norbs)
        except IOError:
            print "File:" + "\"" + fname + "\"" +  " doesn't exist!"


