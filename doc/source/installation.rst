.. _installation:

Installation of WannierTools (Linux or Mac)
===========================================

Prerequisites
-------------

You need to install the following (mandatory) packages:

* Fortran compiler (Gfortran or ifort)
* MPICH version higher than 2.1.5 
* Lapack and Blas library, (Intel MKL recommended)
* Arpack-ng

Compilation
-----------

First Check out the repository by ::

  git clone https://github.com/quanshengwu/wannier_tools.git

Or download the .zip file directly from https://github.com/quanshengwu/wannier_tools, then uncompress it

Then Go into wannier_tools/src directory, Choose and Edit Makefile, Change the blas library " libs= " to your lapack+blas library

At present, we prepared 3 typical Makefiles, which are squential+gfotran, sequential+ifort and mpi+ifort. 

For the mpi compiler, you should switch on the compile flag "-DMPI", see Makefile.intel-mpi

After the compliation, the binary 'wann_tools' is copied to wannier_tools/bin/, you can put this path to the system PATH with ::

  export PATH=/where/you/downloaded/wannier_tools/bin:$PATH

to the ``.bashrc`` file in your home directory.

.. _usage:

Usage
-----
Now you can enjoy your exploration for topological materials with WannierTools.

There are two files you have to prepare, 

1. wt.in. All the control and user specified parameters are inclued in this file.
2. wannier90_hr.dat. Tight binding model constructed by Wannier90 or written in the format as wannier90_hr.dat.

After the preparation of these two files, you can just run wann_tools in the same folder ::

   wt.x &

or in multi-cores  ::

   mpirun -np 4 wt.x &

The output information during the running are written in WT.out. 

.. _usefultools: 

Plotting tools
---------------------

1. `gnuplot <http://gnuplot.sourceforge.net>`_
2. `xmgrace <http://plasma-gate.weizmann.ac.il/Grace/>`_
3. `xcrysden <http://www.xcrysden.org>`_
4. `matlab <http://mathworks.com>`_


