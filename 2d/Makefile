obj =  module.o aux.o symmetry.o readHmnR.o inverse.o\
       eigen.o ham_qlayer2qlayer.o \
       readinput.o surfgreen.o surfstat.o \
		 mat_mul.o ham_ribbon.o ek_ribbon.o ham_bulk.o ek_bulk.o \
		 wanniercenter.o fermisurface.o dos.o \
       main.o

# compiler
f90  = mpif90 # -check all -pg -traceback

#FLAGS =  -O3 -nogen-interface  -warn all 
flag = -O3 -nogen-interface #-warn all 

# blas and lapack libraries
libs = -L/opt/intel/mkl/lib/ \
		-lmkl_intel_lp64 -lmkl_sequential \
		-lmkl_core -liomp5
 
main :  $(obj)
	$(f90) $(obj) -o wann_tools_2d $(libs) 
	cp -f wann_tools_2d ../bin

.SUFFIXES: .o .f90

.f90.o :
	$(f90) -c $(flag) $(includes) $*.f90

clean :
	rm -f *.o *.mod *~ wann_tools
