obj =  module.o readinput.o generate_crystal_structure.o \
       generate_hr.o main.o

# compiler
#f90  = gfortran  -cpp
#fcheck =   -ffree-line-length-0 -g -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan

f90=ifort -fpp
fcheck = -nogen-interface  -warn all  -check all
flag = -O3 ${fcheck} #-nogen-interface #  -warn all 

# blas and lapack libraries
libs = -L/opt/intel/mkl/lib/ \
		-lmkl_intel_lp64 -lmkl_sequential \
		-lmkl_core        
 
main :  $(obj)
	$(f90) $(obj) -o tgtbgen $(libs) 
	cp tgtbgen  ../../bin

.SUFFIXES: .o .f90

.f90.o :
	$(f90) -c $(flag) $(includes) $*.f90

clean :
	rm -f *.o *.mod *~ tgtbgen
