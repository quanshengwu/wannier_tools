obj =  module.o eigen.o readinput.o gen_hk.o \
       ek_bulk.o main.o

# compiler
#f90  = gfortran  -cpp
#fcheck =   -ffree-line-length-0 -g -Wextra  -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
#flag = -O3 ${fcheck} #-nogen-interface #  -warn all 

#f90  = mpif90  -fpp -DMPI
#f90  = gfortran -cpp
f90  = ifort -fpp

FLAGS =  -O3 #-nogen-interface  
flag = ${FLAGS}  # -check all -traceback -pg

# blas and lapack libraries
 libs = -L/opt/intel/mkl/lib/ \
		-lmkl_intel_lp64 -lmkl_sequential \
		-lmkl_core        

#libs = -L/usr/lib/ -llapack  

main :  $(obj)
	$(f90) $(obj) -o tg_kpgen $(libs) 
	cp tg_kpgen  ../../bin

.SUFFIXES: .o .f90

.f90.o :
	$(f90) -c $(flag) $(includes) $*.f90

clean :
	rm -f *.o *.mod *~ tg_kpgen
