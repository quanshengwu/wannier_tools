obj =  module.o main.o readHmnR.o inverse.o\
       eigen.o ham_qlayer2qlayer.o \
		 ham_slab.o ek_slab.o \
       readinput.o fermisurface.o surfgreen.o surfstat.o \
		 mat_mul.o ham_bulk.o ek_bulk.o \
       fermiarc.o spintexture.o


f90  = mpif90 # -check all

#FLAGS =  -O3 -nogen-interface  -warn all 
flag = -O3 -nogen-interface  -warn all 

libs = -L/opt/intel/mkl/lib/ \
		-lmkl_intel_lp64 -lmkl_sequential \
		-lmkl_core -liomp5
 
main :  $(obj)
	$(f90) $(obj) -o wann_surf_soc $(libs) 
	cp -f wann_surf_soc ../bin

.SUFFIXES: .o .f90

.f90.o :
	$(f90) -c $(flag) $(includes) $*.f90

clean :
	rm -f *.o *.mod *~ wann_surf_soc
