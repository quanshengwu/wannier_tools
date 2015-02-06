obj =  module.o main.o readinput.o readHmnR.o parse.o

f90  = mpif90 # -check all

#FLAGS =  -O3 -nogen-interface  -warn all 
flag = -O3 -nogen-interface  -warn all 

libs = -L/opt/intel/mkl/lib/ \
		-lmkl_intel_lp64 -lmkl_intel_thread \
		-lmkl_core -liomp5
 
main :  $(obj)
	$(f90) $(obj) -o get_hopping  $(libs) 

.SUFFIXES: .o .f90

.f90.o :
	$(f90) -c $(flag) $(includes) $*.f90

clean :
	rm -f *.o *.mod *~ get_hopping
