obj = module.o addsoc.o amoeba.o calfun.o eigen.o fitsoc.o ham_bulk.o \
		readHmnR.o readinput.o main.o  

f90  = mpif90 # -check all -traceback -p -g

#FLAGS =  -O3 -nogen-interface  -warn all 
flag = -O3 -nogen-interface  -warn all 

libs = -L/opt/intel/mkl/lib/ \
		-lmkl_intel_lp64 -lmkl_intel_thread \
		-lmkl_core -liomp5
 
main :  $(obj)
	$(f90) $(obj) -o fit_soc  $(libs) 
	cp -f  fit_soc ../bin/

.SUFFIXES: .o .f90

.f90.o :
	$(f90) -c $(flag) $(includes) $*.f90

clean :
	rm -f *.o *.mod *~ fit_soc 
