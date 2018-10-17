
include ../make.inc

.SUFFIXES: .o .f90

MOD =  module.o
OBJ =  aux.o symmetry.o readHmnR.o inverse.o proteus.o \
			 eigen.o ham_qlayer2qlayer.o psi.o \
			 ham_slab.o ek_slab.o ek_bulk_polar.o ek_bulk2D.o \
			 readinput.o fermisurface.o surfgreen.o surfstat.o \
			 mat_mul.o ham_ribbon.o ek_ribbon.o ham_bulk.o ek_bulk.o \
			 fermiarc.o spintexture.o berrycurvature.o \
			 wanniercenter.o dos.o  orbital_momenta.o \
			 landau_level.o berry.o wanniercenter_adaptive.o \
			 effective_mass.o findnodes.o sigma_AHC.o \
			 main.o

default : wt.x

$(OBJ) : $(MOD)

wt.x : $(MOD) $(OBJ)
	@echo "Linking $@"
	@$(F90) $(LFLAG) $(OBJ) $(MOD) -o $@ $(LIBS)
	@cp -f $@ ../bin

.f90.o :
	@echo "Compiling $@"
	@$(CPP) $(CPPFLAG) $*.f90 > $*_cpp.f90
	@$(F90) $(FFLAG) $(IFLAG) -c $*_cpp.f90 -o $*.o

.PHONY : clean
clean :
	@echo "Cleaning..."
	@rm -f *.o *.mod *_cpp.f90 wt.x
