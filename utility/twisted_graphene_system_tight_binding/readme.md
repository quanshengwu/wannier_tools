A Fortran code to generate the tight-binding model TG\_hr.dat for twisted graphene system. TG\_hr.dat is stored as a sparse formated hr.dat.

1. Complilation

```
make &
```

2. prepare system.in file like this

```
! an example for Twisted A-AB system
&PARAMETERS
number_layers = 3    ! number of layers 
twisted_index_m= 1   ! twisted index m, the twist angle theta= acos((3d0*m*m + 3d0*m + 0.5d0)/(3d0*m*m + 3d0*m + 1d0));
twisted_angle_array_input = 0 1 1  ! twisted angle array, unit is theta; Number_of_layers numbers
stacking_sequences_input = "A" "A"  "B" ! AB stacking sequences, only three values "A", "B", "C"; Number_of_layers numbers
use_poscar = F    ! use POSCAR or not, "=F" we will generate POSCAR, "=T" we will use the provided POSCAR to generate hr.dat.
hr_generate = T    ! Generate hr.dat or not,  "=T" we will generate hr.dat
gen_sparse_hr = T   ! Choose sparse or dense stored hr.dat,  "=T" we will generate hr.dat in sparse format
hr_cutoff=0.00010  ! set HmnR=0 if HmnR<hr_cutoff in unit of eV
vpppi=-2.81       ! pi bond of p orbital, in unit of eV
iR_cut = 1   ! R is in [-iR_cut, -iR_cut+1, ..., iR_cut]
/
```

3. run it with 

```
./tgtbgen
```

4. Then you can run WannierTools with the generated wt.in and TG\_hr.dat
