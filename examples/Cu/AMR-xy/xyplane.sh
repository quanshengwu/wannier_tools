#!/bin/bash

# alpha is the angle between the magnetic field and the z' axis in the z'-b plane, where z' axis is
# perpendicular to a and b axis.

for ((iphi=0; iphi<=12; iphi++))
do

theta=90
phi=`echo "$iphi*15"|bc`
dir='Btheta'$theta'Bphi'$phi
echo $theta $phi $dir
mkdir $dir

cat >$dir/wt.in <<EOF
&TB_FILE
Hrfile = 'wannier90_hr.dat_nsymm48'
/


&CONTROL
Boltz_OHE_calc        = T
Symmetry_Import_calc = T ! please set it to be true for magnetoresistance calculation
/

&SYSTEM
SOC = 0                ! without soc : SOC=0; with soc : SOC=1
E_FERMI = 7.7083       ! e-fermi
Btheta= $theta, Bphi= $phi    ! magnetic field direction, Btheta is the angle with z axial, Bphi is the angle with respect to x axial in the x-y plane
NumOccupied = 6        ! set it anyway even don't use it.
/

&PARAMETERS
OmegaNum = 1        ! omega number
OmegaMin = 0     ! energy interval
OmegaMax = 0     ! energy interval E_i= OmegaMin+ (OmegaMax-OmegaMin)/(OmegaNum-1)*(i-1)
EF_integral_range = 0.05  ! in eV, a broadening factor to choose the k points for integration
Nk1 =41            ! Kmesh(1) for KCUBE_BULK
Nk2 =41            ! Kmesh(2) for KCUBE_BULK
Nk3 =41            ! Kmesh(3) for KCUBE_BULK
BTauNum= 101        ! Number of B*tau we calculate
BTauMax = 10.0      ! The maximum B*tau, starting from Btau=0.
Tmin = 30           ! Temperature in Kelvin
Tmax = 120          ! Temperature in Kelvin
NumT = 4           ! number temperature we calculate. T_i=Tmin+(Tmax-Tmin)*(i-1)/(NumT-1)
Nslice_BTau_Max = 20000 ! increase this number if negative magnetoresistance occurs, default =5000 
RKF45_PERIODIC_LEVEL = 1 !
/

LATTICE
Angstrom
   0.0000000   1.8075000   1.8075000
   1.8075000   0.0000000   1.8075000
   1.8075000   1.8075000   0.0000000

ATOM_POSITIONS
1                               ! number of atoms for projectors
Cartisen                          ! Direct or Cartisen coordinate
Cu    0.000000      0.000000      0.000000

PROJECTORS
9            ! number of projectors
Cu s s s  s s dxy dyz dzx dx2-y2 dz2

SURFACE            ! should be given even don't use
 1  0  0
 0  1  0

SELECTEDBANDS
1
6

KCUBE_BULK
 0.00  0.00  0.00   ! Original point for 3D k plane
 1.00  0.00  0.00   ! The first vector to define 3d k space plane
 0.00  1.00  0.00   ! The second vector to define 3d k space plane
 0.00  0.00  1.00   ! The third vector to define 3d k cube
EOF

cat>$dir/wt-theta.sh<<EOF2
#!/bin/bash
#SBATCH -J $phi
#SBATCH -p wzhctdnormal
#SBATCH -N 1
#SBATCH --ntasks-per-node=2
#SBATCH -o out
#SBATCH -e error

module purge
module load compiler/intel/2021.3.0
module load mpi/intelmpi/2021.3.0

export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi.so
export PATH=~/wt2024/bin:$PATH

srun --mpi=pmi2 wt.x
EOF2

cp wannier90_hr.dat_nsymm48 $dir/
cd $dir
sbatch wt-theta.sh
cd ..
done

