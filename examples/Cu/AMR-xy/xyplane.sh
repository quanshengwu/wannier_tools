#!/bin/bash

# alpha is the angle between the magnetic field and the z' axis in the z'-b plane, where z' axis is
# perpendicular to a and b axis.

for ((iphi=0; iphi<=36; iphi++))
do

theta=90
phi=`echo "$iphi*5"|bc`
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
EF_broadening = 0.05  ! in eV, a broadening factor to choose the k points for integration
Nk1 =81            ! Kmesh(1) for KCUBE_BULK
Nk2 =81            ! Kmesh(2) for KCUBE_BULK
Nk3 =81            ! Kmesh(3) for KCUBE_BULK
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
#!/bin/bash -l

##SBATCH --exclusive
#SBATCH --account=hmt03
#SBATCH --time=72:00:00
##SBATCH --exclude=hpcc[154,155,156,114]
#SBATCH --nodes=1
##SBATCH --gres=gpu:4
#SBATCH --partition=long
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=56
#SBATCH --job-name=vasp_run
#SBATCH --output=./log
#SBATCH --error=./errormsg
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MV2_ENABLE_AFFINITY=0
echo "The current job ID is $SLURM_JOB_ID"
echo "Running on $SLURM_JOB_NUM_NODES nodes:"
echo $SLURM_JOB_NODELIST
echo "Using $SLURM_NTASKS_PER_NODE tasks per node"
echo "A total of $SLURM_NTASKS tasks is used"

ulimit -s unlimited
ulimit -c unlimited
module load cuda11.8
module load oneapi22.3
module load nvhpc/22.11

mpirun /home/liuzh/Wanniertools/wannier_tools/bin/wt.x

echo work done
EOF2

cp wannier90_hr.dat_nsymm48 $dir/
cd $dir
sbatch wt-theta.sh
cd ..
done

