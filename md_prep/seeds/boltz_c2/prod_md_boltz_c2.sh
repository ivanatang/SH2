#!/bin/bash

#SBATCH --job-name=boltz_c2_prod
#SBATCH --output=output_prod_boltz_c2_%j.out
#SBATCH --error=error_prod_boltz_c2_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate SH2

DIR=/projects/ivta1597/SH2/md_prep
MDP=$DIR/mdp
TOP=$DIR/seeds/boltz_c2/gromacs/boltz_c2_dodecahedron.top
SEQ_DIR=/scratch/alpine/ivta1597/SH2/seeds/boltz_c2

cd "$SEQ_DIR"
mkdir -p prod_md_48x1
cd prod_md_48x1
gmx grompp -f $MDP/prod_md_200ns.mdp -c $SEQ_DIR/NPT/npt.gro -t $SEQ_DIR/NPT/npt.cpt -p $TOP -o prod_md_200ns.tpr

# Same benchmark-determined settings as prod_md_SH2.sh: 48-rank real MPI,
# ntomp=1, auto DD/PME (measured 260.36 +/- 2.54 ns/day for this system size
# on this partition -- see prod_md_SH2.sh's comment for the full benchmark
# rationale). 200 ns at that rate is ~18.4 h of compute, well past this
# script's 6h allocation -- this launches the run and checkpoints; resubmit
# resume_prod_boltz_c2.sh (unmodified, repeatedly) to continue to completion.
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_200ns -ntomp $SLURM_CPUS_PER_TASK
