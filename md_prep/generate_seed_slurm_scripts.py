"""
Generate per-seed EM / equilibration / 200 ns production (+ resume) SLURM
scripts for the 8 AF3+Boltz-2 cluster-medoid seeds, following the exact
structure/conventions of the existing single-system scripts
(em_SH2.sh, equil_SH2.sh, prod_md_SH2.sh, xtnd_prod_SH2.sh) -- same SBATCH
directives, module loads, conda env, and gmx invocation pattern, just
pointed at each seed's own parameterized system
(md_prep/seeds/<name>/gromacs/<name>_dodecahedron.gro/.top) and its own
scratch output tree (kept separate per seed so 8 concurrent/sequential jobs
don't collide).

Production target: 200 ns (mdp/prod_md_200ns.mdp, nsteps=100000000). At the
benchmarked 260.36 ns/day (see prod_md_SH2.sh), that's ~18.4 h of compute,
so prod_md_{seed}.sh requests 22h (~20% margin) and should complete the
full 200 ns in one submission under normal conditions. resume_prod_{seed}.sh
is still generated as a fallback (matching xtnd_prod_SH2.sh's -s/-cpi/-append
pattern) in case a run gets interrupted before finishing.

IMPORTANT: unlike prod_md_SH2.sh/xtnd_prod_SH2.sh, neither prod nor resume
script here uses --exclusive -- sbatch rejects it on this partition/qos
("Error 22: The oversubscribe Slurm directive cannot be used with the
provided partition"). An earlier attempt misdiagnosed this as a walltime
issue (24h + --exclusive was rejected, 6h + --exclusive had "worked" before,
so walltime looked like the variable) before confirming --exclusive itself
is the blocker, independent of walltime. Since xtnd_prod_SH2.sh's own
comment documents a ~0.35 ns/day node-sharing slowdown without --exclusive
on the original single-system runs, throughput here isn't guaranteed to
match the 260.36 ns/day benchmark -- use check_prod_performance.sh once
jobs are running to confirm actual throughput rather than assuming it.

Usage: python3 md_prep/generate_seed_slurm_scripts.py
"""
import os

BASE = "/Users/ivanatang/Developer/SH2"
PROJ_DIR = "/projects/ivta1597/SH2/md_prep"   # cluster path, matches existing scripts
SCRATCH_ROOT = "/scratch/alpine/ivta1597/SH2/seeds"  # per-seed, unlike the single-system SEQ_DIR

SEEDS = ["af3_c1", "af3_c2", "af3_c3", "af3_c4", "boltz_c1", "boltz_c2", "boltz_c3", "boltz_c4"]

COMMON_HEADER = """#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --output=output_{tag}_%j.out
#SBATCH --error=error_{tag}_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time={time}
#SBATCH --nodes=1
{extra_sbatch}#SBATCH --qos=cpu-normal
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

DIR={proj_dir}
MDP=$DIR/mdp
TOP=$DIR/seeds/{seed}/gromacs/{seed}_dodecahedron.top
SEQ_DIR={scratch_root}/{seed}
"""

EM_BODY = """
mkdir -p "$SEQ_DIR/EM"
cd "$SEQ_DIR/EM"
gmx grompp -f $MDP/em.mdp \\
  -c $DIR/seeds/{seed}/gromacs/{seed}_dodecahedron.gro \\
  -p $TOP \\
  -o em.tpr
gmx mdrun -deffnm em
"""

EQUIL_BODY = """
# NVT
cd "$SEQ_DIR"
mkdir -p NVT
cd NVT
gmx grompp -f $MDP/nvt.mdp -c $SEQ_DIR/EM/em.gro -r $SEQ_DIR/EM/em.gro -p $TOP -o nvt.tpr
gmx mdrun -deffnm nvt

# NPT
cd "$SEQ_DIR"
mkdir -p NPT
cd NPT
gmx grompp -f $MDP/npt.mdp -c $SEQ_DIR/NVT/nvt.gro -t $SEQ_DIR/NVT/nvt.cpt -r $SEQ_DIR/NVT/nvt.gro -p $TOP -o npt.tpr
gmx mdrun -deffnm npt
"""

PROD_BODY = """
cd "$SEQ_DIR"
mkdir -p prod_md_48x1
cd prod_md_48x1
gmx grompp -f $MDP/prod_md_200ns.mdp -c $SEQ_DIR/NPT/npt.gro -t $SEQ_DIR/NPT/npt.cpt -p $TOP -o prod_md_200ns.tpr

# Same benchmark-determined settings as prod_md_SH2.sh: 48-rank real MPI,
# ntomp=1, auto DD/PME (measured 260.36 +/- 2.54 ns/day for this system size
# on this partition -- see prod_md_SH2.sh's comment for the full benchmark
# rationale). NOTE: unlike prod_md_SH2.sh/xtnd_prod_SH2.sh, --exclusive is
# NOT used here -- sbatch rejects it on this partition/qos ("Error 22: The
# oversubscribe Slurm directive cannot be used with the provided partition").
# xtnd_prod_SH2.sh's own comment says --exclusive was added after a run
# without it measured ~0.35 ns/day from node-sharing DD load imbalance, so
# throughput here isn't guaranteed to match the 260.36 ns/day benchmark --
# check it with check_prod_performance.sh once running, don't just assume.
mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_200ns -ntomp $SLURM_CPUS_PER_TASK
"""

RESUME_HEADER = """#!/bin/bash

#SBATCH --job-name={job_name}
#SBATCH --output=output_{tag}_%j.out
#SBATCH --error=error_{tag}_%j.err
#SBATCH --account=ucb351_asc4
#SBATCH --partition=acpu
#SBATCH --time=22:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --qos=cpu-normal
#SBATCH --mail-user=ivana.tang@colorado.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#
# Resumes an interrupted prod_md_{seed}.sh run from its checkpoint, continuing
# toward the SAME original target (prod_md_200ns.mdp's nsteps=100000000, i.e.
# 200 ns -- NOT extending past 200 ns, just finishing the run that ran out of
# wall time). Mirrors xtnd_prod_SH2.sh's -s/-cpi/-append pattern, but WITHOUT
# --exclusive (matching prod_md_{seed}.sh) -- this partition/qos rejects it
# at sbatch time ("Error 22: The oversubscribe Slurm directive cannot be
# used with the provided partition"), unlike the original single-system
# scripts where it worked. xtnd_prod_SH2.sh's own comment documents a
# ~0.35 ns/day node-sharing slowdown without --exclusive, so check actual
# throughput with check_prod_performance.sh rather than assuming it's fine.
#
# Usage: sbatch resume_prod_{seed}.sh
# (resubmit again, unmodified, if it runs out of wall time again before 200 ns)

export TMPDIR=$SLURM_SCRATCH
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module purge
module load gcc
module load openmpi
module load anaconda
module load gromacs

conda activate SH2

SEQ_DIR={scratch_root}/{seed}

cd "$SEQ_DIR/prod_md_48x1"

if [ ! -f prod_md_200ns.cpt ]; then
    echo "ERROR: no checkpoint found at $SEQ_DIR/prod_md_48x1/prod_md_200ns.cpt -- nothing to resume from." >&2
    exit 1
fi

mpirun -np $SLURM_NTASKS gmx_mpi mdrun -deffnm prod_md_200ns \\
  -s prod_md_200ns.tpr -cpi prod_md_200ns.cpt -append -ntomp $SLURM_CPUS_PER_TASK
"""


def write(path, content):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(content)
    os.chmod(path, 0o755)


for seed in SEEDS:
    seed_dir = f"{BASE}/md_prep/seeds/{seed}"

    em_extra = "#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=12\n#SBATCH --constraint=ib\n"
    em = COMMON_HEADER.format(
        job_name=f"em_{seed}", tag=f"em_{seed}", time="00:30:00",
        extra_sbatch=em_extra, proj_dir=PROJ_DIR, seed=seed, scratch_root=SCRATCH_ROOT,
    ) + EM_BODY.format(seed=seed)
    write(f"{seed_dir}/em_{seed}.sh", em)

    eq_extra = "#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=24\n#SBATCH --constraint=ib\n"
    eq = COMMON_HEADER.format(
        job_name=f"eq_{seed}", tag=f"eq_{seed}", time="01:00:00",
        extra_sbatch=eq_extra, proj_dir=PROJ_DIR, seed=seed, scratch_root=SCRATCH_ROOT,
    ) + EQUIL_BODY
    write(f"{seed_dir}/equil_{seed}.sh", eq)

    prod_extra = "#SBATCH --ntasks=48\n#SBATCH --cpus-per-task=1\n"
    prod = COMMON_HEADER.format(
        job_name=f"{seed}_prod", tag=f"prod_{seed}", time="22:00:00",
        extra_sbatch=prod_extra, proj_dir=PROJ_DIR, seed=seed, scratch_root=SCRATCH_ROOT,
    ) + PROD_BODY.format(seed=seed)
    write(f"{seed_dir}/prod_md_{seed}.sh", prod)

    resume = RESUME_HEADER.format(job_name=f"{seed}_prod_xtnd", tag=f"resume_{seed}", seed=seed, scratch_root=SCRATCH_ROOT)
    write(f"{seed_dir}/resume_prod_{seed}.sh", resume)

    print(f"wrote scripts for {seed}")

print("done")
