#!/bin/bash
#SBATCH --account=mXXXX
#SBATCH --job-name=VASP_CBaMnBr4_pol2
#SBATCH --qos=regular
#SBATCH --constraint=gpu
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --gpus=4
#SBATCH --time=1:00:00
#SBATCH --array=0-10

# Perlumtter has 4 GPUs per node

# OpenMP Options
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Patch unreliable Perlmutter interconnect settings
export FI_CXI_DEFAULT_TX_SIZE=16384
export FI_CXI_RDZV_THRESHOLD=65536

# VASP Module
# Use ISYM=0 for magnetic systems if using VASP <= 6.3.2 on Perlmutter
module load vasp/6.4.3-gpu

# tee VASP stdout to file but still track if VASP failed.
set -o pipefail

echo "Doing step $SLURM_ARRAY_TASK_ID"
cd step_$SLURM_ARRAY_TASK_ID || exit 1

echo "Doing magnetization calculation at `date`"

cd 1.mag || exit 1

# Print the info from INCAR to stdout, minus comments and empty lines.
# Keeps a record if tweaking INCAR to coax convergence.
sed '/^#/d; /^$/d' INCAR

srun --nodes=$SLURM_JOB_NUM_NODES --gpus=$SLURM_GPUS --ntasks=$SLURM_GPUS --cpus-per-task=32 --ntasks-per-node=4 --gpus-per-task=1 --cpu-bind=cores --gpu-bind=none --unbuffered vasp_std | tee -a stdout.log

if (($?)); then
	echo "Error in VASP at `date`"
	exit 1
fi

echo "Zipping files at `date`"

sed -n 's/^ *TITEL *=//p' POTCAR > POTCAR.spec
rm CHG PCDAT REPORT XDATCAR CONTCAR
zip -T CHGCAR.zip CHGCAR
zip -T Main.zip INCAR KPOINTS POSCAR OUTCAR
zip -mT Main.zip POTCAR.spec
zip -mT Data.zip DOSCAR EIGENVAL IBZKPT OSZICAR OUTCAR
gzip -f vasprun.xml
zip --junk-paths Main.zip $0
zip -T Main.zip stdout.log

cp WAVECAR ../2.pol

cd .. || exit 1

echo "Doing polarization calculation at `date`"

cd 2.pol || exit 1

# Print the info from INCAR to stdout, minus comments and empty lines.
# Keeps a record if tweaking INCAR to coax convergence.
sed '/^#/d; /^$/d' INCAR

srun --nodes=$SLURM_JOB_NUM_NODES --gpus=$SLURM_GPUS --ntasks=$SLURM_GPUS --cpus-per-task=32 --ntasks-per-node=4 --gpus-per-task=1 --cpu-bind=cores --gpu-bind=none --unbuffered vasp_std | tee -a stdout.log

if (($?)); then
	echo "Error in VASP at `date`"
	exit 1
fi

echo "Zipping files at `date`"

sed -n 's/^ *TITEL *=//p' POTCAR > POTCAR.spec
rm CHG PCDAT REPORT XDATCAR CONTCAR
zip -T CHGCAR.zip CHGCAR
zip -T Main.zip INCAR KPOINTS POSCAR OUTCAR
zip -mT Main.zip POTCAR.spec
zip -mT Data.zip DOSCAR EIGENVAL IBZKPT OSZICAR
gzip -f OUTCAR
gzip -f vasprun.xml
zip --junk-paths Main.zip $0
zip -T Main.zip stdout.log

echo "Finished VASP at `date`"
