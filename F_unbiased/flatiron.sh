#!/bin/bash
#SBATCH --job-name="wt-act"
#SBATCH --output=stdout_file
#SBATCH --error=stderr_file
#SBATCH --partition=ccb
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
##SBATCH --mem=128GB
#SBATCH --export=ALL
#SBATCH -t 7-00:00:00



#####################################################
echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURMTMPDIR="$SLURMTMPDIR

echo "working directory = "$SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

#####################################################
# configure environment
module load gcc
module load openmpi
#module load mpi/gcc_openmpi
gmxdir="/mnt/home/kpatil/src_gromacs_2018.6/bin"
source $gmxdir/GMXRC
gmx_mpi grompp -f md.mdp -c md.part0001.gro -r md.part0001.gro -p system.top -o md.part0002.tpr
gmx_mpi mdrun -v -deffnm md.part0002
