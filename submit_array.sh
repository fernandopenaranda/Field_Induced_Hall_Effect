#!/bin/bash
#SBATCH --job-name=eval_subcubes
#SBATCH --array=1-64
#SBATCH --output=job_%A_%a.out
#SBATCH --error=slurm-%A.%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=23:00:00
JULIA_SCRIPT="/Users/fernandopenaranda/.julia/packages/Optics_in_the_length_gauge/47RmS/src/cluster/cluster_evaluate_subcube.jl"
dirj=$1
dirE=$2
dirB=$3
T=$4
evals=$5
omega_switch=$6
ps_switch=$7
qm_switch=$8
fermi_surface=$9
epsilon=${10}
which_mm=${11}
integration_method=${12}
t1=${13}
t2=${14}
t3=${15}
mumin=${16}
mumax=${17}
mupoints=${18}

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Array length: $SLURM_ARRAY_TASK_MAX"

CMD="/scratch/ferpe/julia-1.9.4/bin/julia --compiled-modules=no $JULIA_SCRIPT \
$SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX $SLURM_ARRAY_JOB_ID \
$dirj $dirE $dirB $T $evals \
$omega_switch $ps_switch $qm_switch $fermi_surface \
$epsilon $which_mm $integration_method \
$t1 $t2 $t3 \
$mumin $mumax $mupoints"

if [ "${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
    printf "JOBID=%s CMD=%s\n" "$SLURM_ARRAY_JOB_ID" "$CMD" >> "julia_command_interpolations.txt"
fi

exec $CMD

