
"""
this function determines BZ cubic partitions by
splitting the BZ parametrized as a cube in u1, u2, u3, 
with BZ -> k = u1 *b1 + u2*b2 + u3*b3
"""
function split_cube(divisions::Int)
    # Determine number of divisions along each axis
    if divisions <= 0
        N = 1
    else
        N = floor(Int, divisions^(1/3))  
    end
    x = range(-0.5, 0.5, length=N+1)
    subcubes = []
    x = range(-0.5, 0.5, length=N+1)  # +1 because these are boundaries
    subcubes = []
    for i in 1:N
        for j in 1:N
            for k in 1:N
                x0, x1 = x[i], x[i+1]
                y0, y1 = x[j], x[j+1]
                z0, z1 = x[k], x[k+1]
                push!(subcubes, ((x0, x1), (y0, y1), (z0, z1)))
            end
        end
    end
    return subcubes
end

""" creates the .sh file for slurm computations. The number of jobs
it runs coincide with the number of divisions of the intergration BZ cube
This .sh file will call to `julia_script` with the preset directives.
"""
function create_bashfile(divisions)
    julia_script = dirname(pathof(Field_Induced_Hall_Effect)) * "/cluster/cluster_evaluate_subcube.jl"  # script each array job runs
    subcubes_file = pwd() * "/subcubes.jls"        # File that stores the domain integration partitions
    output_script = pwd() * "/submit_array.sh"     # SLURM submission script

    subcubes = split_cube(divisions) # Compute partitions
    total_subcubes = length(subcubes)
    println("Total subcubes (SLURM array size): ", subcubes)
    
    serialize(subcubes_file, subcubes) # Save partitions
    println("Saved ", length(subcubes), " subcubes to $subcubes_file")

    # ===== Write SLURM script =====
    # open(output_script, "w") do io
    #     println(io, """#!/bin/bash
    # #SBATCH --job-name=eval_subcubes
    # #SBATCH --array=1-$total_subcubes
    # #SBATCH --output=job_%A_%a.out
    # #SBATCH --cpus-per-task=1
    # #SBATCH --mem-per-cpu=2G
    # #SBATCH --time=02:00:00 
    # #SBATCH --error=slurm-%A.%a.err 
    # # Run the Julia script for this subcube
    # CMD="/scratch/ferpe/julia-1.9.4/bin/julia --compiled-modules=no $julia_script \$SLURM_ARRAY_TASK_ID \$SLURM_ARRAY_TASK_MAX \$SLURM_ARRAY_JOB_ID"
    
    # if [ "\${SLURM_ARRAY_TASK_ID:-0}" -eq 1 ]; then
    #     printf "JOBID=%s CMD=%s\n" "\$SLURM_ARRAY_JOB_ID" "\$CMD" >> "julia_command_interpolations.txt"
    # fi
    # exec \$CMD
    # """)
    open(output_script, "w") do io
        println(io, "#!/bin/bash")
        println(io, "#SBATCH --job-name=eval_subcubes")
        println(io, "#SBATCH --array=1-$total_subcubes")
        println(io, "#SBATCH --output=slurm-%A_%a.out")
        println(io, "#SBATCH --error=slurm-%A.%a.err")
        println(io, "#SBATCH --cpus-per-task=1")
        println(io, "#SBATCH --mem-per-cpu=8G")
        println(io, "#SBATCH --time=23:00:00")
        println(io, "JULIA_SCRIPT=\"$julia_script\"")

        println(io, raw"""
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
        printf "JOBID=%s CMD=%s\n" "$SLURM_ARRAY_JOB_ID" "$CMD" >> "julia_commands.txt"
    fi
    
    exec $CMD
    """)
    end
    println("SLURM submission script written to $output_script")
end