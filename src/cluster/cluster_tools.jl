using JLD2
using FilePathsBase, Field_Induced_Hall_Effect
function slurm_conductivities(; dirj = :x, dirE = :y, dirB = :x, T = 5, evals = 1000, omega_switch = true, 
    ps_switch = true, qm_switch = true, fermi_surface = true, epsilon = 1e-5, which_mm = :orbital, integration_method = :montecarlo,
    t1 = 1, t2 = 0.5, t3 = 0.1, mumin = -3, mumax = 3, mupoints = 1, dryrun=false)
    script = pwd() * "/submit_array.sh"
    cmd = `sbatch --wait $script \
        $dirj $dirE $dirB $T $evals \
        $omega_switch $ps_switch $qm_switch $fermi_surface \
        $epsilon $which_mm $integration_method \
        $t1 $t2 $t3 \
        $mumin $mumax $mupoints`
    dryrun && return cmd
    run(cmd)
end
#_________________________________________________________________________________________
# reshape data after calculation

#using LMC, JLD2 # needs to load LMC structures to approately read the preset files

""" PID folder finder: e.g. PID = 0000001, 
    `find_folder(PID)`
 find the absolute path to PID folder regardless the data
 store phase diagrams, drude, qah, LMCS responses... """
find_folder(PID) = find_folder(string(PID))
function find_folder(target::AbstractString)
    root = pwd()
    for (dirpath, dirnames, _) in walkdir(root)
        if target in dirnames
            return joinpath(dirpath, target)
        end
    end
    return nothing
end



postprocessing(PID::Number) = postprocessing(string(PID))
function postprocessing(PID::String)
    destination =  homedir() * "/Projects/Field_Induced_Hall_Effect/" * PID 
    calcfile = destination * "/merged_calculation.jld"
    presfile = destination * "/merged_presets.jld"
    pid_folder =  find_folder(PID)
    subfolders = filter(isdir, joinpath.(pid_folder, readdir(pid_folder)))
    first_vector = nothing
    summed_vector = nothing
    for folder in subfolders
        file = joinpath(folder, "calculation.jld")
        if isfile(file)
            data = load(file)
            v1 = data["muvec"]
            v2 = data["sijks"]
            if first_vector === nothing
                first_vector = v1
                summed_vector = copy(v2)
            else
                summed_vector .+= v2
            end
        end
    end

    save(calcfile,
        "muvec" => first_vector,
        "sijks" => summed_vector
    )

    cp(pid_folder * "/1/presets.jld", presfile)

    println("Saved summed result to $destination")

end