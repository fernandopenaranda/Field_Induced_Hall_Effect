DEPOT_PATH[:] .= ["/scratch/ferpe/julia_depot"]
using Optics_in_the_length_gauge, CSV, JLD2, Serialization, Field_Induced_Hall_Effect
job_id = parse(Int, ARGS[1]) # Get the task ID from SLURM
jobs_num = parse(Int, ARGS[2])
PID = ARGS[3]
dirj =  Symbol(ARGS[4])
dirE =  Symbol(ARGS[5])
dirB =  Symbol(ARGS[6])
T = parse(Float64, ARGS[7]) 
evals = parse(Int, ARGS[8])
omega_switch = parse(Bool, ARGS[9])
ps_switch = parse(Bool, ARGS[10])
qm_switch = parse(Bool, ARGS[11])
fermi_surface = parse(Bool, ARGS[12])
epsilon = parse(Float64, ARGS[13])
which_mm = Symbol(ARGS[14])
integration_method = Symbol(ARGS[15])
t1 = parse(Float64, ARGS[16])
t2 = parse(Float64, ARGS[17])
t3 = parse(Float64, ARGS[18])
mumin = parse(Float64, ARGS[19])
mumax = parse(Float64, ARGS[20])
mupoints= parse(Int, ARGS[21])

print("Starting...")
if mupoints == 0 || mupoints == 1
    muvec = [mumin]
else
    muvec = collect(range(mumin, mumax, length=mupoints))
end
a = 1
c = 1

subcubes_file = pwd() * "/subcubes.jls"
subcubes = deserialize(subcubes_file)
my_subcube = subcubes[job_id]
botbounds = [my_subcube[1][1], my_subcube[2][1], my_subcube[3][1]]
topbounds = [my_subcube[1][2], my_subcube[2][2], my_subcube[3][2]]

p = modelpresets(;t1 =t1 ,t2 =t2, t3=t3, a= a, c=c)
keyws = (
    dirj=dirj, 
    dirE=dirE, 
    dirB=dirB, 
    T = T, 
    evals = evals, 
    Ω_MM_switch = omega_switch, 
    PS_switch = ps_switch, 
    QM_switch = qm_switch, 
    fermi_surface = fermi_surface,
    ϵ = epsilon, 
    which_mm = which_mm, 
    integration_method = integration_method, 
    botbounds = botbounds, 
    topbounds = topbounds)

sijks = Float64[]
for mu in muvec  
    sijk_pres = sigma_abc_TaRh2B2_wrapper(TaRh2B2Params(p, μ =mu, t1 =t1 ,t2 =t2, t3=t3); keyws...) #computing struct for each mu
    push!(sijks,quantum_contribution(sijk_pres)) # compute
end
sijk_pres = sigma_abc_TaRh2B2_wrapper(TaRh2B2Params(p, t1 =t1 ,t2 =t2, t3=t3); keyws...) 

#_________________________________________________________________________________________
#store
#_________________________________________________________________________________________
data_folder = pwd() * "/Data/" * string(PID) * "/" * string(job_id)
mkpath(data_folder)
@save data_folder * "/presets.jld" sijk_pres
@save data_folder * "/calculation.jld" muvec sijks

str = pwd() * "/slurm-" * string(PID) * "." * string(job_id)
isfile(str * ".out") && mv(str * ".out", data_folder * "/output.out")
isfile(str * ".err") && mv(str * ".err", data_folder * "/error.err")