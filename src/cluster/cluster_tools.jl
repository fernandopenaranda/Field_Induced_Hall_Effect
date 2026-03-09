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
