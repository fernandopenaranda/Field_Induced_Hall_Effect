
#=
here we build the appropiate structs in Optics_in_the_length_gauge to pass the model
=#
modelpresets(; μ =0, t1 =1 ,t2 =1, t3=1, a= 1, c=1) = TaRh2B2Params(μ, t1,t2,t3,a,c)
@with_kw struct TaRh2B2Params
    μ::Float64
    t1::Float64
    t2::Float64
    t3::Float64
    a::Float64
    c::Float64
end

"""
structure wrapper botbounds and topbounds in [-0.5,0.5]
"""
function sigma_abc_TaRh2B2_wrapper(p; dirj=:x, dirE=:y, dirB=:z, T = 1, τ = 200, evals = 100, 
        omega_MM_switch = true, PS_switch = true, QM_switch = true, fermi_surface = false,
        which_mm = :orbital, epsilon = 1e-5, integration_method = :hcubature,
        botbounds = [-0.5,-0.5,-0.5], topbounds = [0.5,0.5,0.5])
    a = p.a
    c = p.c
    unit_convention_two_packages_t = 1e-15
    τ *= unit_convention_two_packages_t
    #hamiltonians and derivatives
    h(q) = hamiltonian_TaRh2B2(p.μ,p.t1,p.t2,p.t3,q,p.a,p.c)
    dhx(q) = d_hamiltonian_TaRh2B2(p.t1,p.t2,p.t3, q, :x, p.a, p.c)
    dhy(q) = d_hamiltonian_TaRh2B2(p.t1,p.t2,p.t3, q, :y, p.a, p.c)
    dhz(q) = d_hamiltonian_TaRh2B2(p.t1,p.t2,p.t3, q, :z, p.a, p.c)
    dh(q) = [dhx(q), dhy(q), dhz(q)]
    didjh(q,i,j) = d_d_hamiltonian_TaRh2B2(p.t1,p.t2,p.t3,q,i,j,p.a,p.c) 
    ddh(q) = [[didjh(q,:x,:x), didjh(q,:x,:y), didjh(q,:x,:z)],
            [didjh(q,:y,:x), didjh(q,:y,:y), didjh(q,:y,:z)], 
            [didjh(q,:z,:x), didjh(q,:z,:y), didjh(q,:z,:z)]]

    # integral presets
    Rs = (a .* [1.0, 0, 0], a .* [-1/2, √3/2, 0], [0, 0.0, c])
    Gs = dualbasis(Rs)
    computation = Transport_computation_3d_presets(botbounds,topbounds, evals, integration_method)
    return Quantum_correction_σijk_antisym(a, dirj, dirE, dirB, h, dh, ddh, Gs, τ, T, computation,
         which_mm, omega_MM_switch, PS_switch, QM_switch, fermi_surface, epsilon)
    # old bounds _________________________________________________________________________________________
    # cell = wignerseitz(Gs) 
    # v = cartesianize!(cell).verts #vertices of the BZ
    # largest_component(v, i) = maximum(x -> x[i], v)
    # smallest_component(v,i) = minimum(x -> x[i], v)
    # lmin = smallest_component(v,1)
    # lmax = largest_component(v,1)
    # lx = lmax-lmin
    # botbounds =1 .* [a* lmin, smallest_component(v,2), c*smallest_component(v,3)] 
    # topbounds =1 .*[a*(lmin+3/4*lx), largest_component(v,2), c*largest_component(v,3)]
    # a*(lmin+3/4*lx) this factor ensures that the square integral is limited to the BZ (without FS duplicities)
    # need a prefactor due to periodicity
end

function qc_sweep(μlist, p, cp_pres)
    s = []
    for μ in μlist
        h(q) = hamiltonian_TaRh2B2(μ,p.t1,p.t2,p.t3,q,p.a,p.c)
        new_comp_pres = Quantum_correction_σijk_antisym(cp_pres.a0, cp_pres.dirJ, 
            cp_pres.dirE, cp_pres.dirB, h, cp_pres.nabla_h, cp_pres.nabla_nabla_h, 
            cp_pres.gs, cp_pres.τ, cp_pres.T, cp_pres.computation, cp_pres.which_mm)
        append!(s, quantum_contribution(new_comp_pres))
    end
    return μlist ,s
end

