# UNITS E -> eV, time -> s, length -> Å, temperature -> K
""" model for TaRh2B2 (3D)
three bands three hoppings """
hamiltonian_TaRh2B2(μ, ts::Array, q, a =1 ,c = 1) = hamiltonian_TaRh2B2(μ, ts[1],ts[2],ts[3], q, a, c)
function hamiltonian_TaRh2B2(μ, t1::Number,t2,t3,q, a =1 ,c = 1)
    ts = [t1,t2,t3]
    δs = site_distances(a, c)
    phase(i) = 1im*(δs[i]'* q)
    st1 = sum([ts[i] * exp(phase(i)) for i in 1:3])
    st2 = sum([ts[mod1(i+1,3)] * exp(-phase(i)) for i in 1:3])
    st3 = sum([ts[mod1(i+2,3)] * exp(phase(i)) for i in 1:3])
    return [-μ st1 st2; conj(st1) -μ st3; conj(st2) conj(st3) -μ]
end

""" first derivative ∂dir H """
d_hamiltonian_TaRh2B2(ts::Array,q, dir, a =1 ,c = 1) = d_hamiltonian_TaRh2B2(ts[1],ts[2],ts[3],q, dir, a, c)
function d_hamiltonian_TaRh2B2(t1::Number,t2,t3,q, dir, a =1 ,c = 1)
    ts = [t1,t2,t3]
    δs = site_distances(a, c)
    phase(i) = 1im*(δs[i]'* q)
    st1 = sum([ts[i]* δs[i][symb_to_ind(dir)] * exp(phase(i)) for i in 1:3])
    st2 = sum([ts[mod1(i+1,3)]* δs[i][symb_to_ind(dir)] * exp(-phase(i)) for i in 1:3])
    st3 = sum([ts[mod1(i+2,3)]* δs[i][symb_to_ind(dir)] * exp(phase(i)) for i in 1:3])
    return [0 1im * st1 -1im * st2; -1im * conj(st1) 0 1im * st3; 1im * conj(st2) -1im * conj(st3) 0]
end

""" second derivative ∂dirout ∂dirin H """
d_d_hamiltonian_TaRh2B2(ts::Array,q, dir_out, dir_in, a =1 ,c = 1) = 
    d_d_hamiltonian_TaRh2B2(ts[1],ts[2],ts[3],q, dir_out, dir_in, a, c)
function d_d_hamiltonian_TaRh2B2(t1::Number,t2,t3, q, dir_out, dir_in, a =1 ,c = 1)
    ts = [t1,t2,t3]
    δs = site_distances(a, c)
    phase(i) = 1im*(δs[i]'* q)
    st1 = sum([ts[i]* δs[i][symb_to_ind(dir_out)] * δs[i][symb_to_ind(dir_in)] * exp(phase(i)) for i in 1:3])
    st2 = sum([ts[mod1(i+1,3)]* δs[i][symb_to_ind(dir_out)] * δs[i][symb_to_ind(dir_in)] * exp(-phase(i)) for i in 1:3])
    st3 = sum([ts[mod1(i+2,3)]* δs[i][symb_to_ind(dir_out)] * δs[i][symb_to_ind(dir_in)] * exp(phase(i)) for i in 1:3])
    return -1 .* [0 st1 st2; conj(st1) 0  st3; conj(st2) conj(st3) 0]
end

site_distances(a = 1, c = 1) = [[0, a/√3, c/3], [-a/2, -a/(2√3), c/3], [a/2, -a/(2√3), c/3]]

function symb_to_ind(dir)
    if dir == :x
        1
    elseif dir == :y
        2
    elseif dir == :z 
        3
    else 
        0
        throw(ArgumentError("dir not in [:x,:y,:z]"))
    end
end


function bands(i,q,μ,t1,t2,t3)
    # mat = hamiltonian_TaRh2B2(a,c,t1,t2,t3,q)
    mat = hamiltonian_TaRh2B2(μ,t1,t2,t3,q)
    ϵs, ψs = eigen(mat)   
    return ϵs[i]
end

bands_TaRh2B2(μ, ts) = bands_TaRh2B2(μ, ts[1],ts[2],ts[3])
function bands_TaRh2B2(μ, t1,t2,t3; a = 1, c = 1)
    Rs = (a .* [1.0, 0, 0], a .* [-1/2, √3/2, 0], [0, 0.0, c])  # Lattice vectors
    sgnum = 144                                         # point group (144 = p3_1)
    N = 200    #kpoints
    fig = Figure(); ax = Axis(fig[1,1], ylabel = "E [eV]", title = "t1 = $(t1) eV, t2 = $(t2) eV, t3 = $(t3) eV")
    for i in 1:3
        dispersion(q) = bands(i,q,μ, t1,t2,t3)
        BZpaths.plot_observable_in_kpath!(ax, dispersion, Rs, sgnum, N)#, high_sym_line high_sym_line = [:Γ, :M, :K, :Γ, :A, :L, :H, :A] )
    end
    fig
end