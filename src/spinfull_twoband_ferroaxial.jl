

function ferroaxial_ham3d_spinfull(k,  μ = 0, Δ = 1, t = 1, tp = 0.5, tp_z = 0.5, tc = 0.5, tλ = 1, tsx = 1, tsy = 1)
    a = 1.0
    a1, a2, a3, δ1, δ2, δ3 = lattice_vectors(a)
    ds = [ δ1, δ2, δ3 ]

    sx = [0 1; 1 0]
    sy = [0 -im; im 0]
    h = ferroaxial_ham3d(k, μ, Δ, t, tp, tp_z, tc)

    return kron(I(2), h) .+ kane_mele_soc(k, tλ, ds)  .+ tsx * kron(I(2), sx) + tsy * kron(I(2), sy)
end

function d_ferroaxial_ham3d_spinfull(t,tp,Δ,tpz,tc, tλ, q, dir)
    a = 1.0
    a1, a2, a3, δ1, δ2, δ3 = lattice_vectors(a)
    ds = [ δ1, δ2, δ3 ]
    dh = d_ferroaxial_ham3d(t,tp,Δ,tpz,tc, q, dir) 
    return kron([1 0; 0 1], dh) + d_kane_mele_soc(q, tλ, ds, dir) 
end

function d2_ferroaxial_ham3d_spinfull(t,tp,tpz,tc, tλ, q, dir1, dir2)
    a = 1.0
    a1, a2, a3, δ1, δ2, δ3 = lattice_vectors(a)
    ds = [ δ1, δ2, δ3 ]
    dh = d2_ferroaxial_ham3d(t, tp, tpz,tc, q, dir1, dir2) 
    return kron([1 0; 0 1], dh) + d2_kane_mele_soc(q, tλ, ds, dir1, dir2) 
end

kane_mele_soc(k, tλ, ds) = 
    sigmazsz_directproduct(kane_mele_soc(k, tλ, ds[1], ds[2], ds[3]))

d_kane_mele_soc(k, tλ, ds, dir) = 
    sigmazsz_directproduct(d_kane_mele_soc(k, tλ, ds[1], ds[2], ds[3], dir))

d2_kane_mele_soc(k, tλ, ds, dir1, dir2) = 
    sigmazsz_directproduct(d2_kane_mele_soc(k, tλ, ds[1], ds[2], ds[3], dir1, dir2))

kane_mele_soc(k, tλ, d1, d2, d3) = 
    2*tλ*sum([sin(dot(k, δ)) for δ in [d1-d2, d2-d3, d3-d1]])
    
d_kane_mele_soc(k, tλ, d1, d2, d3, dir) = 
    2*tλ*sum([δ[symb_to_ind(dir)]*cos(dot(k, δ)) for δ in [d1-d2, d2-d3, d3-d1]])
    
d2_kane_mele_soc(k, tλ, d1, d2, d3, dir1, dir2) = 
    2*tλ*sum([-δ[symb_to_ind(dir1)]*δ[symb_to_ind(dir2)]*sin(dot(k, δ)) for δ in [d1-d2, d2-d3, d3-d1]])
    
sigmazsz_directproduct(soc) = soc * kron([1 0; 0 -1], [1 0; 0 -1])

# function d_ferroaxial_ham3d_spinfull()

# end

# function d2_ferroaxial_ham3d_spinfull()
    
# end

