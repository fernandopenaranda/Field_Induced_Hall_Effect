
function spinfull_ferroaxial(k,  μ = 0, Δ = 1, t = 1, tp = 0.5, tp_z = 0.5, tc = 0.5, tλ = 1)
    h = ferroaxial_ham3d(k, μ, Δ, t, tp, tp_z, tc)
    return kron(h, [1 0; 0 1]) .+ kane_mele_soc(k, tλ, ds) 
end

"""
    kane_mele_term(k, λ)

Kane-Mele spin-orbit coupling: next-nearest neighbor term with spin structure.
Structure: I₃ ⊗ [i*λ*g(k)*σ_z]

where g(k) = (nnn_cw - nnn_ccw) with alternating signs on clockwise/counterclockwise loops.

Returns 6×6 matrix acting in full spin-orbital space.
"""

function kane_mele_soc(k, tλ, ds) 
    soc = kane_mele_soc(k, tλ, ds[1], ds[2], ds[3])
    sublatsoc = kron([1 0; 0 -1], soc)
    return kron([1 0; 0 -1], sublatsoc)
end

kane_mele_soc(k, λ, d1, d2, d3) = 
     2*tλ*sum([sin(dot(k, δ)) for δ in [d1-d2, d2-d3, d3-d1]])

