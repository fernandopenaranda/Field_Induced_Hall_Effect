
#a = 1, L in Å, E in eV# Δ = sublattice mass 
# Model of graphene with three neighbours with opposite signs and sublatice masses
# breaking Mx and My and enabling a ferroaxial response
function ferroaxial_ham2d(k, μ = 0, Δ = 1, t = 1, tp = 0.5)
    a = 1.0
    a1 = a * [3/2,  sqrt(3)/2]
    a2 = a * [3/2, -sqrt(3)/2]
    δ1 = a * [0.0, 1/sqrt(3)]
    δ2 = a * [-1/2, -1/(2sqrt(3))]
    δ3 = a * [ 1/2, -1/(2sqrt(3))]
    ds = [δ1, δ2, δ3]
    fk = f(k, ds)
    gk = g(k, ds)
    return [Δ-μ  t*fk+tp*gk; conj(t*fk+tp*gk) -Δ-μ]
end


""" unbounded in the z direction"""
function ferroaxial_ham3d(k, μ = 0, Δ = 1, t = 1, tp = 0.5)
    a = 1.0
    a1 = a * [3/2,  sqrt(3)/2, 0]
    a2 = a * [3/2, -sqrt(3)/2, 0]
    a3 = a * [0, 0, 1]
    δ1 = a * [0.0, 1/sqrt(3), 0]
    δ2 = a * [-1/2, -1/(2sqrt(3)), 0]
    δ3 = a * [ 1/2, -1/(2sqrt(3)), 0]
    ds = [δ1, δ2, δ3]
    fk = f(k, ds)
    gk = g(k, ds)
    fzk = fz(k, a3)
    return [Δ-μ  t*(fzk+fk)+tp*(gk); conj(t*(fzk+fk)+tp*gk) -Δ-μ]
end

fz(k, a3) = cis(dot(k, a3))
f(k, deltas) = sum(exp(im * dot(k, δ)) for δ in deltas)
g(k, ds) = g(k, ds[1], ds[2], ds[3])
g(k,d1,d2,d3) = cis(dot(k,2d1-d2)) + cis(dot(k,2d2-d3)) + cis(dot(k,2d3-d1)) - cis(dot(k,2d1-d3)) - cis(dot(k,2d2-d1)) -cis(dot(k,2d3-d2 ))

function d_ferroaxial_ham2d(t,tp,Δ,q, dir)
    a = 1.0
    a1 = a * [3/2,  sqrt(3)/2]
    a2 = a * [3/2, -sqrt(3)/2]
    δ1 = a * [0.0, 1/sqrt(3)]
    δ2 = a * [-1/2, -1/(2sqrt(3))]
    δ3 = a * [ 1/2, -1/(2sqrt(3))]
    ds = [δ1, δ2, δ3]
    dfk = df(q, ds, dir)
    dgk = dg(q, ds, dir)
    return [0 t*dfk + tp*dgk;conj(t*dfk + tp*dgk) 0]
end

function d_ferroaxial_ham3d(t,tp,Δ,q, dir)
    a = 1.0
    a1 = a * [3/2,  sqrt(3)/2, 0]
    a2 = a * [3/2, -sqrt(3)/2, 0]
    a3 = a * [0, 0, 1]
    δ1 = a * [0.0, 1/sqrt(3), 0]
    δ2 = a * [-1/2, -1/(2sqrt(3)), 0]
    δ3 = a * [ 1/2, -1/(2sqrt(3)), 0]
    ds = [δ1, δ2, δ3]
    dfk = df(q, ds, dir)
    dfzk = dfz(q, a3, dir)
    dgk = dg(q, ds, dir)
    return [0 t*(dfk+dfzk) + tp*dgk; conj(t*(dfk+dfzk) + tp*dgk) 0]
end

dfz(k, a3 , dir) = ifelse(dir != :z, 0, im*a3[3]*cis(dot(k, a3)))
df(k, deltas, a) = sum(im*δ[symb_to_ind(a)]*exp(im * dot(k, δ)) for δ in deltas)
dg(k, ds,a) = dg(k, ds[1], ds[2], ds[3], a)
dg(k,d1,d2,d3,a) = im*(2d1-d2)[symb_to_ind(a)]*cis(dot(k,2d1-d2)) + 
                    im*(2d2-d3)[symb_to_ind(a)]*cis(dot(k,2d2-d3)) + 
                    im*(2d3-d1)[symb_to_ind(a)]*cis(dot(k,2d3-d1)) -
                    im*(2d1-d3)[symb_to_ind(a)]*cis(dot(k,2d1-d3)) -
                    im*(2d2-d1)[symb_to_ind(a)]*cis(dot(k,2d2-d1)) -
                    im*(2d3-d2)[symb_to_ind(a)]*cis(dot(k,2d3-d2))

#second derivative
function d2_ferroaxial_ham2d(t, tp, q, dir1, dir2)
    a = 1.0
    δ1 = a * [0.0, 1/sqrt(3)]
    δ2 = a * [-1/2, -1/(2sqrt(3))]
    δ3 = a * [ 1/2, -1/(2sqrt(3))]
    ds = [δ1, δ2, δ3]
    d2fk = d2f(q, ds, dir1, dir2)
    d2gk = d2g(q, ds, dir1, dir2)
    off = t*d2fk + tp*d2gk
    [0 off; conj(off) 0]
end

function d2_ferroaxial_ham3d(t, tp, q, dir1, dir2)
    a = 1.0
    δ1 = a * [0.0, 1/sqrt(3),0.0]
    δ2 = a * [-1/2, -1/(2sqrt(3)),0.0]
    δ3 = a * [ 1/2, -1/(2sqrt(3)),0.0]
    a3 = a * [0, 0, 1]
    ds = [δ1, δ2, δ3]
    d2fk = d2f(q, ds, dir1, dir2)
    d2fzk = d2fz(q, a3, dir1, dir2)
    d2gk = d2g(q, ds, dir1, dir2)
    off = t*(d2fk+d2fzk) + tp*d2gk
    [0 off; conj(off) 0]
end

function d2fz(k, a3, dir1, dir2)
    i = symb_to_ind(dir1)
    j = symb_to_ind(dir2)
    if dir1 == dir2 && dir1 == :z
        return -a3[j]^2*cis(dot(k,a3))
    else 
        return 0
    end
end

function d2f(k, deltas, dir1, dir2)
    i = symb_to_ind(dir1)
    j = symb_to_ind(dir2)
    sum(-δ[i]*δ[j]*cis(dot(k,δ)) for δ in deltas)
end

function d2g(k, d1, d2, d3, dir1, dir2)
    i = symb_to_ind(dir1)
    j = symb_to_ind(dir2)
    r1 = 2*d1 - d2
    r2 = 2*d2 - d3
    r3 = 2*d3 - d1
    r4 = 2*d1 - d3
    r5 = 2*d2 - d1
    r6 = 2*d3 - d2
    return -(r1[i]*r1[j]*cis(dot(k,r1)) +
        r2[i]*r2[j]*cis(dot(k,r2)) +
        r3[i]*r3[j]*cis(dot(k,r3)) -
        r4[i]*r4[j]*cis(dot(k,r4)) -
        r5[i]*r5[j]*cis(dot(k,r5)) -
        r6[i]*r6[j]*cis(dot(k,r6)))
end

d2g(k, ds, dir1, dir2) = d2g(k, ds[1], ds[2], ds[3], dir1, dir2)


