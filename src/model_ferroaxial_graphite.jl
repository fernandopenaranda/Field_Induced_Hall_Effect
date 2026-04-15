
function lattice_vectors(a)
    a1 = a/2 * [3, √3, 0] # Γ - M line along the x axis 
    a2 = a/2 * [3, -√3, 0]
    a3 = a * [0, 0, 1]
    δ1 = a/2 * [1, √3, 0]
    δ2 = a/2 * [1, -√3, 0]
    δ3 = a * [-1, 0, 0]
    return a1, a2, a3, δ1, δ2, δ3
end

""" unbounded in the z direction. With opposite third neighbour hoppings (opposite) for  the AA and BB hoppings in z
tp, tp_z is a first order"""
function ferroaxial_ham3d(k, μ = 0, Δ = 1, t = 1, tp = 0.5, tp_z = 0.5, tc = 0.5)
    a = 1.0
    a1, a2, a3, δ1, δ2, δ3 = lattice_vectors(a)
    ds = [δ1, δ2, δ3]
    fk = f(k, ds) # firs neighbour hoppings
    gk = g(k, ds) # third neightbour hoppings (opposite signs see matricial form)
    mk = m(k, ds, a3) # z- AB hoppings
    nnk = nn(k, ds) # second neighbour hoppings 
    # println(gk)
    # throw(ArgumentError(""))
    return [Δ-μ+tc*nnk  t*fk+tp*gk+tp_z*mk; conj(t*fk+tp*gk+tp_z*mk) -Δ-μ + tc*nnk ]
end
nn(k, deltas)= nn(k, deltas[1], deltas[2], deltas[3])
nn(k, d1,d2,d3) = sum([2 * cos(dot(k, δ)) for δ in [d1-d2,d2-d3,d3-d1]])

# nn(k, d1, d2, d3) = cis(dot(k,d1-d2)) + cis(dot(k,d2-d3)) + cis(dot(k,d3-d1)) + cis(dot(k,-d1+d2)) + cis(dot(k,-d2+d3)) + cis(dot(k,-d3+d1))
fz(k, a3) = cis(dot(k, a3))
f(k, deltas) = sum(exp(im * dot(k, δ)) for δ in deltas)
m(k, ds, a3) = sum(exp(im * dot(k, δ+a3)) for δ in ds) + sum(exp(im * dot(k, δ-a3)) for δ in ds) # same signs for AB and BA
g(k, ds) = g(k, ds[1], ds[2], ds[3])
g(k,d1,d2,d3) = cis(dot(k,2d1-d2)) + cis(dot(k,2d2-d3)) + cis(dot(k,2d3-d1)) - cis(dot(k,2d1-d3)) - cis(dot(k,2d2-d1)) -cis(dot(k,2d3-d2 ))

function d_ferroaxial_ham3d(t,tp,Δ,tpz,tc, q, dir) 
    a = 1.0
    a1, a2, a3, δ1, δ2, δ3 = lattice_vectors(a)
    ds = [δ1, δ2, δ3]
    dfk = df(q, ds, dir)
    dmk = dm(q, ds, a3, dir)
    dgk = dg(q, ds, dir)
    dnnk = dnn(q, ds, dir)
    return [tc*dnnk t*dfk+tp*dgk+tpz*dmk; conj(t*dfk+tp*dgk+tpz*dmk) tc*dnnk]
end

dm(k, deltas, a3, a) = sum(im*(δ+a3)[symb_to_ind(a)]*exp(im * dot(k, (δ+a3))) for δ in deltas) + 
    sum(im*(δ-a3)[symb_to_ind(a)]*exp(im * dot(k, (δ-a3))) for δ in deltas)  # same signs for AB and BA
dfz(k, a3 , dir) = ifelse(dir != :z, 0, im*a3[3]*cis(dot(k, a3)))
df(k, deltas, a) = sum(im*δ[symb_to_ind(a)]*exp(im * dot(k, δ)) for δ in deltas)
dg(k, ds,a) = dg(k, ds[1], ds[2], ds[3], a)
dg(k,d1,d2,d3,a) = im*(2d1-d2)[symb_to_ind(a)]*cis(dot(k,2d1-d2)) + 
                    im*(2d2-d3)[symb_to_ind(a)]*cis(dot(k,2d2-d3)) + 
                    im*(2d3-d1)[symb_to_ind(a)]*cis(dot(k,2d3-d1)) -
                    im*(2d1-d3)[symb_to_ind(a)]*cis(dot(k,2d1-d3)) -
                    im*(2d2-d1)[symb_to_ind(a)]*cis(dot(k,2d2-d1)) -
                    im*(2d3-d2)[symb_to_ind(a)]*cis(dot(k,2d3-d2))


dnn(k,deltas,a) = dnn(k,deltas[1], deltas[2], deltas[3],a)
dnn(k, d1, d2, d3, a) = sum([-2*δ[symb_to_ind(a)]* sin(dot(k, δ)) for δ in [d1-d2,d2-d3,d3-d1]])

d2nn(k,deltas,dir1,dir2) = d2nn(k,deltas[1],deltas[2],deltas[3],dir1,dir2)

function d2nn(k, d1,d2,d3, dir1, dir2)
    i = symb_to_ind(dir1)
    j = symb_to_ind(dir2)
    sum([-2*δ[i]*δ[j]* cos(dot(k, δ)) for δ in [d1-d2,d2-d3,d3-d1]])
end


function d2_ferroaxial_ham3d(t, tp, tpz,tc, q, dir1, dir2) # model with the AA hoppings in z
    a = 1.0
    a1, a2, a3, δ1, δ2, δ3 = lattice_vectors(a)
    ds = [δ1, δ2, δ3]
    d2fk = d2f(q, ds, dir1, dir2)
    d2gk = d2g(q, ds, dir1, dir2)
    d2mk = d2m(q, ds,a3, dir1, dir2)
    off = t*(d2fk) + tp*d2gk + tpz * d2mk
    d2nnk= d2nn(q, ds, dir1, dir2)
    return [tc*d2nnk off; conj(off)  tc*d2nnk]
end

function symb_to_ind(dir)
    if dir == :x
        1
    elseif dir == :y
        2
    elseif dir == :z 
        3
    end
end

function d2f(k, deltas, dir1, dir2)
    i = symb_to_ind(dir1)
    j = symb_to_ind(dir2)
    sum(-δ[i]*δ[j]*cis(dot(k,δ)) for δ in deltas)
end

function d2m(k, deltas, a3, dir1, dir2)
    i = symb_to_ind(dir1)
    j = symb_to_ind(dir2)
    sum(-(δ+a3)[i]*(δ+a3)[j]*cis(dot(k,δ+a3)) for δ in deltas) + 
        sum(-(δ-a3)[i]*(δ-a3)[j]*cis(dot(k,δ-a3)) for δ in deltas)
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

function d2fz(k, a3, dir1, dir2)
    i = symb_to_ind(dir1)
    j = symb_to_ind(dir2)
    if dir1 == dir2 && dir1 == :z
        return -a3[j]^2*cis(dot(k,a3))
    else 
        return 0
    end
end