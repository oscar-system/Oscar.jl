################################################################################
#
#  Tropicalization of zero-dimensional ideals
#
#  WARNING: assumes without test that `I` is zero-dimensional
#
#  WARNING: currently p-adic valuation only
#
################################################################################

function tropical_variety_zerodimensional_eigenvalue(I::MPolyIdeal,nu::TropicalSemiringMap{QQField,QQFieldElem,<:Union{typeof(min),typeof(max)}})
    k,(a,_) = number_field(I)
    zk = maximal_order(k)
    p = uniformizer(nu)
    lp = [x[1] for x = prime_decomposition(zk,p)]
    ma = representation_matrix(a)
    mb = representation_matrix(k(lp[1].gen_two*lp[2].gen_two^2))
    @assert iszero(ma*mb - mb*ma)
    Qp = padic_field(p; precision=10)
    TropVDict = simultaneous_diagonalization([map_entries(Qp, ma),map_entries(Qp, mb)])

    TropVPoints = collect(values(TropVDict))
    TropVPointsUnique = unique(TropVPointsMults)
    Sigma = polyhedral_complex(IncidenceMatrix([[i] for i in 1:length(TropVPointsUnique)]), TropVPointsUnique)
    TropVMults = [ZZ(length(findall(isequal(p),TropVPoints))) for p in TropVPointsUnique]
    TropV = tropical_variety(Sigma,TropVMults)
    set_attribute!(TropV,:algebraic_points,collect(keys(TropVDict)))
    return TropV
end


###
# Code by Claus:
###
function slope_eigenspace(M::MatElem{T}) where T <: Hecke.NonArchLocalFieldElem
    f = charpoly(M)
    lf = Hecke.slope_factorization(f)
    # @req all(x->x==1, values(lf))

    se = Dict{typeof(f), typeof(M)}()
    k = base_ring(M)
    zk = maximal_order(k)

    for f = keys(lf)
        se[f] = kernel(f(M), side = :right) #hopefully, this is in rref
    end
    @assert sum(ncols(x) for x = values(se)) == nrows(M)
    return se
end

function _intersect(M::MatElem{T}, N::MatElem{T}) where T <: Hecke.FieldElem
    k = base_ring(M)
    I = [M N]
    PR = maximum(precision, I)
    pr = minimum(precision, I)
    if pr != PR
        for i = eachindex(I)
            I[i] = setprecision(I[i], pr)
        end
    end

    v = kernel(I, side = :right) #precision issues...
    l = M*v[1:ncols(M), 1:ncols(v)]
    return transpose(rref(transpose(l))[2])
end

function valuation_of_roots(f::PolyRingElem{<:Hecke.NonArchLocalFieldElem})
    iszero(f) && error("polynomial must not be zero")
    return (valuation(constant_coefficient(f)) - valuation(leading_coefficient(f)))//degree(f)
end

function simultaneous_diagonalization(v::Vector{<:MatElem{T}}) where T <: Hecke.NonArchLocalFieldElem

    k = base_ring(v[1])
    @assert all(x->base_ring(x) == k, v)
    n = nrows(v[1])
    @assert all(x->ncols(x) == nrows(x) == n, v)

    vv = map(slope_eigenspace, v)

    d = Dict(v => [valuation_of_roots(k)] for (k,v) = vv[1])
    @assert sum(ncols(x) for x = keys(d)) == n
    for i=2:length(vv)
        dd = typeof(d)()
        for (mat, pol_vec) = d
            for (p, m) = vv[i]
                j = _intersect(mat, m)
                if ncols(j) > 0
                    dd[j] = push!(copy(pol_vec), valuation_of_roots(p))
                end
            end
        end
        d = dd
        @assert sum(ncols(x) for x = keys(d)) == n
    end

    return d
end
