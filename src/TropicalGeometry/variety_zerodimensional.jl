################################################################################
#
#  Tropicalization of zero-dimensional ideals
#
#  WARNING: assumes without test that `I` is zero-dimensional
#
#  WARNING: currently p-adic valuation only
#
################################################################################

@doc raw"""
    Oscar.tropical_variety_zerodimensional(I::MPolyIdeal,nu::TropicalSemiringMap{QQField,ZZRingElem,<:Union{typeof(min),typeof(max)}}; precision::Int=64)

Internal function for computing zero-dimensional tropical varieties over p-adic numbers via
finite precision Eigenvalue computation.  Assumes without test that `I` is zero-dimensional.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> R,(x1,x2,x3) = polynomial_ring(QQ,3);

julia> I = ideal([28*x3^2 - 1*x3 - 1,
                  2*x2 - x3,
                  2*x1 - x2]);

julia> nu = tropical_semiring_map(QQ,2)
Map into Min tropical semiring encoding the 2-adic valuation on Rational field

julia> TropI = Oscar.tropical_variety_zero_dimensional(I,nu)
Min tropical variety

julia> vertices(TropI)
2-element SubObjectIterator{PointVector{QQFieldElem}}:
 [-2, -1, 0]
 [-4, -3, -2]

julia> nu = tropical_semiring_map(QQ,3,max)
Map into Max tropical semiring encoding the 3-adic valuation on Rational field

julia> TropI = Oscar.tropical_variety_zero_dimensional(I,nu)
Max tropical variety

julia> vertices(TropI)
1-element SubObjectIterator{PointVector{QQFieldElem}}:
 [0, 0, 0]

```
"""
function tropical_variety_zero_dimensional(I::MPolyIdeal,nu::TropicalSemiringMap{QQField,ZZRingElem,<:Union{typeof(min),typeof(max)}}; precision::Int=64)
    # Construct the representation matrices of the multiplications by xi in K[x]/I
    _,x = number_field(I)
    mx = representation_matrix.(x)

    # Compute their simultaneous diagonalization numerically
    Qp = padic_field(uniformizer(nu), precision=precision)
    TropVDict = Oscar.simultaneous_diagonalization(map_entries.(Ref(Qp), mx))

    # Construct their tropical variety as a polyhedral complex consisting only of vertices
    # and a list of multiplicities
    TropVPoints = convention(nu)==min ? collect(values(TropVDict)) : -collect(values(TropVDict))
    TropVPointsUnique = unique(TropVPoints)
    Sigma = polyhedral_complex(IncidenceMatrix([[i] for i in 1:length(TropVPointsUnique)]), TropVPointsUnique)
    TropVMults = [ZZ(length(findall(isequal(p),TropVPoints))) for p in TropVPointsUnique]
    TropV = tropical_variety(Sigma,TropVMults,convention(nu))
    set_attribute!(TropV,:algebraic_points,collect(keys(TropVDict)))
    return TropV
end


###
# Code by Claus:
###
function slope_eigenspace(M::MatElem{T}) where T <: Hecke.NonArchLocalFieldElem
    f = charpoly(M)
    lf = Hecke.slope_factorization(f)
    # @req all(==(1), values(lf))

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
