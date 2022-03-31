###
# Tropical linear spaces in Oscar
# ===============================
###

export TropicalLinearSpace


###
# 1. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention, affecting initial ideals, Pluecker vector, etc.
# EMB = true or false:
#   embedded or abstract tropical linear space
#   embedded tropical linear space = weighted polyhedral complex in euclidean space
#   abstract tropical linear space = weighted hypergraph with enumerated vertices
###

@attributes mutable struct TropicalLinearSpace{M,T} <: TropicalVarietySupertype{M,T}
    # GapTV::GapObj
    polyhedralComplex::PolyhedralComplex
    function TropicalLinearSpace{M,T}(Sigma::PolyhedralComplex) where {M,T}
        return new{M,T}(Sigma)
    end
end



###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalLinearSpace(I::MPolyIdeal, val::TropicalSemiringMap)

Constructs a tropical linear space from a degree 1 polynomial ideal `I` and a map to the tropical semiring `val`.

# Examples

```jldoctest
julia> R,(x1,x2,x3,x4,x5,x6) = PolynomialRing(ZZ,6)
(Multivariate Polynomial(Multivariate Polynomial Ring in 6 variables x1, x2, x3, x4, ..., x6 over Integer Ring, fmpz_mpoly[x1, x2, x3, x4, x5, x6])

julia> I = ideal(R,[-x1+x3+x4,-x2+x3+x5,-x1+x2+x6])
ideal(-x1 + x3 + x4, -x2 + x3 + x5, -x1 + x2 + x6)

julia> val = TropicalSemiringMap(ZZ)
The trivial valuation on Integer Ring

julia> TropicalLinearSpace(I,val)
TropicalLinearSpace{min, true}(A polyhedral complex in ambient dimension 6, #undef)
```
"""
function TropicalLinearSpace(I::MPolyIdeal, val)
    R = base_ring(I)
    n = ngens(R)
    g = gens(I)
    l = length(g)
    id = [ v for v in identity_matrix(ZZ,n)]
    A = [ [ coeff(p,Vector{Int}(v)) for v in [id[i,:] for i in 1:nrows(id)] ] for p in gens(I) ]
    B = kernel_basis(matrix(QQ,A))
    return TropicalLinearSpace(matrix(B),val)
end

@doc Markdown.doc"""
    TropicalLinearSpace(plv::Vector{TropicalSemiringElem}, rank::IntegerUnion, nElements::IntegerUnion)

Constructs a tropical linear space from a (tropical) Pluecker vector `plv`, rank `rank`, and size of the ground set `nElements`.

#Examples
```jldoctest
julia> R = TropicalSemiring(min);

julia> plv = [R(e) for e in [2,1,1,0,0,zero(R)]];

julia> L = TropicalLinearSpace(plv, 2, 4)
TropicalLinearSpace{min, true}(A polyhedral complex in ambient dimension 4, #undef)

julia> f_vector(L)
2-element Vector{Int64}:
 1
 3
```
"""
function TropicalLinearSpace_impl(plv, rank, nElements, M)
    Zero = zero(TropicalSemiring(M))
    indexSet = findall(i->i!=Zero, plv)
    bases = [ sort(Hecke.subsets(Vector{Int}(0:nElements-1), rank))[i] for i in indexSet ]
    val = Polymake.matroid.ValuatedMatroid{M}(BASES = bases, N_ELEMENTS = nElements,VALUATION_ON_BASES = [plv[i].data for i in indexSet])
    #return Polymake.tropical.linear_space{min}(val)
    P = Polymake.tropical.linear_space{M}(val)
    P = PolyhedralComplex{fmpq}(P)
    return TropicalLinearSpace{M,true}(P)
end
TropicalLinearSpace(plv::Vector{TropicalSemiringElem{typeof(min)}},rank::IntegerUnion, nElements::IntegerUnion) =
TropicalLinearSpace_impl(plv, rank, nElements, min)
TropicalLinearSpace(plv::Vector{TropicalSemiringElem{typeof(max)}},rank::IntegerUnion, nElements::IntegerUnion) =
TropicalLinearSpace_impl(plv, rank, nElements, max)


@doc Markdown.doc"""
    TropicalLinearSpace(M::MatElem,val::TropicalSemiringMap)

Constructs a tropical linear space from a matrix `M` and a map to the tropical semiring `val`.

# Examples
```jldoctest
julia> Kt, t = RationalFunctionField(QQ,"t");

julia> val = TropicalSemiringMap(Kt,t);

julia> A = matrix(Kt,[[t,4*t,0,2],[1,4,1,t^2]]);

julia> TropicalLinearSpace(A, val)
TropicalLinearSpace{min, true}(A polyhedral complex in ambient dimension 4, #undef)

julia> p = 3;

julia> val = TropicalSemiringMap(QQ, p);

julia> A = matrix(QQ, [[3,7,5,1], [9,7,1,2]])
[3   7   5   1]
[9   7   1   2]

julia> TropicalLinearSpace(A,val)
TropicalLinearSpace{min, true}(A polyhedral complex in ambient dimension 4, #undef)
```
"""
function TropicalLinearSpace(M::MatElem, val)
  plv = [val(p) for p in Nemo.minors(M, min(nrows(M), ncols(M)))]
  rk = rank(M)
  nelement = max(nrows(M), ncols(M))
  return TropicalLinearSpace(plv, rk, nelement)
end


###
# 3. Basic properties
# -------------------
###
