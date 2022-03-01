###
# Tropical linear spaces in Oscar
# ===============================
###



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
export TropicalLinearSpace



###
# 2. Basic constructors
# ---------------------
###

@doc Markdown.doc"""
    TropicalLinearSpace(ideal::MPolyIdeal{fmpq_poly})

Construct a tropical linear space from a degree 1 polynomial ideal

# Examples
julia> R,(x_0,x_1,x_2,x_3,x_4,x_5)=PolynomialRing(ZZ,["x_0","x_1","x_2","x_3","x_4","x_5"])
julia> I=ideal(R,[-x_0+x_2+x_3,-x_1+x_2+x_4,-x_0+x_1+x_5])
julia> val = ValuationMap(QQ)
julia> TropicalLinearSpace(I,val)
"""
function TropicalLinearSpace(I::MPolyIdeal{fmpz_mpoly}, val)
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
    TropicalLinearSpace(plv::Vector)

Construct a tropical linear space from its Pluecker vector

# Examples
julia> R = tropical_semiring(min);
julia> plv = [R(e) for e in [2,1,1,0,0,zero(R)]];
julia> L = TropicalLinearSpace(plv, 2, 4)
julia> f_vector(L)

# Examples
julia> Kt, t = RationalFunctionField(QQ,"t");
julia> val = TropicalSemiringMap(Kt,t);
julia> A = matrix(Kt,[[t,4t,0,2],[1,4,1,t^2]])
julia> TropicalLinearSpace(A, val)
"""
function TropicalLinearSpace_impl(plv, rank, nElements, M)
    Zero = zero(tropical_semiring(M))
    indexSet = findall(i->i!=Zero, plv)
    bases = [ sort(Hecke.subsets(Vector{Int}(0:nElements-1), rank))[i] for i in indexSet ]
    val = Polymake.matroid.ValuatedMatroid{M}(BASES = bases, N_ELEMENTS = nElements,VALUATION_ON_BASES = [plv[i].data for i in indexSet])
    #return Polymake.tropical.linear_space{min}(val)
    P = Polymake.tropical.linear_space{M}(val)
    P = PolyhedralComplex(P)
    return TropicalLinearSpace{M,true}(P)
end

TropicalLinearSpace(plv::Vector{TropicalSemiringElem{typeof(min)}},rank::IntegerUnion, nElements::IntegerUnion) =
TropicalLinearSpace_impl(plv, rank, nElements, min)

TropicalLinearSpace(plv::Vector{TropicalSemiringElem{typeof(max)}},rank::IntegerUnion, nElements::IntegerUnion) =
TropicalLinearSpace_impl(plv, rank, nElements, max)

#needs Oscar type as entry

#TODO requires a fix of TropicalSemiringMap
function TropicalLinearSpace(tropicalmatrix::MatElem, val)
  plv = [val(p) for p in Nemo.minors(tropicalmatrix, min( nrows(tropicalmatrix), ncols(tropicalmatrix)) )]
  rk = rank(tropicalmatrix)
  nelement = max( nrows(tropicalmatrix), ncols(tropicalmatrix))
  println(typeof(plv))
  return TropicalLinearSpace(plv, rk, nelement)
end

function TropicalLinearSpace(tropicalmatrix::Matrix{Int})
  #which valuation?
  return TropicalLinearSpace(matrix(ZZ, tropicalmatrix))
end

function TropicalLinearSpace(tropicalmatrix::Matrix{Union{fmpq, fmpz}})
  #which valuation?
  return TropicalLinearSpace(matrix(base_ring(tropicalmatrix), tropicalmatrix))
end

# -------------------
###
