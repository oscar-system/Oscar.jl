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
"""
function TropicalLinearSpace(plv::Vector{Union{Rational,PosInf,Oscar.IntegerUnion}},nElement::Int,rank::Int; M::Union{typeof(min),typeof(max)}=min)
    indexSet = findall(i->i!=inf, plv)
    bases = [ sort(Hecke.subsets(Vector{Int}(0:nElement-1),rank))[i] for i in indexSet ]
    val = Polymake.matroid.ValuatedMatroid{M}(BASES = bases, N_ELEMENTS = nElement,VALUATION_ON_BASES = Vector{Rational}(plv[indexSet]) )
    #return Polymake.tropical.linear_space{min}(val)
    P = Polymake.tropical.linear_space{M}(val)
    P = PolyhedralComplex(P)
    return TropicalLinearSpace{M,true}(P)
end


@doc Markdown.doc"""
    TropicalLinearSpace(plv::Vector)

Construct a tropical linear space from its Pluecker vector

# Examples
"""
function TropicalLinearSpace(plv::Vector{Union{Rational,PosInf,Oscar.IntegerUnion}},nElement::Int,rank::Int; M::Union{typeof(min),typeof(max)}=min)
    indexSet = findall(i->i!=inf, plv)
    bases = [ sort(Hecke.subsets(Vector{Int}(0:nElement-1),rank))[i] for i in indexSet ]
    val = Polymake.matroid.ValuatedMatroid{M}(BASES = bases, N_ELEMENTS = nElement,VALUATION_ON_BASES = Vector{Rational}(plv[indexSet]) )
    #return Polymake.tropical.linear_space{min}(val)
    P = Polymake.tropical.linear_space{M}(val)
    P = PolyhedralComplex(P)
    return TropicalLinearSpace{M,true}(P)
end

#needs Oscar type as entry
function TropicalLinearSpace(tropicalmatrix::MatElem)
  plv = Nemo.minors(tropicalmatrix, min(size(tropicalmatrix)[1], size(tropicalmatrix)[2])) 
  rk = rank(tropicalmatrix)
  nelement = size(tropicalmatrix)[2]
  return (plv,nelement, rk)
  ### TropicalLinearSpace(plv, nelement, rk)
end

function TropicalLinearSpace(tropicalmatrix::Matrix{Int})
  return TropicalLinearSpace(matrix(ZZ, tropicalmatrix))
end


function TropicalLinearSpace(tropicalmatrix::Matrix{fmpq})
  return TropicalLinearSpace(matrix(base_ring(tropicalmatrix), tropicalmatrix))
end


function TropicalLinearSpace(tropicalmatrix::Matrix{fmpz})
  return TropicalLinearSpace(matrix(base_ring(tropicalmatrix), tropicalmatrix))
end

function TropicalLinearSpace(tropicalmatrix::Matrix{Rational})
  return TropicalLinearSpace(matrix(QQ, tropicalmatrix))  
end


###
# 3. Basic properties
# -------------------
###
