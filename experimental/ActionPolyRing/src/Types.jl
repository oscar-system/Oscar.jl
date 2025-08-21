###############################################################################
#
#  Action polynomial rings 
#
###############################################################################

abstract type ActionPolyRing{T} <: Ring end 

abstract type ActionPolyRingElem{T} <: RingElem end

### Difference ###
mutable struct DifferencePolyRing{T} <: ActionPolyRing{T}
  upoly_ring::AbstractAlgebra.Generic.UniversalPolyRing{T}
  elementary_symbols::Vector{Symbol}
  ndiffs::Int
  are_perms_up_to_date::Bool
  jet_to_var::Any #Always of type Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}
  var_to_jet::Any #Always of type Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}
  jet_to_upoly_idx::Dict{Tuple{Int, Vector{Int}}, Int}
  ranking::Any #Alyways of type DifferenceRanking{T}
  permutation::Vector{Int}

  function DifferencePolyRing{T}(R::Ring, n_elementary_symbols::Int, ndiffs::Int) where {T}
    @req n_elementary_symbols >= 1 "The number of elementary symbols must be positive"
    elementary_symbols = map(x -> Symbol('u', x), 1:n_elementary_symbols)
    return DifferencePolyRing{T}(R, elementary_symbols, ndiffs)
  end
 
  function DifferencePolyRing{T}(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int) where {T}
    @req !is_empty(elementary_symbols) "The number of elementary symbols must be positive"
    @req ndiffs >= 1 "The number of endomorphisms must be positive"
    upoly_ring = universal_polynomial_ring(R; cached = false)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
    
    return new{T}(upoly_ring, elementary_symbols, ndiffs, false, jet_to_var, var_to_jet, jet_to_upoly_idx)
  end

end

mutable struct DifferencePolyRingElem{T} <: ActionPolyRingElem{T}
  p::AbstractAlgebra.Generic.UniversalPolyRingElem{T}
  parent::DifferencePolyRing{T}
  is_perm_up_to_date::Bool
  permutation::Vector{Int}

  DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}) where {T} = new{T}(zero(dpr.upoly_ring), dpr, true, Int[])

  function DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}, upre::AbstractAlgebra.Generic.UniversalPolyRingElem{T}) where {T}
    @req dpr.upoly_ring === parent(upre) "The parent does not match"
    new{T}(upre, dpr, false, zeros(Int, length(upre)))
  end

  function DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}, mpre::MPolyRingElem{T}) where {T}
    upr = dpr.upoly_ring
    @req upr.mpoly_ring === parent(mpre) "The parent does not match"
    new{T}(upr(collect(coefficients(mpre)), collect(exponents(mpre))), dpr, false, zeros(Int, length(mpre)))
  end

end

### Differential ###
mutable struct DifferentialPolyRing{T} <: ActionPolyRing{T}
  upoly_ring::AbstractAlgebra.Generic.UniversalPolyRing{T}
  elementary_symbols::Vector{Symbol}
  ndiffs::Int
  are_perms_up_to_date::Bool
  jet_to_var::Any #Always of type Dict{Tuple{Int, Vector{Int}}, DifferentialPolyRingElem{T}}
  var_to_jet::Any #Always of type Dict{DifferentialPolyRingElem{T}, Tuple{Int, Vector{Int}}}
  jet_to_upoly_idx::Dict{Tuple{Int, Vector{Int}}, Int}
  ranking::Any #Alyways of type DifferentialRanking{T}
  permutation::Vector{Int}

  function DifferentialPolyRing{T}(R::Ring, n_elementary_symbols::Int, ndiffs::Int) where {T}
    @req n_elementary_symbols >= 1 "The number of elementary symbols must be positive"
    elementary_symbols = map(x -> Symbol('u', x), 1:n_elementary_symbols)
    return DifferentialPolyRing{T}(R, elementary_symbols, ndiffs)
  end
 
  function DifferentialPolyRing{T}(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int) where {T}
    @req !is_empty(elementary_symbols) "The number of elementary symbols must be positive"
    @req ndiffs >= 1 "The number of endomorphisms must be nonnegative"
    upoly_ring = universal_polynomial_ring(R; cached = false)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferentialPolyRingElem{T}}()
    var_to_jet = Dict{DifferentialPolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
    
    return new{T}(upoly_ring, elementary_symbols, ndiffs, false, jet_to_var, var_to_jet, jet_to_upoly_idx)
  end

end

mutable struct DifferentialPolyRingElem{T} <: ActionPolyRingElem{T}
  p::AbstractAlgebra.Generic.UniversalPolyRingElem{T}
  parent::DifferentialPolyRing{T}
  is_perm_up_to_date::Bool
  permutation::Vector{Int}

  DifferentialPolyRingElem{T}(dpr::DifferentialPolyRing{T}) where {T} = new{T}(zero(dpr.upoly_ring), dpr, true, Int[])

  function DifferentialPolyRingElem{T}(dpr::DifferentialPolyRing{T}, upre::AbstractAlgebra.Generic.UniversalPolyRingElem{T}) where {T}
    @req dpr.upoly_ring === parent(upre) "The parent does not match"
    new{T}(upre, dpr, false, zeros(Int, length(upre)))
  end

  function DifferentialPolyRingElem{T}(dpr::DifferentialPolyRing{T}, mpre::MPolyRingElem{T}) where {T}
    upr = dpr.upoly_ring
    @req upr.mpoly_ring === parent(mpre) "The parent does not match"
    new{T}(upr(collect(coefficients(mpre)), collect(exponents(mpre))), dpr, false, zeros(Int, length(mpre)))
  end

end

###############################################################################
#
#  Iterator types 
#
###############################################################################

struct ActionPolyCoeffs{T, PolyT<:ActionPolyRingElem{T}}
   poly::PolyT
end

struct ActionPolyExponentVectors{T, PolyT<:ActionPolyRingElem{T}}
   poly::PolyT
end

struct ActionPolyTerms{T, PolyT<:ActionPolyRingElem{T}}
   poly::PolyT
end

struct ActionPolyMonomials{T, PolyT<:ActionPolyRingElem{T}}
   poly::PolyT
end

###############################################################################
#
#  Rankings 
#
###############################################################################

# For a ranking of jet variables, e.g. u1[123], u3[041], ... one needs:
# - An ordering to compare the elementary symbols u_1, ..., u_n
# - An ordering of the multiindices (this is like a monomial_ordering)
# - A decision on which of the two has priority, e.g. position-over-term (pot)
#   or term-over-position (top). Here, the position of the jet variable corresponding
#   to (i, [a_1, ..., a_n]) is the first coordinate, i.e. i. More generally, this decision means
#   choosing an ordered partition of the set {1, ...,m}, pot and top being the choices of the two trivial
#   partitions

abstract type ActionRanking end

mutable struct ActionPolyRingRanking{PolyT} <: ActionRanking where {T, PolyT <: ActionPolyRing{T}}
  ring::PolyT
  partition::Vector{Vector{Int}}
  index_ordering_matrix::ZZMatrix
  riquier_matrix::ZZMatrix

  function ActionPolyRingRanking{PolyT}(dpr::PolyT, partition::Vector{Vector{Int}}, index_ordering_matrix::ZZMatrix) where {T, PolyT <: ActionPolyRing{T}}
    return new{PolyT}(dpr, partition, index_ordering_matrix)
  end
  
end
