
#######################################
#
#  Action polynomial rings 
#
#######################################

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

  function DifferencePolyRing{T}(R::Ring, nelementary_symbols::Int, ndiffs::Int, cached::Bool, internal_ordering::Symbol) where {T}
    @req nelementary_symbols >= 0 "The number of elementary symbols must be nonnegative"
    @req ndiffs >= 0 "The number of endomorphisms must be nonnegative"
    elementary_symbols = map(x -> Symbol("u" * string(x)), 1:nelementary_symbols)
    upoly_ring = universal_polynomial_ring(R; cached = cached, internal_ordering = internal_ordering)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
    
    return new{T}(upoly_ring, elementary_symbols, ndiffs, false, jet_to_var, var_to_jet, jet_to_upoly_idx)
  end
 
  function DifferencePolyRing{T}(R::Ring, elementary_symbols::Vector{Symbol}, ndiffs::Int, cached::Bool, internal_ordering::Symbol) where {T}
    @req ndiffs >= 0 "The number of endomorphisms must be nonnegative"
    upoly_ring = universal_polynomial_ring(R; cached = cached, internal_ordering = internal_ordering)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
    
    return new{T}(upoly_ring, elementary_symbols, ndiffs, false, jet_to_var, var_to_jet, jet_to_upoly_idx)
  end

end

mutable struct DifferencePolyRingElem{T} <: ActionPolyRingElem{T}
  p::AbstractAlgebra.Generic.UniversalPolyRingElem{T}
  parent::DifferencePolyRing{T}
  initial::Any #Always of type DifferencePolyRingElem{T}
  permutation::Vector{Int}

  DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}) where {T} = new{T}(zero(dpr.upoly_ring), dpr)

  function DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}, upre::AbstractAlgebra.Generic.UniversalPolyRingElem{T}) where {T}
    @req dpr.upoly_ring === parent(upre) "The parent does not match"
    new{T}(upre, dpr)
  end

  function DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}, mpre::MPolyRingElem{T}) where {T}
    @req dpr.upoly_ring.mpoly_ring === parent(mpre) "The parent does not match"
    new{T}(mpre, dpr)
  end

end

#######################################
#
#  Rankings 
#
#######################################

# For a ranking of jet variables, e.g. u1[123], u3[041], ... one needs:
# - An ordering to compare the elementary symbols u_1, ..., u_n
# - An ordering of the multiindices (this is like a monomial_ordering)
# - A decision on which of the two has priority, e.g. position-over-term (pot)
#   or term-over-position (top). Here, the position of the jet variable corresponding
#   to (i, [a_1, ..., a_n]) is the first coordinate, i.e. i. More generally, this decision means
#   choosing an ordered partition of the set {1, ...,m}, pot and top being the choices of the two trivial
#   partitions

abstract type Ranking end

mutable struct DifferenceRanking{T} <: Ranking where {T}
  ring::DifferencePolyRing{T}
  partition::Vector{Vector{Int}}
  index_ordering_matrix::ZZMatrix
  riquier_matrix::ZZMatrix

  function DifferenceRanking{T}(dpr::DifferencePolyRing{T}, partition::Vector{Vector{Int}}, index_ordering_matrix::ZZMatrix) where {T}  
    return new{T}(dpr, partition, index_ordering_matrix)
  end
  
end

