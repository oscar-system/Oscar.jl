###############################################################################
#
#  Action maps 
#
###############################################################################

abstract type ActionMap{D <: Ring} <: Map{D, D, Any, Any} end

# =============================================================================
# Shifts
# =============================================================================

abstract type ActionShift{D <: Ring} <: ActionMap{D} end

struct TrivialActionShift{D <: Ring} <: ActionShift{D}
  base_ring::D
end

struct NontrivialActionShift{D <: Ring} <: ActionShift{D}
  underlying_map::Map{D, D}
    
  function NontrivialActionShift{D}(m::Map{D, D}) where {D <: Ring}
    @req domain(m) === codomain(m) "The domain and codomain of a shift operator must coincide"
    return new{D}(m)
  end
end

# =============================================================================
# Derivations
# =============================================================================

abstract type ActionDerivation{D <: Ring} <: ActionMap{D} end

struct TrivialActionDerivation{D <: Ring} <: ActionDerivation{D}
  base_ring::D
end

struct NontrivialActionDerivation{D <: Ring} <: ActionDerivation{D}
  underlying_map::Map{D, D}
    
  function NontrivialActionDerivation{D}(m::Map{D, D}) where {D <: Ring}
    @req domain(m) === codomain(m) "The domain and codomain of a derivation must coincide"
    return new{D}(m)
  end
end

###############################################################################
#
#  Action polynomial rings 
#
###############################################################################

abstract type ActionPolyRing{T} <: Ring end
abstract type ActionPolyRingElem{T} <: RingElem end

 # To be implemented by subtypes: 
 # Mandatory:
 #  Getters:
 #   __are_perms_up_to_date(R::MyActionPolyRing) -> Bool
 #   __is_perm_up_to_date(x::MyActionPolyRing) -> Bool
 #   __jtu_idx(R::MyActionPolyRing) -> Dict{Tuple{Int, Vector{Int}}, Int}
 #   __jtv(R::MyActionPolyRing) -> Dict{Tuple{Int, Vector{Int}}, MyActionPolyRingElem{T}} 
 #   __perm_for_sort(R::MyActionPolyRing) -> Vector{Int}
 #   __perm_for_sort_poly(x::MyActionPolyRingElem) -> Vector{Int}
 #   __vtj(R::MyActionPolyRing) -> Dict{MyActionPolyRingElem{T}, Tuple{Int, Vector{Int}}}
 #   base_ring(R::MyActionPolyRing) -> AbstractAlgebra.UniversalPolyRing{T}
 #   data(x::MyActionPolyRingElem) -> AbstractAlgebra.UniversalPolyRingElem{T}
 #   elem_type(::Type{MyActionPolyRing{T}}) = MyActionPolyRingElem{T} 
 #   action_indeterminates(R::MyActionPolyRing) -> Vector{Symbols}
 #   action_maps(R::MyActionPolyRing) -> Vector{ActionMap}
 #   parent(x::MyActionPolyRingElem{T}) -> MyActionPolyRing{T} 
 #   parent_type(::Type{MyActionPolyRingElem{T}}) = MyActionPolyRing{T} 
 #   ranking(R::MyActionPolyRing{T}) -> ActionPolyRingRanking{MyActionPolyRing{T}}
 #
 #  Setters:
 #   function __set_are_perms_up_to_date!(R::MyActionPolyRing, update::Bool)
 #     R.are_perms_up_to_date = update
 #   end
 #
 #   function __set_is_perm_up_to_date!(x::MyActionPolyRingElem, update::Bool)
 #     x.is_perm_up_to_date = update
 #   end
 #
 #   function __set_perm_for_sort!(R::MyActionPolyRing)
 #     R.permutation = sortperm(R.(gens(base_ring(R))); rev = true)
 #     __set_are_perms_up_to_date!(R, true)
 #   end
 #
 #   function __set_perm_for_sort_poly!(x::MyActionPolyRingElem)
 #     exps = collect(exponents(data(x)))
 #     n = length(exps)
 #     if n <= 1
 #       x.permutation = collect(1:n)
 #       return __set_is_perm_up_to_date!(x, true)
 #     end
 #
 #     perm = (parent(x).permutation)[1:min(end, length(exps[1]))] #trim unused indices (avoids padding with zeros)
 #
 #     x.permutation = sortperm(exps; lt=__my_lt_for_vec(perm), rev=true)
 #     __set_is_perm_up_to_date!(x, true)
 #   end  
 #  
 #  Constructors:
 #    function (R::MyActionPolyRing{T})(x::MyActionPolyRingElem{T}) where {T}
 #      @req parent(x) === R "Wrong parent"
 #      return a
 #    end
 #

### Difference ###
mutable struct DifferencePolyRing{T} <: ActionPolyRing{T}
  upoly_ring::AbstractAlgebra.UniversalPolyRing{T}
  action_indeterminates::Vector{Symbol}
  action_maps::Any #Always of type Vector{Union{TrivialActionShift{parent_type(T)}, NontrivialActionShift{parent_type(T)}}}
  are_perms_up_to_date::Bool
  jet_to_var::Any #Always of type Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}
  var_to_jet::Any #Always of type Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}
  jet_to_upoly_idx::Dict{Tuple{Int, Vector{Int}}, Int}
  ranking::Any #Alyways of type ActionPolyRanking{T, DifferencePolyRing{T}}
  permutation::Vector{Int}

  function DifferencePolyRing{T}(R::D, n_action_indeterminates::Int, action_maps::Vector{<:Union{TrivialActionShift{D}, NontrivialActionShift{D}}}) where {D <: Ring, T}
    @req n_action_indeterminates >= 1 "The number of action indeterminates must be positive"
    action_indeterminates = map(x -> Symbol('u', x), 1:n_action_indeterminates)

    return DifferencePolyRing{T}(R, action_indeterminates, action_maps)
  end
 
  function DifferencePolyRing{T}(R::D, action_indeterminates::Vector{Symbol}, action_maps::Vector{<:Union{TrivialActionShift{D}, NontrivialActionShift{D}}}) where {D <: Ring, T}
    @req !is_empty(action_indeterminates) "The number of action indeterminates must be positive"
    @req !is_empty(action_maps) "The number of shift operators must be positive"
    upoly_ring = universal_polynomial_ring(R; cached = false)
    
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferencePolyRingElem{T}}()
    var_to_jet = Dict{DifferencePolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
    
    return new{T}(upoly_ring, action_indeterminates, action_maps, false, jet_to_var, var_to_jet, jet_to_upoly_idx)
  end

end

mutable struct DifferencePolyRingElem{T} <: ActionPolyRingElem{T}
  p::AbstractAlgebra.UniversalPolyRingElem{T}
  parent::DifferencePolyRing{T}
  is_perm_up_to_date::Bool
  permutation::Vector{Int}

  DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}) where {T} = new{T}(zero(dpr.upoly_ring), dpr, true, Int[])

  function DifferencePolyRingElem{T}(dpr::DifferencePolyRing{T}, upre::AbstractAlgebra.UniversalPolyRingElem{T}) where {T}
    @req dpr.upoly_ring === parent(upre) "The parent does not match"
    new{T}(upre, dpr, false, zeros(Int, length(upre)))
  end

end

### Differential ###
mutable struct DifferentialPolyRing{T} <: ActionPolyRing{T}
  upoly_ring::AbstractAlgebra.UniversalPolyRing{T}
  action_indeterminates::Vector{Symbol}
  action_maps::Any #Always of type Vector{Union{TrivialActionDerivation{parent_type(T)}, NontrivialActionDerivation{parent_type(T)}}}
  are_perms_up_to_date::Bool
  jet_to_var::Any #Always of type Dict{Tuple{Int, Vector{Int}}, DifferentialPolyRingElem{T}}
  var_to_jet::Any #Always of type Dict{DifferentialPolyRingElem{T}, Tuple{Int, Vector{Int}}}
  jet_to_upoly_idx::Dict{Tuple{Int, Vector{Int}}, Int}
  ranking::Any #Always of type ActionPolyRanking{T, DifferentialPolyRing{T}}
  permutation::Vector{Int}

  function DifferentialPolyRing{T}(R::D, n_action_indeterminates::Int, action_maps::Vector{<:Union{TrivialActionDerivation{D}, NontrivialActionDerivation{D}}}) where {D <: Ring, T}
    @req n_action_indeterminates >= 1 "The number of action indeterminates must be positive"
    action_indeterminates = map(x -> Symbol('u', x), 1:n_action_indeterminates)

    return DifferentialPolyRing{T}(R, action_indeterminates, action_maps)
  end
 
  function DifferentialPolyRing{T}(R::D, action_indeterminates::Vector{Symbol}, action_maps::Vector{<:Union{TrivialActionDerivation{D}, NontrivialActionDerivation{D}}}) where {D <: Ring, T}
    @req !is_empty(action_indeterminates) "The number of action indeterminates must be positive"
    @req !is_empty(action_maps) "The number of derivations must be positive"
    upoly_ring = universal_polynomial_ring(R; cached = false)
        
    jet_to_var = Dict{Tuple{Int, Vector{Int}}, DifferentialPolyRingElem{T}}()
    var_to_jet = Dict{DifferentialPolyRingElem{T}, Tuple{Int, Vector{Int}}}()
    jet_to_upoly_idx = Dict{Tuple{Int, Vector{Int}}, Int}()
        
    return new{T}(upoly_ring, action_indeterminates, action_maps, false, jet_to_var, var_to_jet, jet_to_upoly_idx)
  end

end

mutable struct DifferentialPolyRingElem{T} <: ActionPolyRingElem{T}
  p::AbstractAlgebra.UniversalPolyRingElem{T}
  parent::DifferentialPolyRing{T}
  is_perm_up_to_date::Bool
  permutation::Vector{Int}

  DifferentialPolyRingElem{T}(dpr::DifferentialPolyRing{T}) where {T} = new{T}(zero(dpr.upoly_ring), dpr, true, Int[])

  function DifferentialPolyRingElem{T}(dpr::DifferentialPolyRing{T}, upre::AbstractAlgebra.UniversalPolyRingElem{T}) where {T}
    @req dpr.upoly_ring === parent(upre) "The parent does not match"
    new{T}(upre, dpr, false, zeros(Int, length(upre)))
  end

end

###############################################################################
#
#  Iterator types 
#
###############################################################################

struct ActionPolyCoeffs{PolyT<:ActionPolyRingElem} poly::PolyT end
struct ActionPolyExponentVectors{PolyT<:ActionPolyRingElem} poly::PolyT end
struct ActionPolyTerms{PolyT<:ActionPolyRingElem} poly::PolyT end
struct ActionPolyMonomials{PolyT<:ActionPolyRingElem} poly::PolyT end

###############################################################################
#
#  Rankings 
#
###############################################################################

# For a ranking of jet variables, e.g. u1[123], u3[041], ... one needs:
# - An ordering to compare the action indeterminates u_1, ..., u_n
# - An ordering of the multiindices (this is like a monomial_ordering)
# - A decision on which of the two has priority, e.g. position-over-term (pot)
#   or term-over-position (top). Here, the position of the jet variable corresponding
#   to (i, [a_1, ..., a_n]) is the first coordinate, i.e. i. More generally, this decision means
#   choosing an ordered partition of the set {1, ...,m}, pot and top being the choices of the two trivial
#   partitions

mutable struct ActionPolyRingRanking{PolyT <: ActionPolyRing}
  ring::PolyT
  partition::Vector{Vector{Int}}
  index_ordering_matrix::ZZMatrix
  riquier_matrix::ZZMatrix

  function ActionPolyRingRanking{PolyT}(S::PolyT, partition::Vector{Vector{Int}}, index_ordering_matrix::ZZMatrix) where {T, PolyT <: ActionPolyRing{T}}
    return new{PolyT}(S, partition, index_ordering_matrix)
  end
  
end

