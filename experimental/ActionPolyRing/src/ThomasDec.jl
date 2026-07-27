#=
@doc raw"""
    leader_shift_for_partial_reduction(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the jet difference vector `k` such that `apply_action(q, k)` has the same leader
as the highest-ranking unreduced jet variable in `p`. If `p` is already partially reduced 
with respect to `q`, return `nothing`. Note that the notion of an 'unreduced' jet variable
depends on the type of `p` (and `q`):

- If the inputs are of type `DifferentialPolyRingElem`, then the occurrence of the jet variable in
  question is sufficient.
- If the inputs are of type `DifferencePolyRingElem`, then the degree of `p` in this variable must be
  greater or equal to the degree of `q` in its leader.
"""=#
function _leader_shift_for_partial_reduction(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_constant(q) "Cannot reduce with respect to a constant polynomial"
  is_constant(p) && return nothing
  
  apr = parent(q)
  ld_q = leader(q)
  ld_i, ld_idx = __vtj(apr)[ld_q]
  
  for var in vars(p)
    var_i, var_idx = __vtj(apr)[var]
    
    if var_i == ld_i
      is_derivative = true
      for k in 1:length(ld_idx)
        if var_idx[k] < ld_idx[k]
          is_derivative = false
          break
        end
      end 
    
      if is_derivative && var_idx != ld_idx
        if __is_proper_shift_reducible(p, q, var)
          return var_idx .- ld_idx 
        end
      end 
    end 
  end 
  
  return nothing
end

__is_proper_shift_reducible(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: DifferencePolyRingElem} = degree(p, var) >= degree(q, leader(q))
__is_proper_shift_reducible(p::PolyT, q::PolyT, var::PolyT) where {PolyT <: DifferentialPolyRingElem} = true

@doc raw"""
    is_partially_reduced(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return `true` if the polynomial `p` is partially reduced with respect to the non-constant
polynomial `q`. 
"""
is_partially_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem} = isnothing(_leader_shift_for_partial_reduction(p, q))

@doc raw"""
    partially_reduce(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the partial remainder of `p` with respect to the non-constant polynomial `q`, by performing successive pseudo-divisions
of `p` by proper derivatives (or shifts) of `q`. This means that the returned polynomial has strictly smaller degree in each jet
variable that is a proper derivative (or shift) of the leader of `q`.
"""
function partially_reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_constant(q) "Cannot reduce with respect to a constant polynomial"
    
  p_red = p
  shift = _leader_shift_for_partial_reduction(p_red, q)
    
  while shift !== nothing
    p_red = pseudorem(p_red, apply_action(q, shift))
    shift = _leader_shift_for_partial_reduction(p_red, q)
  end
    
  return p_red
end

@doc raw"""
    is_reduced(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return `true` if the polynomial `p` is fully reduced with respect to the non-constant
polynomial `q`. This means `p` is partially reduced with respect to `q`, and the degree 
of `p` in the leader of `q` is strictly less than the degree of `q` in its leader.
"""
function is_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_constant(q) "Cannot reduce with respect to a constant polynomial"
  
  is_constant(p) && return true
  is_partially_reduced(p, q) || return false
  
  ld_q = leader(q)
  return degree(p, ld_q) < degree(q, ld_q)
end

@doc raw"""
    reduce(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the full remainder of `p` with respect to the non-constant polynomial `q`. This remainder is
reduced with respect to `q` in the sense that the degree of `p` in each derivative (or shift) of the
leader of `q` is strictly smaller than the degree of the respective derivative (or shift) of `q` in that
same jet variable.
"""
function reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_constant(q) "Cannot reduce with respect to a constant polynomial"
  return pseudorem(partially_reduce(p,q), q)
end
