
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
  @req !is_constant(q) "Cannot compute leader shift with respect to a constant polynomial"
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

Return `true` if the action polynomial `p` is partially reduced with respect to the nonzero action polynomial `q`.
"""
function is_partially_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  @req !is_zero(q) "Cannot partially reduce with respect to the zero polynomial"
  is_constant(q) && return is_zero(p)
  return isnothing(_leader_shift_for_partial_reduction(p, q))
end

@doc raw"""
    is_partially_reduced(p::ActionPolyRingElem, S::Vector{<:ActionPolyRingElem})

Return `true` if `p` is partially reduced with respect to every element in the set `S`.
"""
is_partially_reduced(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem} = all(q -> is_partially_reduced(p, q), S)

@doc raw"""
    partially_reduce(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the partial remainder of `p` with respect to the non-zero polynomial `q`, by performing successive pseudo-divisions
of `p` by proper derivatives (or shifts) of `q`. This means that the returned polynomial has strictly smaller degree in each jet
variable that is a proper derivative (or shift) of the leader of `q`. If `q` is a non-zero constant then the zero polynomial is
returned.
"""
function partially_reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_zero(q) "Cannot partially reduce with respect to the zero polynomial"
  is_constant(q) && return zero(p)  

  p_red = p
  shift = _leader_shift_for_partial_reduction(p_red, q)
    
  while shift !== nothing
    p_red = pseudorem(p_red, apply_action(q, shift))
    shift = _leader_shift_for_partial_reduction(p_red, q)
  end
    
  return p_red
end

@doc raw"""
    partially_reduce(p::ActionPolyRingElem, S::Vector{ActionPolyRingElem})

Partially reduce the action polynomial `p` with respect to the vector `S`. This is done by pre-sorting `S`
with respect to Ritt ordering, filtering out zero polynomials and then performing top-down partial
reductions of `p` by the remaining elements of `S` until no further reductions are possible.
"""
function partially_reduce(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
  S = filter(!is_zero, S)
  any(is_constant, S) && return zero(p)

  sorted_S = sort(S, lt=ritt_is_less)

  res = p
  changed = true
  while changed && !iszero(res)
    changed = false
    for q in Iterators.reverse(sorted_S)
      if !is_partially_reduced(res, q)
        res = partially_reduce(res, q)
        changed = true
      end
    end
  end

  return res
end

@doc raw"""
    is_reduced(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return `true` if the polynomial `p` is fully reduced with respect to the non-constant
polynomial `q`. This means `p` is partially reduced with respect to `q`, and the degree 
of `p` in the leader of `q` is strictly less than the degree of `q` in its leader.
"""
function is_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_zero(q) "Cannot reduce with respect to the zero polynomial"
  is_constant(q) && return is_zero(p)
  
  is_constant(p) && return true
  is_partially_reduced(p, q) || return false
  
  ld_q = leader(q)
  return degree(p, ld_q) < degree(q, ld_q)
end

@doc raw"""
    is_reduced(p::ActionPolyRingElem, S::Vector{<:ActionPolyRingElem})

Return `true` if `p` is fully reduced with respect to every element in the set `S`.
"""
is_reduced(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem} = all(q -> is_reduced(p, q), S)

@doc raw"""
    reduce(p::ActionPolyRingElem, q::ActionPolyRingElem)

Return the full remainder of `p` with respect to the nonzero polynomial `q`. This remainder is
reduced with respect to `q` in the sense that the degree of `p` in each derivative (or shift) of the
leader of `q` is strictly smaller than the degree of the respective derivative (or shift) of `q` in that
same jet variable. If `q` is a nonzero constant then the zero polynomial is returned.
"""
function reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
  check_parent(p, q)
  @req !is_zero(q) "Cannot reduce with respect to the zero polynomial"
  is_constant(q) && return zero(p)
  return pseudorem(partially_reduce(p,q), q)
end

@doc raw"""
    reduce(p::ActionPolyRingElem, S::Vector{ActionPolyRingElem})

Reduce the action polynomial `p` with respect to the vector `S`. This is done by pre-sorting `S`
with respect to Ritt ordering, filtering out zero polynomials and then performing top-down
reductions of `p` by the elements of `S` until no further reductions are possible.
"""
function reduce(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
  S = filter(!is_zero, S)
  any(is_constant, S) && return zero(parent(p))

  sorted_S = sort(S, lt=ritt_is_less)

  res = p
  changed = true
  while changed && !iszero(res)
    changed = false
    for q in Iterators.reverse(sorted_S)
      if !is_reduced(res, q)
        res = reduce(res, q)
        changed = true
      end
    end
  end

  return res
end

@doc raw"""
    is_autoreduced(S::Vector{<:ActionPolyRingElem})

Return `true` if `S` is an autoreduced set. A set is autoreduced if every polynomial 
in the set is fully reduced with respect to all other polynomials in the set, and 
the set is sorted by Ritt ordering.
"""
function is_autoreduced(S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
  any(is_constant, S) && return (length(S) == 1 && is_constant(S[1]))
  
  for i in 1:length(S)
    for j in 1:length(S)
      i == j && continue
      !is_reduced(S[i], S[j]) && return false
    end
  end
  
  return issorted(S, lt=ritt_is_less)
end

@doc raw"""
    autoreduce(S::Vector{ActionPolyRingElem})

Compute an autoreduced set from the vector of action polynomials `S`. If at any point in the reduction
process a nonzero constant is discovered, the vector containing just this constant is returned.
"""
function autoreduce(S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
  # S will serve as the working list throughout this algorithm
  S = filter(!iszero, S)
  sort!(S, lt=ritt_is_less)
  
  A = PolyT[] # The "successively" built up autoreduced set 
  while !isempty(S)
    if !isempty(A) && is_constant(A[1])
      return [A[1]]
    elseif !isempty(S) && is_constant(S[1])
      return [S[1]]
    end
    
    p = popfirst!(S)
    original_p = p
    
    rank_dropped = false
    for q in Iterators.reverse(A)
      p = reduce(p, q) 
      
      if iszero(p)
        break
      end
      
      # Reduce by coefficients in the base ring 
      c_gcd = foldl(gcd, coefficients(p))
      if !isone(c_gcd)
        p = divexact(p, c_gcd)
      end
      
      # Detect rank drops of p
      if ritt_is_less(p, original_p)
        rank_dropped = true
        break 
      end
    end
    
    # Nothing to do here
    if iszero(p)
      continue 
    end
    
    # Reinsert both p and all element of A greater than p back into the working list S. 
    if rank_dropped
      # inserting elts of A
      while !isempty(A) && ritt_is_less(p, A[end])
        q_popped = pop!(A)
        idx = searchsortedfirst(S, q_popped, lt=ritt_is_less)
        insert!(S, idx, q_popped)
      end
      
      # inserting p
      idx = searchsortedfirst(S, p, lt=ritt_is_less)
      insert!(S, idx, p)
      
      continue
    end
    
    # No rank drop ever occurred, so appending p to A will give an autoreduced set
    push!(A, p)
  end
  
  return A
end
