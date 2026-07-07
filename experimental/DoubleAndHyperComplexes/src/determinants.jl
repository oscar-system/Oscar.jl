@doc raw"""
    det(
      c::AbsHyperComplex{T};
      lower_bound::Union{Int, Nothing}=has_lower_bound(c, 1) ? lower_bound(c, 1) : nothing,
      upper_bound::Union{Int, Nothing}=has_upper_bound(c, 1) ? upper_bound(c, 1) : nothing,
      direction::Symbol=:from_left_to_right
    ) where {T<:OFPModule}

Compute the determinant of a one-dimensional complex as outlined for instance in [GelfandKapranovZelevinsky94](@cite). The keyword `direction` can be set to either `:from_left_to_right` to start from an upper bound on the non-zero modules, or `:from_right_to_left` to start from a lower bound. All modules in the complex need to be free and defined over a polynomial ring.

!!! note
    It is assumed, but not checked that there exists a single range of integers `r = a:b` such that the modules `c[i]` are non-zero for `i` in `r` and zero otherwise. 
"""
function det(
    c::AbsHyperComplex{T};
    lower_bound::Union{Int, Nothing}=has_lower_bound(c, 1) ? lower_bound(c, 1) : nothing,
    upper_bound::Union{Int, Nothing}=has_upper_bound(c, 1) ? upper_bound(c, 1) : nothing,
    direction::Symbol= isnothing(lower_bound) ? :from_right_to_left : :from_left_to_right, 
    reduction_map=nothing
  ) where {T <: OFPModule}
  @assert is_one(dim(c)) "complex must be 1-dimensional"
  # We have a `direction` of the complex, depending on whether this is a chain, or a 
  # cochain complex. In addition, we have a direction for the computations to be 
  # carried out, i.e. with which morphism to start and in which direction to proceed. 
  # While the first is predetermined by the input, there might be preferable choices 
  # for the latter. 
  if direction == :from_left_to_right
    @assert !isnothing(lower_bound) "lower bound needs to be provided"
    return _det(c, Val{direction}, Val{Oscar.direction(c, 1)}, lower_bound::Int; reduction_map)
  elseif direction == :from_right_to_left
    @assert !isnothing(upper_bound) "upper bound needs to be provided"
    return _det(c, Val{direction}, Val{Oscar.direction(c, 1)}, upper_bound::Int; reduction_map)
  end
end

function _det(c::AbsHyperComplex, ::Type{Val{:from_left_to_right}}, 
    ::Type{Val{:cochain}}, ind::Int; reduction_map=nothing
  )
  while !can_compute_index(c, ind) || is_zero(c[ind]) 
    ind += 1
  end
  c[ind]::FreeMod
  R = base_ring(c[ind])::MPolyRing{<:RingElem}
  result = one(R)
  r = ngens(c[ind]) # the rank of the current map
  I = collect(1:r)
  while can_compute_map(c, ind) && !is_zero(map(c, ind))
    A = matrix(map(c, ind))
    m = nrows(A)
    n = ncols(A)
    J0 = first(combinations(n, r)) # Initialize variable
    p = one(R)
    # TODO: Maybe we can speed things up here by making choices of "small" determinants?
    for J in combinations(n, r)
      p = det(A[I, data(J)])
      is_zero(p) && continue
      J0 = J
      break
    end

    result = is_odd(ind) ? result*p : result//p
    
    I = Int[i for i in 1:n if !(i in data(J0))]
    r = n - r
    ind += 1
  end
  return result
end

function _det(c::AbsHyperComplex, ::Type{Val{:from_left_to_right}}, 
    ::Type{Val{:chain}}, ind::Int;
    reduction_map=nothing
  )
  while !can_compute_index(c, ind) || is_zero(c[ind]) 
    ind += 1
  end
  c[ind]::FreeMod
  R = base_ring(c[ind])::MPolyRing{<:RingElem}
  result = one(R)
  r = ngens(c[ind]) # the rank of the current map
  I = collect(1:r)
  while can_compute_map(c, ind+1) && !is_zero(map(c, ind+1))
    A = matrix(map(c, ind+1))
    m = nrows(A)
    n = ncols(A)
    A_subs = A[:, I]
    non_zero_rows = [i for i in 1:nrows(A_subs) if !is_zero(A[i, :])]
    A_ess = A_subs[non_zero_rows, :]

    if isnothing(reduction_map)
      J0 = Int[]
      p = one(R)
      # TODO: Maybe we can speed things up here by making choices of "small" determinants?

      for J in combinations(length(non_zero_rows), r)
        #@show data(J)
        #p = det(A[data(J), I])
        B = A_ess[data(J), :]
        any(iszero(B[:, k]) for k in 1:ncols(B)) && continue
        p = det(B)
        is_zero(p) && continue
        J0 = non_zero_rows[data(J)]
        break
      end

      # update variables
      result = is_even(ind) ? result*p : result//p
      I = Int[i for i in 1:m if !(i in J0)]
    else
      A_red = map_entries(reduction_map, A_ess)::MatrixElem{<:RingElem}
      rk, M = rref(transpose(A_red))
      @assert r == rk "reduction map did not recover the anticipated rank ($r anticipated vs. $rk real)"
      pivots = Int[findfirst(!is_zero(a) for a in M[i, :]) for i in 1:nrows(M) if !is_zero(M[i, :])]
      @assert length(pivots) == r "collection of pivot indices failed"
      @assert !is_zero(det(M[1:r, pivots]))
      @assert !is_zero(det(A_ess[pivots, :]))
      J = non_zero_rows[pivots]
      p = det(A[J, I])
      @assert !is_zero(p)

      # update variables
      result = is_even(ind) ? result*p : result//p
      I = Int[i for i in 1:m if !(i in J)]
    end
    
    r = m - r
    ind += 1
  end
  return result
end

function _det(c::AbsHyperComplex, ::Type{Val{:from_right_to_left}}, 
    ::Type{Val{:cochain}}, ind::Int;
    reduction_map=nothing
  )
  while !can_compute_index(c, ind) || is_zero(c[ind]) 
    ind -= 1
  end
  c[ind]::FreeMod
  R = base_ring(c[ind])::MPolyRing{<:RingElem}
  result = one(R)
  r = ngens(c[ind]) # the rank of the current map
  I = collect(1:r)
  while can_compute_map(c, ind-1) && !is_zero(map(c, ind-1))
    A = matrix(map(c, ind-1))
    m = nrows(A)
    n = ncols(A)
    J0 = first(combinations(m, r)) # Initialize variable
    p = one(R)
    # TODO: Maybe we can speed things up here by making choices of "small" determinants?
    for J in combinations(m, r)
      p = det(A[data(J), I])
      is_zero(p) && continue
      J0 = J
      break
    end

    result = is_even(ind) ? result*p : result//p
    
    I = Int[i for i in 1:m if !(i in data(J0))]
    r = m - r
    ind -= 1
  end
  return result
end

function _det(c::AbsHyperComplex, ::Type{Val{:from_right_to_left}}, 
    ::Type{Val{:chain}}, ind::Int;
    reduction_map=nothing
  )
  while !can_compute_index(c, ind) || is_zero(c[ind]) 
    ind -= 1
  end
  c[ind]::FreeMod
  R = base_ring(c[ind])
  result = one(R)
  r = ngens(c[ind]) # the rank of the current map
  I = collect(1:r)
  while can_compute_map(c, ind) && !is_zero(map(c, ind))
    A = matrix(map(c, ind))
    m = nrows(A)
    n = ncols(A)
    J0 = first(combinations(n, r)) # Initialize variable
    p = one(R)
    # TODO: Maybe we can speed things up here by making choices of "small" determinants?
    for J in combinations(n, r)
      p = det(A[I, data(J)])
      is_zero(p) && continue
      J0 = J
      break
    end

    result = is_odd(ind) ? result*p : result//p
    
    I = Int[i for i in 1:n if !(i in data(J0))]
    r = n - r
    ind -= 1
  end
  return result
end

