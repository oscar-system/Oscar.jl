@doc raw"""
    det(
      c::AbsHyperComplex{T};
      lower_bound::Union{Int, Nothing}=has_lower_bound(c, 1) ? lower_bound(c, 1) : nothing,
      upper_bound::Union{Int, Nothing}=has_upper_bound(c, 1) ? upper_bound(c, 1) : nothing,
      direction::Symbol=:from_left_to_right
    ) where {T<:ModuleFP}

Compute the determinant of a one-dimensional complex as outlined for instance in [GelfandKapranovZelevinsky94](@cite). The keyword `direction` can be set to either `:from_left_to_right` to start from an upper bound on the non-zero modules, or `:from_right_to_left` to start from a lower bound. All modules in the complex need to be free and defined over a polynomial ring.
"""
function det(
    c::AbsHyperComplex{T};
    lower_bound::Union{Int, Nothing}=has_lower_bound(c, 1) ? lower_bound(c, 1) : nothing,
    upper_bound::Union{Int, Nothing}=has_upper_bound(c, 1) ? upper_bound(c, 1) : nothing,
    direction::Symbol=:from_left_to_right
  ) where {RET<:MPolyRingElem, T<:ModuleFP{RET}}
  @assert is_one(dim(c)) "complex must be 1-dimensional"
  dir = Oscar.direction(c, 1) == :chain ? -1 : 1
  inc = direction == :from_left_to_right ? 1 : -1
  ind = 0 # Initialize variable
  if direction == :from_left_to_right
    @assert !isnothing(lower_bound) "lower bound needs to be provided"
    ind = lower_bound::Int
  elseif direction == :from_right_to_left
    @assert !isnothing(upper_bound) "upper bound needs to be provided"
    ind = upper_bound::Int
  end

  while !can_compute_index(c, (ind,)) || is_zero(c[ind])
    ind += inc
  end
  r = ngens(c[ind])
  R = base_ring(c[ind])
  F = fraction_field(R)

  result = one(F)
  r = ngens(c[ind])
  I = collect(1:r)
  @assert c[ind] isa FreeMod "all modules must be free"
  while can_compute_index(c, (ind+inc,)) && !is_zero(c[ind+inc])
    @assert c[ind+inc] isa FreeMod "all modules must be free"
    i = ind
    if Oscar.direction(c, 1) == :chain && direction == :from_left_to_right
      i = ind + 1
    elseif Oscar.direction(c, 1) == :cochain && direction == :from_right_to_left
      i = ind - 1
    end
    A = matrix(map(c, i))
    if Oscar.direction(c, 1) == :chain && direction == :from_left_to_right
      A = transpose(A)
    elseif Oscar.direction(c, 1) == :cochain && direction == :from_right_to_left
      A = transpose(A)
    end
    m = ngens(c[ind])
    n = ngens(c[ind+inc])
    J0 = first(combinations(n, r)) # Initialize variable
    p = one(R)
    # TODO: Maybe we can speed things up here by making choices of "small" determinants?
    for J in combinations(n, r)
      p = det(A[I, data(J)])
      is_zero(p) && continue
      J0 = J
    end
    result = (is_even(ind) && direction == :from_left_to_right) || (is_odd(ind) && direction == :from_right_to_left) ? result*p : result//p
    I = Int[i for i in 1:n if !(i in data(J0))]
    r = n - r
    ind += inc
  end
  return result
end

