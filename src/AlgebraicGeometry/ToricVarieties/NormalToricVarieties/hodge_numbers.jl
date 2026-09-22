############################
# Hodge numbers
############################

function _hodge_number_matrix(d::Int)
  H = Matrix{ZZRingElem}(undef, d + 1, d + 1)
  for p in 1:(d + 1)
    for q in 1:(d + 1)
      p != q && (H[p, q] = ZZ(0))
    end
  end
  return H
end

@doc raw"""
    hodge_number(v::NormalToricVarietyType, p::Int, q::Int)

For a complete and simplicial toric variety `v`, return the  rational Hodge
number ``h^{p,q}`` of the variety.

# Examples
```jldoctest
julia> X = weighted_projective_space(NormalToricVariety, [1, 1, 2]);

julia> (hodge_number(X, 1, 1), hodge_number(X, 1, 0))
(1, 0)
```
"""
function hodge_number(v::NormalToricVarietyType, p::Int, q::Int)
  @req is_complete(v) && is_simplicial(v) "Hodge numbers are currently supported only for complete and simplicial toric varieties"
  d = dim(v)
  !(0 <= p <= d && 0 <= q <= d) && return ZZ(0)
  H = get_attribute!(() -> _hodge_number_matrix(d), v, :hodge_numbers)::Matrix{ZZRingElem}
  if !isassigned(H, p + 1, q + 1)
    H[p + 1, q + 1] = betti_number(v, 2 * p)
  end
  return deepcopy(H[p + 1, q + 1])
end

@doc raw"""
    hodge_numbers(v::NormalToricVarietyType) -> ZZMatrix

For a complete and simplicial toric variety `v`, return the matrix `H`
of rational Hodge numbers of `v`. The entry `H[p + 1, q + 1]` is ``h^{p,q}``.

Use [`print_hodge_diamond`](@ref) for diamond-shaped printing.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2);

julia> H = hodge_numbers(P2)
[1   0   0]
[0   1   0]
[0   0   1]
```
"""
function hodge_numbers(v::NormalToricVarietyType)
  @req is_complete(v) && is_simplicial(v) "Hodge numbers are currently available only for complete and simplicial toric varieties"
  d = dim(v)
  H = if has_attribute(v, :hodge_numbers)
    get_attribute(v, :hodge_numbers)::Matrix{ZZRingElem}
  else
    _hodge_number_matrix(d)
  end
  for p in 0:d
    if !isassigned(H, p + 1, p + 1)
      H[p + 1, p + 1] = betti_number(v, 2 * p)
    end
  end
  set_attribute!(v, :hodge_numbers, H)
  return matrix(ZZ, deepcopy(H))
end

@doc raw"""
    print_hodge_diamond(v::NormalToricVarietyType)
    print_hodge_diamond(io::IO, v::NormalToricVarietyType)

Print the Hodge numbers of `v` in diamond form.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2);

julia> print_hodge_diamond(P2)
  1
 0 0
0 1 0
 0 0
  1
```
"""
function print_hodge_diamond(io::IO, v::NormalToricVarietyType)
  H = hodge_numbers(v)
  nrows = size(H, 1)
  entry_width = maximum(textwidth(string(x)) for x in H)
  row_width = nrows * entry_width + nrows - 1
  for degree in 0:(2 * nrows - 2)
    p_min = max(0, degree - nrows + 1)
    p_max = min(degree, nrows - 1)
    row = join(
      (lpad(string(H[p + 1, degree - p + 1]), entry_width) for p in p_max:-1:p_min),
      " ",
    )
    println(io, " "^((row_width - textwidth(row)) >> 1), row)
  end
end

print_hodge_diamond(v::NormalToricVarietyType) = print_hodge_diamond(stdout, v)
