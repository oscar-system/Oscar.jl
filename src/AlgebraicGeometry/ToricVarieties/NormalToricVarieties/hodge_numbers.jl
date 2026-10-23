############################
# Hodge numbers
############################

@doc raw"""
    hodge_number(v::NormalToricVarietyType, p::Int, q::Int)

Return the rational Hodge number $h^{p,q}$ of a complete and simplicial
toric variety `v`.

Both `p` and `q` must lie between zero and `dim(v)`, inclusive.

By Theorem 9.3.2 of [CLS11](@cite), the Hodge numbers of such a variety
vanish off the diagonal: $h^{p,q}(X_\Sigma) = 0$ for $p \neq q$. Consequently,
the diagonal Hodge numbers agree with the even Betti numbers,
\[
h^{k,k}(X_\Sigma) = b_{2k}(X_\Sigma) \, .
\]

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
  @req 0 <= p <= d && 0 <= q <= d "Hodge number indices must lie between zero and the dimension of the variety"
  p != q && return ZZ(0)
  return betti_number(v, 2 * p)
end

@doc raw"""
    hodge_numbers(v::NormalToricVarietyType) -> ZZMatrix

For a complete and simplicial toric variety `v`, return the matrix `H`
of rational Hodge numbers of `v`. The entry `H[p + 1, q + 1]` is ``h^{p,q}``.

Use [`print_hodge_diamond`](@ref) for diamond-shaped printing.
The returned matrix is cached; use `deepcopy` before modifying it.

Note that by Theorem 9.3.2 of [CLS11](@cite), the Hodge numbers of such a
toric variety vanish off the diagonal: $h^{p,q}(X_\Sigma) = 0$ for
$p \neq q$. Consequently, the diagonal Hodge numbers agree with the even
Betti numbers,
\[
h^{k,k}(X_\Sigma) = b_{2k}(X_\Sigma) \, .
\]

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2);

julia> H = hodge_numbers(P2)
[1   0   0]
[0   1   0]
[0   0   1]
```
"""
@attr ZZMatrix function hodge_numbers(v::NormalToricVarietyType)
  @req is_complete(v) && is_simplicial(v) "Hodge numbers are currently available only for complete and simplicial toric varieties"
  d = dim(v)
  H = zero_matrix(ZZ, d + 1, d + 1)
  for p in 0:d
    H[p + 1, p + 1] = hodge_number(v, p, p)
  end
  return H
end

@doc raw"""
    print_hodge_diamond(v::NormalToricVarietyType)
    print_hodge_diamond(io::IO, v::NormalToricVarietyType)

Print the Hodge numbers of a toric variety in diamond shape.

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
