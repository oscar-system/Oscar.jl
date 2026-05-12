@doc raw"""
    convolution(a::Vector{T}, b::Vector{T}) where T

Return the convolution of a and b, that is the vector `c`
such that `c[i]` is the sum of `a[k] * b[i-k+1]` for all `k`.

This algorithm is optimized for cases where `a` and `b` do not contain
a particularly large number of entries. For very large `a` and `b`,
an algorithm based on the fast Fourier transform is faster.

# Examples

```jldoctest
julia> convolution([1, 2, 3, 4], [5, 6, 7])
6-element Vector{Int64}:
  5
 16
 34
 52
 45
 28
```
"""
function convolution(a::Vector{T}, b::Vector{T}) where T
  (isempty(a) || isempty(b)) && return T[]
  return [sum(a[k]*b[i-k+1] for k in max(1, i - length(b) + 1):min(i, length(a))) for i in 1:length(a) + length(b) - 1]
end
