################################################################################
#
#  Constructors and basic functionality
#
################################################################################

data(C::WeakComposition) = C.c

@doc raw"""
    weak_composition(parts::Vector{T}; check::Bool = true) where T <: IntegerUnion

Return the weak composition given by the integer sequence `parts` as an object of type
`WeakComposition{T}`.

If `check` is `true` (default), it is checked whether the given sequence defines
a weak composition, that is, whether all elements of `parts` are non-negative.

# Examples
```jldoctest
julia> W = weak_composition([6, 0, 2, 3]) # the weak composition 6, 0, 2, 3 of 11
[6, 0, 2, 3]

julia> W = weak_composition(Int8[6, 0, 2, 3]) # save the elements in 8-bit integers
Int8[6, 0, 2, 3]
```
"""
function weak_composition(parts::Vector{T}; check::Bool = true) where {T <: IntegerUnion}
  if check
    @req all(>=(0), parts) "The integers must be non-negative"
  end
  return WeakComposition{T}(parts)
end

function Base.show(io::IO, ::MIME"text/plain", C::WeakComposition)
  c = data(C)
  if isempty(c)
    print(io, "Empty weak composition")
    return
  end
  print(io, c)
end

################################################################################
#
#  Array-like functionality
#
################################################################################

function Base.size(C::WeakComposition)
  return size(data(C))
end

function Base.length(C::WeakComposition)
  return length(data(C))
end

function Base.getindex(C::WeakComposition, i::IntegerUnion)
  return getindex(data(C), Int(i))
end

function Base.setindex!(C::WeakComposition, x::IntegerUnion, i::IntegerUnion)
  return setindex!(data(C), x, Int(i))
end

function Base.copy(C::WeakComposition)
  return weak_composition(copy(data(C)), check = false)
end

################################################################################
#
#  Generating weak compositions
#
################################################################################

@doc raw"""
    weak_compositions(n::IntegerUnion, k::IntegerUnion; inplace::Bool=false)

Return an iterator over all weak compositions of a non-negative integer `n` into
`k` parts, produced in lexicographically *descending* order.
Using a smaller integer type for `n` (e.g. `Int8`) may increase performance.

If `inplace` is `true`, the elements of the iterator may share their memory. This
means that an element returned by the iterator may be overwritten 'in place' in
the next iteration step. This may result in significantly fewer memory allocations.

By a weak composition of `n` into `k` parts we mean a sequence of `k` non-negative
integers whose sum is `n`.

# Examples
```jldoctest
julia> W = weak_compositions(3, 2)
Iterator over the weak compositions of 3 into 2 parts

julia> length(W)
4

julia> collect(W)
4-element Vector{WeakComposition{Int64}}:
 [3, 0]
 [2, 1]
 [1, 2]
 [0, 3]
```
"""
function weak_compositions(n::IntegerUnion, k::IntegerUnion; inplace::Bool=false)
  kk = Int(k)
  return WeakCompositions(n, kk, inplace)
end

# I have no idea what to call these getter functions
base(W::WeakCompositions) = W.n
parts(W::WeakCompositions) = W.k

Base.eltype(::Type{WeakCompositions{T}}) where T = WeakComposition{T}

@doc raw"""
    number_of_weak_compositions(n::IntegerUnion, k::IntegerUnion)

Return the number of weak compositions of the non-negative integer `n` into
`k >= 0` parts.
If `n < 0` or `k < 0`, return `0`.
"""
function number_of_weak_compositions(n::IntegerUnion, k::IntegerUnion)
  if n < 0 || k < 0
    return ZZ(0)
  end
  return binomial(ZZ(n) + ZZ(k) - 1, ZZ(n))
end

Base.length(W::WeakCompositions) = BigInt(number_of_weak_compositions(base(W), parts(W)))

function Base.iterate(W::WeakCompositions{T}, state::Nothing = nothing) where T
  n = base(W)
  k = parts(W)
  if n == 0
    s = zeros(T, k)
    c = W.inplace ? s : copy(s)
    return weak_composition(c, check = false), s
  end
  if k == 0
    return nothing
  end

  s = zeros(T, k)
  s[1] = n
  c = W.inplace ? s : copy(s)
  return weak_composition(c, check = false), s
end

@inline function Base.iterate(W::WeakCompositions{T}, s::Vector{T}) where T
  n = base(W)
  k = parts(W)
  if k == 0 || s[k] == n
    return nothing
  end
  if s[k - 1] == 1 && s[k] == n - 1
    s[k - 1] = 0
    s[k] = n
    c = W.inplace ? s : copy(s)
    return weak_composition(c, check = false), s
  end

  i = findlast(!iszero, view(s, 1:(k - 1)))
  @assert !isnothing(i)

  s[i] -= 1
  s[i + 1] = s[k] + 1
  if i + 1 != k
    s[k] = 0
  end
  c = W.inplace ? s : copy(s)
  return weak_composition(c, check = false), s
end

function Base.show(io::IO, W::WeakCompositions)
    print(pretty(io), "Iterator over the weak compositions of $(base(W)) into ", ItemQuantity(parts(W), "part"))
end
