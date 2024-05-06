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
    @req all(x -> x >= 0, parts) "The integers must be non-negative"
  end
  return WeakComposition{T}(parts)
end

function Base.show(io::IO, ::MIME"text/plain", C::WeakComposition)
  print(io, data(C))
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
    weak_compositions(n::IntegerUnion, k::IntegerUnion)

Return an iterator over all weak compositions of a non-negative integer `n` into
`k` parts, produced in lexicographically *descending* order.
Using a smaller integer type for `n` (e.g. `Int8`) may increase performance.

By a weak composition of `n` into `k` parts we mean a sequence of `k` non-negative
integers whose sum is `n`.

# Examples
```jldoctest
julia> W = weak_compositions(3, 2)
Iterator over the weak compositions of 3 into 2 parts

julia> length(W)
4

julia> collect(W)
4-element Vector{Oscar.WeakComposition{Int64}}:
 [3, 0]
 [2, 1]
 [1, 2]
 [0, 3]
```
"""
function weak_compositions(n::IntegerUnion, k::IntegerUnion)
  kk = Int(k)
  return WeakCompositions(n, kk)
end

# I have no idea what to call these getter functions
base(W::WeakCompositions) = W.n
parts(W::WeakCompositions) = W.k

Base.eltype(W::WeakCompositions{T}) where T = WeakComposition{T}

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
    if k > 0
      s[k] = n + 1
    end
    return (weak_composition(zeros(T, k), check = false), s)
  end
  if k == 0
    return nothing
  end

  c = zeros(T, k)
  c[1] = n
  s = zeros(T, k)
  if isone(k)
    s[1] = n + 1
    return (weak_composition(c, check = false), s)
  end

  s[1] = n - 1
  s[2] = 1
  return (weak_composition(c, check = false), s)
end

function Base.iterate(W::WeakCompositions{T}, s::Vector{T}) where T
  n = base(W)
  k = parts(W)
  if k == 0 || s[k] == n + 1
    return nothing
  end
  c = copy(s)
  if s[k] == n
    s[k] += 1
    return (weak_composition(c, check = false), s)
  end

  for i = k - 1:-1:1
    if !iszero(s[i])
      s[i] -= 1
      if i + 1 == k
        s[k] += 1
      else
        s[i + 1] = 1
        if !iszero(s[k])
          s[i + 1] += s[k]
          s[k] = 0
        end
      end
      return (weak_composition(c, check = false), s)
    end
  end
end

function Base.show(io::IO, W::WeakCompositions)
  if get(io, :supercompact, false)
    print(io, "Iterator")
  else
    io = pretty(io)
    print(io, "Iterator over the weak compositions of $(base(W)) into ", ItemQuantity(parts(W), "part"))
  end
end
