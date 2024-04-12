################################################################################
# Compositions of an integer.
#
# Copyright (C) 2021 Ulrich Thiel, ulthiel.com/math
#
# Originally taken from the JuLie [repository](https://github.com/ulthiel/JuLie)
################################################################################

################################################################################
#
#  Constructors and basic functionality
#
################################################################################

data(C::Composition) = C.c

@doc raw"""
    composition(parts::Vector{T}; check::Bool = true) where T <: IntegerUnion

Return the composition given by the integer sequence `parts` as an object of
type `Composition{T}`.

If `check` is `true` (default), it is checked whether the given sequence defines
a composition, that is, whether all elements of `parts` are positive.

# Examples
```jldoctest
julia> C = composition([6, 1, 2, 3]) # the composition 6, 1, 2, 3 of 12
[6, 1, 2, 3]

julia> C = composition(Int8[6, 1, 2, 3]) # save the elements in 8-bit integers
Int8[6, 1, 2, 3]
```
"""
function composition(parts::Vector{T}; check::Bool = true) where {T <: IntegerUnion}
  if check
    @req all(x -> x > 0, parts) "The integers must be positive"
  end
  return Composition{T}(parts)
end

function Base.show(io::IO, ::MIME"text/plain", C::Composition)
  print(io, data(C))
end

################################################################################
#
#  Array-like functionality
#
################################################################################

function Base.size(C::Composition)
  return size(data(C))
end

function Base.length(C::Composition)
  return length(data(C))
end

function Base.getindex(C::Composition, i::IntegerUnion)
  return getindex(data(C), Int(i))
end

function Base.setindex!(C::Composition, x::IntegerUnion, i::IntegerUnion)
  return setindex!(data(C), x, Int(i))
end

function Base.copy(C::Composition{T}) where T<:IntegerUnion
  return Composition{T}(copy(data(C)))
end

################################################################################
#
#  Number of compositions
#
################################################################################

@doc raw"""
    number_of_compositions(n::IntegerUnion, k::IntegerUnion)

Return the number of compositions of the non-negative integer `n` into `k >= 0`
parts.
If `n < 0` or `k < 0`, return `0`.
"""
function number_of_compositions(n::IntegerUnion, k::IntegerUnion)
  if n < 0 || k < 0
    return ZZ(0)
  end

  if n == 0
    if k == 0
      return ZZ(1)
    else
      return ZZ(0)
    end
  else
    return binomial(ZZ(n - 1), ZZ(k - 1))
  end
end

@doc raw"""
    number_of_compositions(n::IntegerUnion)

Return the number of compositions of the non-negative integer `n`.
For `n < 0`, return `0`.
"""
function number_of_compositions(n::IntegerUnion)
  if n < 0
    return ZZ(0)
  elseif n == 0
    return ZZ(1)
  else
    return ZZ(2)^(n - 1)
  end
end

################################################################################
#
#  Generating compositions with fixed number of parts
#
################################################################################

@doc raw"""
    compositions(n::IntegerUnion, k::IntegerUnion)

Return an iterator over all compositions of a non-negative integer `n` into
`k` parts, produced in lexicographically *descending* order.

By a composition of `n` into `k` parts we mean a sequence of `k` positive integers
whose sum is `n`.

# Examples
```jldoctest
julia> C = compositions(4, 2)
Iterator over the compositions of 4 into 2 parts

julia> collect(C)
3-element Vector{Composition{Int64}}:
[1, 3]
[2, 2]
[3, 1]
```
"""
function compositions(n::IntegerUnion, k::IntegerUnion)
  kk = Int(k)
  return CompositionsFixedNumParts(n, kk)
end

base(C::CompositionsFixedNumParts) = C.n
parts(C::CompositionsFixedNumParts) = C.k

Base.eltype(C::CompositionsFixedNumParts{T}) where T = Composition{T}

function Base.show(io::IO, C::CompositionsFixedNumParts)
  if get(io, :supercompact, false)
    print(io, "Iterator")
  else
    io = pretty(io)
    print(io, "Iterator over the compositions of $(base(C)) into ", ItemQuantity(parts(C), "part"))
  end
end

Base.length(C::CompositionsFixedNumParts) = BigInt(number_of_compositions(base(C), parts(C)))

function Base.iterate(C::CompositionsFixedNumParts{T}, state::Union{Nothing,Vector{T}}=nothing) where T
  w = iterate(C.weak_comp_iter, state)
  if w === nothing
    return nothing
  end

  # w[1] is a weak composition of n - k into k parts, we get a composition of n
  # into k parts by adding 1 to every entry.
  c = data(w[1])
  s = w[2]
  for i in 1:length(c)
    c[i] += 1
  end

  return composition(c, check=false), s
end

################################################################################
#
# Generating all compositions
#
################################################################################

@doc raw"""
    compositions(n::IntegerUnion)

Return an iterator over all compositions of a non-negative integer `n`.

By a composition of `n` we mean a sequence of positive integers whose sum is `n`.

# Examples
```jldoctest
julia> C = compositions(4)
Iterator over the compositions of 4

julia> collect(C)
8-element Vector{Composition{Int64}}:
 [4]
 [3, 1]
 [2, 2]
 [1, 3]
 [2, 1, 1]
 [1, 2, 1]
 [1, 1, 2]
 [1, 1, 1, 1]
```
"""
function compositions(n::IntegerUnion)
  return Compositions(n)
end

base(C::Compositions) = C.n

Base.eltype(C::Compositions{T}) where T = Composition{T}

function Base.show(io::IO, C::Compositions)
  if get(io, :supercompact, false)
    print(io, "Iterator")
  else
    io = pretty(io)
    print(io, "Iterator over the compositions of $(base(C))")
  end
end

Base.length(C::Compositions) = BigInt(number_of_compositions(base(C)))

function iterate(C::Compositions, state::Nothing = nothing)
  n = base(C)

  if n == 0
    # The only composition is the empty one
    Ck = compositions(n, 0)
  else
    # There are no compositions of n > 0 into k == 0 parts,
    # so we can start at 1
    Ck = compositions(n, 1)
  end
  next = iterate(Ck)
  return next[1], (Ck, next[2])
end

function iterate(C::Compositions{T}, state::Tuple{CompositionsFixedNumParts{T}, Vector{T}}) where T <: IntegerUnion
  n = base(C)
  Ck = state[1]
  Ck_state = state[2]

  next = iterate(Ck, Ck_state)

  if next !== nothing
    return next[1], (Ck, next[2])
  else
    k = parts(Ck)
    if k == n
      return nothing
    else
      k += 1
      Ck = compositions(n, k)
      next = iterate(Ck)
      return next[1], (Ck, next[2])
    end
  end
end
