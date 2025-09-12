# This file is based on an implementation from CoxeterGroups.jl by Ulrich Thiel (@ulthiel), Cameron Braunstein (@CameronBraunstein),
# Joel Gibson (University of Sydney, @joelgibson), and Tom Schmit (@schto223)

###############################################################################
#
#   Weyl Groups
#
###############################################################################

@doc raw"""
    weyl_group(cartan_matrix::ZZMatrix; check::Bool=true) -> WeylGroup
    weyl_group(cartan_matrix::Matrix{<:IntegerUnion}; check::Bool=true) -> WeylGroup

Construct the Weyl group defined by the given (generalized) Cartan matrix.

If `check=true` the function will verify that `cartan_matrix` is indeed a generalized Cartan matrix.
"""
function weyl_group(cartan_matrix::ZZMatrix; kwargs...)
  return weyl_group(root_system(cartan_matrix; kwargs...))
end

function weyl_group(cartan_matrix::Matrix{<:IntegerUnion}; kwargs...)
  return weyl_group(root_system(cartan_matrix; kwargs...))
end

@doc raw"""
    weyl_group(fam::Symbol, rk::Int) -> WeylGroup

Construct the Weyl group of the given type.

The input must be a valid Cartan type, see [`is_cartan_type(::Symbol, ::Int)`](@ref).

# Examples
```jldoctest
julia> weyl_group(:A, 2)
Weyl group
  of root system of rank 2
    of type A2
```
"""
function weyl_group(fam::Symbol, rk::Int)
  return weyl_group(root_system(fam, rk))
end

@doc raw"""
    weyl_group(type::Vector{Tuple{Symbol,Int}}) -> WeylGroup
    weyl_group(type::Tuple{Symbol,Int}...) -> WeylGroup

Construct the Weyl group of the given type.

Each element of `type` must be a valid Cartan type, see [`is_cartan_type(::Symbol, ::Int)`](@ref).
The vararg version needs at least one element.

# Examples
```jldoctest
julia> weyl_group([(:G, 2), (:D, 4)])
Weyl group
  of root system of rank 6
    of type G2 x D4
```
"""
function weyl_group(type::Vector{Tuple{Symbol,Int}})
  return weyl_group(root_system(type))
end

function weyl_group(type1::Tuple{Symbol,Int}, type::Tuple{Symbol,Int}...)
  return weyl_group(root_system([type1; type...]))
end

@doc raw"""
    (W::WeylGroup)(word::Vector{<:Integer}) -> WeylGroupElem

Construct a Weyl group element from the given word.

The word must be a list of integers, where each integer is the index of a simple reflection.

If the word is known to be in short lex normal form, the normalization can be skipped by setting `normalize=false`.

# Examples
```jldoctest
julia> W = weyl_group(:A, 2);

julia> x = W([1, 2, 1])
s1 * s2 * s1
```
"""
function (W::WeylGroup)(word::Vector{<:Integer}; normalize::Bool=true)
  return WeylGroupElem(W, word; normalize)
end

@doc raw"""
    (W::WeylGroup)(sylls::AbstractVector{<:Pair{<:Integer,<:IntegerUnion}}; normalize::Bool=true) -> WeylGroupElem

Construct a Weyl group element from the given syllables.

The syllables must be a list of pairs, where each entry is the index of a simple reflection and its exponent.

If the syllables describe a word in short lex normal form, the normalization can be skipped by setting `normalize=false`.

# Examples
```jldoctest
julia> W = weyl_group(:A, 2);

julia> x = W([1 => 1, 2 => 1])
s1 * s2
```
"""
function (W::WeylGroup)(
  sylls::AbstractVector{<:Pair{<:Integer,<:IntegerUnion}}; normalize::Bool=true
)
  res = [pair.first for pair in sylls if isodd(pair.second)]
  return WeylGroupElem(W, res; normalize)
end

function Base.IteratorSize(::Type{WeylGroup})
  return Base.SizeUnknown()
end

function Base.eltype(::Type{WeylGroup})
  return WeylGroupElem
end

function Base.iterate(W::WeylGroup)
  @req is_finite(W) "currently only finite Weyl groups are supported"
  state = (weyl_vector(root_system(W)), one(W))
  return one(W), state
end

function Base.iterate(W::WeylGroup, state::WeylIteratorNoCopyState)
  state = _iterate_nocopy(state)
  if isnothing(state)
    return nothing
  end

  return deepcopy(state[2]), state
end

@doc raw"""
    is_finite(W::WeylGroup) -> Bool

Return whether `W` is finite.
"""
function is_finite(W::WeylGroup)
  return W.finite
end

@doc raw"""
    one(W::WeylGroup) -> WeylGroupElem

Return the identity element of `W`.
"""
function Base.one(W::WeylGroup)
  return W(UInt8[]; normalize=false)
end

function Base.show(io::IO, mime::MIME"text/plain", W::WeylGroup)
  @show_name(io, W)
  @show_special(io, mime, W)
  io = pretty(io)
  println(io, LowercaseOff(), "Weyl group")
  print(io, Indent(), "of ", Lowercase())
  show(io, mime, root_system(W))
  print(io, Dedent())
end

function Base.show(io::IO, W::WeylGroup)
  @show_name(io, W)
  @show_special(io, W)
  io = pretty(io)
  if is_terse(io)
    print(io, LowercaseOff(), "Weyl group")
  else
    print(io, LowercaseOff(), "Weyl group of ", Lowercase(), root_system(W))
  end
end

function elem_type(::Type{WeylGroup})
  return WeylGroupElem
end

@doc raw"""
    gen(W::WeylGroup, i::Int) -> WeylGroupElem

Return the `i`-th simple reflection (with respect to the underlying root system) of `W`.

This is a more efficient version for `gens(W)[i]`.
"""
function gen(W::WeylGroup, i::Integer)
  @req 1 <= i <= ngens(W) "invalid index"
  return W(UInt8[i]; normalize=false)
end

@doc raw"""
    gens(W::WeylGroup) -> WeylGroupElem

Return the simple reflections (with respect to the underlying root system) of `W`.

See also: [`gen(::WeylGroup, ::Int)`](@ref).
"""
function gens(W::WeylGroup)
  return [gen(W, i) for i in 1:ngens(W)]
end

@doc raw"""
    longest_element(W::WeylGroup) -> WeylGroupElem

Return the unique longest element of `W`.
This only exists if `W` is finite.
"""
function longest_element(W::WeylGroup)
  @req is_finite(W) "Weyl group is not finite"

  _, w0 = conjugate_dominant_weight_with_elem(-weyl_vector(root_system(W)))
  return w0
end

@doc raw"""
    number_of_generators(W::WeylGroup) -> Int

Return the number of generators of the `W`, i.e. the rank of the underlying root system.
"""
function number_of_generators(W::WeylGroup)
  return rank(root_system(W))
end

@doc raw"""
    order(W::WeylGroup) -> ZZRingELem
    order(::Type{T}, W::WeylGroup) where {T} -> T

Return the order of `W`.

If `W` is infinite, an `InfiniteOrderError` exception will be thrown.
Use [`is_finite(::WeylGroup)`](@ref) to check this prior to calling this function if in doubt.
"""
function order(::Type{T}, W::WeylGroup) where {T}
  if !is_finite(W)
    throw(InfiniteOrderError(W))
  end

  ord = T(1)
  for (fam, rk) in root_system_type(root_system(W))
    if fam == :A
      ord *= T(factorial(rk + 1))
    elseif fam == :B || fam == :C
      ord *= T(2^rk * factorial(rk))
    elseif fam == :D
      ord *= T(2^(rk - 1) * factorial(rk))
    elseif fam == :E
      if rk == 6
        ord *= T(51840)
      elseif rk == 7
        ord *= T(2903040)
      elseif rk == 8
        ord *= T(696729600)
      end
    elseif fam == :F
      ord *= T(1152)
    else
      ord *= T(12)
    end
  end

  return ord
end

@doc raw"""
    cartan_matrix(W::WeylGroup) -> RootSystem

Return the Cartan matrix of `W`.
"""
function cartan_matrix(W::WeylGroup)
  return cartan_matrix(root_system(W))
end

@doc raw"""
    root_system(W::WeylGroup) -> RootSystem

Return the underlying root system of `W`.
"""
function root_system(W::WeylGroup)
  return W.root_system
end

@doc raw"""
    torsion_subgroup(W::WeylGroup) -> WeylGroup, Map

Return the torsion subgroup of `W`, i.e. the subgroup of all elements of `W` with finite order.  

For a finite Weyl group, every element has finite order, so the torsion subgroup
is the entire group `W`.

If `W` is infinite, an `NotImplementedError` exception will be thrown.

# Examples
```jldoctest
julia> W = weyl_group(:A, 2);

julia> torsion_subgroup(W)
(Weyl group of root system of type A2, Identity map of W)
```
"""
function torsion_subgroup(W::WeylGroup)
  is_finite(W) || throw(
    NotImplementedError(
      :torsion_subgroup, "torsion_subgroup only implemented for finite Weyl groups"
    ),
  )
  return (W, id_hom(W))
end

###############################################################################
# Weyl group elements

function Base.:*(x::WeylGroupElem, y::WeylGroupElem)
  @req parent(x) === parent(y) "parent mismatch"

  p = deepcopy(x)
  for s in word(y)
    rmul!(p, s)
  end
  return p
end

function Base.:*(::WeylGroupElem, ::Union{RootSpaceElem,WeightLatticeElem})
  error("OSCAR only supports the right action of Weyl groups")
end

@doc raw"""
    *(r::RootSpaceElem, x::WeylGroupElem) -> RootSpaceElem
    *(w::WeightLatticeElem, x::WeylGroupElem) -> WeightLatticeElem

Return the result of acting with `x` **from the right** on `r` or `w`.

See also: [`*(::WeylGroupElem, ::Union{RootSpaceElem,WeightLatticeElem})`](@ref).
"""
function Base.:*(rw::Union{RootSpaceElem,WeightLatticeElem}, x::WeylGroupElem)
  @req root_system(parent(x)) === root_system(rw) "Incompatible root systems"

  rw2 = deepcopy(rw)
  for s in word(x)
    reflect!(rw2, Int(s))
  end

  return rw2
end

# to be removed once GroupCore is supported
function Base.:(^)(x::WeylGroupElem, n::Int)
  if n == 0
    return one(parent(x))
  elseif n < 0
    return inv(x)^(-n)
  end

  px = deepcopy(x)
  for _ in 2:n
    for s in word(x)
      rmul!(px, s)
    end
  end

  return px
end

@doc raw"""
    <(x::WeylGroupElem, y::WeylGroupElem) -> Bool

Return whether `x` is smaller than `y` with respect to the Bruhat order,
i.e., whether some (not necessarily connected) subexpression of a reduced
decomposition of `y`, is a reduced decomposition of `x`.
"""
function Base.:(<)(x::WeylGroupElem, y::WeylGroupElem)
  @req parent(x) === parent(y) "$x, $y must belong to the same Weyl group"

  if length(x) >= length(y)
    return false
  elseif isone(x)
    return true
  end

  tx = deepcopy(x)
  for i in length(y):-1:1
    b, j, _ = explain_rmul(tx, y[i])
    if !b
      deleteat!(word(tx), j)
      if isone(tx)
        return true
      end
    end

    if length(tx) > i - 1
      return false
    end
  end

  return false
end

function Base.:(==)(x::WeylGroupElem, y::WeylGroupElem)
  return parent(x) === parent(y) && word(x) == word(y)
end

function Base.deepcopy_internal(x::WeylGroupElem, dict::IdDict)
  return get!(dict, x) do
    parent(x)(deepcopy_internal(word(x), dict); normalize=false)
  end
end

@doc raw"""
    getindex(x::WeylGroupElem, i::Int) -> UInt8

Return the index of simple reflection at the `i`-th position in the normal form of `x`.
"""
function Base.getindex(x::WeylGroupElem, i::Int)
  return word(x)[i]
end

function Base.hash(x::WeylGroupElem, h::UInt)
  b = 0x80f0abce1c544784 % UInt
  h = hash(parent(x), h)
  h = hash(word(x), h)

  return xor(h, b)
end

function Base.inv(x::WeylGroupElem)
  y = parent(x)(sizehint!(UInt8[], length(x)); normalize=false)
  for s in Iterators.reverse(word(x))
    rmul!(y, s)
  end
  return y
end

@doc raw"""
    isone(x::WeylGroupElem) -> Bool

Return whether `x` is the identity element of its parent.
"""
function Base.isone(x::WeylGroupElem)
  return isempty(word(x))
end

@doc raw"""
    length(x::WeylGroupElem) -> Int

Return the length of `x`.
"""
function Base.length(x::WeylGroupElem)
  return length(word(x))
end

@doc raw"""
    parent(x::WeylGroupElem) -> WeylGroup

Return the Weyl group that `x` is an element of.
"""
function Base.parent(x::WeylGroupElem)
  return x.parent
end

@doc raw"""
    is_finite_order(x::WeylGroupElem) -> Bool

Return whether `x` is of finite order, i.e. there is a natural number `n` such that `is_one(x^n)`.
"""
function is_finite_order(x::WeylGroupElem)
  is_finite(parent(x)) && return true
  length(x) <= 1 && return true
  _, hom = geometric_representation(parent(x))
  return is_finite_order(hom(x))
end

@doc raw"""
    order(x::WeylGroupElem) -> ZZRingELem
    order(::Type{T}, x::WeylGroupElem) where {T} -> T

Return the order of `x`, i.e. the smallest natural number `n` such that `is_one(x^n)`.

If `x` is of infinite order, an `InfiniteOrderError` exception will be thrown.
Use [`is_finite_order(::WeylGroupElem)`](@ref) to check this prior to calling this function if in doubt.
"""
function order(::Type{T}, x::WeylGroupElem) where {T}
  is_one(x) && return T(1)
  length(x) == 1 && return T(2)
  _, hom = geometric_representation(parent(x))
  try
    return order(T, hom(x))
  catch e
    if e isa InfiniteOrderError
      throw(InfiniteOrderError(x))
    end
    rethrow()
  end
end

@doc raw"""
    rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{WeylGroup})

Return a pseudo-random element of the Weyl group. The elements are not uniformly distributed.
"""
function Base.rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{WeylGroup})
  W = rs[]
  if is_finite(W)
    return W(Int.(Random.randsubseq(rng, word(longest_element(W)), 2 / 3)))
  else
    m = 2^rand(1:10)
    n = rand(0:m)
    return W(rand(1:ngens(W), n))
  end
end

function Base.show(io::IO, x::WeylGroupElem)
  @show_name(io, x)
  @show_special_elem(io, x)
  if length(word(x)) == 0
    print(io, "id")
  else
    print(io, join(Iterators.map(i -> "s$i", word(x)), " * "))
  end
end

@doc raw"""
    rmul(x::WeylGroupElem, i::Integer) -> WeylGroupElem

Return the result of multiplying `x` from the right by the `i`-th simple reflection.
"""
function rmul(x::WeylGroupElem, i::Integer)
  return rmul!(deepcopy(x), i)
end

@doc raw"""
    rmul!(x::WeylGroupElem, i::Integer) -> WeylGroupElem

Multiply `x` in-place from the right by the `i`-th simple reflection, and return the result.
"""
function rmul!(x::WeylGroupElem, i::Integer)
  b, j, r = explain_rmul(x, i)
  if b
    insert!(word(x), j, r)
  else
    deleteat!(word(x), j)
  end

  return x
end

# explains what multiplication of s_i from the right will do.
# Return a tuple where the first entry is true/false, depending on whether an insertion or deletion will happen,
# the second entry is the position, and the third is the simple root.
function explain_rmul(x::WeylGroupElem, i::Integer)
  @req 1 <= i <= rank(root_system(parent(x))) "Invalid generator"

  insert_index = length(x) + 1
  insert_letter = UInt8(i)

  root = insert_letter
  for s in length(x):-1:1
    if x[s] == root
      return false, s, x[s]
    end

    root = parent(x).refl[Int(x[s]), Int(root)]
    if iszero(root)
      # r is no longer a minimal root, meaning we found the best insertion point
      return true, insert_index, insert_letter
    end

    # check if we have a better insertion point now. Since word[i] is a simple
    # root, if root < word[i] it must be simple.
    if root < x[s]
      insert_index = s
      insert_letter = UInt8(root)
    end
  end

  return true, insert_index, insert_letter
end

function parent_type(::Type{WeylGroupElem})
  return WeylGroup
end

@doc raw"""
    reduced_expressions(x::WeylGroupElem; up_to_commutation::Bool=false)

Return an iterator over all reduced expressions of `x` as `Vector{UInt8}`.

If `up_to_commutation` is `true`, the iterator will not return an expression that only
differs from a previous one by a swap of two adjacent commuting simple reflections.

The type of the iterator and the order of the elements are considered
implementation details and should not be relied upon.

# Examples
```jldoctest
julia> W = weyl_group(:A, 3);

julia> x = W([1,2,3,1]);

julia> collect(reduced_expressions(x))
3-element Vector{Vector{UInt8}}:
 [0x01, 0x02, 0x01, 0x03]
 [0x01, 0x02, 0x03, 0x01]
 [0x02, 0x01, 0x02, 0x03]

julia> collect(reduced_expressions(x; up_to_commutation=true))
2-element Vector{Vector{UInt8}}:
 [0x01, 0x02, 0x01, 0x03]
 [0x02, 0x01, 0x02, 0x03]
```
The second expression of the first iterator is not contained in the second iterator
because it only differs from the first expression by a swap of two the two commuting simple reflections `s1` and `s3`.
"""
function reduced_expressions(x::WeylGroupElem; up_to_commutation::Bool=false)
  return ReducedExpressionIterator(x, up_to_commutation)
end

@doc raw"""
    word(x::WeylGroupElem) -> Vector{UInt8}

Return `x` as a list of indices of simple reflections, in reduced form.

This function is right inverse to calling `(W::WeylGroup)(word::Vector{<:Integer})`.

# Examples

```jldoctest
julia> W = weyl_group(:A, 2);

julia> x = longest_element(W)
s1 * s2 * s1

julia> word(x)
3-element Vector{UInt8}:
 0x01
 0x02
 0x01
```
"""
function word(x::WeylGroupElem)
  return x.word
end

@doc raw"""
    letters(x::WeylGroupElem) -> Vector{UInt8}

Return `x` as a list of indices of simple reflections, in reduced form.

This function is right inverse to calling `(W::WeylGroup)(word::Vector{<:Integer})`.

# Examples

```jldoctest
julia> W = weyl_group(:A, 2);

julia> x = longest_element(W)
s1 * s2 * s1

julia> letters(x)
3-element Vector{Int64}:
 1
 2
 1
```
"""
function letters(x::WeylGroupElem)
  return Int.(x.word)
end

@doc raw"""
    syllables(x::WeylGroupElem) -> Vector{Pair{UInt8, <:IntegerUnion}}

Return `x` as a list of pairs, each entry corresponding to indices of simple reflections and its exponent.

This function is right inverse to calling `(W::WeylGroup)(sylls::Vector{Pair{UInt8, <:IntegerUnion}})`.

# Examples

```jldoctest
julia> W = weyl_group(:A, 2);

julia> x = longest_element(W)
s1 * s2 * s1

julia> syllables(x)
3-element Vector{Pair{Int64, ZZRingElem}}:
 1 => 1
 2 => 1
 1 => 1
```
"""
function syllables(x::WeylGroupElem)
  return Pair{Int,ZZRingElem}[i => 1 for i in x.word]
end

@doc raw"""
    map_word(w::WeylGroupElem, genimgs::Vector; genimgs_inv::Vector = genimgs, init = nothing)

If `init` is `nothing` and `word(w) = [`$i_1$`, ..., `$i_n$`]`,
then return the product $R_1 R_2 \cdots R_n$ with $R_j =$ `genimgs[`$i_j$`]`.
Otherwise return the product $xR_1 R_2 \cdots R_n$ where $x =$ `init`.

The length of `genimgs` must be equal to the rank of the parent of `w`.
If `w` is the trivial element, then `init` is returned if it is different
from `nothing`, and otherwise `one(genimgs[1])` is returned if `genimgs` is non-empty.
If `w` is trivial, `init` is nothing and `genimgs` is empty, an error occurs.

See also: [`map_word(::Union{FPGroupElem, SubFPGroupElem}, ::Vector)`](@ref),
[`map_word(::Union{PcGroupElem, SubPcGroupElem}, ::Vector)`](@ref).
Note that `map_word(::WeylGroupElem)` accepts the `genimgs_inv` keyword argument
for consistency with other `map_word` methods, but ignores it because the
generators of a Weyl group are always self-inverse.

# Examples
```jldoctest
julia> W = weyl_group(:B, 3); imgs = [2, 3, 5];

julia> map_word(one(W), imgs)
1

julia> map_word(W([1]), imgs)
2

julia> map_word(W([1]), imgs; init=7)
14

julia> map_word(W([1,2,1]), imgs)
12

julia> map_word(W([2,1,2]), imgs) # W([2,1,2]) == W([1,2,1])
12

julia> map_word(W([3, 2, 1, 3, 2, 3]), imgs)
2250
```
"""
function map_word(
  w::WeylGroupElem, genimgs::Vector; genimgs_inv::Vector=genimgs, init=nothing
)
  @req length(genimgs) == number_of_generators(parent(w)) begin
    "Length of vector of images does not equal rank of Weyl group"
  end
  return map_word(Int.(word(w)), genimgs; init=init)
end

@doc raw"""
    reflection(beta::RootSpaceElem) -> WeylGroupElem

Return the Weyl group element corresponding to the reflection at the hyperplane orthogonal to the root `beta`.
"""
function reflection(beta::RootSpaceElem)
  R = root_system(beta)
  W = weyl_group(R)
  rk = number_of_simple_roots(R)

  b, index_of_beta = is_root_with_index(beta)
  @req b "Not a root"
  # for a negative root we want the index of the corresponding positive root
  if index_of_beta > number_of_positive_roots(R)
    index_of_beta -= number_of_positive_roots(R)
  end

  found_simple_root = index_of_beta <= rk
  current_index = index_of_beta
  list_of_indices = Int[]
  while !found_simple_root
    for j in 1:rk
      next_index = W.refl[j, current_index]
      if !iszero(next_index) &&
        next_index < current_index
        current_index = next_index
        if current_index <= rk
          found_simple_root = true
        end
        push!(list_of_indices, j)
        break
      end
    end
  end
  return W([list_of_indices; current_index; reverse(list_of_indices)])
end

###############################################################################
# ReducedExpressionIterator

function Base.IteratorSize(::Type{ReducedExpressionIterator})
  return Base.SizeUnknown()
end

function Base.eltype(::Type{ReducedExpressionIterator})
  return Vector{UInt8}
end

function Base.iterate(iter::ReducedExpressionIterator)
  w = deepcopy(word(iter.el))
  return w, w
end

function Base.iterate(iter::ReducedExpressionIterator, word::Vector{UInt8})
  isempty(word) && return nothing

  rk = rank(root_system(parent(iter.el)))

  # we need to copy word; iterate behaves differently when length is (not) known
  next = deepcopy(word)
  weight = reflect!(weyl_vector(root_system(parent(iter.el))), Int(next[end]))

  i = length(next)
  s = rk + 1
  while true
    # search for new simple reflection to add to the word
    while s <= rk && is_positive_entry(weight, s)
      s += 1
    end

    if s == rk + 1
      i -= 1
      if i == 0
        return nothing
      elseif i == length(next)
        return next, next
      end

      # revert last reflection and continue with next one
      s = Int(next[i])
      reflect!(weight, s)
      s += 1
    else
      if iter.up_to_commutation &&
        i > 1 &&
        s < next[i - 1] &&
        is_zero_entry(cartan_matrix(parent(iter.el)), s, Int(next[i - 1]))
        s += 1
        continue
      end

      next[i] = UInt8(s)
      reflect!(weight, s)
      i += 1
      s = 1
    end
  end
end

###############################################################################
# WeylIteratorNoCopy

# Iterates over all weights in the Weyl group orbit of the dominant weight `weight`,
# or analogously over all elements in the quotient W/W_P
# The iterator returns a tuple (wt, x), such that wt*x == iter.weight;
# this choice is made to align with conjugate_dominant_weight_with_elem

function Base.IteratorSize(::Type{WeylIteratorNoCopy})
  return Base.SizeUnknown()
end

function Base.eltype(::Type{WeylIteratorNoCopy})
  return WeylIteratorNoCopyState
end

function Base.iterate(iter::WeylIteratorNoCopy)
  state = deepcopy(iter.weight), one(iter.weyl_group)
  return state, state
end

function Base.iterate(iter::WeylIteratorNoCopy, state::WeylIteratorNoCopyState)
  state = _iterate_nocopy(state)
  if isnothing(state)
    return nothing
  end
  return state, state
end

# based on [Ste01], 4.C and 4.D
function _iterate_nocopy(state::WeylIteratorNoCopyState)
  wt, path = state[1], word(state[2])

  ai = isempty(path) ? UInt8(0) : first(path)
  # compute next descendant index
  di = UInt8(0)
  while true
    di = next_descendant_index(Int(ai), Int(di), wt)
    if !iszero(di)
      break
    elseif isempty(path)
      return nothing
    elseif iszero(di)
      reflect!(wt, Int(ai))
      di = popfirst!(path)
      ai = isempty(path) ? UInt8(0) : first(path)
    end
  end

  pushfirst!(path, di)
  reflect!(wt, Int(di))
  return state
end

# based on [Ste01], 4.D
function next_descendant_index(ai::Int, di::Int, wt::WeightLatticeElem)
  if iszero(ai)
    for j in (di + 1):rank(parent(wt))
      if !is_zero_entry(wt, j)
        return j
      end
    end
    return 0
  end

  for j in (di + 1):(ai - 1)
    if !is_zero_entry(wt, j)
      return j
    end
  end

  for j in (max(ai, di) + 1):rank(parent(wt))
    if is_zero_entry(cartan_matrix(root_system(wt)), ai, j)
      continue
    end

    ok = true
    reflect!(wt, j)
    for k in ai:(j - 1)
      if is_negative_entry(wt, k)
        ok = false
        break
      end
    end
    reflect!(wt, j)
    if ok
      return j
    end
  end

  return 0
end

###############################################################################
# WeylOrbitIterator

@doc raw"""
    weyl_orbit(wt::WeightLatticeElem)

Return an iterator over the Weyl group orbit at the weight `wt`.

The type of the iterator and the order of the elements are considered
implementation details and should not be relied upon.

!!! note
    This function currently only supports finite Weyl groups.
"""
function weyl_orbit(wt::WeightLatticeElem)
  @req is_finite(weyl_group(root_system(parent(wt)))) "currently only finite Weyl groups are supported"
  return WeylOrbitIterator(wt)
end

@doc raw"""
    weyl_orbit(R::RootSystem, vec::Vector{<:Integer})

Shorthand for `weyl_orbit(WeightLatticeElem(R, vec))`.
"""
function weyl_orbit(R::RootSystem, vec::Vector{<:Integer})
  return weyl_orbit(WeightLatticeElem(R, vec))
end

@doc raw"""
    weyl_orbit(W::WeylGroup, vec::Vector{<:Integer})

Shorthand for `weyl_orbit(root_system(R), vec)`.
"""
function weyl_orbit(W::WeylGroup, vec::Vector{<:Integer})
  return weyl_orbit(root_system(W), vec)
end

function Base.IteratorSize(::Type{WeylOrbitIterator})
  return Base.IteratorSize(WeylIteratorNoCopy)
end

function Base.eltype(::Type{WeylOrbitIterator})
  return WeightLatticeElem
end

function Base.iterate(iter::WeylOrbitIterator)
  (wt, _), data = iterate(iter.nocopy)
  return deepcopy(wt), data
end

function Base.iterate(iter::WeylOrbitIterator, state::WeylIteratorNoCopyState)
  it = iterate(iter.nocopy, state)
  if isnothing(it)
    return nothing
  end

  (wt, _), state = it
  return deepcopy(wt), state
end

###############################################################################
# Misc other things

@doc raw"""
    geometric_representation(W::WeylGroup) -> ZZMatrixGroup, Map{WeylGroup, ZZMatrixGroup}

Return the geometric representation `G` of the Weyl group `W`,
together with the isomorphism `hom` from `W` to `G`.

This representation is defined by `coefficients(a) * hom(x) == coefficients(a * x)`
for all `x` in `W` and `a` a simple root of `root_system(W)`. By linear extension,
this also holds for all elements `a` in the root space of `root_system(W)`.

See [Hum90; Sect. 5.3](@cite) for more details.

# Examples
```jldoctest
julia> R = root_system(:B, 3); W = weyl_group(R);

julia> r = positive_root(R, 8)
a_1 + a_2 + 2*a_3

julia> x = W([1, 2, 1, 3])
s1 * s2 * s1 * s3

julia> G, hom = geometric_representation(W)
(Matrix group of degree 3 over ZZ, Map: W -> G)

julia> coefficients(r) * hom(x)
[1   1   0]

julia> coefficients(r * x)
[1   1   0]
```
"""
@attr Tuple{ZZMatrixGroup,MapFromFunc{WeylGroup,ZZMatrixGroup}} function geometric_representation(
  W::WeylGroup
)
  gcm = cartan_matrix(W)
  n = ngens(W)
  imgs = ZZMatrix[]
  for j in 1:n
    g = identity_matrix(ZZ, n)
    GC.@preserve g gcm begin
      for i in 1:n
        # a_i * s_j = a_i - gcm_{j,i} * a_j
        # efficient version of: g[i, j] -= gcm[j, i]
        sub!(mat_entry_ptr(g, i, j), mat_entry_ptr(gcm, j, i))
      end
    end
    push!(imgs, g)
  end

  G = matrix_group(imgs)

  iso = function (w::WeylGroupElem)
    map_word(w, gens(G); init=one(G))
  end

  return G, MapFromFunc(W, G, iso)
end

@doc raw"""
    dual_geometric_representation(W::WeylGroup) -> ZZMatrixGroup, Map{WeylGroup, ZZMatrixGroup}

Return the dual geometric representation `G` of the Weyl group `W`,
together with the isomorphism `hom` from `W` to `G`.

This representation is defined by `coefficients(w) * hom(x) == coefficients(w * x)`
for all `x` in `W` and `w` a fundamental weight of `root_system(W)`. By linear extension,
this also holds for all elements `w` in `weight_lattice(root_system(W))`.

# Examples
```jldoctest
julia> R = root_system(:B, 3); W = weyl_group(R);

julia> w = WeightLatticeElem(R, [1, 4, -3])
w_1 + 4*w_2 - 3*w_3

julia> x = W([1, 2, 1, 3])
s1 * s2 * s1 * s3

julia> G, hom = dual_geometric_representation(W)
(Matrix group of degree 3 over ZZ, Map: W -> G)

julia> coefficients(w) * hom(x)
[-4   6   -7]

julia> coefficients(w * x)
[-4   6   -7]
```
"""
@attr Tuple{ZZMatrixGroup,MapFromFunc{WeylGroup,ZZMatrixGroup}} function dual_geometric_representation(
  W::WeylGroup
)
  gcm = cartan_matrix(W)
  n = ngens(W)
  imgs = ZZMatrix[]
  for i in 1:n
    g = identity_matrix(ZZ, n)
    GC.@preserve g gcm begin
      for j in 1:n
        # efficient version of: g[i, j] -= gcm[j, i]
        sub!(mat_entry_ptr(g, i, j), mat_entry_ptr(gcm, j, i))
      end
    end
    push!(imgs, g)
  end

  G = matrix_group(imgs)

  iso = function (w::WeylGroupElem)
    map_word(w, gens(G); init=one(G))
  end

  return G, MapFromFunc(W, G, iso)
end
