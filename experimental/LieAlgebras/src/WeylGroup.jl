# This file is based on an implementation from CoxeterGroups.jl by Ulrich Thiel (@ulthiel), Cameron Braunstein (@CameronBraunstein),
# Joel Gibson (University of Sydney, @joelgibson), and Tom Schmit (@schto223)

###############################################################################
#
#   Weyl Groups
#
###############################################################################

@doc raw"""
    weyl_group(cartan_matrix::ZZMatrix) -> WeylGroup

Construct the Weyl group defined by the given (generalized) Cartan matrix.
"""
function weyl_group(cartan_matrix::ZZMatrix)
  return weyl_group(root_system(cartan_matrix))
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
"""
function (W::WeylGroup)(word::Vector{<:Integer}; normalize::Bool=true)
  return WeylGroupElem(W, word; normalize)
end

function Base.IteratorSize(::Type{WeylGroup})
  return Base.SizeUnknown()
end

function Base.eltype(::Type{WeylGroup})
  return WeylGroupElem
end

function Base.iterate(W::WeylGroup)
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

function coxeter_matrix(W::WeylGroup)
  return cartan_to_coxeter_matrix(cartan_matrix(root_system(W)))
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
    root_system(W::WeylGroup) -> RootSystem

Return the underlying root system of `W`.
"""
function root_system(W::WeylGroup)
  return W.root_system
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
  if haskey(dict, x)
    return dict[x]
  end

  y = parent(x)(deepcopy_internal(word(x), dict); normalize=false)
  dict[x] = y
  return y
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
    rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{WeylGroup})

Return a random element of the Weyl group. The elements are not uniformly distributed.
"""
function Base.rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{WeylGroup})
  W = rs[]
  return W(Int.(Random.randsubseq(rng, word(longest_element(W)), 2 / 3)))
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
    reduced_expressions(x::WeylGroupElem; up_to_commutation::Bool=false) -> ReducedExpressionIterator

Return an iterator over all reduced expressions of `x`.

If `up_to_commutation` is `true`, the iterator will not return an expression that only
differs from a previous one by a swap of two adjacent commuting simple reflections.

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
"""
function word(x::WeylGroupElem)
  return x.word
end

@doc raw"""
    fp_group(W::WeylGroup) -> FPGroup

Construct a group of type `FPGroup` that is isomorphic to `W`.

The `FPGroup` will be the quotient of a free group with the same rank as `W`,
where we have the natural 1-to-1 correspondence of generators, modulo the Coxeter relations of `W`.

Also see: [`isomorphism(::Type{FPGroup}, ::WeylGroup)`](@ref).
"""
function fp_group(W::WeylGroup; set_properties::Bool=true)
  return codomain(isomorphism(FPGroup, W; set_properties))
end

@doc raw"""
    isomorphism(::Type{FPGroup}, W::WeylGroup) -> Map{WeylGroup, FPGroup}

Construct an isomorphism between `W` and a group of type `FPGroup`.

The properties of the codomain group and the isomorphism are described in [`fp_group(::WeylGroup)`](@ref).
"""
function isomorphism(::Type{FPGroup}, W::WeylGroup; set_properties::Bool=true)
  R = root_system(W)
  F = free_group(rank(R))

  gcm = cartan_matrix(R)
  rels = [
    (gen(F, i) * gen(F, j))^coxeter_matrix_entry_from_cartan_matrix(gcm, i, j) for
    i in 1:rank(R) for j in i:rank(R)
  ]

  G, _ = quo(F, rels)

  if set_properties
    set_is_finite(G, is_finite(W))
    is_finite(W) && set_order(G, order(W))
  end

  iso = function (w::WeylGroupElem)
    return G([i => 1 for i in word(w)])
  end

  isoinv = function (g::FPGroupElem)
    return W(abs.(letters(g))) # TODO: check if normalize=false can be added here (probably not)
  end

  return MapFromFunc(W, G, iso, isoinv)
end

function permutation_group(W::WeylGroup; set_properties::Bool=true)
  return codomain(isomorphism(PermGroup, W; set_properties))
end

function isomorphism(::Type{PermGroup}, W::WeylGroup; set_properties::Bool=true)
  @req is_finite(W) "Weyl group is not finite"
  R = root_system(W)
  type, ordering = root_system_type_with_ordering(R)

  if length(type) != 1
    error("Not implemented (yet)")
  end
  if !issorted(ordering)
    error("Not implemented (yet)")
  end
  coxeter_type, n = only(type)
  if coxeter_type == :A
    G = symmetric_group(n + 1)

    iso = function (w::WeylGroupElem)
      reduce(*, [cperm(G, [i, i + 1]) for i in word(w)]; init=cperm(G))
    end

    isoinv = function (p::PermGroupElem)
      word = UInt8[]
      for cycle in cycles(p)
        transpositions = [
          sort([c, cycle[i + 1]]) for (i, c) in enumerate(cycle) if i < length(cycle)
        ]
        for t in transpositions
          word = reduce(
            vcat,
            [
              [i for i in t[1]:(t[2] - 1)],
              [i for i in reverse(t[1]:(t[2] - 2))],
              word,
            ],
          )
        end
      end
      return W(word) # TODO: check if normalize=false can be added here
    end
  else
    error("Not implemented (yet)")
  end

  if set_properties
    set_order(G, order(W))
  end

  return MapFromFunc(W, G, iso, isoinv)
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
    while s <= rk && weight.vec[s] > 0
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
        is_zero_entry(cartan_matrix(root_system(parent(iter.el))), s, Int(next[i - 1]))
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
      if !iszero(wt[j])
        return j
      end
    end
    return 0
  end

  for j in (di + 1):(ai - 1)
    if !iszero(wt[j])
      return j
    end
  end

  for j in (max(ai, di) + 1):rank(parent(wt))
    if is_zero_entry(cartan_matrix(root_system(wt)), ai, j)
      continue
    end

    ok = true
    for k in ai:(j - 1)
      if reflect(wt, j)[k] < 0
        ok = false
        break
      end
    end
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
"""
function weyl_orbit(wt::WeightLatticeElem)
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
