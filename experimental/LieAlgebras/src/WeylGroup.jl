# This file is based on an implementation from CoxeterGroups.jl by Ulrich Thiel (@ulthiel), Cameron Braunstein (@CameronBraunstein),
# Joel Gibson (University of Sydney, @joelgibson), and Tom Schmit (@schto223)

###############################################################################
#
#   Weyl Groups
#
###############################################################################

struct WeylGroup <: AbstractAlgebra.Group
  finite::Bool              # finite indicates whether the Weyl group is finite
  refl::Matrix{UInt}        # see positive_roots_and_reflections
  root_system::RootSystem   # root_system is the RootSystem from which the Weyl group was constructed

  function WeylGroup(finite::Bool, refl::Matrix{UInt}, root_system::RootSystem)
    return new(finite, refl, root_system)
  end
end

struct WeylGroupElem <: AbstractAlgebra.GroupElem
  parent::WeylGroup     # parent group
  word::Vector{UInt8}   # short revlex normal form of the word

  function WeylGroupElem(W::WeylGroup, word::Vector{<:Integer}; normalize::Bool=true)
    if !normalize
      if word isa Vector{UInt8}
        return new(W, word)
      else
        return new(W, UInt8.(word))
      end
    end

    @req all(1 <= i <= ngens(W) for i in word) "word contains invalid generators"
    x = new(W, sizehint!(UInt8[], length(word)))
    for s in Iterators.reverse(word)
      lmul!(x, s)
    end

    return x
  end
end

const WeylIteratorNoCopyState = Tuple{WeightLatticeElem,WeylGroupElem}

@doc raw"""
    weyl_group(cartan_matrix::ZZMatrix) -> WeylGroup

Returns the Weyl group defined by a generalized Cartan matrix `cartan_matrix`.
"""
function weyl_group(cartan_matrix::ZZMatrix)
  return weyl_group(root_system(cartan_matrix))
end

@doc raw"""
    weyl_group(fam::Symbol, rk::Int) -> WeylGroup

Returns the Weyl group of the given type. See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.
"""
function weyl_group(fam::Symbol, rk::Int)
  return weyl_group(root_system(fam, rk))
end

@doc raw"""
    weyl_group(type::Vector{Tuple{Symbol,Int}}) -> WeylGroup

Returns the Weyl group of the given type. See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.
"""
function weyl_group(type::Vector{Tuple{Symbol,Int}})
  return weyl_group(root_system(type))
end

@doc raw"""
    weyl_group(type::Tuple{Symbol,Int}...) -> WeylGroup

Returns the Weyl group of the given type. See `cartan_matrix(fam::Symbol, rk::Int)` for allowed combinations.
"""
function weyl_group(type::Tuple{Symbol,Int}...)
  return weyl_group(root_system(collect(type)))
end

@doc raw"""
    (W::WeylGroup)(word::Vector{Int}) -> WeylGroupElem
"""
function (W::WeylGroup)(word::Vector{<:Integer}; normalize::Bool=true)
  return weyl_group_elem(W, word; normalize=normalize)
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
    isfinite(W::WeylGroup) -> Bool
"""
function is_finite(W::WeylGroup)
  return W.finite
end

@doc raw"""
    one(W::WeylGroup) -> WeylGroupElem
"""
function Base.one(W::WeylGroup)
  return W(UInt8[]; normalize=false)
end

function Base.show(io::IO, W::WeylGroup)
  print(io, "Weyl group for $(W.root_system)")
end

function coxeter_matrix(W::WeylGroup)
  return cartan_to_coxeter_matrix(cartan_matrix(root_system(W)))
end

function elem_type(::Type{WeylGroup})
  return WeylGroupElem
end

@doc raw"""
    gen(W::WeylGroup, i::Int) -> WeylGroupElem

Returns the `i`th simple reflection (with respect to the underlying root system) of `W`.
"""
function gen(W::WeylGroup, i::Integer)
  @req 1 <= i <= ngens(W) "invalid index"
  return W(UInt8[i]; normalize=false)
end

@doc raw"""
    gens(W::WeylGroup) -> WeylGroupElem

Returns the simple reflections (with respect to the underlying root system) of `W`.
"""
function gens(W::WeylGroup)
  return [gen(W, i) for i in 1:ngens(W)]
end

@doc raw"""
    longest_element(W::WeylGroup) -> WeylGroupElem

Returns the unique longest element of `W`.
"""
function longest_element(W::WeylGroup)
  @req is_finite(W) "$W is not finite"

  _, w0 = conjugate_dominant_weight_with_elem(-weyl_vector(root_system(W)))
  return w0
end

@doc raw"""
    number_of_generators(W::WeylGroup) -> Int

Returns the number of generators of the `W`, i.e. the rank of the underyling root system.
"""
function number_of_generators(W::WeylGroup)
  return rank(root_system(W))
end

@doc raw"""
    order(::Type{T}, W::WeylGroup) where {T} -> T

Returns the order of `W`.
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

Returns the underlying root system of `W`.
"""
function root_system(W::WeylGroup)
  return W.root_system
end

###############################################################################
# Weyl group elements

function weyl_group_elem(R::RootSystem, word::Vector{<:Integer}; normalize::Bool=true)
  return WeylGroupElem(weyl_group(R), word; normalize=normalize)
end

function weyl_group_elem(W::WeylGroup, word::Vector{<:Integer}; normalize::Bool=true)
  return WeylGroupElem(W, word; normalize=normalize)
end

function Base.:(*)(x::WeylGroupElem, y::WeylGroupElem)
  @req x.parent === y.parent "$x, $y must belong to the same Weyl group"

  p = deepcopy(y)
  for s in Iterators.reverse(word(x))
    lmul!(p, s)
  end
  return p
end

function Base.:(*)(x::WeylGroupElem, w::WeightLatticeElem)
  @req root_system(parent(x)) === root_system(w) "Incompatible root systems"

  w2 = deepcopy(w)
  for s in Iterators.reverse(x.word)
    reflect!(w2, Int(s))
  end

  return w2
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
    for s in Iterators.reverse(x.word)
      lmul!(px, s)
    end
  end

  return px
end

@doc raw"""
    Base.:(<)(x::WeylGroupElem, y::WeylGroupElem)

Returns whether `x` is smaller than `y` with respect to the Bruhat order.
"""
function Base.:(<)(x::WeylGroupElem, y::WeylGroupElem)
  @req parent(x) === parent(y) "$x, $y must belong to the same Weyl group"

  if length(x) >= length(y)
    return false
  end

  # x < y in the Bruhat order, iff some (not necessarily connected) subexpression
  # of a reduced decomposition of y, is a reduced decomposition of x
  j = length(x)
  for i in length(y):-1:1
    if word(y)[i] == word(x)[j]
      j -= 1
      if j == 0
        return true
      end
    end
  end

  return false
end

function Base.:(==)(x::WeylGroupElem, y::WeylGroupElem)
  return x.parent === y.parent && x.word == y.word
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

Returns the index of simple reflection at the `i`th position in the normal form of `x`.
"""
function Base.getindex(x::WeylGroupElem, i::Int)
  return word(x)[i]
end

function Base.hash(x::WeylGroupElem, h::UInt)
  b = 0x80f0abce1c544784 % UInt
  h = hash(x.parent, h)
  h = hash(x.word, h)

  return xor(h, b)
end

@doc raw"""
    inv(x::WeylGroupElem) -> WeylGroupElem

Returns the inverse of `x`.
"""
function Base.inv(x::WeylGroupElem)
  y = parent(x)(sizehint!(UInt8[], length(x)); normalize=false)
  for s in word(x)
    lmul!(y, s)
  end
  return y
end

@doc raw"""
    isone(x::WeylGroupElem) -> Bool

Returns whether `x` is the identity.
"""
function Base.isone(x::WeylGroupElem)
  return isempty(word(x))
end

@doc raw"""
    length(x::WeylGroupElem) -> Int

Returns the length of `x`.
"""
function Base.length(x::WeylGroupElem)
  return length(word(x))
end

@doc raw"""
    parent(x::WeylGroupElem) -> WeylGroup

Returns the Weyl group that `x` is an element of.
"""
function Base.parent(x::WeylGroupElem)
  return x.parent
end

@doc raw"""
    rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{WeylGroup})

Returns a random element of the Weyl group. The elements are not uniformally distributed.
"""
function Base.rand(rng::Random.AbstractRNG, rs::Random.SamplerTrivial{WeylGroup})
  W = rs[]
  return W(Int.(Random.randsubseq(rng, word(longest_element(W)), 2 / 3)))
end

function Base.show(io::IO, x::WeylGroupElem)
  if length(x.word) == 0
    print(io, "id")
  else
    print(io, join(Iterators.map(i -> "s$i", x.word), " * "))
  end
end

@doc raw"""
    lmul(x::WeylGroupElem, i::Integer) -> WeylGroupElem

Returns the result of multiplying `x` from the left by the `i`th simple reflection.
"""
function lmul(x::WeylGroupElem, i::Integer)
  return lmul!(deepcopy(x), i)
end

@doc raw"""
    lmul!(x::WeylGroupElem, i::Integer) -> WeylGroupElem

Returns the result of multiplying `x` in place from the left by the `i`th simple reflection.
"""
function lmul!(x::WeylGroupElem, i::Integer)
  @req 1 <= i <= rank(root_system(parent(x))) "Invalid generator"

  insert_index = 1
  insert_letter = UInt8(i)

  root = insert_letter
  for s in 1:length(x)
    if x[s] == root
      deleteat!(word(x), s)
      return x
    end

    root = parent(x).refl[Int(x[s]), Int(root)]
    if iszero(root)
      # r is no longer a minimal root, meaning we found the best insertion point
      break
    end

    # check if we have a better insertion point now. Since word[i] is a simple
    # root, if root < word[i] it must be simple.
    if root < x[s]
      insert_index = s + 1
      insert_letter = UInt8(root)
    end
  end

  insert!(word(x), insert_index, insert_letter)
  return x
end

function parent_type(::Type{WeylGroupElem})
  return WeylGroup
end

# rename to reduced decompositions ?
function reduced_expressions(x::WeylGroupElem; up_to_commutation::Bool=false)
  return ReducedExpressionIterator(x, up_to_commutation)
end

@doc raw"""
    word(x::WeylGroupElem) -> Vector{UInt8}
"""
function word(x::WeylGroupElem)
  return x.word
end

###############################################################################
# ReducedExpressionIterator

struct ReducedExpressionIterator
  el::WeylGroupElem         # the Weyl group element for which we a searching reduced expressions
  #letters::Vector{UInt8}   # letters are the simple reflections occuring in one (hence any) reduced expression of el
  up_to_commutation::Bool   # if true and say s1 and s3 commute, we only list s3*s1 and not s1*s3
end

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
  rk = rank(root_system(parent(iter.el)))

  # we need to copy word; iterate behaves differently when length is (not) known
  next = deepcopy(word)
  weight = reflect!(weyl_vector(root_system(parent(iter.el))), Int(next[1]))

  i = 1
  s = rk + 1
  while true
    # search for new simple reflection to add to the word
    while s <= rk && weight.vec[s] > 0
      s += 1
    end

    if s == rk + 1
      i += 1
      if i == length(next) + 1
        return nothing
      elseif i == 1
        return next, next
      end

      # revert last reflection and continue with next one
      s = Int(next[i])
      reflect!(weight, s)
      s += 1
    else
      if iter.up_to_commutation &&
        i < length(word) &&
        s < next[i + 1] &&
        is_zero_entry(cartan_matrix(root_system(parent(iter.el))), s, Int(next[i + 1]))
        s += 1
        continue
      end

      next[i] = UInt8(s)
      reflect!(weight, s)
      i -= 1
      s = 1
    end
  end
end

###############################################################################
# WeylIteratorNoCopy

# Iterates over all weights in the Weyl group orbit of the dominant weight `weight`,
# or analogously over all elements in the quotient W/W_P
# The iterator returns a tuple (wt, x), such that x*wt == iter.weight;
# this choice is made to align with conjugate_dominant_weight_with_elem
struct WeylIteratorNoCopy
  weight::WeightLatticeElem # dominant weight
  weyl_group::WeylGroup

  function WeylIteratorNoCopy(wt::WeightLatticeElem)
    return new(conjugate_dominant_weight(wt), weyl_group(root_system(wt)))
  end
end

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

# based on [Ste01], 4.C and 4.D
function Base.iterate(iter::WeylIteratorNoCopy, state::WeylIteratorNoCopyState)
  state = _iterate_nocopy(state)
  if isnothing(state)
    return nothing
  end
  return state, state
end

function _iterate_nocopy(state::WeylIteratorNoCopyState)
  wt, path = state[1], word(state[2])
  R = root_system(wt)

  ai = isempty(path) ? UInt8(0) : path[end]
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
      di = pop!(path)
      ai = isempty(path) ? UInt8(0) : path[end]
    end
  end

  push!(path, di)
  reflect!(wt, Int(di))
  return state
end

# based on [Ste01], 4.D
function next_descendant_index(ai::Int, di::Int, wt::WeightLatticeElem)
  if iszero(ai)
    for j in (di + 1):rank(root_system(wt))
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

  for j in (max(ai, di) + 1):rank(root_system(wt))
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

struct WeylOrbitIterator
  nocopy::WeylIteratorNoCopy

  function WeylOrbitIterator(wt::WeightLatticeElem)
    return new(WeylIteratorNoCopy(wt))
  end
end

@doc raw"""
    weyl_orbit(wt::WeightLatticeElem)

Returns an iterator over the Weyl group orbit at the weight `wt`.
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
  # wt is already a copy, so here we don't need to make one
  return wt, data
end

function Base.iterate(iter::WeylOrbitIterator, state::WeylIteratorNoCopyState)
  it = iterate(iter.nocopy, state)
  if isnothing(it)
    return nothing
  end

  (wt, _), state = it
  return deepcopy(wt), state
end
