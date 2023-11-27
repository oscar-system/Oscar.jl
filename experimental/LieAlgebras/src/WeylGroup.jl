# This file is based on an implementation from CoxeterGroups.jl by Ulrich Thiel (@ulthiel), Cameron Braunstein (@CameronBraunstein),
# Joel Gibson (University of Sydney, @joelgibson), and Tom Schmit (@schto223)

###############################################################################
#
#   Weyl Groups
#
###############################################################################

struct WeylGroup <: CoxeterGroup
  finite::Bool              # finite indicates whether the Weyl group is finite
  refl::ZZMatrix            # see positive_roots_and_reflections
  root_system::RootSystem   # root_system is the RootSystem from which the Weyl group was constructed

  function WeylGroup(finite::Bool, refl::ZZMatrix, root_system::RootSystem)
    return new(finite, refl, root_system)
  end
end

function weyl_group(cartan_matrix::ZZMatrix)
  return weyl_group(root_system(cartan_matrix))
end

function weyl_group(fam::Symbol, rk::Int)
  return weyl_group(root_system(fam, rk))
end

function (W::WeylGroup)(word::Vector{Int})
  return WeylGroupElem(W, word)
end

function Base.IteratorSize(::Type{WeylGroup})
  return Base.SizeUnknown()
end

function Base.isfinite(G::WeylGroup)
  return G.finite
end

function Base.one(W::WeylGroup)
  return WeylGroupElem(W, UInt8[])
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

function gen(W::WeylGroup, i::Int)
  @req 1 <= i <= ngens(W) "invalid index"
  return WeylGroupElem(W, UInt8[i])
end

function gens(W::WeylGroup)
  return [gen(W, i) for i in 1:ngens(W)]
end

function longest_element(W::WeylGroup)
  @req isfinite(W) "$W is not finite"

  rk = rank(root_system(W))
  w = -weyl_vector(root_system(W))

  word = zeros(UInt8, ncols(W.refl))
  s = 1
  i = length(word)
  while s <= rk
    if w[s] < 0
      word[i] = UInt8(s)
      reflect!(w, s)
      s = 1
      i -= 1
    else
      s += 1
    end
  end

  return WeylGroupElem(W, word)
end

function ngens(W::WeylGroup)
  return rank(root_system(W))
end

#function order(G::WeylGroup)
#end

function root_system(W::WeylGroup)
  return W.root_system
end

###############################################################################
# Weyl group elements

struct WeylGroupElem
  parent::WeylGroup     # parent group
  word::Vector{UInt8}   # short revlex normal form of the word

  function WeylGroupElem(W::WeylGroup, word::Vector{UInt8})
    return new(W, word)
  end
end

function WeylGroupElem(W::WeylGroup, word::Vector{Int})
  @req all(1 <= i <= ngens(W) for i in word) "word $word contains invalid generators"

  w = UInt8[]
  for s in Iterators.reverse(word)
    _lmul!(W.refl, w, UInt8(s))
  end

  return WeylGroupElem(W, w)
end

function Base.:(*)(x::WeylGroupElem, y::WeylGroupElem)
  @req x.parent === y.parent "$x, $y must belong to the same Weyl group"

  word = deepcopy(y.word)
  for s in Iterators.reverse(x.word)
    _lmul!(x.parent.refl, word, s)
  end

  return WeylGroupElem(x.parent, word)
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

function Base.:(==)(x::WeylGroupElem, y::WeylGroupElem)
  return x.parent === y.parent && x.word == y.word
end

function Base.deepcopy_internal(x::WeylGroupElem, dict::IdDict)
  if haskey(dict, x)
    return dict[x]
  end

  y = WeylGroupElem(x.parent, deepcopy_internal(x.word, dict))
  dict[x] = y
  return y
end

function Base.hash(x::WeylGroupElem, h::UInt)
  b = 0x80f0abce1c544784 % UInt
  h = hash(x.parent, h)
  h = hash(x.word, h)

  return xor(h, b)
end

function Base.inv(x::WeylGroupElem)
  w = UInt8[]
  sizehint!(w, length(word(x)))
  for s in word(x)
    _lmul!(parent(x).refl, w, s)
  end
  return WeylGroupElem(parent(x), w)
end

function Base.isone(x::WeylGroupElem)
  return isempty(x.word)
end

function Base.length(x::WeylGroupElem)
  return length(x.word)
end

function Base.parent(x::WeylGroupElem)
  return x.parent
end

function Base.show(io::IO, x::WeylGroupElem)
  if length(x.word) == 0
    print(io, "id")
  else
    print(io, join(Iterators.map(i -> "s$i", x.word), " * "))
  end
end

function lmul(x::WeylGroupElem, i::Integer)
  return lmul!(deepcopy(x), i)
end

function lmul!(x::WeylGroupElem, i::Integer)
  @req 1 <= i <= rank(root_system(parent(x))) "Invalid generator"
  _lmul!(parent(x).refl, word(x), UInt8(i))
  return x
end

function parent_type(::Type{WeylGroupElem})
  return WeylGroup
end

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
# Iterators

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

# ----- internal -----

function _lmul!(refl::ZZMatrix, word::Vector{T}, s::T) where {T<:Unsigned}
  insert_index = 1
  insert_letter = s

  root = s
  for i in 1:length(word)
    if word[i] == root
      deleteat!(word, i)
      return nothing
    end

    root = refl[Int(word[i]), Int(root)]
    if root == 0
      # r is no longer a minimal root, meaning we found the best insertion point
      break
    end

    # check if we have a better insertion point now. Since word[i] is a simple
    # root, if root < word[i] it must be simple.
    if root < word[i]
      insert_index = i + 1
      insert_letter = T(root)
    end
  end

  insert!(word, insert_index, insert_letter)
  return nothing
end
