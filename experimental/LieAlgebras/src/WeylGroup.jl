###############################################################################
#
#   Weyl Groups
#
###############################################################################

struct WeylGroup <: CoxeterGroup
  # finite indicates whether the Weyl group is finite
  finite::Bool

  # see positive_roots_and_reflections
  refl::ZZMatrix

  # root_system is the RootSystem from which the Weyl group was constructed
  root_system::RootSystem
end

function weyl_group(gcm::ZZMatrix)
  return weyl_group(root_system(gcm))
end

function weyl_group(fam::Symbol, rk::Int)
  return weyl_group(root_system(fam, rk))
end

function Base.isfinite(G::WeylGroup)
  return G.finite
end

function Base.IteratorSize(::Type{WeylGroup})
  return Base.SizeUnknown()
end

function Base.one(W::WeylGroup)
  return WeylGroupElem(W, [])
end

function Base.show(io::IO, W::WeylGroup)
  print(io, "Weyl group for $(W.root_system)")
end

function coxeter_matrix(W::WeylGroup)
  return cartan_to_coxeter_matrix(W.root_system.cartan_matrix)
end

function elem(W::WeylGroup, word::Vector{Int})
  w = one(W)
  for s in Iterators.reverse(word)
    lmul!(W.refl, w.word, UInt8(s))
  end

  return w
end

function gens(W::WeylGroup)
  return [WeylGroupElem(W, [i]) for i in 1:rank(root_system(W))]
end

#function order(G::WeylGroup)
#end

function longest_element(W::WeylGroup)
  @req W.finite "$W is not finite"

  rk = rank(root_system(W))
  w = -weyl_vector(root_system(W))

  word = zeros(UInt8, ncols(W.refl))
  s = 1
  i = length(word)
  while s <= rk
    if w.vec[s] < 0
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

function root_system(W::WeylGroup)
  return W.root_system
end

###############################################################################
# Weyl group elements

struct WeylGroupElem
  # parent group
  parent::WeylGroup

  # short revlex normal form of the word
  word::Vector{UInt8}
end

function Base.:(*)(x::WeylGroupElem, y::WeylGroupElem)
  @req x.parent === y.parent "$x, $y must belong to the same Weyl group"

  word = copy(y.word)
  for s in Iterators.reverse(x.word)
    lmul!(x.parent.refl, word, s)
  end

  return WeylGroupElem(x.parent, word)
end

function Base.:(*)(x::WeylGroupElem, w::WeightLatticeElem)
  @req x.parent.root_system === w.root_system "the Weyl group of $x and the weight lattice of $w come from different root systems"

  w2 = copy(w)
  for s in Iterators.reverse(x.word)
    reflect!(w2, Int(s))
  end

  return w2
end

function Base.:(==)(x::WeylGroupElem, y::WeylGroupElem)
  return x.parent === y.parent && x.word == y.word
end

function Base.copy(x::WeylGroupElem)
  return WeylGroupElem(x.parent, deepcopy(x.word))
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
  word = UInt8[]
  sizehint!(word, length(x.word))

  for s in x.word
    lmul!(x.parent.refl, word, s)
  end

  return WeylGroupElem(x.parent, word)
end

function Base.isone(x::WeylGroupElem)
  return length(x.word) == 0
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

function Base.length(x::WeylGroupElem)
  return length(x.word)
end

function reduced_expressions(x::WeylGroupElem)
  return ReducedExpressionIterator(x)
end

function lmul!(x::WeylGroupElem, i::Int)
  lmul!(x.parent.refl, x.word, UInt8(i))
  return x
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
  # el is the Weyl group element for which we a searching reduced expressions
  el::WeylGroupElem

  # letters are the simple reflections occuring in one (hence any) reduced expression of el
  #letters::Vector{UInt8}
end

function Base.eltype(::Type{ReducedExpressionIterator})
  return Vector{UInt8}
end

function Base.IteratorSize(::Type{ReducedExpressionIterator})
  return Base.SizeUnknown()
end

function Base.iterate(iter::ReducedExpressionIterator)
  word = copy(iter.el.word)
  return word, word
end

function Base.iterate(iter::ReducedExpressionIterator, word::Vector{UInt8})
  rk = rank(iter.el.parent.root_system)

  # we need to copy word; iterate behaves differently when length is (not) known
  next = copy(word)
  weight = reflect!(weyl_vector(iter.el.parent.root_system), Int(next[1]))

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
      next[i] = UInt8(s)
      reflect!(weight, s)
      i -= 1
      s = 1
    end
  end
end

# ----- internal -----

function lmul!(refl::ZZMatrix, word::Vector{T}, s::T) where {T<:Unsigned}
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
