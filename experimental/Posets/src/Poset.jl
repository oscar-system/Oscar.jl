# Poset

struct Poset
  cov::Matrix{Int} # covering relations of the poset
  rel::BitMatrix # general relations
  set::BitMatrix # indicates whether a general relation was computed
  elems::Vector{Symbol} # symbols to use for the elements
end

function Base.length(P::Poset)
  return length(P.elems)
end

function Base.show(io::IO, P::Poset)
  print(io, "poset with $(length(P)) elements")
end

# PosetElem

struct PosetElem
  i::Int
  parent::Poset
end

function index(x::PosetElem)
  return x.i
end

function parent(x::PosetElem)
  return x.parent
end

function Base.show(io::IO, elem::PosetElem)
  return parent(elem).elems[index(elem)]
end

# MaximalChainsIterator

struct MaximalChainsIterator
  poset::Poset
  inplace::Bool

  function MaximalChainsIterator(P::Poset, inplace::Bool=false)
    return new(P, inplace)
  end
end

# Poset constructors

@doc raw"""
    poset(cov::Matrix{Int}, elems::Vector{<:VarName}) -> Poset

Construct a poset from covering relations `cov`, given in the form of the adjacency matrix
of a Hasse diagram. The covering relations must be given in topological order, i.e. `cov`
must be strictly upper triangular. By default the elements of the poset are named `x_i`.
"""
function poset(cov::Matrix{Int}, elems::Vector{<:VarName}=["x_$i" for i in 1:ncols(cov)])
  @req is_upper_triangular(cov, 1) "matrix must be strictly upper triangular"
  @req nrows(cov) == ncols(cov) "must be a square matrix"
  @req ncols(cov) == length(elems) "size of matrix must match number of elements"

  d = nrows(cov)
  rel = BitMatrix(!iszero(cov[i, j]) for i in 1:d, j in 1:d)
  return Poset(cov, rel, copy(rel), Symbol.(elems))
end

# PosetElem constructors

function (P::Poset)(i::Int)
  @req 1 <= i <= length(P.elems) "index out of range"
  return PosetElem(i, P)
end

function (P::Poset)(elem::VarName)
  i = findfirst(==(Symbol(elem)), P.elems)
  if isnothing(i)
    error("unknown element")
  end

  return Poset(i, P)
end

# MaximalChainsIterator constructors

@doc raw"""
    maximal_chains(P::Poset) -> MaximalChainsIterator

Returns an iterator over the maximal chains of `P`.
"""
function maximal_chains(P::Poset)
  return MaximalChainsIterator(P)
end

# PosetElem functions

function Base.:(<)(x::PosetElem, y::PosetElem)
  @req parent(x) === parent(y) "elements must belong to the same poset"
  ps = parent(x)

  # upper triangularity
  if y.i <= x.i
    return false
  end

  # linearised index
  len = ncols(ps.cov)

  # fast path
  if ps.set[x.i, y.i]
    return ps.rel[x.i, y.i]
  end

  # slow path using covering relations
  q = Int[x.i]

  while !isempty(q)
    @label outer
    n = last(q)

    # because of upper triangularity we only need to go to y.i
    for k in (n + 1):(y.i)
      if !iszero(ps.cov[n, k])
        # set the relation for all previous elements
        for m in q
          ps.set[m, k] = true
          ps.rel[m, k] = true
        end

        # we are done
        if k == y.i
          return ps.rel[x.i, y.i]
        end

        if ps.set[k, y.i]
          # check if can take the fast path
          if ps.rel[k, y.i]
            # set relation for elements in the stack
            for m in q
              ps.set[m, y.i] = true
              ps.rel[m, y.i] = true
            end
            return ps.rel[x.i, y.i]
          else
            # k is not comparable to y
            continue
          end
        end

        # add k to the stack and continue from k
        push!(q, k)
        @goto outer
      end
    end

    # we now know that n is not comparable y
    ps.set[n, y.i] = true
    ps.rel[n, y.i] = false

    # continue with previous element
    pop!(q)
  end

  return false
end

# MaximalChainsIterator implementation

function Base.IteratorSize(::Type{MaximalChainsIterator})
  return Base.SizeUnknown()
end

function Base.eltype(::Type{MaximalChainsIterator})
  return Vector{Int}
end

function Base.iterate(iter::MaximalChainsIterator, chain::Vector{Int}=[1, 1])
  j = 0
  while true
    s = pop!(chain) + 1
    if isempty(chain)
      return nothing
    end

    j = findnext(!=(0), cov[last(chain), :], s)
    if !isnothing(j)
      break
    end
  end

  while !isnothing(j)
    push!(c, j)
    j = findnext(!=(0), cov[last(chain), :], j)
  end

  if iter.inplace
    return chain, chain
  end
  return deepcopy(chain), chain
end
