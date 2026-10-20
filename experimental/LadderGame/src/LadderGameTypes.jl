# These will be used to keep track of subgroup transversals
# and related coset membership testing
# (following the description of "test functions"
# in the literature)

# Required interface for tf::AbstractTestFunction:
# elements(tf): returns a vector of distinct coset representatives representing the transversal
# coset_rep(tf, x): returns an object from elements(tf) representing the coset containing x
# domain(tf): the domain of coset_rep
# _ambient_group(tf): returns from coset_rep are reparented to here (if necessary)
# set_ambient(tf, G): change the ambient group to G
# sibling_lt(tf, x, y): returns whether sibling elements x,y in elements(tf)
# satisfy the relation x<y with respect to some predetermined total ordering
# (important: this will only be meaningful for x,y having a common parent at the previous step of the subgroup ladder;
# the total ordering on elements(tf) will be determined recursively utilizing this previous step,
# and this `sibling_lt` method will be used as a tie-breaker.)
# A user can provide a custom `sibling_lt` by setting a `:sibling_order` attribute on the corresponding `LadderStep`.


abstract type AbstractTestFunction end


mutable struct GAPTransversalTest <: AbstractTestFunction
  T::Oscar.SubgroupTransversal{PermGroup, PermGroup, PermGroupElem}
  Tvec::Vector{PermGroupElem}
  domain::PermGroup
  ambient::PermGroup
  __index_cache::Dict{PermGroupElem,Int}
end

function GAPTransversalTest(G::PermGroup, H::PermGroup; ambient::PermGroup = G)
  T = right_transversal(G, H)
  Tvec = PermGroupElem[ reparent(ambient, t) for t in T ]
  return GAPTransversalTest(T, Tvec, G, ambient)
end

_coset_rep(tf::GAPTransversalTest, x::PermGroupElem) = tf.Tvec[index_of_coset(Int, tf.T, reparent(tf.domain, x))]



mutable struct FunctionTest{F} <: AbstractTestFunction
  f::F   # maps to coset rep
  Tvec::Vector{PermGroupElem}
  domain::PermGroup
  ambient::PermGroup
  __index_cache::Dict{PermGroupElem,Int}
end

function FunctionTest(Tvec::Vector{PermGroupElem}, f::F; ambient::PermGroup = parent(first(Tvec))) where F
  tf = FunctionTest{F}(f, copy(Tvec), parent(first(Tvec)), ambient)
  set_ambient!(tf, ambient)
  return tf
end

_coset_rep(tf::FunctionTest, x::PermGroupElem) = reparent(tf.ambient, tf.f(reparent(tf.domain, x)))


coset_rep(tf::AbstractTestFunction, x::PermGroupElem) = _coset_rep(tf, x)::PermGroupElem
Oscar.elements(tf::AbstractTestFunction) = tf.Tvec::Vector{PermGroupElem}
Oscar.domain(tf::AbstractTestFunction) = tf.domain::PermGroup
_ambient_group(tf::AbstractTestFunction) = tf.ambient::PermGroup
function set_ambient!(tf::AbstractTestFunction, G::PermGroup)
  map!(x -> reparent(G, x), tf.Tvec, tf.Tvec)
  tf.ambient = G
  return tf
end

function _index_lookup(tf::AbstractTestFunction)
  if isempty(tf.__index_cache)
    for (i, a) in enumerate(elements(tf))
      tf.__index_cache[a] = i
    end
  end
  return tf.__index_cache
end


function sibling_lt(tf::AbstractTestFunction, x::PermGroupElem, y::PermGroupElem)
  idx = _index_lookup(tf)
  return idx[x] < idx[y]
end

function Base.show(io::IO, tf::GAPTransversalTest)
  if is_terse(io)
    print(io, "Coset-sifting test function")
  else
    print(io, "Coset-sifting test function with ", ItemQuantity(length(elements(tf)), "element"))
  end
end

function Base.show(io::IO, tf::FunctionTest)
  if is_terse(io)
    print(io, "Combinatorial test function")
  else
    print(io, "Combinatorial test function with ", ItemQuantity(length(elements(tf)), "element"))
  end
end


mutable struct TransversalChain
  C::Vector{Tuple{PermGroup, Vector{PermGroupElem}}}
end

@inline data(T::TransversalChain) = T.C
@inline Base.size(T::TransversalChain) = size(data(T))
@inline Base.length(T::TransversalChain) = length(data(T))
@inline Base.firstindex(T::TransversalChain) = 1
@inline Base.lastindex(T::TransversalChain) = length(T)
@inline Base.getindex(T::TransversalChain, i::IntegerUnion) = getindex(data(T), Int(i))


function TransversalChain(U::PermGroup)
  U0 = U
  C = Tuple{PermGroup, Vector{PermGroupElem}}[]
  while length(U0) != 1
    s = first(moved_points(U0))
    U = stabilizer(U0, s)[1]
    T = right_transversal(U0, U)
    push!(C, (U0, T))
    U0 = U
  end
  push!(C, (U0, [one(U0)]))

  return TransversalChain(C)
end

get_transversal_chain(U::PermGroup) = TransversalChain(U)

function Base.show(io::IO, T::TransversalChain)
  if is_terse(io)
    print(io, "Transversal chain")
  else
    print(io, "Transversal chain of length ", length(T))
  end
end


@attributes mutable struct LadderStep
  is_up_step::Bool
  Aprev::PermGroup
  A::PermGroup
  test::AbstractTestFunction
  DSt::Dict{PermGroupElem, TransversalChain}
  Im::Dict{PermGroupElem, Tuple{PermGroupElem, PermGroupElem}}
  F

  function LadderStep(Aprev::PermGroup, A::PermGroup)
    if is_subset(Aprev, A)
      u = true
    elseif is_subset(A, Aprev)
      u = false
    else
      error("Aprev and A must satisfy an inclusion relation")
    end
    return new(u, Aprev, A)
  end

  function LadderStep(Aprev::PermGroup, A::PermGroup, test::AbstractTestFunction)
    if is_subset(Aprev, A)
      u = true
    elseif is_subset(A, Aprev)
      u = false
    else
      error("Aprev and A must satisfy an inclusion relation")
    end
    return new(u, Aprev, A, test)
  end
end

function Base.getproperty(S::LadderStep, s::Symbol)
  if s === :F
    return getfield(S, :F)::SubArray{LadderStep, 1, SubgroupLadder, Tuple{UnitRange{Int64}}, false}
  else
    return getfield(S, s)
  end
end

@inline _get_previous_steps(S::LadderStep) =  S.F

function sibling_lt(S::LadderStep, x::PermGroupElem, y::PermGroupElem)
  f = get_attribute(S, :sibling_order, nothing)
  f !== nothing && return f(x, y)
  return sibling_lt(S.test, x, y)   # fall back to the test function's own default
end


function Base.show(io::IO, ::MIME"text/plain", S::LadderStep)
  @show_name(io, S)
  @show_special(io, MIME"text/plain"(), S)
  dir = S.is_up_step ? "up-step" : "down-step"
  println(io, "Ladder $dir")
  io = pretty(io)
  println(io, Indent(), Lowercase(), S.Aprev)
  print(io, "-> ", Lowercase(), S.A)
  print(io, Dedent())
end

function Base.show(io::IO, S::LadderStep)
  @show_name(io, S)
  @show_special(io, S)
  if is_terse(io)
    print(io, "Ladder step")
  else
    dir = S.is_up_step ? "up" : "down"
    print(io, "Ladder $dir-step ")
    print(terse(pretty(io)), Lowercase(), S.Aprev, " -> ", Lowercase(), S.A)
  end
end



# SubgroupLadder stores Vector{LadderStep}
@attributes mutable struct SubgroupLadder <: AbstractVector{LadderStep}
  S::Vector{LadderStep}
  G::PermGroup

  function SubgroupLadder(S::Vector{LadderStep})
    @req !isempty(S) "Ladder must be nonempty"
    L = new(S, S[1].A)
    for (i, s) in enumerate(L)
      s.F = view(L, 1:i)
    end
    return L
  end
end


# When extending a subgroup ladder with a "down" step, need to ensure that we still have a strong subgroup ladder. (See Koch Theorem 6.5.2)
function _is_valid_strong_down_extension(Aprev::PermGroup, A::PermGroup, C::PermGroup)
  D = _my_intersect(A, C)
  return index(Aprev, C) == index(A, D), D
end

function extend_ladder_down!(S::Vector{LadderStep}, A::PermGroup, C::PermGroup)
  Aprev = S[end].A
  @req is_subset(A, Aprev) "A must be a subgroup of the previous group"
  is_strong, D = _is_valid_strong_down_extension(Aprev, A, C)
  @req is_strong "This does not give a valid strong subgroup ladder"
  push!(S, LadderStep(Aprev, A))
  return D
end

function extend_ladder_up!(S::Vector{LadderStep}, A::PermGroup, C::PermGroup)
  Aprev = S[end].A
  @req is_subset(Aprev, A) "A must contain the previous group"
  push!(S, LadderStep(Aprev, A))
  return C
end

function extend_ladder!(S::Vector{LadderStep}, A::PermGroup, C::PermGroup)
  Aprev = S[end].A
  is_subset(A, Aprev) && return extend_ladder_down!(S, A, C)
  is_subset(Aprev, A) && return extend_ladder_up!(S, A, C)
  error("A and the previous group must satisfy an inclusion relation")
end

function SubgroupLadder(S::Vector{PermGroup})
  @req !isempty(S) "Ladder must be nonempty"
  C = first(S)
  LS = LadderStep[LadderStep(C, C)]
  for i in 2:length(S)
    C = extend_ladder!(LS, S[i], C)
  end
  return SubgroupLadder(LS)
end

_ambient_group(L::SubgroupLadder) = L.G::PermGroup
Base.parent(S::LadderStep) = S.F.parent::SubgroupLadder
_ambient_group(S::LadderStep) = _ambient_group(parent(S))

@inline data(L::SubgroupLadder) = L.S::Vector{LadderStep}

function Base.show(io::IO, ::MIME"text/plain", L::SubgroupLadder)
  @show_name(io, L)
  @show_special(io, MIME"text/plain"(), L)
  println(io, "Subgroup ladder of length ", length(L))
  io = pretty(io)
  print(io, Indent())
  println(io, "from ", Lowercase(), _ambient_group(L))
  print(io, "to ", Lowercase(), L[end].A)
  print(io, Dedent())
end

function Base.show(io::IO, L::SubgroupLadder)
  @show_name(io, L)
  @show_special(io, L)
  if is_terse(io)
    print(io, "Subgroup ladder")
  else
    print(io, "Subgroup ladder of length ", length(L), " from ")
    print(terse(pretty(io)), Lowercase(), _ambient_group(L))
    print(io, " to ")
    print(terse(pretty(io)), Lowercase(), L[end].A)
  end
end

################################################################################
#
#  Array-like functionality
#
################################################################################

@inline Base.size(L::SubgroupLadder) =  size(data(L))
@inline Base.length(L::SubgroupLadder) = length(data(L))
@inline Base.firstindex(L::SubgroupLadder) = 1
@inline Base.lastindex(L::SubgroupLadder) = length(L)
@inline Base.getindex(L::SubgroupLadder, i::IntegerUnion) = getindex(data(L), Int(i))
@inline Base.copy(L::SubgroupLadder) = SubgroupLadder(copy(data(L)))
