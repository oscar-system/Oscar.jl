###############################################################################
#
# hacks related to graded rings
#

# FIXME this way `simplify` can be used when the ring is not a quotient
# it will remove all the higher degree / codimension stuff
function Oscar.simplify(x::MPolyDecRingElem)
  R = parent(x)
  n = get_attribute(R, :abstract_variety_dim)
  n === nothing && return x
  return sum(x[0:n])
end

# FIXME this only treat the case when the grading is ZZ
# same for `total_degree` and `getindex` below
function gradings(R::MPolyDecRingOrQuo)
  [Int(x.coeff[1]) for x in (R isa MPolyQuoRing ? base_ring(R).d : R.d)]
end

function _total_degree(x::MPolyDecRingElem)
  R = parent(x)
  x = simplify(x)
  x == 0 && return 0
  d = R isa MPolyDecRing ? R.d : R.R.d
  f = R isa MPolyDecRing ? x.f : x.f.f
  max([Int(sum(d .* degrees(t))[1]) for t in terms(f)]...)
end

function Base.getindex(x::MPolyDecRingOrQuoElem, d::GrpAbFinGenElem)
  parent(x)(homogeneous_component(x, d))
end

function Base.getindex(x::MPolyDecRingOrQuoElem, d::Int)
  R = parent(x)
  D = grading_group(R)
  getindex(x, D([d]))
end

function Base.getindex(x::MPolyDecRingOrQuoElem, degs::Vector{GrpAbFinGenElem})
  R = parent(x)
  comps = homogeneous_components(x)
  ans = typeof(x)[]
  for d in degs
    push!(ans, d in keys(comps) ? R(comps[d]) : R())
  end
  ans
end

function Base.getindex(x::MPolyDecRingOrQuoElem, I::UnitRange)
  R = parent(x)
  D = grading_group(R)
  return getindex(x, [D([n]) for n in I])
end

###############################################################################
#
# Combinatoric functions
# 

# partitions of n with at most k numbers each ≤ m
function partitions(n::Int, k::Int=n, m::Int=-1)
  ans = Partition[]
  (n > 0 && k == 0) && return ans
  if m < 0 m = n end
  n <= m && push!(ans, Partition(n > 0 ? [n] : Int[]))
  for v in Int(min(n-1,m)):-1:1
    for p in partitions(n-v, k-1, v)
      push!(ans, Partition(pushfirst!(collect(p), v)))
    end
  end
  ans
end

# make combinations work for arrays
# function combinations(I::UnitRange, k::Int) combinations(collect(I), k) end
# function combinations(l::Vector, k::Int)
#   [[l[i] for i in c] for c in combinations(length(l), k)]
# end

###############################################################################
#
# pretty printing
#

# generate a list of symbols [x₁,…,xₙ] using LaTeX / unicode for IJulia / REPL
import AbstractAlgebra.Generic: subscriptify
function _parse_symbol(symbol::String, I::UnitRange)
  isdefined(Main, :IJulia) && Main.IJulia.inited && return [symbol*"_{$i}" for i in I]
  [symbol*subscriptify(i) for i in I]
end
function _parse_symbol(symbol::String, n::Int, I::UnitRange)
  isdefined(Main, :IJulia) && Main.IJulia.inited && return [symbol*"_{$n,$i}" for i in I]
  [symbol*subscriptify(n)*","*subscriptify(i) for i in I]
end

