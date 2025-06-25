###############################################################################
#
# hacks related to graded rings
#

# FIXME this way `simplify` can be used when the ring is not a quotient
# it will remove all the higher degree / codimension stuff
function Oscar.simplify(x::MPolyDecRingElem)
  R = parent(x)
  n = get_attribute(R, :abstract_variety_dim)::Union{Nothing, Int}
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

function Base.getindex(x::MPolyDecRingOrQuoElem, d::FinGenAbGroupElem)
  parent(x)(homogeneous_component(x, d))
end

function Base.getindex(x::MPolyDecRingOrQuoElem, d::Int)
  R = parent(x)
  D = grading_group(R)
  getindex(x, D([d]))
end

function Base.getindex(x::MPolyDecRingOrQuoElem, degs::Vector{FinGenAbGroupElem})
  R = parent(x)
  comps = homogeneous_components(x)
  ans = typeof(x)[]
  for d in degs
    push!(ans, d in keys(comps) ? R(comps[d]) : R())
  end
  ans
end

function Base.getindex(x::MPolyDecRingOrQuoElem, I::AbstractUnitRange)
  R = parent(x)
  D = grading_group(R)
  return getindex(x, [D([n]) for n in I])
end
