# Small helper function to connect to vectors of indices
function _get_edges_with_singularity(a::Vector{Int}, b::Vector{Int}, si::Int)
  na = length(a)
  nb = length(b)
  nafter = na - si
  nbefore = si - 1
  result = Vector{Int}[]
  if nbefore > 0
    for i in 1:nbefore
      push!(result, [i, i])
    end
  end
  if nafter > 0
    for i in 1:nafter
      push!(result, [na - i + 1, nb - i + 1])
    end
  end
  for i in si:(length(b) - nafter)
    push!(result, [si, i])
  end
  return result
end
_get_edges_with_singularity(na::Int, nb::Int, si::Int) =
  _get_edges_with_singularity(collect(1:na), collect(1:nb), si)

# Simply encoding a 2D-Point
struct _Point
  xcoord::QQFieldElem
  ycoord::QQFieldElem
end
Base.:*(M::QQMatrix, pt::_Point) = _Point(M * [pt.xcoord, pt.ycoord]...)
Base.:*(i::Int, pt::_Point) = _Point(i*pt.xcoord, i*pt.ycoord)
Base.:*(r::QQFieldElem, pt::_Point) = _Point(r*pt.xcoord, r*pt.ycoord)
Base.:*(r::Rational, pt::_Point) = _Point(r*pt.xcoord, r*pt.ycoord)
Base.:+(pt0::_Point, pt1::_Point) = _Point(pt0.xcoord+pt1.xcoord, pt0.ycoord+pt1.ycoord)
Base.:-(pt0::_Point, pt1::_Point) = _Point(pt0.xcoord-pt1.xcoord, pt0.ycoord-pt1.ycoord)

################################################################################
################################################################################
##
## _IsotopyGraph
## 
################################################################################
################################################################################

mutable struct _IsotopyGraph
  G::Graph{Undirected}
  node2pair::Vector{Vector{Int}}
  pair2node::Dict{Vector{Int},Int}
  singularNodes::Vector{Int}
  ytangentNodes::Vector{Int}
  node2coordinates::Vector{_Point}
  edge4sequences::Vector{Vector{Int}}
  nosequenceedges::Vector{Vector{Int}}
end

function _IsotopyGraph()
  G = Graph{Undirected}(0)
  node2pair = Vector{Vector{Int}}()
  pair2node = Dict{Vector{Int},Int}()
  singularNodes = Vector{Int}()
  ytangentNodes = Vector{Int}()
  node2coordinates = Vector{_Point}()
  edge4sequences = Vector{Vector{Int}}()
  nosequenceedges = Vector{Vector{Int}}()
  return _IsotopyGraph(
    G,
    node2pair,
    pair2node,
    singularNodes,
    ytangentNodes,
    node2coordinates,
    edge4sequences,
    nosequenceedges,
  )
end

function add_node!(IG::_IsotopyGraph, p::Vector{Int}, coords::_Point; type::Symbol=:default)
  if !(p in keys(IG.pair2node))
    add_vertex!(IG.G)
    newindex = nv(IG.G)
    IG.pair2node[p] = newindex
    push!(IG.node2pair, p)
    push!(IG.node2coordinates, coords)
    if type == :singularity
      push!(IG.singularNodes, newindex)
    elseif type == :ytangent
      push!(IG.ytangentNodes, newindex)
    end
  end
end

function get_scale(IG::_IsotopyGraph, random_transform)
  if length(IG.node2coordinates) == 0
    return 1
  end
  pt = IG.node2coordinates[1]
  (xpt, ypt) = random_transform * [pt.xcoord, pt.ycoord]
  (xmin, xmax) = (xpt, xpt)
  (ymin, ymax) = (ypt, ypt)
  for pt in IG.node2coordinates
    (xpt, ypt) = random_transform * [pt.xcoord, pt.ycoord]
    xmin = min(xmin, xpt)
    xmax = max(xmax, xpt)
    ymin = min(ymin, ypt)
    ymax = max(ymax, ypt)
  end
  scale = max(xmax - xmin, ymax - ymin)
  if scale != 0
    return 100/scale
  else
    return 1
  end
end

function add_edge_sequence!(IG::_IsotopyGraph, seq::Vector{Vector{Int}})
  for p in seq
    @assert p in keys(IG.pair2node) "$p does not specify a node."
  end
  seqi = [IG.pair2node[p] for p in seq]
  push!(IG.edge4sequences, seqi)
  for i in 1:(length(seqi) - 1)
    Oscar.add_edge!(IG.G, seqi[i], seqi[i + 1])
  end
end

function add_edge!(
  IG::_IsotopyGraph, p0::Vector{Int}, p1::Vector{Int}; nosequence::Bool=false
)
  @assert p0 in keys(IG.pair2node) "$p0 does not specify a node."
  @assert p1 in keys(IG.pair2node) "$p1 does not specify a node."
  Oscar.add_edge!(IG.G, IG.pair2node[p0], IG.pair2node[p1])
  if nosequence
    push!(IG.nosequenceedges, [IG.pair2node[p0], IG.pair2node[p1]])
  end
end

function _assemble_isotopy_graph!(
  IG::_IsotopyGraph,
  f,
  singcoords,
  stypes,
  svecs,
  sindices,
  random_transform,
  selected_precision,
  root_finder,
)
  # println("singcoords: ",singcoords)
  xcoords = [s[1] for s in singcoords]
  xmin = minimum(xcoords)
  xmax = maximum(xcoords)

  # Collect singularities
  for i in 1:length(singcoords)
    ptp = [3 * i, sindices[i]]
    (xcoord, ycoord) = singcoords[i]
    pt = _Point(xcoord, ycoord)
    add_node!(IG, ptp, random_transform * pt; type=stypes[i])
  end

  # We need this to intersect with parallels to y-axis
  Rxy = parent(f)
  (x, y) = gens(Rxy)
  K = base_ring(Rxy)
  Ry, t = polynomial_ring(K, [:y])
  projy = hom(Rxy, Ry, [0, t[1]])

  # Collect points in between "singular" places by intersecting with parallels
  # to the y-axis
  for i in 1:(length(svecs) - 1)
    x0 = xcoords[i]
    x3 = xcoords[i + 1]
    dist = x3 - x0
    x1 = x0 + dist / 3
    x2 = x3 - dist / 3
    y1 = root_finder(projy(Oscar.evaluate(f, [Rxy(x1), y])); selected_precision)
    y2 = root_finder(projy(Oscar.evaluate(f, [Rxy(x2), y])); selected_precision)
    # Apparently HomotopyContinuation.jl will sometimes produce weird results
    if length(y1) == length(y2)
      start_edges = _get_edges_with_singularity(length(svecs[i]), length(y1), sindices[i])
      end_edges = _get_edges_with_singularity(
        length(svecs[i + 1]), length(y2), sindices[i + 1]
      )
      for (a, b) in start_edges
        ee = findfirst(e -> e[2] == b, end_edges)
        pt0p = Int[3 * i, a]
        pt1p = Int[3 * i + 1, b]
        pt2p = Int[3 * i + 2, b]
        pt3p = Int[3 * i + 3, end_edges[ee][1]]
        pt0 = _Point(x0, svecs[i][a])
        pt1 = _Point(x1, y1[b])
        pt2 = _Point(x2, y2[b])
        pt3 = _Point(x3, svecs[i + 1][end_edges[ee][1]])
        add_node!(IG, pt0p, random_transform * pt0)
        add_node!(IG, pt1p, random_transform * pt1)
        add_node!(IG, pt2p, random_transform * pt2)
        add_node!(IG, pt3p, random_transform * pt3)
        add_edge_sequence!(IG, [pt0p, pt1p, pt2p, pt3p])
      end
    end
  end

  # Collect some points before the first and after the last one.
  dist = (xmax - xmin) / 20
  # Loops are very similar, to reduce duplication:
  # - coeff means direction, 1 for before, -1 for after
  # - pos means position of singularity, 1 for first, ... for last
  for (coeff, pos) in ((1, 1), (-1, length(svecs)))
    x1 = xcoords[pos]
    x0 = x1 - coeff * dist
    y0 = root_finder(projy(Oscar.evaluate(f, [Rxy(x0), y])); selected_precision)
    y1 = svecs[pos]
    edges = _get_edges_with_singularity(length(y1), length(y0), sindices[pos])
    for (a, b) in edges
      pt0p = Int[3 * pos - coeff, b]
      pt1p = Int[3 * pos, a]
      pt0 = _Point(x0, y0[b])
      pt1 = _Point(x1, y1[a])
      add_node!(IG, pt0p, random_transform * pt0)
      add_node!(IG, pt1p, random_transform * pt1)
      add_edge!(IG, pt0p, pt1p; nosequence=true)
    end
  end
end
