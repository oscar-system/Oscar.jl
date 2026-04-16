# B = transpose(matrix(QQ, ([1 0 1; 0 1 1; 0 0 1])))
# ptmp = fc((inv(B)*gens(Rh))...)
# IG,_ = Oscar._compute_isotopy_graph(zchart(ptmp), Oscar.get_random_transform_matrices()[1], 3; solver_precision=1024, ybox_tolerance=512)

function analyse_component_position(vpos::Vector{Vector{Union{Int, QQFieldElem}}}, xfilter::Dict{Int, Vector{Int}}, c::Vector{Int})
  height = 0
  top = vpos[c[1]][2]
  bottom = vpos[c[1]][2]
  for e in c
    pt = vpos[e]
    xslice = intersect(xfilter[pt[1]], c)
    yvals = [vpos[n][2] for n in xslice]
    ymax = maximum(yvals)
    ymin = minimum(yvals)
    height = max(height, ymax-ymin)
    top = max(top, ymax)
    bottom = min(bottom, ymin)
  end
  return height, top, bottom
end

function align_component!(vpos::Vector{Vector{Union{Int, QQFieldElem}}}, xfilter::Dict{Int, Vector{Int}}, c::Vector{Int})
  xs = unique([vpos[e][1] for e in c])
  h, t, b = analyse_component_position(vpos, xfilter, c)
  targetb = 0
  for x in xs
    xslice = intersect(xfilter[x], c)
    yvals = [vpos[n][2] for n in xslice]
    ymax = maximum(yvals)
    ymin = minimum(yvals)
    if ymax != ymin
      targetb = max(targetb, ymin)
    end
  end
  targett = targetb + h
  for x in xs
    xslice = intersect(xfilter[x], c)
    yvals = [vpos[n][2] for n in xslice]
    ymax = maximum(yvals)
    ymin = minimum(yvals)
    if ymax != ymin
      for n in xfilter[x]
        pt = vpos[n]
        if pt[2] >= ymin
          if pt[2] >= ymax
            vpos[n][2] = targett + pt[2] - ymax
          else
            vpos[n][2] = targetb + pt[2] - ymin
          end
        end
      end
    end
  end
end

function align_nodes_vertically(IG::_IsotopyGraph)
  xfilter = Dict{Int, Vector{Int}}()
  for i in 1:length(IG.node2pair)
    v = IG.node2pair[i]
    if !haskey(xfilter, v[1])
      xfilter[v[1]] = [i]
    else
      push!(xfilter[v[1]], i)
    end
  end

  aligned = Vector{Vector{Union{Int, QQFieldElem}}}(IG.node2pair)
  CC = connected_components(IG.G)
  node2component = zeros(Int, nv(IG.G))
  for i in 1:length(CC)
    for e in CC[i]
      node2component[e] = i
    end
  end

  queue = collect(1:length(CC))
  donenodes = Int[]
  currenty = 0
  while length(queue) > 0
    for i in queue
      h, t, b = analyse_component_position(aligned, xfilter, CC[i])
      if h == t-b && b == currenty
        align_component!(aligned, xfilter, CC[i])
        queue = Base.setdiff(queue, i)
        append!(donenodes, CC[i])
      end
    end
    ynodes = filter(i->aligned[i][2]==currenty, 1:nv(IG.G))
    ynodes = Base.setdiff(ynodes, donenodes)
    for xn in ynodes
      x = aligned[xn][1]
      for n in Base.setdiff(xfilter[x], donenodes)
        if aligned[n][2] >= currenty
          aligned[n][2] += 1
        end
      end
    end
    currenty += 1
  end

  for n in vcat(IG.singularNodes, IG.ytangentNodes)
    nbs = neighbors(IG.G, n)
    yvals = [aligned[nb][2] for nb in nbs]
    aligned[n][2] = QQ(1,2)*(maximum(yvals) + minimum(yvals))
  end
  return aligned
end

function draw_graph_tikz_aligned(IG::_IsotopyGraph, io, marked_components)
  scale = QQ(1)
  node2pair = align_nodes_vertically(IG)
  for e in edges(IG.G)
    a = src(e)
    b = dst(e)
    type = :normal
    if a in IG.ytangentNodes
      type = :ytangent
    elseif a in IG.singularNodes
      type = :singularity
    end
    if b in IG.ytangentNodes
      type = :ytangent
      a,b = b,a
    elseif b in IG.singularNodes
      type = :singularity
      a,b = b,a
    end
    pta = _vec2point(node2pair[a])
    ptb = _vec2point(node2pair[b])
    if marked_components !== nothing
      for (mc, c) in marked_components
        a = IG.node2coordinates[src(e)]
        b = IG.node2coordinates[dst(e)]
        mca = (1 + Float64(mc(a.xcoord, a.ycoord))) - 1
        mcb = (1 + Float64(mc(b.xcoord, b.ycoord))) - 1
        if iszero(mca) && iszero(mcb)
          _draw_edge_tikz_aligned(io, pta, ptb, type, scale; color=c, opts="line width=5pt, opacity=.3")
        end
      end
    end
    _draw_edge_tikz_aligned(io, pta, ptb, type, scale)
  end
  for p in IG.singularNodes
    pt = _vec2point(node2pair[p])
    _draw_pt_tikz(io, pt, scale; color="red", size=4)
  end
  for p in IG.ytangentNodes
    pt = _vec2point(node2pair[p])
    _draw_pt_tikz(io, pt, scale; color="blue", size=4)
  end
end

function _draw_edge_tikz_aligned(
  io::IO, pt0::_Point, pt1::_Point, type::Symbol, scale::QQFieldElem; color::String="black", opts=nothing
)
  pt0str = _pt_2_string(pt0, scale)
  pt1str = _pt_2_string(pt1, scale)
  tostr = "--"
  if type != :normal
    tostr = get_to_str(type, pt0, pt1)
  end
  if opts === nothing
    write(io, "\\draw[$color] $pt0str $tostr $pt1str;\n")
  else
    write(io, "\\draw[$color, $opts] $pt0str $tostr $pt1str;\n")
  end
end

function get_to_str(type::Symbol, pta::_Point, ptb::_Point)
  result = ""
  let inang, outang
    if pta.xcoord < ptb.xcoord
      inang = 180
    else
      inang = 0
    end
    if type == :singularity
      if pta.xcoord < ptb.xcoord
        outang = 0
      else
        outang = 180
      end
    else
      if pta.ycoord < ptb.ycoord
        outang = 270
      else
        outang = 90
      end
    end
    result = "to[out=$outang, in=$inang]"
  end
  return result
end
