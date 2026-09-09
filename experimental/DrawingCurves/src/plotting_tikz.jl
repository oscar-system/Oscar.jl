function _pt_2_string(pt::_Point, scale::QQFieldElem)
  (a, b) = scale * [pt.xcoord, pt.ycoord]
  return "($(round(Float64(a); digits=6)), $(round(Float64(b); digits=6)))"
end

function _draw_pt_tikz(
  io::IO, pt::_Point, scale::QQFieldElem; color::String="black", size=1
)
  ptstr = _pt_2_string(pt, scale)
  sizef = Float64((scale * size)//2)
  write(io, "\\fill[$color] $ptstr circle ($(sizef)pt);\n")
end

function _draw_edge_tikz(
  io::IO, pt0::_Point, pt1::_Point, scale::QQFieldElem; color::String="black", opts=nothing
)
  pt0str = _pt_2_string(pt0, scale)
  pt1str = _pt_2_string(pt1, scale)
  if opts === nothing
    write(io, "\\draw[$color] $pt0str -- $pt1str;\n")
  else
    write(io, "\\draw[$color, $opts] $pt0str -- $pt1str;\n")
  end
end

function _draw_edge_sequence_tikz(
  io::IO, pts::Vector{_Point}, scale::QQFieldElem; color::String="black"
)
  for i in 1:(length(pts) - 1)
    _draw_edge_tikz(io, pts[i], pts[i + 1], scale; color)
  end
end

function _vecint2point(v::Vector{Int})
  @req length(v) == 2 "Not a two dimensional point."
  return _Point(QQ(1, 3) * QQ(v[1]), QQ(v[2]))
end
function draw_graph_tikz(IG::_IsotopyGraph, io, marked_components)
  scale = QQ(1)
  if marked_components !== nothing
    for (mc, c) in marked_components
      for e in edges(IG.G)
        a = IG.node2coordinates[src(e)]
        b = IG.node2coordinates[dst(e)]
        mca = (1 + Float64(mc(a.xcoord, a.ycoord))) - 1
        mcb = (1 + Float64(mc(b.xcoord, b.ycoord))) - 1
        if iszero(mca) && iszero(mcb)
          pt0 = _vecint2point(IG.node2pair[src(e)])
          pt1 = _vecint2point(IG.node2pair[dst(e)])
          _draw_edge_tikz(io, pt0, pt1, scale; color=c, opts="line width=5pt, opacity=.3")
        end
      end
    end
  end
  for e in IG.nosequenceedges
    pt0 = _vecint2point(IG.node2pair[e[1]])
    pt1 = _vecint2point(IG.node2pair[e[2]])
    _draw_edge_tikz(io, pt0, pt1, scale; color="black")
  end
  for s in IG.edge4sequences
    pts = [_vecint2point(IG.node2pair[si]) for si in s]
    _draw_edge_sequence_tikz(io, pts, scale; color="black")
  end
  for p in IG.singularNodes
    pt = _vecint2point(IG.node2pair[p])
    _draw_pt_tikz(io, pt, scale; color="red", size=4)
  end
  for p in IG.ytangentNodes
    pt = _vecint2point(IG.node2pair[p])
    _draw_pt_tikz(io, pt, scale; color="blue", size=4)
  end
end
draw_graph_tikz(IG::_IsotopyGraph, io::IO) = draw_graph_tikz(IG, io, nothing)

function draw_curve_tikz(IG::_IsotopyGraph, scale, io; custom_edge_plot=nothing)
  for e in IG.nosequenceedges
    pt0 = IG.node2coordinates[e[1]]
    pt1 = IG.node2coordinates[e[2]]
    _draw_edge_tikz(io, pt0, pt1, scale; color="gray")
  end
  for s in IG.edge4sequences
    pts = [IG.node2coordinates[si] for si in s]
    _draw_edge_sequence_tikz(io, pts, scale; color="purple")
    if custom_edge_plot !== nothing
      custom_edge_plot(io, pts, scale; color="blue")
    end
  end
  for p in IG.singularNodes
    pt = IG.node2coordinates[p]
    _draw_pt_tikz(io, pt, scale; color="red", size=2)
  end
  for p in IG.ytangentNodes
    pt = IG.node2coordinates[p]
    _draw_pt_tikz(io, pt, scale; color="blue", size=2)
  end
end
