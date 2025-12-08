function _pt_2_string(pt::_Point, scale)
   (a, b) = scale*[pt.xcoord, pt.ycoord]
   return "($(Float64(a)), $(Float64(b)))";
end

function _draw_pt_tikz(io, pt::_Point, scale; color::String="black", size=1)
   ptstr = _pt_2_string(pt, scale)
   sizef = Float64((scale*size)//2)
   write(io, "\\fill[$color] $ptstr circle ($(sizef)pt);\n")
end

function _draw_edge_tikz(io, pt0::_Point, pt1::_Point, scale; color::String="black")
   pt0str = _pt_2_string(pt0, scale)
   pt1str = _pt_2_string(pt1, scale)
   write(io, "\\draw[$color] $pt0str -- $pt1str;\n")
end


function _draw_edge_sequence_tikz(io, pts::Vector{_Point}, scale; color::String="black")
   for i in 1:length(pts)-1
      _draw_edge_tikz(io, pts[i], pts[i+1], scale; color)
   end
end

function draw_graph_tikz(IG::_IsotopyGraph, io)
   for e in edges(IG.G)
      pt0 = IG.node2pair[src(e)]
      pt1 = IG.node2pair[dst(e)]
      pt0str = "($(pt0[1]/3), $(pt0[2]))"
      pt1str = "($(pt1[1]/3), $(pt1[2]))"
      write(io, "\\draw[black] $pt0str -- $pt1str;\n")
   end
   for p in IG.singularNodes
      pt = IG.node2pair[p]
      ptstr = "($(pt[1]/3), $(pt[2]))"
      write(io, "\\fill[red] $ptstr circle (2pt);\n")
   end
   for p in IG.ytangentNodes
      pt = IG.node2pair[p]
      ptstr = "($(pt[1]/3), $(pt[2]))"
      write(io, "\\fill[yellow] $ptstr circle (2pt);\n")
   end
end

function draw_curve_tikz(IG::_IsotopyGraph, scale, io; custom_edge_plot=nothing)
   for e in IG.nosequenceedges
      pt0 = IG.node2coordinates[e[1]]
      pt1 = IG.node2coordinates[e[2]]
      _draw_edge_tikz(io, pt0, pt1, scale; color="gray")
   end
   for s in IG.edge4sequences
      pts = [IG.node2coordinates[si] for si in s]
      _draw_edge_sequence_tikz(io, pts, scale; color="purple");
      if custom_edge_plot !== nothing
        custom_edge_plot(io, pts, scale; color="blue");
      end
   end
   for p in IG.singularNodes
      pt = IG.node2coordinates[p]
      _draw_pt_tikz(io, pt, scale; color="red", size=2);
   end
   for p in IG.ytangentNodes
      pt = IG.node2coordinates[p]
      _draw_pt_tikz(io, pt, scale; color="yellow", size=2);
   end
end

