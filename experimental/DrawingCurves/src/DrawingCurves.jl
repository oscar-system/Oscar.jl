
include("IsotopyGraph.jl")
include("plotting_tikz.jl")
include("backend_msolve.jl")
include("patchwork.jl")

################################################################################
################################################################################
##
## Main interface
## 
################################################################################
################################################################################

isotopy_graph_from_curve_inner(
  IG::_IsotopyGraph, f_in::QQMPolyRingElem, random_transform::QQMatrix,
  selected_precision::Int,
) = isotopy_graph_from_msolve(IG, f_in, random_transform, selected_precision)::Bool

function _compute_isotopy_graph(f_in; selected_precision::Int=128)
  # Oscar.@check !isfile(filename) "Output file exists"
  IG = _IsotopyGraph()
  # Small rotation
  base_transform = matrix(QQ, [2052//2055 -111//2055; 111//2055 2052//2055])
  # Small shearing
  random_transform = matrix(QQ, base_transform)
  # base_transform *= matrix(QQ, [1 1//16; 0 1])
  result = isotopy_graph_from_curve_inner(IG, f_in, random_transform, selected_precision)
  println("result is $result")
  while !result
    # println("Turning!")
    random_transform *= base_transform
    result = isotopy_graph_from_curve_inner(IG, f_in, random_transform, selected_precision)
    println("result is $result")
  end
  scale = get_scale(IG, random_transform)
  return IG, scale
end

@doc raw"""
    draw_curve_tikz(f_in; filename::String="curve.tikz", selected_precision::Int=128, graph::Bool=false, custom_edge_plot=nothing)

Takes a polynomial in two variables and constructs a plot of the resulting real
algebraic curve in TikZ.
"""
function draw_curve_tikz(
  f_in;
  filename::String="curve.tikz",
  selected_precision::Int=128,
  graph::Bool=false,
  custom_edge_plot=nothing,
)
  IG, scale = _compute_isotopy_graph(f_in; selected_precision)
  io = open(filename, "w")
  if graph
    draw_graph_tikz(IG, io)
  else
    draw_curve_tikz(IG, scale, io; custom_edge_plot)
  end
  close(io)
end

export draw_curve_tikz
