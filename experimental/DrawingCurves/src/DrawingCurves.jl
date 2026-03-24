
include("IsotopyGraph.jl")
include("plotting_tikz.jl")
include("backend_msolve.jl")
include("patchwork.jl")
include("bezier.jl")

################################################################################
################################################################################
##
## Main interface
##
################################################################################
################################################################################

isotopy_graph_from_curve_inner(
  IG::_IsotopyGraph, f_in::QQMPolyRingElem, random_transform::QQMatrix,
  solver_precision::Int, ybox_tolerance::Int,
) = isotopy_graph_from_msolve(
  IG, f_in, random_transform, solver_precision, ybox_tolerance
)::Bool

function _compute_isotopy_graph(
  f_in, transform::MatrixElem, ntries::Int; solver_precision::Int=128,
  ybox_tolerance::Int=32,
)
  IG = _IsotopyGraph()
  random_transform = transform
  success = isotopy_graph_from_curve_inner(
    IG, f_in, random_transform, solver_precision, ybox_tolerance
  )
  @vprintln :DrawingCurves 2 "Was isotopy graph generated successfully? $success"
  counter = 0
  while counter < ntries && !success
    random_transform *= transform
    success = isotopy_graph_from_curve_inner(
      IG, f_in, random_transform, solver_precision, ybox_tolerance
    )
    @vprintln :DrawingCurves 2 "Was isotopy graph generated successfully? $success"
    counter += 1
  end
  scale = get_scale(IG, random_transform)
  return IG, scale, success
end

# For some background on this, see <https://xkcd.com/221/> 
function get_random_transform_matrices()
  random_transform_list = [
    matrix(QQ, [2052//2055 -111//2055; 111//2055 2052//2055]),
    matrix(QQ, [1 1//16; 0 1]),
  ]
  return random_transform_list
end

@doc raw"""
    draw_curve_tikz(filename::String, f_in; solver_precision::Int=128, graph::Bool=false, custom_edge_plot=nothing)

Takes a polynomial in two variables and constructs a plot of the resulting real
algebraic curve in TikZ.

The algorithm is based on the critical points of $f$, points where both $f$ and
its $y$-derivative vanish. If there are no such points, e.g. if $f$ is strictly
positive everywhere or consists of a bunch of parallel lines, currently an
`NotImplementedError` is thrown.

The red points are the singularities of $f$ and the blue points are the places
where the tangent is parallel to the $y$-axis. However, a random transformation
may have been applied beforehand, so the blue points may not be critical points
in the original coordinate system.

Keyword arguments:
  - `transform::Union{MatrixElem, AbstractMatrix}`: Provide a matrix such that
    transforming the curve with said transformation results in a curve with
    generic set of critical points, i.e. distinct $x$-coordinates.
    The default is to loop over a hard-coded list of transformation matrices.
  - `ntries::Int`: Apply `transform` repeatedly to reach generic critical
    points, but at most `ntries` times.
  - `graph::Bool=false`: Only draw isotopy graph. 

# Examples
```jldoctest
julia> fn = tempname();

julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> f = (x^2+y^2-1)*((x-1)^2+y^2-1)
x^4 - 2*x^3 + 2*x^2*y^2 - x^2 - 2*x*y^2 + 2*x + y^4 - y^2

julia> b = draw_curve_tikz(fn, f)
true
```
"""
function draw_curve_tikz(
  filename::String,
  f_in;
  overwrite::Bool=false,
  kwargs...,
)
  println("overwrite: $overwrite")
  if !overwrite
    @req !isfile(filename) "Output file exists"
  end
  io = open(filename, "w")
  result = draw_curve_tikz(io, f_in; kwargs...)
  close(io)
  return result
end
function draw_curve_tikz(
  io::IO,
  f_in;
  transform::Union{MatrixElem,AbstractMatrix,Nothing}=nothing,
  ntries::Int=1,
  solver_precision::Int=128,
  ybox_tolerance::Int=32,
  graph::Bool=false,
  custom_edge_plot=_draw_edge_sequence_bernstein,
  marked_components=nothing,
)
  if transform === nothing
    success = false
    for R in get_random_transform_matrices()
      success = _draw_curve_tikz(
        io, f_in, R, 30, solver_precision, ybox_tolerance, graph, custom_edge_plot,
        marked_components,
      )
      if success
        break
      end
    end
    return success
  else
    @req ntries > 0 "Number of tries needs to be positive"
    return _draw_curve_tikz(
      io, f_in, transform, ntries, solver_precision, ybox_tolerance, graph, custom_edge_plot
    )
  end
end
function _draw_curve_tikz(
  io::IO,
  f_in,
  transform::Union{MatrixElem,AbstractMatrix},
  ntries::Int,
  solver_precision::Int,
  ybox_tolerance::Int,
  graph::Bool,
  custom_edge_plot,
  marked_components,
)
  T = matrix(base_ring(f_in), transform)
  IG, scale, success = _compute_isotopy_graph(
    f_in, T, ntries; solver_precision, ybox_tolerance
  )
  if success
    if graph
      draw_graph_tikz(IG, io, marked_components)
    else
      draw_curve_tikz(IG, scale, io; custom_edge_plot=custom_edge_plot(T))
    end
  else
    # Some error?
  end
  return success
end

export draw_curve_tikz
