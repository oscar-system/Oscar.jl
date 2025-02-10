const visual_supported_types = Union{PolyhedralObjectUnion,Graph,SimplicialComplex}

@doc raw"""
    visualize(P::Union{Polyhedron{T}, Cone{T}, PolyhedralFan{T}, PolyhedralComplex{T}, SubdivisionOfPoints{T}}; kwargs...) where T<:Union{FieldElem, Float64}

Visualize a polyhedral object of dimension at most four (in 3-space).
In dimensions up to 3 a usual embedding is shown.
Four-dimensional polytopes are visualized as a Schlegel diagram, which is a projection onto one of the facets; e.g., see Chapter 5 of [Zie95](@cite).

In higher dimensions there is no standard method; use projections to lower dimensions or try ideas from [GJRW10](@cite).

# Extended help

# Keyword Arguments

## Colors

Colors can be given as
- a literal `String`, e.g. `"green"`.
- a `String` of the format `"r g b"` where $r, g, b \in {0, \dots, 255}$ are integers corresponding to the R/G/B values of the color.
- a `String` of the format `"r g b"` where $r, g, b \in [0, 1]$ are decimal values corresponding to the R/G/B values of the color.

Possible arguments are:

- `FacetColor`: Filling color of the polygons.
- `EdgeColor`: Color of the boundary lines.
- `PointColor`/`VertexColor`: Color of the spheres or rectangles representing the points.

## Scaling and other gradient properties

These arguments can be given as a floating point number:

- `FacetTransparency`: Transparency factor of the polygons between 0 (opaque) and 1 (completely translucent).
- `EdgeThickness`: Scaling factor for the thickness of the boundary lines.
- `PointThickness`/VertexThickness`: Scaling factor for the size of the spheres or rectangles representing the points.

## Camera

These arguments can be given as a 3-element vector over floating point numbers:

- `ViewPoint`: Position of the camera.
- `ViewDirection`: Direction of the camera.

## Appearance and Texts

These arguments can be given as a string:

- `FacetStyle`: If set to `"hidden"`, the inner area of the polygons are not rendered at all.
- `FacetLabels`: If set to `"hidden"`, the facet labels are not displayed (in the most cases this is the default behavior). TODO
- `EdgeStyle`: If set to `"hidden"`, the boundary lines are not rendered.
- `Name`: The name of this visual object in the drawing.
- `PointLabels`/`VertexLabels`: If set to `"hidden"`, no point labels are displayed.
- `PointStyle`/`VertexStyle`: If set to `"hidden"`, neither point nor its label is rendered.
- `LabelAlignment`: Defines the alignment of the vertex labels: `"left"`, `"right"` or `"center"`.
"""
function visualize(
  P::Union{
    Polyhedron{<:Union{Float64,FieldElem}},
    Cone{<:Union{Float64,FieldElem}},
    PolyhedralFan{<:Union{Float64,FieldElem}},
    PolyhedralComplex{<:Union{Float64,FieldElem}},
    SubdivisionOfPoints{<:Union{Float64,FieldElem}},
    Graph,
    SimplicialComplex,
  };
  kwargs...,
)
  pmo = _prepare_visualization(P)
  Polymake.visual(pmo; kwargs...)
end

@doc raw"""
    visualize(P::Vector; kwargs...)

Visualize a vector of polyhedral objects `P`, i.e., all polyhedral objects in `P` in one visualization.  See [`visualize`](@ref Oscar.visualize(::Union{SimplicialComplex, Cone{<:Union{Float64, FieldElem}}, Graph, PolyhedralComplex{<:Union{Float64, FieldElem}}, PolyhedralFan{<:Union{Float64, FieldElem}}, Polyhedron, SubdivisionOfPoints{<:Union{Float64, FieldElem}}})) for which polyhedral objects and keyword arguments are possible.

# Examples
```julia
julia> p = simplex(3)
Polytope in ambient dimension 3

julia> P = [p,p+[3,0,0],p+[0,3,0],p+[0,0,3]]
4-element Vector{Polyhedron{QQFieldElem}}:
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3
 Polytope in ambient dimension 3

julia> visualize(P)

```
"""
function visualize(P::Vector; kwargs::Dict=Dict{Int,Nothing}())
  P = map(
    p -> begin
      @req p isa visual_supported_types "Can not visualize objects of type $(typeof(P))"
      _prepare_visualization(p)
    end, P)
  vis = [
    Polymake.visual(
      Polymake.Visual, P[i]; get(kwargs, i, Vector{Nothing}(undef, 0))...
    ) for i in 1:length(P)
  ]
  if isdefined(Main, :IJulia) && Main.IJulia.inited
    # this will return a visual object,
    # the visualization is then triggered by the show method
    return Polymake.call_function(:common, :compose, vis...)
  else
    # this will call visual in void context and trigger the background viewer
    Polymake.call_function(Nothing, :common, :compose, vis...)
    return nothing
  end
end

function _prepare_visualization(
  P::Union{
    Polyhedron{<:Union{Float64,FieldElem}},
    Cone{<:Union{Float64,FieldElem}},
    PolyhedralFan{<:Union{Float64,FieldElem}},
    PolyhedralComplex{<:Union{Float64,FieldElem}},
  },
)
  d = ambient_dim(P)
  b = P isa Polyhedron
  if b && d == 4
    @req is_fulldimensional(P) "Can only visualize full-dimensional $(typeof(P)) of ambient dimension $d"
  else
    @req d < 4 "Can not visualize $(typeof(P)) of ambient dimension $d. Supported range: 1 <= d <= $(3 + b)"
  end
  # polymake will by default use 0:n-1 as ray labels so we assign labels
  # starting from 1 here if there are no labels yet
  # (note: labels are mutable, i.e. they can be changed again later)
  if !Polymake.exists(pm_object(P), "RAY_LABELS")
    pm_object(P).RAY_LABELS = string.(1:(Oscar.pm_object(P).N_RAYS))
  end
  return pm_object(P)
end

function _prepare_visualization(P::SubdivisionOfPoints{<:Union{Float64,FieldElem}})
  d = ambient_dim(P)
  @req d <= 3 "Can not visualize $(typeof(P)) of ambient dimension $d. Supported range: 1 <= d <= 3"
  # polymake will by default use 0:n-1 as labels so we assign labels
  # starting from 1 here if there are no labels yet
  # (note: labels are mutable, i.e. they can be changed again later)
  if !Polymake.exists(pm_object(P), "POINT_LABELS")
    pm_object(P).POINT_LABELS = string.(1:(Oscar.pm_object(P).N_POINTS))
  end
  return pm_object(P)
end

function _prepare_visualization(sc::SimplicialComplex)
  return pm_object(sc)
end

function _prepare_visualization(
  g::Graph{T}
) where {T<:Union{Polymake.Directed,Polymake.Undirected}}
  bg = Polymake.graph.Graph{T}(; ADJACENCY=pm_object(g))
  bg.NODE_LABELS = string.(1:n_vertices(g))
  return bg
end
