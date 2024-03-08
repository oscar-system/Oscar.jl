@doc raw"""
    visualize(P::Union{Polyhedron{T}, Cone{T}, PolyhedralFan{T}, PolyhedralComplex{T}, SubdivisionOfPoints{T}}) where T<:Union{FieldElem, Float64}

Visualize a polyhedral object of dimension at most four (in 3-space).
In dimensions up to 3 a usual embedding is shown.
Four-dimensional polytopes are visualized as a Schlegel diagram, which is a projection onto one of the facets; e.g., see Chapter 5 of [Zie95](@cite).

In higher dimensions there is no standard method; use projections to lower dimensions or try ideas from [GJRW10](@cite).
"""
function visualize(P::Union{Polyhedron{T}, Cone{T}, PolyhedralFan{T}, PolyhedralComplex{T}}) where T<:Union{Float64, FieldElem}
  d = ambient_dim(P)
  b = P isa Polyhedron
  if b && d == 4
    @req is_fulldimensional(P) "Can only visualize full-dimensional $(typeof(P)) of ambient dimension $d"
  else
    @req d < 4  "Can not visualize $(typeof(P)) of ambient dimension $d. Supported range: 1 <= d <= $(3 + b)"
  end
  # polymake will by default use 0:n-1 as ray labels so we assign labels
  # starting from 1 here if there are no labels yet
  # (note: labels are mutable, i.e. they can be changed again later)
  if !Polymake.exists(pm_object(P), "RAY_LABELS")
    pm_object(P).RAY_LABELS = string.(1:Oscar.pm_object(P).N_RAYS)
  end
  pmo = pm_object(P)
  Polymake.visual(pmo)
end

function visualize(P::SubdivisionOfPoints{T}) where T<:Union{FieldElem, Float64}
  d = ambient_dim(P)
  @req d <= 3 "Can not visualize $(typeof(P)) of ambient dimension $d. Supported range: 1 <= d <= 3"
  # polymake will by default use 0:n-1 as labels so we assign labels
  # starting from 1 here if there are no labels yet
  # (note: labels are mutable, i.e. they can be changed again later)
  if !Polymake.exists(pm_object(P), "POINT_LABELS")
    pm_object(P).POINT_LABELS = string.(1:Oscar.pm_object(P).N_POINTS)
  end
  pmo = pm_object(P)
  Polymake.visual(pmo)
end
