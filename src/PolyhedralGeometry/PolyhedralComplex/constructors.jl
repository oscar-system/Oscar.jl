###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

struct PolyhedralComplex{T} <: PolyhedralObject{T}
  pm_complex::Polymake.BigObject
  parent_field::Field

  PolyhedralComplex{T}(pm::Polymake.BigObject, p::Field) where {T<:scalar_types} =
    new{T}(pm, p)
  PolyhedralComplex{QQFieldElem}(pm::Polymake.BigObject) = new{QQFieldElem}(pm, QQ)
end

function polyhedral_complex(p::Polymake.BigObject)
  T, f = _detect_scalar_and_field(PolyhedralComplex, p)
  return PolyhedralComplex{T}(p, f)
end

pm_object(pc::PolyhedralComplex) = pc.pm_complex

@doc raw"""
    polyhedral_complex(::T, polyhedra, vr, far_vertices, L) where T<:scalar_types

# Arguments
- `T`: `Type` or parent `Field` of scalar to use, defaults to `QQFieldElem`.
- `polyhedra::IncidenceMatrix`: An incidence matrix; there is a 1 at position
  (i,j) if the ith polytope contains point j and 0 otherwise.
- `vr::AbstractCollection[PointVector]`: The points whose convex hulls make up
  the polyhedral complex. This matrix also contains the far vertices.
- `far_vertices::Vector{Int}`: Vector containing the indices of the rows
  corresponding to the far vertices in `vr`.
- `L::AbstractCollection[RayVector]`: Generators of the lineality space of the
  polyhedral complex.

A polyhedral complex formed from points, rays, and lineality combined into
polyhedra indicated by an incidence matrix, where the columns represent the
points and the rows represent the polyhedra.

# Examples
```jldoctest
julia> IM = incidence_matrix([[1,2,3],[1,3,4]])
2×4 IncidenceMatrix
 [1, 2, 3]
 [1, 3, 4]

julia> vr = [0 0; 1 0; 1 1; 0 1]
4×2 Matrix{Int64}:
 0  0
 1  0
 1  1
 0  1

julia> PC = polyhedral_complex(IM, vr)
Polyhedral complex in ambient dimension 2
```

Polyhedral complex with rays and lineality:
```jldoctest
julia> VR = [0 0 0; 1 0 0; 0 1 0; -1 0 0];

julia> IM = incidence_matrix([[1,2,3],[1,3,4]]);

julia> far_vertices = [2,3,4];

julia> L = [0 0 1];

julia> PC = polyhedral_complex(IM, VR, far_vertices, L)
Polyhedral complex in ambient dimension 3

julia> lineality_dim(PC)
1
```
"""
function polyhedral_complex(
  f::scalar_type_or_field,
  polyhedra::IncidenceMatrix,
  vr::AbstractCollection[PointVector],
  far_vertices::Union{Vector{Int},Nothing}=nothing,
  L::Union{AbstractCollection[RayVector],Nothing}=nothing;
  non_redundant::Bool=false,
)
  parent_field, scalar_type = _determine_parent_and_scalar(f, vr, L)
  points = homogenized_matrix(parent_field, vr, 1)
  LM = if isnothing(L) || isempty(L)
    zero_matrix(parent_field, 0, size(points, 2))
  else
    homogenized_matrix(parent_field, L, 0)
  end

  # Rays and Points are homogenized and combined and
  # If some vertices are far vertices, give them a leading 0
  if !isnothing(far_vertices)
    points[far_vertices, 1] .= 0
  end

  if non_redundant
    return PolyhedralComplex{scalar_type}(
      Polymake.fan.PolyhedralComplex{_scalar_type_to_polymake(scalar_type)}(;
        VERTICES=points, LINEALITY_SPACE=LM, MAXIMAL_CONES=polyhedra
      ),
      parent_field,
    )
  else
    return PolyhedralComplex{scalar_type}(
      Polymake.fan.PolyhedralComplex{_scalar_type_to_polymake(scalar_type)}(;
        POINTS=points, INPUT_LINEALITY=LM, INPUT_CONES=polyhedra
      ),
      parent_field,
    )
  end
end
# default scalar type: guess from input, fallback QQ
polyhedral_complex(
  polyhedra::IncidenceMatrix,
  vr::AbstractCollection[PointVector],
  far_vertices::Union{Vector{Int},Nothing}=nothing,
  L::Union{AbstractCollection[RayVector],Nothing}=nothing;
  non_redundant::Bool=false,
) = polyhedral_complex(
  _guess_fieldelem_type(vr, L),
  polyhedra,
  vr,
  far_vertices,
  L;
  non_redundant=non_redundant,
)

function polyhedral_complex(
  f::scalar_type_or_field,
  v::AbstractCollection[PointVector],
  vi::IncidenceMatrix,
  r::AbstractCollection[RayVector],
  ri::IncidenceMatrix,
  L::Union{AbstractCollection[RayVector],Nothing}=nothing;
  non_redundant::Bool=false,
)
  vr = [unhomogenized_matrix(v); unhomogenized_matrix(r)]
  far_vertices = collect((size(v, 1) + 1):size(vr, 1))
  polyhedra = hcat(vi, ri)
  return polyhedral_complex(f, polyhedra, vr, far_vertices, L; non_redundant=non_redundant)
end

@doc raw"""
    polyhedral_complex(polytopes::AbstractVector{Polyhedron{T}}) where T<:scalar_types

Assemble a polyhedral complex from a non-empty list of polyhedra.
"""
function polyhedral_complex(
  polytopes::AbstractVector{Polyhedron{T}}; non_redundant::Bool=false
) where {T<:scalar_types}
  @req length(polytopes) > 0 "list of polytopes must be non-empty"
  if non_redundant
    pmcplx = Polymake.fan.complex_from_polytopes(pm_object.(polytopes)...)
  else
    pmcplx = Polymake.fan.check_complex_objects(pm_object.(polytopes)...)
  end
  return PolyhedralComplex{T}(pmcplx, coefficient_field(iterate(polytopes)[1]))
end

# shortcut for maximal polytopes of an existing complex
function polyhedral_complex(iter::SubObjectIterator{Polyhedron{T}}) where {T<:scalar_types}
  if iter.Acc == _maximal_polyhedron && iter.Obj isa PolyhedralComplex
    return deepcopy(iter.Obj)
  end
  return polyhedral_complex(collect(iter))
end

polyhedral_complex(P::Polyhedron{T}) where {T<:scalar_types} = polyhedral_complex([P])

@doc raw"""
    polyhedral_complex(F::PolyhedralFan)

Turn a polyhedral fan into a polyhedral complex.
"""
function polyhedral_complex(F::PolyhedralFan{T}) where {T<:scalar_types}
  pmo_in = pm_object(F)
  pmo_out = Polymake.fan.PolyhedralComplex{_scalar_type_to_polymake(T)}()
  for prop in ("RAYS", "INPUT_RAYS", "FACET_NORMALS")
    if Polymake.exists(pmo_in, prop)
      Polymake.take(pmo_out, prop, embed_at_height_one(Polymake.give(pmo_in, prop), true))
    end
  end
  for prop in ("INPUT_LINEALITY", "LINEALITY_SPACE", "LINEAR_SPAN_NORMALS")
    if Polymake.exists(pmo_in, prop)
      Polymake.take(pmo_out, prop, embed_at_height_one(Polymake.give(pmo_in, prop), false))
    end
  end
  for prop in ("MAXIMAL_CONES", "INPUT_CONES", "MAXIMAL_CONES_FACETS")
    if Polymake.exists(pmo_in, prop)
      inprop = Polymake.give(pmo_in, prop)
      outprop = if inprop isa Polymake.IncidenceMatrix
        Polymake.IncidenceMatrix(Polymake.nrows(inprop), Polymake.ncols(inprop) + 1)
      else
        Polymake.spzeros(Int, Polymake.nrows(inprop), Polymake.ncols(inprop) + 1)
      end
      outprop[1:end, 2:end] = inprop
      for i in 1:Polymake.nrows(inprop)
        outprop[i, 1] = 1
      end
      Polymake.take(pmo_out, prop, outprop)
    end
  end
  if Polymake.exists(pmo_in, "MAXIMAL_CONES_LINEAR_SPAN_NORMALS")
    Polymake.take(
      pmo_out,
      "MAXIMAL_CONES_LINEAR_SPAN_NORMALS",
      Polymake.give(pmo_in, "MAXIMAL_CONES_LINEAR_SPAN_NORMALS"),
    )
  end
  return PolyhedralComplex{T}(pmo_out, coefficient_field(F))
end

###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################
function Base.show(io::IO, PC::PolyhedralComplex{T}) where {T<:scalar_types}
  try
    ad = ambient_dim(PC)
    print(io, "Polyhedral complex in ambient dimension $(ad)")
    T != QQFieldElem && print(io, " with $T type coefficients")
  catch e
    print(io, "Polyhedral complex without ambient dimension")
  end
end
