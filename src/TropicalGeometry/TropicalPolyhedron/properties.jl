@doc raw"""
    vertices([as::Type{T} = PointVector,] P::TropicalPolyhedron)

Returns an iterator over the vertices of `P`, that is, a minimal generating set `V`
such that `P` equals `tropical_convex_hull(V)`.

# Examples
The following computes retrieves the vertices of a tropical polytope with redundant generating set:
```jldoctest
julia> P = tropical_convex_hull(min, QQ[0 -1 0; 0 -4 -1; 0 -3 0]);

julia> vertices(P)
2-element SubObjectIterator{PointVector{TropicalSemiringElem{typeof(min)}}}:
 [(0), (-1), (0)]
 [(0), (-4), (-1)]
```
"""
function vertices(as::Type{PointVector{T}}, P::TropicalPolyhedron) where {T<:TropicalSemiringElem}
  n = n_vertices(P)

  return SubObjectIterator{as}(P, _vertex_of_polyhedron, n)
end
vertices(P::TropicalPolyhedron{M}) where {M<:MinOrMax} = vertices(PointVector{TropicalSemiringElem{M}}, P)
n_vertices(P::TropicalPolyhedron) = size(pm_object(P).VERTICES, 1)

function _vertex_of_polyhedron(::Type{PointVector{TropicalSemiringElem{M}}}, P::TropicalPolyhedron, i::Int) where {M<:MinOrMax}
  T = tropical_semiring(convention(P))

  return point_vector(
    T,
    T.(P.pm_tpolytope.VERTICES[i,:])
  )
end


function pseudovertices(as::Type{PointVector{T}}, P::TropicalPolyhedron) where {T<:TropicalSemiringElem}
  CT = pm_object(P).PSEUDOVERTEX_COARSE_COVECTORS
  ind = findall(1:size(CT, 1)) do i
    all(!=(0), CT[i,:])
  end
  n = length(ind)
  TT = tropical_semiring(convention(P))

  return SubObjectIterator{as}(
    P,
    (_,_,i)->point_vector(TT,TT.(pm_object(P).PSEUDOVERTICES[ind[i],2:end])),
    n
  )
end

function pseudovertices(as::Type{PointVector{T}}, P::TropicalPointConfiguration) where {T<:TropicalSemiringElem}
  n = n_pseudovertices(P)

  return SubObjectIterator{as}(P, _pseudovertices, n)
end
pseudovertices(P::Union{TropicalPolyhedron{M},TropicalPointConfiguration{M}}) where {M<:MinOrMax} = pseudovertices(PointVector{TropicalSemiringElem{M}}, P)
n_pseudovertices(P::Union{TropicalPolyhedron,TropicalPointConfiguration}) = pm_object(P).PSEUDOVERTICES |> size |> first

function _pseudovertices(::Type{PointVector{TropicalSemiringElem{M}}}, P::Union{TropicalPolyhedron,TropicalPointConfiguration}, i::Int) where {M<:MinOrMax}
  T = tropical_semiring(convention(P))

  return point_vector(
    T,
    T.(pm_object(P).PSEUDOVERTICES[i,2:end])
  )
end

function points(as::Type{PointVector{T}}, P::TropicalPointConfiguration) where {T<:TropicalSemiringElem}
  n = n_points(P)

  return SubObjectIterator{as}(P, _points, n)
end
points(P::TropicalPointConfiguration{M}) where {M<:MinOrMax} = points(PointVector{TropicalSemiringElem{M}}, P)
n_points(P::TropicalPointConfiguration) = pm_object(P).POINTS |> size |> first

function _points(::Type{PointVector{TropicalSemiringElem{M}}}, P::TropicalPointConfiguration, i::Int) where {M<:MinOrMax}
  T = tropical_semiring(convention(P))

  return point_vector(
    T,
    T.(P.pm_tpolytope.POINTS[i,:])
  )
end

dim(P::TropicalPolyhedron) = covector_decomposition(P) |> dim
ambient_dim(P::Union{TropicalPolyhedron,TropicalPointConfiguration}) = pm_object(P).PROJECTIVE_AMBIENT_DIM

function maximal_covectors(as::Type{IncidenceMatrix}, P::TropicalPointConfiguration)
  n = pm_object(P).MAXIMAL_COVECTORS |> length

  return SubObjectIterator{as}(P, _maximal_covectors, n)
end
maximal_covectors(P::TropicalPointConfiguration) = maximal_covectors(IncidenceMatrix, P)

function _maximal_covectors(::Type{IncidenceMatrix}, P::TropicalPointConfiguration, i::Int)
  return P.pm_tpolytope.MAXIMAL_COVECTORS[i]
end

function maximal_covectors(as::Type{IncidenceMatrix}, P::TropicalPolyhedron)
  n = pm_object(P).POLYTOPE_MAXIMAL_COVECTORS |> length

  return SubObjectIterator{as}(P, _maximal_covectors, n)
end
maximal_covectors(P::TropicalPolyhedron) = maximal_covectors(IncidenceMatrix, P)

function _maximal_covectors(::Type{IncidenceMatrix}, P::TropicalPolyhedron, i::Int)
  return P.pm_tpolytope.POLYTOPE_MAXIMAL_COVECTORS[i]
end

@doc raw"""
    covector_decomposition(P::TropicalPointConfiguration)

Get the covector decomposition of the tropical projective torus induced by `P` as a `PolyhedralComplex`.
"""
function covector_decomposition(P::TropicalPointConfiguration)
  return Polymake.tropical.torus_subdivision_as_complex(P.pm_tpolytope) |> polyhedral_complex
end

@doc raw"""
    covector_decomposition(P::TropicalPolyhedon)

Get the covector decomposition of the tropical polytope `P` as a `PolyhedralComplex`.
"""
function covector_decomposition(P::TropicalPolyhedron; dehomogenize_by=1)
  pv = pm_object(P).PSEUDOVERTICES
  cov = pm_object(P).POLYTOPE_MAXIMAL_COVECTOR_CELLS
  ct = pm_object(P).PSEUDOVERTEX_COARSE_COVECTORS
  ind = findall(1:size(ct, 1)) do i
    all(!=(0), ct[i,:])
  end

  if !isnothing(dehomogenize_by)
    coords = filter(!=(dehomogenize_by+1), 1:size(pv, 2))
    for i in 1:size(pv, 1)
      pv[i, 2:end] .-= pv[i, dehomogenize_by+1]
    end
    return Polymake.fan.PolyhedralComplex(VERTICES=pv[ind,coords],MAXIMAL_POLYTOPES=cov[:,ind]) |> polyhedral_complex
  else
    return Polymake.fan.PolyhedralComplex(VERTICES=pv[ind,:],MAXIMAL_POLYTOPES=cov[:,ind]) |> polyhedral_complex
  end
  ## The following is sufficient for Polymake v4.14 and above
  #if !isnothing(dehomogenize_by)
  #  return Polymake.tropical.polytope_subdivision_as_complex(P.pm_tpolytope, dehomogenize_by) |> polyhedral_complex
  #else
  # pv = pm_object(P).PSEUDOVERTICES
  # cov = pm_object(P).POLYTOPE_MAXIMAL_COVECTOR_CELLS
  # ct = pm_object(P).PSEUDOVERTEX_COARSE_COVECTORS
  # ind = findall(1:size(ct, 1)) do i
  #   all(!=(0), ct[i,:])
  # end
  # return Polymake.fan.PolyhedralComplex(VERTICES=pv[ind,:],MAXIMAL_POLYTOPES=cov[:,ind]) |> polyhedral_complex
  #end
end

@doc raw"""
    is_bounded(P::TropicalPolyhedron)

Checks whether `P` is bounded in the tropical projective torus.
"""
function is_bounded(P::TropicalPolyhedron)
  return all(!iszero, pm_object(P).VERTICES)
end

function Base.show(io::IO, P::TropicalPolyhedron{M}) where {M<:MinOrMax}
  d = ambient_dim(P)
  if M == typeof(min)
    print(io, "Min tropical polyhedron in tropical projective torus of dimension $d")
  else
    print(io, "Max tropical polyhedron in tropical projective torus of dimension $d")
  end
end
function Base.show(io::IO, P::TropicalPointConfiguration{M}) where {M<:MinOrMax}
  d = ambient_dim(P)
  if M == typeof(min)
    print(io, "Min tropical point configuration in tropical projective torus of dimension $d")
  else
    print(io, "Max tropical point configuration in tropical projective torus of dimension $d")
  end
end
