function vertices(as::Type{PointVector{T}}, P::TropicalPolyhedron) where {T<:TropicalSemiringElem}
  n = n_vertices(P)

  return SubObjectIterator{as}(P, _vertex_of_polyhedron, n)
end
vertices(P::TropicalPolyhedron{M}) where {M<:MinOrMax} = vertices(PointVector{TropicalSemiringElem{M}}, P)

function _vertex_of_polyhedron(::Type{PointVector{TropicalSemiringElem{M}}}, P::TropicalPolyhedron, i::Int) where {M<:MinOrMax}
  T = tropical_semiring(convention(P))

  return point_vector(
    T,
    T.(P.pm_tpolytope.POINTS[i,:])
  )
end

n_vertices(P::TropicalPolyhedron) = size(pm_object(P).VERTICES, 1)

function pseudovertices(as::Type{PointVector{T}}, P::TropicalPolyhedron) where {T<:TropicalSemiringElem}
  n = pm_object(P).PSEUDOVERTICES |> size |> first

  return SubObjectIterator{as}(P, _pseudovertices, n)
end
pseudovertices(P::TropicalPolyhedron{M}) where {M<:MinOrMax} = pseudovertices(PointVector{TropicalSemiringElem{M}}, P)

function _pseudovertices(::Type{PointVector{TropicalSemiringElem{M}}}, P::TropicalPolyhedron, i::Int) where {M<:MinOrMax}
  T = tropical_semiring(convention(P))

  return point_vector(
    T,
    T.(pm_object(P).PSEUDOVERTICES[i,2:end])
  )
end

n_pseudovertices(P::TropicalPolyhedron) = size(pm_object(P).PSEUDOVERTICES, 1)

dim(P::TropicalPolyhedron) = polytope_covector_decomposition(P) |> dim
ambient_dim(P::TropicalPolyhedron) = pm_object(P).PROJECTIVE_AMBIENT_DIM

function torus_covector_decomposition(P::TropicalPolyhedron)
  return Polymake.tropical.torus_subdivision_as_complex(P.pm_tpolytope) |> polyhedral_complex
end

function polytope_covector_decomposition(P::TropicalPolyhedron)
  return Polymake.tropical.polytope_subdivision_as_complex(P.pm_tpolytope) |> polyhedral_complex
end

function Base.show(io::IO, P::TropicalPolyhedron{T}) where {T<:MinOrMax}
  d = ambient_dim(P)
  if T == typeof(min)
    print(io, "Min tropical polyhedron in tropical projective torus of dimension $d")
  else
    print(io, "Max tropical polyhedron in tropical projective torus of dimension $d")
  end
end
