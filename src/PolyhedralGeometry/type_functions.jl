################################################################################
######## Scalar types
################################################################################

function detect_scalar_type(
  n::Type{T}, p::Polymake.BigObject
) where {T<:Union{Polyhedron,Cone,PolyhedralFan,SubdivisionOfPoints,PolyhedralComplex}}
  typename = Polymake.bigobject_eltype(p)
  return typename == "OscarNumber" ? nothing : scalar_type_to_oscar[typename]
end

scalar_type(
  ::Union{Polyhedron{T},Cone{T},Hyperplane{T},Halfspace{T}}
) where {T<:scalar_types} = T

################################################################################
######## SubObjectIterator
################################################################################

# Matrices with rational only elements
for (sym, name) in (
  ("linear_inequality_matrix", "Linear Inequality Matrix"),
  ("affine_inequality_matrix", "Affine Inequality Matrix"),
  ("linear_equation_matrix", "Linear Equation Matrix"),
  ("affine_equation_matrix", "Affine Equation Matrix"),
)
  M = Symbol(sym)
  _M = Symbol("_", sym)
  @eval begin
    $M(
      iter::SubObjectIterator{
        <:Union{
          Halfspace{QQFieldElem},
          Hyperplane{QQFieldElem},
          Polyhedron{QQFieldElem},
          Cone{QQFieldElem},
          Pair{Matrix{QQFieldElem},QQFieldElem},
        },
      },
    ) = matrix(QQ, Matrix{QQFieldElem}($_M(Val(iter.Acc), iter.Obj; iter.options...)))
    $M(
      iter::SubObjectIterator{
        <:Union{Halfspace{T},Hyperplane{T},Polyhedron{T},Cone{T},Pair{Matrix{T},T}}
      },
    ) where {T<:scalar_types} =
      matrix(coefficient_field(iter.Obj), $_M(Val(iter.Acc), iter.Obj; iter.options...))
    $_M(::Any, ::PolyhedralObject) =
      throw(ArgumentError(string($name, " not defined in this context.")))
  end
end

function halfspace_matrix_pair(
  iter::SubObjectIterator{
    <:Union{Halfspace{T},Hyperplane{T},Polyhedron{T},Cone{T},Pair{Matrix{T},T}}
  },
) where {T<:scalar_types}
  try
    f = coefficient_field(iter.Obj)
    h = affine_matrix_for_polymake(iter)
    return (A=matrix(f, h[:, 2:end]), b=[f(x) for x in -h[:, 1]])
  catch e
    throw(ArgumentError("Halfspace-Matrix-Pair not defined in this context."))
  end
end

for fun in (cones, faces, facets, maximal_cones, maximal_polyhedra, rays, vertices)
  F = Symbol(fun)
  @eval $F(::Type{IncidenceMatrix}, x...) = IncidenceMatrix($F(x...))
end
