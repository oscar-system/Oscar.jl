## Iterator over FiniteRationalPointSet
function Base.length(PS::FiniteRationalPointSet)
  return length(PS.coordinates_list)
end

function Base.iterate(PS::Oscar.FiniteRationalPointSet,i::Int=1)
  i >length(PS) && return nothing
  return PS.codomain(PS.coordinates_list[i]),i+1
end

# constructors for FiniteRationalPointSet
function finite_point_set(X::Union{AffineAlgebraicSet{T},
                                  AffineVariety{T}},
                         v::Vector{Vector{<:FieldElem}}) where T <: Field
  return FiniteRationalPointSet(base_scheme(X),X,v)
end

function finite_point_set(X::Union{ProjectiveAlgebraicSet{T},
                                  ProjectiveVariety{T}},
                         v::Vector{Vector{<:FieldElem}}) where T<: Field
  return FiniteRationalPointSet(base_scheme(X),X,v)
end