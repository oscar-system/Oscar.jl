function Base.length(PS::Oscar.FiniteRationalPointSet)
  return length(PS.coordinates_list)
end

function Base.iterate(PS::Oscar.FiniteRationalPointSet,i::Int=1)
  i >length(PS) && return nothing
  return PS.codomain(PS.coordinate_list[i]),i+1
end

export finite_pointset

function finite_pointset(X::Union{AffineAlgebraicSet{T},
                                  AffineVariety{T}},
                         v::Vector{Vector{<:FieldElem}}) where T <: Field
  return Oscar.FiniteRationalPointSet(base_scheme(X),X,v)
end

function finite_pointset(X::Union{ProjectiveAlgebraicSet{T},
                                  ProjectiveVariety{T}},
                         v::Vector{Vector{<:FieldElem}}) where T<: Field
  return Oscar.FiniteRationalPointSet(base_scheme(X),X,v)
end