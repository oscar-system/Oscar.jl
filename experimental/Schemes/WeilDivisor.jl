export WeilDivisor
export scheme_type, ideal_sheaf_type, coefficient_ring_type, coefficient_type
export scheme, components, coefficient_dict, coefficient_ring

@Markdown.doc """
    WeilDivisor{
      CoveredSchemeType<:CoveredScheme, 
      IdealSheafType<:IdealSheaf, 
      CoefficientRingType<:AbstractAlgebra.Ring, 
      CoefficientRingElemType<:AbstractAlgebra.RingElem
    }

Describes a Weil divisor on a covered scheme ``X`` of type 
`CoveredSchemeType` as a formal linear combination of 
ideal sheaves of type `IdealSheafType` with coefficients 
in a ring ``R`` of type `CoefficientRingType`.
"""
@attributes mutable struct WeilDivisor{
    CoveredSchemeType<:CoveredScheme, 
    IdealSheafType<:IdealSheaf, 
    CoefficientRingType<:AbstractAlgebra.Ring, 
    CoefficientRingElemType<:AbstractAlgebra.RingElem
   }
  X::CoveredSchemeType # the parent
  R::CoefficientRingType # the ring of coefficients
  coefficients::Dict{IdealSheafType, CoefficientRingElemType} # the formal linear combination

  function WeilDivisor(
      X::CoveredSchemeType, 
      R::CoefficientRingType, 
      coefficients::Dict{IdealSheafType, CoefficientRingElemType}
    ) where {CoveredSchemeType, IdealSheafType, CoefficientRingType, CoefficientRingElemType}
    # TODO: Do we want to require that the different effective divisors 
    # have the same underlying covering? Probably not.
    for D in keys(coefficients)
      scheme(D) == X || error("component of divisor does not lay in the given scheme")
      parent(coefficients[D]) == R || error("coefficients do not lay in the given ring")
    end
    return new{CoveredSchemeType, IdealSheafType, CoefficientRingType, CoefficientRingElemType}(X, R, coefficients)
  end
end

### type getters 
scheme_type(D::WeilDivisor{S, T, U, V}) where{S, T, U, V} = S
scheme_type(::Type{WeilDivisor{S, T, U, V}}) where{S, T, U, V} = S
ideal_sheaf_type(D::WeilDivisor{S, T, U, V}) where{S, T, U, V} = T
ideal_sheaf_type(::Type{WeilDivisor{S, T, U, V}}) where{S, T, U, V} = T
coefficient_ring_type(D::WeilDivisor{S, T, U, V}) where{S, T, U, V} = U
coefficient_ring_type(::Type{WeilDivisor{S, T, U, V}}) where{S, T, U, V} = U
coefficient_type(D::WeilDivisor{S, T, U, V}) where{S, T, U, V} = V
coefficient_type(::Type{WeilDivisor{S, T, U, V}}) where{S, T, U, V} = V

### getter methods
scheme(D::WeilDivisor) = D.X
getindex(D::WeilDivisor, I::IdealSheaf) = (D.coefficients)[I]
components(D::WeilDivisor) = [ Z for Z in keys(D.coefficients)]
coefficient_dict(D::WeilDivisor) = D.coefficients
coefficient_ring(D::WeilDivisor) = D.R

set_name!(X::WeilDivisor, name::String) = set_attribute!(X, :name, name)
name_of(X::WeilDivisor) = get_attribute(X, :name)::String
has_name(X::WeilDivisor) = has_attribute(X, :name)

function setindex!(D::WeilDivisor, c::RET, I::IdealSheaf) where {RET<:AbstractAlgebra.RingElem}
  parent(c) == coefficient_ring(D) || error("coefficient does not belong to the correct ring")
  coefficient_dict(D)[I] = c
end

function WeilDivisor(X::CoveredScheme, R::RingType) where {RingType<:AbstractAlgebra.Ring}
  D = Dict{ideal_sheaf_type(X), elem_type(R)}()
  return WeilDivisor(X, R, D)
end

function WeilDivisor(I::IdealSheaf, R::RT) where {RT<:AbstractAlgebra.Ring} 
  D = WeilDivisor(scheme(I), R)
  D[I] = one(R)
  return D
end

function WeilDivisor(I::IdealSheaf)
  D = WeilDivisor(scheme(I), ZZ)
  D[I] = one(ZZ)
  return D
end

### copy constructor
# Only copies the coefficient dictionary; everything else has parent character.
function WeilDivisor(D::WeilDivisor) 
  return WeilDivisor(scheme(D), coefficient_ring(D), copy(coefficient_dict(D)))
end

function Base.show(io::IO, D::WeilDivisor)
  if has_name(D)
    print(io, name_of(D))
    return
  end
  if length(components(D)) == 0
    println(io, "the zero Weil divisor on $(scheme(D))")
    return
  end
  println(io, "Weil divisor on $(scheme(D)) given as the formal sum")
  comp = ["$(D[I]) ⋅ $(I)" for I in components(D)]
  out_str = comp[1]
  for c in comp[2:end]
    out_str = out_str* " + " * c
  end
  println(io, out_str)
end

function +(D::T, E::T) where {T<:WeilDivisor}
  X = scheme(D)
  X == scheme(E) || error("divisors do not live on the same scheme")
  R = coefficient_ring(D)
  R == coefficient_ring(E) || error("coefficient rings do not coincide")
  C = WeilDivisor(D)
  for I in components(E)
    if haskey(coefficient_dict(C), I)
      c = C[I] + E[I]
      if iszero(c) 
        delete!(coefficient_dict(C), I)
      else 
        C[I] = c
      end
    else
      coefficient_dict(C)[I] = E[I]
    end
  end
  return C
end

function -(D::T) where {T<:WeilDivisor}
  E = WeilDivisor(D)
  for I in components(E)
    E[I] = -E[I]
  end
  return E
end

-(D::T, E::T) where {T<:WeilDivisor} = D + (-E)

function *(a::RET, E::WeilDivisor) where {RET<:AbstractAlgebra.RingElem}
  c = coefficient_ring(E)()
  D = WeilDivisor(E)
  parent(a) == coefficient_ring(E) ? (c = a)::elem_type(coefficient_ring(E)) : c = coefficient_ring(E)(a)
  for I in components(D)
    D[I] = c*D[I]
  end
  return D
end

*(a::Int, E::WeilDivisor) = coefficient_ring(E)(a)*E
*(a::Integer, E::WeilDivisor) = coefficient_ring(E)(a)*E

+(D::WeilDivisor, I::IdealSheaf) = D + WeilDivisor(I)

function intersection(D::T, E::T) where {T<:WeilDivisor}
  X = scheme(D)
  X == scheme(E) || error("divisors do not live on the same scheme")
  R = coefficient_ring(D)
  R == coefficient_ring(E) || error("divisors do not have the same coefficient ring")
  # prepare a copy of the divisors
  D_copy = WeilDivisor(X, R)
  E_copy = WeilDivisor(X, R)
  # check whether a common refinement of the covering is necessary
  CD = covering(D)
  CE = covering(E)
  if CD != CE
    CC, f, g = common_refinement(X, CD, CE)
    D_copy = pullback(f, D)
    E_copy = pullback(g, E)
  else
    D_copy = D
    E_copy = E
  end
  # TODO: Work out the intersection
end


