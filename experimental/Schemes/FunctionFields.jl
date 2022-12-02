export function_field, FunctionField
export representative_patch, variety, representative_field

export representative


### Check for irreducibility
@attr Bool function is_irreducible(X::AbsCoveredScheme)
  for U in basic_patches(default_covering(X))
    is_irreducible(U) || return false
  end
  C = default_covering(X)
  # Check that X is connected
  for U in patches(C)
    for V in patches(C)
      A, _ = glueing_domains(C[U, V])
      is_empty(A) && !is_empty(U) && !is_empty(V) && return false
    end
  end
  return true
end


########################################################################
# lower checks for irreducibility                                      #
########################################################################
@attr Bool function is_irreducible(X::AbsSpec{<:Field, <:MPolyRing}) 
  return true
end

@attr Bool function is_irreducible(X::AbsSpec{<:Ring, <:MPolyQuo}) 
  return is_prime(modulus(OO(X)))
end

@attr Bool function is_irreducible(X::AbsSpec{<:Field, <:MPolyLocalizedRing}) 
  return true
end

@attr Bool function is_irreducible(X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing}) 
  return is_prime(modulus(OO(X)))
end


########################################################################
# Check for emptyness                                                  #
########################################################################

@attr Bool function Base.isempty(U::SpecOpen)
  return all(isempty, affine_patches(U))
end

########################################################################
# Methods for VarietyFunctionField                                     #
########################################################################
### essential getters
representative_patch(KK::VarietyFunctionField) = KK.U
variety(KK::VarietyFunctionField) = KK.X
# provide one alias for the user's convenience
scheme(KK::VarietyFunctionField) = variety(KK)
coefficient_ring(KK::VarietyFunctionField) = KK.kk
representative_field(KK::VarietyFunctionField) = KK.KK

### user facing constructors 
@doc Markdown.doc"""
    function_field(X::CoveredScheme)

Return the function field of the irreducible variety `X`.

Internally, a rational function is represented by an element in the field of
fractions of the `ambient_coordinate_ring` of the `representative_patch`.
"""
function_field(X::CoveredScheme; 
               check::Bool=true, 
               representative_patch::AbsSpec=default_covering(X)[1]
              ) = VarietyFunctionField(X,check=check,representative_patch=representative_patch)


@doc Markdown.doc"""
    FunctionField(X::CoveredScheme; kw...)

Return the function field of the irreducible variety `X`.

See [`function_field(X::CoveredScheme)`](@ref).
"""
FunctionField(X::CoveredScheme; kw... ) = function_field(X; kw...)

########################################################################
# Methods for VarietyFunctionFieldElem                                 #
########################################################################
### essential getters 
numerator(f::VarietyFunctionFieldElem) = numerator(f.f)
denominator(f::VarietyFunctionFieldElem) = denominator(f.f)
parent(f::VarietyFunctionFieldElem) = f.KK
representative(f::VarietyFunctionFieldElem) = f.f

### constructors
one(KK::VarietyFunctionField) = VarietyFunctionFieldElem(KK, one(ambient_coordinate_ring(representative_patch(KK))),
                                                         one(ambient_coordinate_ring(representative_patch(KK)))
                                                        )
zero(KK::VarietyFunctionField) = VarietyFunctionFieldElem(KK, zero(ambient_coordinate_ring(representative_patch(KK))),
                                                          one(ambient_coordinate_ring(representative_patch(KK)))
                                                         )

### arithmetic 
function +(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(representative(a) + representative(b), check=false)
end

function -(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(representative(a) - representative(b), check=false)
end

function -(a::T) where {T<:VarietyFunctionFieldElem}
  return (parent(a))((-1)*representative(a), check=false)
end

function *(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === parent(b) || error("the arguments do not have the same parent ring")
  return (parent(a))(representative(a) * representative(b), check=false)
end

function Base.:(//)(a::Integer, b::T) where {T<:VarietyFunctionFieldElem}
  return (parent(b))(a//representative(b))
end

function Base.:(//)(a::fmpz, b::T) where {T<:VarietyFunctionFieldElem}
  return (parent(b))(a//representative(b))
end

function Base.:(//)(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  return (parent(a))(representative(a) // representative(b), check=false)
end

function ==(a::T, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) == parent(b) || error("the arguments do not have the same parent ring")
  KK = parent(a)
  U = representative_patch(KK)
  return iszero(OO(U)(numerator(a)*denominator(b) - numerator(b)*denominator(a)))
end

# We need to manually split this into three methods, because 
# otherwise it seems that Julia can not dispatch this function.
function ^(a::VarietyFunctionFieldElem, i::Int64)
  return parent(a)(representative(a)^i, check=false)
end
function ^(a::VarietyFunctionFieldElem, i::Integer)
  return parent(a)(representative(a)^i, check=false)
end
function ^(a::VarietyFunctionFieldElem, i::fmpz)
  return parent(a)(representative(a)^i, check=false)
end

# Multiplication with elements of the coefficient ring.
# Eventually this needs to be restricted further in case 
# of ambiguities.
function *(a::RingElem, b::T) where {T<:VarietyFunctionFieldElem}
  parent(a) === coefficient_ring(parent(b)) || error("the arguments do not have the same parent ring")
  return (parent(b))(a * representative(b), check=false)
end
function *(b::T, a::RingElem) where {T<:VarietyFunctionFieldElem}
  return a*b
end

# try to avoid a groebner basis computation
iszero(a::VarietyFunctionFieldElem) = iszero(representative(a)) || iszero(OO(representative_patch(parent(a)))(numerator(a)))
isone(a::VarietyFunctionFieldElem) = isone(representative(a)) || iszero(OO(representative_patch(parent(a)))(numerator(a) - denominator(a)))
isunit(a::VarietyFunctionFieldElem) = !iszero(representative(a))

########################################################################
# Conversion of rational functions on arbitrary patches                #
########################################################################

function (KK::VarietyFunctionField)(a::MPolyQuoElem, b::MPolyQuoElem; check::Bool=true)
  return KK(lift(a), lift(b), check=check)
end

function (KK::VarietyFunctionField)(a::MPolyQuoLocalizedRingElem, 
                                    b::MPolyQuoLocalizedRingElem; 
                                    check::Bool=true)
  return KK(lifted_numerator(a)*lifted_denominator(b), lifted_denominator(a)*lifted_numerator(b), check=check)
end

function (KK::VarietyFunctionField)(a::MPolyLocalizedRingElem, 
                                    b::MPolyLocalizedRingElem; 
                                    check::Bool=true)
  return KK(numerator(a)*denominator(b), denominator(a)*numerator(b), check=check)
end

function (KK::VarietyFunctionField)(a::MPolyElem, b::MPolyElem; check::Bool=true)
  R = parent(a)
  R === parent(b) || error("rings are not compatible")
  R === ambient_coordinate_ring(representative_patch(KK)) && return VarietyFunctionFieldElem(KK, a, b)
  
  # otherwise check whether we can find the ring of h among the affine patches
  R in [ambient_coordinate_ring(V) for V in patches(default_covering(variety(KK)))] || error("ring does not belong to any of the affine charts")
  # allocate a variable for the patch in which a and be are living
  V = representative_patch(KK)
  X = variety(KK)
  C = default_covering(X)
  for i in 1:npatches(C)
    if ambient_coordinate_ring(C[i]) === R
      V = C[i]
      break
    end
  end

  # convert it 
  U = representative_patch(KK)
  h_generic = move_representative(a, b, V, U, C)
  return VarietyFunctionFieldElem(KK, numerator(h_generic),
                                  denominator(h_generic)
                                  )
end

@Markdown.doc """
    function move_representative(
        a::MPolyElem, b::MPolyElem,
        V::AbsSpec, U::AbsSpec,
        C::Covering
      )

Given a fraction ``a/b ∈ Quot(P)`` with ``P = 𝕜[x]`` the `ambient_coordinate_ring`
of an affine patch ``V`` in a covering `C`, move that fraction to 
one in ``Quot(P')`` where ``P'`` is the ambient coordinate ring of another patch ``U``.

**Note:** This is only guaranteed to work for irreducible schemes! 
"""
function move_representative(
    a::MPolyElem, b::MPolyElem,
    V::AbsSpec, U::AbsSpec,
    C::Covering
  )
  G = C[U, V]
  f, _ = glueing_morphisms(G)
  A, B = glueing_domains(G)
  pba = pullback(f)(OO(B)(a))
  pbb = pullback(f)(OO(B)(b))
  iszero(pbb) && error("pullback of denominator is zero")
  # in the next line, A is either a SpecOpen or a PrincipalOpenSubset
  h_generic = generic_fraction(pba, A)//generic_fraction(pbb, A)
  return h_generic
end

function (KK::VarietyFunctionField)(h::AbstractAlgebra.Generic.Frac; check::Bool=true)
  return KK(numerator(h), denominator(h), check=check)
end

function Base.show(io::IO, f::VarietyFunctionFieldElem)
  print(io, representative(f))
end

### given the fraction field of the `ambient_coordinate_ring` in one
# affine chart, return a representative of `f` in that field
function (K::AbstractAlgebra.Generic.FracField)(f::VarietyFunctionFieldElem)
  R = base_ring(K)
  V = representative_patch(parent(f))
  C = default_covering(variety(parent(f)))
  for U in patches(C)
    if ambient_coordinate_ring(U) == R
      V = U
      break
    end
  end
  f_mov = move_representative(numerator(f), denominator(f), 
                              representative_patch(parent(f)),
                              V, C
                             )
  return K(numerator(f_mov), denominator(f_mov))
end

function getindex(f::VarietyFunctionFieldElem, V::AbsSpec)
  C = default_covering(variety(parent(f))) 
  V in C || error("patch not found")
  return move_representative(numerator(f), denominator(f), 
                             representative_patch(parent(f)),
                             V, C
                            )
end

########################################################################
# Implementation of the rest of the interfaces                         #
########################################################################

function elem_type(::Type{T}) where {FracFieldType, 
                                     T<:VarietyFunctionField{<:Field, 
                                                             <:FracFieldType
                                                            }
                                    }
  return VarietyFunctionFieldElem{elem_type(FracFieldType), T}
end

function parent_type(::Type{T}) where {ParentType, T<:VarietyFunctionFieldElem{<:Any, ParentType}}
  return ParentType
end

base_ring(KK::VarietyFunctionField) = base_ring(representative_field(KK))
base_ring(a::VarietyFunctionFieldElem) = base_ring(parent(a))
is_domain_type(::Type{T}) where {T<:VarietyFunctionFieldElem} = true
is_exact_type(::Type{T}) where {T<:VarietyFunctionFieldElem} = true

function Base.hash(f::VarietyFunctionFieldElem, h::UInt)
  r = 57103
  return xor(r, hash(representative(f), h))
end

function deepcopy_internal(f::VarietyFunctionFieldElem, dict::IdDict)
  return parent(f)(deepcopy_internal(numerator(representative(f)), dict), 
                   deepcopy_internal(denominator(representative(f)), dict), 
                   check=false)
end

(KK::VarietyFunctionField)() = zero(KK)
(KK::VarietyFunctionField)(a::Integer) = KK(base_ring(KK)(a), one(base_ring(KK)), check=false)
(KK::VarietyFunctionField)(f::VarietyFunctionFieldElem) = (parent(f) == KK ? f : error("element does not belong to the given field"))
(KK::VarietyFunctionField)(a::MPolyElem) = KK(a, one(parent(a)), check=false)
canonical_unit(f::VarietyFunctionFieldElem) = f # part of the ring interface that becomes trivial for fields

function Base.show(io::IO, KK::VarietyFunctionField)
  print(io, "function field of $(variety(KK))")
end

function divexact(f::VarietyFunctionFieldElem, 
    g::VarietyFunctionFieldElem;
    check::Bool=true
  )
  return f//g
end
inv(f::VarietyFunctionFieldElem) = parent(f)(denominator(representative(f)),
                                      numerator(representative(f)),
                                      check=false
                                     )

AbstractAlgebra.promote_rule(::Type{T}, ::Type{S}) where {T<:VarietyFunctionFieldElem, S<:Integer} = T
AbstractAlgebra.promote_rule(::Type{T}, ::Type{S}) where {T<:VarietyFunctionFieldElem, S<:AbstractAlgebra.Generic.Frac} = T

AbstractAlgebra.promote_rule(::Type{T}, ::Type{T}) where {T<:VarietyFunctionFieldElem} = T

function AbstractAlgebra.promote_rule(::Type{FFET}, ::Type{U}) where {T, FFET<:VarietyFunctionFieldElem{T}, U<:RingElement}
  promote_rule(T, U) == T ? FFET : Union{}
end 

function AbstractAlgebra.promote_rule(::Type{FFET}, ::Type{U}) where {T, FFET<:VarietyFunctionFieldElem{AbstractAlgebra.Generic.Frac{T}}, U<:RingElement}
  promote_rule(T, U) == T ? FFET : Union{}
end 

characteristic(KK::VarietyFunctionField) = characteristic(base_ring(variety(KK)))

@Markdown.doc """
    is_regular(f::VarietyFunctionFieldElem, U::Scheme)

Return whether ``f ∈ K(X)`` restricts to a regular function 
on an open subset ``U ⊂ X``.

**Note:** ``U`` must either be an `affine_chart` of ``X`` or 
its `ambient_scheme` must be an `affine_chart`.
"""
function is_regular(f::VarietyFunctionFieldElem, U::Scheme)
  error("method not implemented")
end

function is_regular(f::VarietyFunctionFieldElem, U::AbsSpec)
  p = numerator(f[U])
  q = denominator(f[U])
  return _is_regular_fraction(OO(U), p, q)
end

function is_regular(f::VarietyFunctionFieldElem, U::PrincipalOpenSubset)
  Y = ambient_scheme(U)
  ff = f[Y]
  p = numerator(ff)
  q = denominator(ff)
  return _is_regular_fraction(OO(U), p, q)
end

function is_regular(f::VarietyFunctionFieldElem, W::SpecOpen)
  return all(U->is_regular(f, U), affine_patches(W))
end

