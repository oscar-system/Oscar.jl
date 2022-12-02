export set_name!


####################################################################################
# (1) Equality
####################################################################################

# TODO: Add further cross-type comparison methods as needed.

function ==(X::T, Y::T) where {T<:AbsSpec{<:Ring, <:MPolyRing}}
  return OO(X) === OO(Y)
end


function ==(X::AbsSpec, Y::AbsSpec)
  X === Y && return true
  return issubset(X, Y) && issubset(Y, X)
end


function ==(X::AbsSpec, Y::EmptyScheme)
  return issubset(X, Y)
end


==(X::EmptyScheme, Y::AbsSpec) = (Y == X)



########################################################
# (2) Display
########################################################

function Base.show(io::IO, X::AbsSpec)
  if has_attribute(X, :name)
    print(io, name(X))
    return
  end
  print(io, "Spec of $(OO(X))")
end



########################################################
# (3) Check for zero divisors in rings
########################################################

@Markdown.doc """
    is_non_zero_divisor(f::RingElem, X::AbsSpec)

Checks if a ring element is a non-zero divisor
in the coordinate ring of an affine scheme.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1, x2, x3) = gens(OO(X))
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> is_non_zero_divisor(x1, X)
true

julia> is_non_zero_divisor(zero(OO(X)), X)
false
```
"""
function is_non_zero_divisor(f::RingElem, X::AbsSpec)
  error("method not implemented for affine schemes of type $(typeof(X))")
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::MPolyQuoElem, X::AbsSpec{<:Ring, <:MPolyQuo})
  R = ambient_coordinate_ring(X)
  I = modulus(OO(X))
  J = ideal(R, lift(f))
  return I == quotient(I, J)
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyLocalizedRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing})
  I = ideal(OO(X), [zero(OO(X))])
  zero_ideal = Oscar.pre_image_ideal(I)
  J = Oscar.pre_image_ideal(ideal(OO(X), [f]))
  Q = quotient(zero_ideal, J)
  return zero_ideal == Q
end

########################################################################
# High level constructors of subschemes                                #
########################################################################

function sub(X::AbsSpec, I::Ideal)
  inc = ClosedEmbedding(X, I) # Will complain if input is not compatible
  return domain(inc), inc
end

########################################################################
# Further intersections                                                #
########################################################################

function intersect(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
  if !(ambient_scheme(U) === ambient_scheme(V))
    return intersect(underlying_scheme(U), underlying_scheme(V))
  end
  return PrincipalOpenSubset(ambient_scheme(U), complement_equation(U)*complement_equation(V))
end

