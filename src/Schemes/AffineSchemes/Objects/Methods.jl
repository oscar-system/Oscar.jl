export set_name


####################################################################################
# (1) Set name for scheme
####################################################################################

# Does this work?
function set_name!(X::AbsSpec, name::String)
  return set_attribute!(X, :name, name)
end



####################################################################################
# (2) Equality
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
# (3) Display
########################################################

function Base.show(io::IO, X::Spec)
  if isdefined(X, :name)
    print(io, name_of(X))
    return
  end
  print(io, "Spec of $(OO(X))")
end



########################################################
# (4) Check for zero divisors in rings
########################################################

function is_non_zero_divisor(f::RingElem, X::AbsSpec{<:Ring, <:MPolyRing})
  return !iszero(OO(X)(f))
end

function is_non_zero_divisor(f::MPolyQuoElem, X::AbsSpec{<:Ring, <:MPolyQuo})
  R = ambient_ring(X)
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
