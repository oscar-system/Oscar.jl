####################################################################################
# (1) Check if a scheme is empty
####################################################################################

@doc raw"""
    is_empty(X::AbsSpec)

Check whether the affine scheme ``X`` is empty.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  with coordinates x1 x2 x3
  over Rational field

julia> isempty(X)
false

julia> is_empty(subscheme(X, one(OO(X))))
true

julia> isempty(EmptyScheme(QQ))
true
```
"""
Base.isempty(X::AbsSpec) = iszero(one(OO(X)))
is_empty(X::EmptyScheme) = true



####################################################################################
# (2) IsSubset for all other schemes
####################################################################################


# (2.0) For empty schemes and whenever issubset cannot be implemented


@doc raw"""
    is_subset(X::AbsSpec, Y::AbsSpec)

Check whether ``X`` is a subset of ``Y`` based on the comparison of their coordinate rings.
See [`inclusion_morphism(::AbsSpec, ::AbsSpec)`](@ref) for the corresponding morphism.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  with coordinates x1 x2 x3
  over Rational field

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate polynomial ring in 3 variables over QQ by ideal(x1*x2)

julia> is_subset(X, Y)
false

julia> is_subset(Y, X)
true
```
"""
function issubset(X::AbsSpec, Y::AbsSpec)
  error("method `issubset(X, Y)` not implemented for `X` of type $(typeof(X)) and `Y` of type $(typeof(Y))")
end

issubset(X::EmptyScheme{BRT}, Y::Scheme{BRT}) where {BRT} = true

function issubset(Y::AbsSpec{BRT, <:Any}, X::EmptyScheme{BRT}) where {BRT}
  return iszero(one(OO(Y)))
end


# (2.1) MPolyRing in first argument

function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyRing}
  return OO(X) === OO(Y)
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing}, 
    Y::AbsSpec{BRT, <:MPolyQuoRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || return false
  return iszero(modulus(OO(Y)))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || return false
  return issubset(inverted_set(OO(Y)), units_of(R))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || return false
  return issubset(inverted_set(OO(Y)), units_of(R)) && iszero(ambient_closure_ideal(Y))
end


# (2.2) MPolyQuoRing in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyQuoRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = ambient_coordinate_ring(Y)
  R === ambient_coordinate_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoRing}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return issubset(ambient_closure_ideal(Y), ambient_closure_ideal(X))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoRing},
    Y::AbsSpec{BRT, <:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return iszero(localized_ring(OO(Y))(modulus(OO(X))))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return issubset(modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end


# (2.3) MPolyLocRing in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyLocRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == ambient_coordinate_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, <:MPolyLocRing},
    Y::AbsSpec{BRT, <:MPolyQuoRing}
  ) where {BRT}
  R = ambient_coordinate_ring(Y)
  R == ambient_coordinate_ring(X) || return false
  return all(x->(iszero(OO(X)(x))), gens(modulus(OO(Y))))
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyLocRing}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX)
end


function issubset(
    X::AbsSpec{BRT, <:MPolyLocRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing}
  ) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX) && iszero(modulus(OO(Y)))
end


# (2.4) MPolyQuoLocRing in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == ambient_coordinate_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocRing},
    Y::AbsSpec{BRT, <:MPolyQuoRing}
  ) where {BRT}
  R = ambient_coordinate_ring(Y)
  R == ambient_coordinate_ring(X) || return false
  L = localized_ring(OO(X))
  return issubset(L(modulus(OO(Y))), modulus(OO(X)))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocRing},
    Y::AbsSpec{BRT, <:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  if !issubset(UY, UX)
    # check whether the inverted elements in Y are units anyway
    for a in denominators(UY)
      is_unit(OO(X)(a)) || return false
    end
  end
  return true  # Spec R/I[S^-1] is a closed subscheme of Spec R[S^-1]
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  if !issubset(UY, UX)
    # check whether the inverted elements in Y are units anyway
    for a in denominators(UY)
      is_unit(OO(X)(a)) || return false
    end
  end
  J = localized_ring(OO(X))(modulus(underlying_quotient(OO(Y))))
  return issubset(J, modulus(OO(X)))
end

function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  issubset(UY,UX) || return false # element of KPointIdeal inverted in UY ?
  J = localized_ring(OO(X))(modulus(underlying_quotient(OO(Y))))
  return issubset(J, modulus(OO(X)))
end


####################################################################################
# (3) Check if a scheme is openly embeeded in another scheme
####################################################################################

#TODO: Add more cross-type methods as needed.

@doc raw"""
    is_open_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is openly embedded in ``Y``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  with coordinates x1 x2 x3
  over Rational field

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate polynomial ring in 3 variables over QQ by ideal(x1*x2)

julia> is_open_embedding(Y, X)
false

julia> Z = hypersurface_complement(X, x1)
Spec of localization of Multivariate polynomial ring in 3 variables over QQ at the powers of QQMPolyRingElem[x1]

julia> is_open_embedding(Z, X)
true
```
"""
function is_open_embedding(X::AbsSpec, Y::AbsSpec)
  return is_open_embedding(standard_spec(X), standard_spec(Y))
end


function is_open_embedding(
    X::Spec{BRT, RT},
    Y::Spec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  issubset(UY, UX) || return false
  J = localized_ring(OO(X))(modulus(underlying_quotient(OO(Y))))
  return modulus(OO(X)) == J 
end


function is_open_embedding(
    X::Spec{BRT, <:MPolyQuoLocRing},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  return OO(Y) == ambient_coordinate_ring(X) && all(iszero, gens(modulus(OO(X))))
end



####################################################################################
# (4) Check if a scheme can be embedded via a closed embeeded in another scheme
####################################################################################

@doc raw"""
    is_closed_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is closed embedded in ``Y``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  with coordinates x1 x2 x3
  over Rational field

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate polynomial ring in 3 variables over QQ by ideal(x1*x2)

julia> is_closed_embedding(Y, X)
true

julia> Z = hypersurface_complement(X, x1)
Spec of localization of Multivariate polynomial ring in 3 variables over QQ at the powers of QQMPolyRingElem[x1]

julia> is_closed_embedding(Z, X)
false
```
"""
function is_closed_embedding(X::AbsSpec, Y::AbsSpec)
  error("`is_closed_embedding(X, Y)` not implemented for X of type $(typeof(X)) and Y of type $(typeof(Y))")
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuoRing},
    Y::AbsSpec{<:Ring, <:MPolyRing}
  )
  return ambient_coordinate_ring(X) === ambient_coordinate_ring(Y)
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyRing},
    Y::AbsSpec{<:Ring, <:MPolyQuoRing}
  )
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return iszero(modulus(OO(Y)))
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuoRing},
    Y::AbsSpec{<:Ring, <:MPolyQuoRing}
  )
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return issubset(modulus(OO(Y)), modulus(OO(X)))
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyRing},
    Y::AbsSpec{<:Ring, <:MPolyRing}
  )
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return true
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyLocRing{<:Any, <:Any, <:Any, <:Any, 
                                            <:MPolyPowersOfElement}},

    Y::AbsSpec{<:Ring, <:MPolyRing}
  )
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  for f in inverted_set(OO(X))
    isunit(OO(Y)(f)) || return false
  end
  return true
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}},

    Y::AbsSpec{<:Ring, <:MPolyQuoRing}
  )
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  for x in inverted_set(OO(X)) 
    isunit(OO(Y)(x)) || return false
  end
  for g in gens(modulus(OO(Y)))
    iszero(OO(X)(g)) || return false
  end
  return true
end


function is_closed_embedding(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  inverted_set(OO(X)) == inverted_set(OO(Y)) || return false
  J = localized_ring(OO(X))(modulus(underlying_quotient(OO(Y))))
  return issubset(J, modulus(OO(X)))
end


function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuoRing},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  OO(Y) === ambient_coordinate_ring(X) || return false
  return true
end


function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuoRing},
    Y::Spec{BRT, <:RT}
  ) where {BRT, RT<:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  all(x->(isunit(OO(X)(x))), denominators(inverted_set(OO(Y)))) || return false
  return issubset(modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end

#############################################################################
# (5) Check, if a Spec is equidimensional
#############################################################################
# TODO: projective schemes, covered schemes

@doc raw"""
   is_equidimensional(X::AbsSpec{<:Field, <:MPolyAnyRing}) 

Check whether the scheme `X` is equidimensional.

Currently this command is available for affine schemes and space germs.

This command relies on [`equidimensional_decomposition_radical`](@ref).

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> I = ideal(R,[(x-y)])
ideal(x - y)

julia> J = ideal(R,[x-1,y-2])
ideal(x - 1, y - 2)

julia> X = Spec(R,I)
Spec of Quotient of Multivariate polynomial ring in 2 variables over QQ by ideal(x - y)

julia> Y = Spec(R,I*J)
Spec of Quotient of Multivariate polynomial ring in 2 variables over QQ by ideal(x^2 - x*y - x + y, x*y - 2*x - y^2 + 2*y)

julia> is_equidimensional(X)
true

julia> is_equidimensional(Y)
false
```
"""
@attr Bool function is_equidimensional(X::AbsSpec{<:Field, <:MPAnyQuoRing})
  I = modulus(OO(X))
# equidimensional decomposition only available for schemes over a field
  P = equidimensional_decomposition_radical(saturated_ideal(I))
  length(P) < 2 && return true
  return false
end

# make is_equidimensional agnostic to quotient
@attr Bool function is_equidimensional(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  return true
end

##############################################################################
# (6) Check, if scheme is reduced
##############################################################################
# TODO: projective schemes

@doc raw"""
   is_reduced(X::AbsSpec{<:Field, <:MPolyAnyRing})

Check whether the affine scheme `X` is reduced.
"""
@attr Bool function is_reduced(X::AbsSpec{<:Field, <:MPAnyQuoRing})
  I = saturated_ideal(modulus(OO(X)))
  return is_reduced(quo(base_ring(I), I)[1])
end

## make is_reduced agnostic to quotient ring
@attr Bool function is_reduced(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  return true
end

@attr Bool function is_geometrically_reduced(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  return true
end

@attr Bool function is_geometrically_reduced(X::AbsSpec{<:Field, <:MPAnyQuoRing})
  F = base_ring(X)
  # is_perfect(F) # needs new AbstractAlgebra version
  if characteristic(F) == 0 || F isa FinField  # F is perfect
    return is_reduced(X)
  end
  throw(NotImplementedError(:is_geometrically_reduced, "currently we can decide this only over a perfect base field"))
end

########################################################################
# (7) Smoothness test based on projective modules.
# The routine checks whether the module for the cotangent sheaf Ω¹(X)
# is locally free over 𝒪(X) and returns `true` if this is the case. 
########################################################################
# TODO: is_regular using Hironaka's criterion

@doc raw"""
    is_smooth(X::AbsSpec{<:Field, <:MPolyAnyRing})

Check whether the scheme `X` is smooth.

Note that smoothness and regularity do not coincide over non-perfect fields. 
Smoothness implies regularity, but regular non-smooth schemes exist.

See also [`singular_locus`](@ref), [`singular_locus_reduced`](@ref).

# Algorithm

This method checks whether the module for the cotangent sheaf Ω¹(X)
is locally free over 𝒪(X).

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> I = ideal(R,[x-y^2])
ideal(x - y^2)

julia> J = ideal(R,[x^2-y^2])
ideal(x^2 - y^2)

julia> X = Spec(R, I)
Spec of Quotient of Multivariate polynomial ring in 2 variables over QQ by ideal(x - y^2)

julia> is_smooth(X)
true

julia> Y = Spec(R, J)
Spec of Quotient of Multivariate polynomial ring in 2 variables over QQ by ideal(x^2 - y^2)

julia> is_smooth(Y)
false

julia> U = MPolyComplementOfKPointIdeal(R,[1,1])
complement of maximal ideal corresponding to point with coordinates QQFieldElem[1, 1]

julia> Z = Spec(R, J, U)
Spec of Localization of Quotient of Multivariate polynomial ring in 2 variables over QQ by ideal(x^2 - y^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates QQFieldElem[1, 1]

julia> is_smooth(Z)
true

```
"""
@attr Bool function is_smooth(X::AbsSpec{<:Field, <:MPolyQuoLocRing})
  R = base_ring(OO(X))
  L = localized_ring(OO(X))
  I = modulus(OO(X))
  f = gens(saturated_ideal(I))
  Df = jacobi_matrix(f)
  A = map_entries(x->OO(X)(x), Df)
  success, _, _ = Oscar._is_projective_without_denominators(A, task=:without_projector)
  return success
end

@attr Bool function is_smooth(X::AbsSpec{<:Field, <:MPolyQuoRing})
  R = base_ring(OO(X))
  I = modulus(OO(X))
  f = gens(I)
  Df = jacobi_matrix(f)
  A = map_entries(x->OO(X)(x), Df)
  success, _, _ = Oscar._is_projective_without_denominators(A, task=:without_projector)
  return success
end

## make is_smooth agnostic to quotient ring
is_smooth(X::AbsSpec{<:Field, <:MPolyRing}) = true
is_smooth(X::AbsSpec{<:Field, <:MPolyLocRing}) = true

###################################################################
# Irreducibility and Integrality                                  #
#    integral iff (irreducible && reduced)                        #
#    integral = OO(X) is integral domain                          #
#    irreducible = nilradical of OO(X) is prime                   #
###################################################################
@doc raw"""
   is_irreducible(X::AbsSpec)

Check whether the affine scheme `X` is irreducible.

!!! note
    Irreducibility is checked over the (computable) base field of the affine scheme as specified upon creation of the ring, not over the algebraic closure thereof.

"""
@attr Bool function is_irreducible(X::AbsSpec{<:Field, <:MPolyAnyRing})
  !is_empty(X) || return false
  !get_attribute(X, :is_integral, false) || return true 
                                           ## integral = irreducible + reduced
  return (length(minimal_primes(saturated_ideal(modulus(OO(X))))) == 1)
end

@doc raw"""
   is_integral(X::AbsSpec)

Check whether the affine scheme `X` is integral, i.e. irreducible and reduced.
"""
@attr Bool function is_integral(X::AbsSpec{<:Field, <:MPolyAnyRing})
  !is_empty(X) || return false
  if has_attribute(X,:is_reduced) && has_attribute(X,:is_irreducible)
     return get_attribute(X,:is_reduced) && get_attribute(X,:is_irreducible)
  end
  return is_prime(modulus(OO(X)))
end

@doc raw"""
    is_geometrically_integral(X::AbsSpec)

Test if ``X/k`` is geometrically integral.

That is if ``X`` is integral when base changed to any field extension of ``k``.
"""
@attr Bool function is_geometrically_integral(X::AbsSpec{<:Field,<:MPolyAnyRing})
  is_integral(X) || return false
  throw(NotImplementedError(:is_geometrically_integral, "absolute primary decomposition is currently only available over the rationals"))
end

@attr Bool function is_geometrically_integral(X::AbsSpec{<:QQField, <:MPolyAnyRing})
  is_integral(X) || return false
  I = ambient_closure_ideal(X)
  AI = absolute_primary_decomposition(I)
  @assert length(AI)==1 # it is prime since X is integral
  return AI[1][4]==1
end

###################################################################
# Connectedness                                                   #
###################################################################
@doc raw"""
   is_connected(X::AbsSpec)

Check whether the affine scheme `X` is connected.
"""
@attr Bool function is_connected(X::AbsSpec)
  error("not implemented yet")
## note for future implementation: expensive property
## 1) do primary decomposition
## 2) check connectedness of lowest two layers of the intersection lattice
  return get_attribute(X,:is_connected)
end
