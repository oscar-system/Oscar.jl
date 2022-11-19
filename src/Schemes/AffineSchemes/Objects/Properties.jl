export is_empty, issubset, is_open_embedding, is_closed_embedding

export is_equidimensional

export is_smooth

####################################################################################
# (1) Check if a scheme is empty
####################################################################################

@Markdown.doc """
    is_empty(X::AbsSpec)

This method returns `true` if the affine scheme ``X`` is empty.
Otherwise, `false` is returned.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

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


@Markdown.doc """
    is_subset(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is a subset of ``Y`` based on the comparison of their coordinate rings.
See [`inclusion_morphism(::AbsSpec, ::AbsSpec)`](@ref) for the corresponding morphism.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)

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
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || return false
  return iszero(modulus(OO(Y)))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || return false
  return issubset(inverted_set(OO(Y)), units_of(R))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_coordinate_ring(Y) || return false
  return issubset(inverted_set(OO(Y)), units_of(R)) && iszero(ambient_closure_ideal(Y))
end


# (2.2) MPolyQuo in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = ambient_coordinate_ring(Y)
  R === ambient_coordinate_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuo}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return issubset(ambient_closure_ideal(Y), ambient_closure_ideal(X))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return iszero(localized_ring(OO(Y))(modulus(OO(X))))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return issubset(modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end


# (2.3) MPolyLocalizedRing in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == ambient_coordinate_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = ambient_coordinate_ring(Y)
  R == ambient_coordinate_ring(X) || return false
  return all(x->(iszero(OO(X)(x))), gens(modulus(OO(Y))))
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyLocalizedRing}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX)
end


function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX) && iszero(modulus(OO(Y)))
end


# (2.4) MPolyQuoLocalizedRing in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == ambient_coordinate_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = ambient_coordinate_ring(Y)
  R == ambient_coordinate_ring(X) || return false
  L = localized_ring(OO(X))
  return issubset(L(modulus(OO(Y))), modulus(OO(X)))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
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
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
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
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyComplementOfKPointIdeal}}) where {BRT}
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

@Markdown.doc """
    is_open_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is openly embedded in ``Y``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)

julia> is_open_embedding(Y, X)
false

julia> Z = hypersurface_complement(X, x1)
Spec of localization of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field at the powers of fmpq_mpoly[x1]

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
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
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
    X::Spec{BRT, <:MPolyQuoLocalizedRing},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  return OO(Y) == ambient_coordinate_ring(X) && all(iszero, gens(modulus(OO(X))))
end



####################################################################################
# (4) Check if a scheme can be embedded via a closed embeeded in another scheme
####################################################################################

@Markdown.doc """
    is_closed_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is closed embedded in ``Y``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X,ideal(R,[x1*x2]))
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1*x2)

julia> is_closed_embedding(Y, X)
true

julia> Z = hypersurface_complement(X, x1)
Spec of localization of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field at the powers of fmpq_mpoly[x1]

julia> is_closed_embedding(Z, X)
false
```
"""
function is_closed_embedding(X::AbsSpec, Y::AbsSpec)
  error("`is_closed_embedding(X, Y)` not implemented for X of type $(typeof(X)) and Y of type $(typeof(Y))")
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuo},
    Y::AbsSpec{<:Ring, <:MPolyRing}
  )
  return ambient_coordinate_ring(X) === ambient_coordinate_ring(Y)
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyRing},
    Y::AbsSpec{<:Ring, <:MPolyQuo}
  )
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  return iszero(modulus(OO(Y)))
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuo},
    Y::AbsSpec{<:Ring, <:MPolyQuo}
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
    X::AbsSpec{<:Ring, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
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
    X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}},

    Y::AbsSpec{<:Ring, <:MPolyQuo}
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
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  inverted_set(OO(X)) == inverted_set(OO(Y)) || return false
  J = localized_ring(OO(X))(modulus(underlying_quotient(OO(Y))))
  return issubset(J, modulus(OO(X)))
end


function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuo},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  OO(Y) === ambient_coordinate_ring(X) || return false
  return true
end


function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuo},
    Y::Spec{BRT, <:RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_coordinate_ring(X)
  R === ambient_coordinate_ring(Y) || return false
  all(x->(isunit(OO(X)(x))), denominators(inverted_set(OO(Y)))) || return false
  return issubset(modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end

#############################################################################
# (5) Check, if a Spec is equidimensional
#############################################################################

@doc Markdown.doc"""
   is_equidimensional(X::AbsSpec{<:Field, <:MPolyAnyRing}) 

Return whether a scheme `X` is equidimensional.

Currently this command is available for affine schemes and space germs.
TODO: projective schemes, covered schemes

This command relies on [`equidimensional_decomposition_radical`](@ref).

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal(R,[(x-y)])
ideal(x - y)

julia> J = ideal(R,[x-1,y-2])
ideal(x - 1, y - 2)

julia> X = Spec(R,I)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x - y)

julia> Y = Spec(R,I*J)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - x*y - x + y, x*y - 2*x - y^2 + 2*y)

julia> is_equidimensional(X)
true

julia> is_equidimensional(Y)
false
'''
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

@doc Markdown.doc"""
   is_reduced(X::AbsSpec{<:Field, <:MPolyAnyRing})

Return the boolean value whether a scheme `X` is reduced.

Currently, this command is available for affine schemes and space germs.
TODO: projective schemes, covered schemes
"""
@attr function is_reduced(X::AbsSpec{<:Ring, <:MPAnyQuoRing})
  I = saturated_ideal(modulus(OO(X)))
  return is_reduced(quo(base_ring(I), I)[1])
end

## make is_reduced agnostic to quotient ring
@attr Bool function is_reduced(X::AbsSpec{<:Ring, <:MPAnyNonQuoRing})
  return true
end

########################################################################
# (7) Smoothness test based on projective modules.
# The routine checks whether the module for the cotangent sheaf Ω¹(X)
# is locally free over 𝒪(X) and returns `true` if this is the case. 
########################################################################

@doc Markdown.doc"""
    is_smooth(X::AbsSpec{<:Field, <:MPolyAnyRing})

Return whether a scheme `X` is smooth.

Currently this command is available for affine schemes and space germs.
TODO: Covered schemes, projective schemes

Note that smoothness and regularity do not coincide over non-perfect fields.
TODO: is_regular using Hironaka's criterion

See also [`singular_locus`](@ref), [`singular_locus_reduced`](@ref).
# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal(R,[x-y^2])
ideal(x - y^2)

julia> J = ideal(R,[x^2-y^2])
ideal(x^2 - y^2)

julia> X = Spec(R, I)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x - y^2)

julia> is_smooth(X)
true

julia> Y = Spec(R, J)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - y^2)

julia> is_smooth(Y)
false

julia> U=MPolyComplementOfKPointIdeal(R,[1,1])
complement of maximal ideal corresponding to point with coordinates fmpq[1, 1]

julia> Z = Spec(R, J, U)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - y^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[1, 1]

julia> is_smooth(Z)
true

```
"""
@attr Bool function is_smooth(X::AbsSpec{<:Field, <:MPolyQuoLocalizedRing})
  R = base_ring(OO(X))
  L = localized_ring(OO(X))
  I = modulus(OO(X))
  f = gens(saturated_ideal(I))
  Df = jacobi_matrix(f)
  A = map_entries(x->OO(X)(x), Df)
  success, _, _ = Oscar._is_projective_without_denominators(A)
  return success
end

@attr Bool function is_smooth(X::AbsSpec{<:Field, <:MPolyQuo})
  R = base_ring(OO(X))
  I = modulus(OO(X))
  f = gens(I)
  Df = jacobi_matrix(f)
  A = map_entries(x->OO(X)(x), Df)
  success, _, _ = Oscar._is_projective_without_denominators(A)
  return success
end

## make is_smooth agnostic to quotient ring
is_smooth(X::AbsSpec{<:Field, <:MPolyRing}) = true
is_smooth(X::AbsSpec{<:Field, <:MPolyLocalizedRing}) = true
