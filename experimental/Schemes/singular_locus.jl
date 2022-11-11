export is_equidimensional
export singular_locus, singular_locus_reduced
export reduced_scheme
export is_smooth


#######################################################################
# Functionality for smoothness/singular locus
#######################################################################

#######################################################################
# MOVE TO Methods
#######################################################################

### TODO: The following two functions (singular_locus, 
###       singular_locus_reduced) need to be made type-sensitive 
###       and reduced=true needs to be set automatically for varieties 
###       as soon as not only schemes, but also varieties as special 
###       schemes have been introduced in OSCAR
### TODO: Make singular locus also available for projective schemes and
###       for covered schemes (the workhorse is already there...).
 
@doc Markdown.doc"""
    singular_locus(X::Scheme) -> Scheme

Return the singular locus of `X`.

For computing the singular locus of the reduced scheme induced by `X`,
please use [`singular_locus_reduced`](@ref).

Currently this command is available for affine schemes and for space germs
(in both settings requiring that these are defined over a field).
TODO: Covered schemes, projective schemes

Over non-perfect fields, this command returns the non-smooth locus and `X`
may still be regular at some points of the returned subscheme.

See also [`is_smooth`](@ref).

# Examples
``` jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"]
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> I = ideal(R, [x^2 - y^2 + z^2])
ideal(x^2 - y^2 + z^2)

julia> A3 = Spec(R)
Spec of Multivariate Polynomial Ring in x, y, z over Rational Field

julia> X = Spec(R,I)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2)

julia> singular_locus(A3)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(1)

julia> singular_locus(X)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2, 2*x, -2*y, 2*z)

julia> U = MPolyComplementOfKPointIdeal(R,[0,0,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> Y = Spec(R,I,U)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

julia> singular_locus(Y)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^2 - y^2 + z^2, x^2 - y^2 + z^2, 2*x, -2*y, 2*z, x^2 - y^2 + z^2, x^2 - y^2 + z^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0, 0]

```
"""
function singular_locus(X::AbsSpec{<:Field, <:MPAnyQuoRing})
  comp = _singular_locus_with_decomposition(X,false)
  if length(comp) == 0 
    return subscheme(X, ideal(OO(X),one(OO(X))))
  end
  R = base_ring(OO(X))
  I= prod([modulus(quotient_ring(OO(Y))) for Y in comp])
  return subscheme(X,I)
end

# make singular_locus agnostic to quotient
function singular_locus(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  return subscheme(X,ideal(OO(X),[one(OO(X))]))
end

@doc Markdown.doc"""
    singular_locus_reduced(X::Scheme) -> Scheme

Return the singular locus of the reduced scheme ``X_{red}`` induced by `X`.

For computing the singular locus of `X` itself, please use 
['singular_locus](@ref)'.

Currently this command is available for affine schemes and for space germs
(in both settings requiring that these are defined over a field).
TODO: Covered schemes, projective schemes

Over non-perfect fields, this command returns the non-smooth locus and
`X_{red}` may still be regular at some points of the returned subscheme.

See also [`is_smooth`](@ref).

# Examples
``` jldoctest
julia> R, (x,y,z) = QQ["x", "y", "z"]
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> I = ideal(R, [(x^2 - y^2 + z^2)^2])
ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4)

julia> X = Spec(R,I)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4)

julia> singular_locus_reduced(X)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4, z, y, x)

julia> singular_locus(X)
Spec of Quotient of Multivariate Polynomial Ring in x, y, z over Rational Field by ideal(x^4 - 2*x^2*y^2 + 2*x^2*z^2 + y^4 - 2*y^2*z^2 + z^4, 4*x^3 - 4*x*y^2 + 4*x*z^2, -4*x^2*y + 4*y^3 - 4*y*z^2, 4*x^2*z - 4*y^2*z + 4*z^3)

```
"""

function singular_locus_reduced(X::AbsSpec{<:Ring, <:MPAnyQuoRing})
  comp =  _singular_locus_with_decomposition(X, true)
  I= ideal(localized_ring(OO(X)),one(localized_ring(OO(X))))
  for Z in comp
    I = intersect(I, modulus(OO(Z)))
  end     
  return subscheme(X,I)
end

# make singular_locus_reduced agnostic to quotient
function singular_locus_reduced(X::AbsSpec{<:Field, <:MPAnyNonQuoRing})
  return subscheme(X,ideal(OO(X),[one(OO(X))]))
end

# internal workhorse, not user-facing
function _singular_locus_with_decomposition(X::AbsSpec{<:Ring, <:MPAnyQuoRing}, reduced::Bool=true)
  I = saturated_ideal(modulus(OO(X)))
  result = typeof(X)[]

  P = []
  if has_attribute(X, :is_equidimensional) && is_equidimensional(X) && !reduced 
    push!(P, I)
  else 
    if reduced
      P = equidimensional_decomposition_radical(I)
    else
      P = equidimensional_decomposition_weak(I)
    end
  end

  if length(P)==1 && !reduced
    d = dim(X)
    R = base_ring(I)
    n = nvars(R) 
    M = _jacobi_matrix_modulus(X)
    minvec = minors(M, n-d)
    J = ideal(R, minvec)
    JX = ideal(OO(X),minvec)
    one(OO(X)) in JX && return result
    return [subscheme(X, J)]
  else
    components = [subscheme(X, J) for J in P]
    for i in 1:length(components)
      for j in (i+1):length(components)
        W = intersect(components[i],components[j])
        if !isempty(W) 
          push!(result, W)
        end
      end
    end
    for Y in components
      result = vcat(result, singular_locus(Y))
    end
  end
  return result
end

##########################################################################
# Basic properties of schemes
##########################################################################

##########################################################################
# MOVE TO: Attributes
##########################################################################
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

#######################################################################
# MOVE TO Methods
#######################################################################
@Markdown.doc """
   reduced_scheme(X::AbsSpec{<:Field, <:MPolyAnyRing})

Return the induced reduced scheme of a `X`.

Currently, this command is available for affine schemes and space germs.
TODO: projective schemes, covered schemes
 
This command relies on [`radical`](@ref).

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> J = ideal(R,[(x-y)^2])
ideal(x^2 - 2*x*y + y^2)

julia> X = Spec(R,J)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - 2*x*y + y^2)

julia> U = MPolyComplementOfKPointIdeal(R,[0,0])
complement of maximal ideal corresponding to point with coordinates fmpq[0, 0]

julia> Y = Spec(R,J,U)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x^2 - 2*x*y + y^2) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0]

julia> reduced_scheme(X)
Spec of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x - y)

julia> reduced_scheme(Y)
Spec of Localization of Quotient of Multivariate Polynomial Ring in x, y over Rational Field by ideal(x - y) at the multiplicative set complement of maximal ideal corresponding to point with coordinates fmpq[0, 0]

```
"""
@attr function reduced_scheme(X::AbsSpec{<:Field, <:MPolyQuoLocalizedRing})
  I = modulus(OO(X))
  J = radical(pre_saturated_ideal(I))
  return Spec(base_ring(J), J, inverted_set(OO(X)))
end

@attr function reduced_scheme(X::AbsSpec{<:Field, <:MPolyQuo})
  J = radical(modulus(OO(X)))
  return Spec(base_ring(J), J)
end

## to make reduced_scheme agnostic for quotient ring
@attr AbsSpec function reduced_scheme(X::AbsSpec{<:Ring, <:MPAnyNonQuoRing})
  return X
end

#######################################################################
# MOVE TO Attributes
#######################################################################
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
# Additional functionality for fractions                               #
########################################################################

########################################################################
# The following functions (derivative, jacobi_matrix) should go to 
# mpoly-localizations.jl and mpolyquo-localizations.jl
# QUESTION:
# documentation for all MPolyAnyRingElem should be happen in to mpoly.jl 
########################################################################
function derivative(f::MPolyLocalizedRingElem, i::Int)
  num = derivative(numerator(f), i)*denominator(f) - derivative(denominator(f), i)*numerator(f)
  den = denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function derivative(f::MPolyQuoLocalizedRingElem, i::Int)
  num = derivative(lifted_numerator(f), i)*lifted_denominator(f) - derivative(lifted_denominator(f), i)*lifted_numerator(f)
  den = lifted_denominator(f)^2
  g = gcd(num, den)
  return parent(f)(divexact(num, g), divexact(den, g), check=false)
end

function jacobi_matrix(f::MPolyLocalizedRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyLocalizedRingElem})
  R = parent(g[1])
  n = nvars(base_ring(R))
  @assert all(x->parent(x) == R, g)
  return matrix(R, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

function jacobi_matrix(f::MPolyQuoLocalizedRingElem)
  L = parent(f)
  n = nvars(base_ring(L))
  return matrix(L, n, 1, [derivative(f, i) for i=1:n])
end

function jacobi_matrix(g::Vector{<:MPolyQuoLocalizedRingElem})
  L = parent(g[1])
  n = nvars(base_ring(L))
  @assert all(x->parent(x) == L, g)
  return matrix(L, n, length(g), [derivative(x, i) for i=1:n for x = g])
end

############################################################################
# jacobi_matrix_modulus should go whereever singular_locus goes, as it is
# an internal helper (less expensive than jacobi_matrix)
############################################################################

## compute some representative of the jacobian matrix of a set of generators 
## of the modulus -- forgetting about the denominators of the generators
## because any contribution of their derivative is killed by the modulus 
function _jacobi_matrix_modulus(X::AbsSpec{<:Ring, <:MPAnyQuoRing})
  g = gens(modulus(quotient_ring(OO(X))))
  L = base_ring(quotient_ring(OO(X)))
  n = nvars(L)
  M = matrix(L, n, length(g),[derivative(f,i) for i=1:n for f in g])
  return M
end

########################################################################
# Smoothness test based on projective modules.
#
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
  f = gens(Oscar.pre_saturated_ideal(I))
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
