########################################################################
# Constructors from ideals                                             #
########################################################################

# TODO: Write one dummy constructor for the documentation with an ideal.
@doc raw"""
    function AffineSchemeOpenSubscheme(X::AbsAffineScheme, I::MPolyAnyIdeal) -> AffineSchemeOpenSubscheme

Return the complement of the zero locus of ``I`` in ``X``.

# Examples
```jldoctest
julia> P, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> A = spec(P)
Spectrum
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field

julia> AffineSchemeOpenSubscheme(A, I)
Open subset
  of affine 3-space
complement to V(x^3 - y^2*z)
```
"""
function AffineSchemeOpenSubscheme(X::AbsAffineScheme, I::MPolyLocalizedIdeal; check::Bool=true)
  base_ring(I) === OO(X) || error("Ideal does not belong to the correct ring")
  g = [numerator(a) for a in gens(I) if !iszero(numerator(a))]
  return AffineSchemeOpenSubscheme(X, g, check=check)
end

function AffineSchemeOpenSubscheme(X::AbsAffineScheme, I::MPolyQuoLocalizedIdeal; check::Bool=true)
  base_ring(I) === OO(X) || error("Ideal does not belong to the correct ring")
  g = [lifted_numerator(a) for a in gens(I) if !iszero(numerator(a))]
  return AffineSchemeOpenSubscheme(X, g, check=check)
end

function AffineSchemeOpenSubscheme(X::AbsAffineScheme, I::MPolyIdeal; check::Bool=true)
  return AffineSchemeOpenSubscheme(X, [g for g in gens(I) if !iszero(OO(X)(g))], check=check)
end

function AffineSchemeOpenSubscheme(X::AbsAffineScheme, I::MPolyQuoIdeal; check::Bool=true)
  return AffineSchemeOpenSubscheme(X, [lift(g) for g in gens(I) if !iszero(OO(X)(g))], check=check)
end

########################################################################
# Constructors from closed subvarieties                                #
########################################################################

@doc raw"""
    complement(X::Scheme, Y::Scheme) -> Scheme

Return the complement ``X \ Y`` of ``Y`` in ``X``.

Since we want the complement `U = X \ Y` to have a well defined scheme
structure,  we require that `Y` is closed in `X`.

# Examples
```jldoctest
julia> P, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> A = spec(P)
Spectrum
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field

julia> Y = spec(P, I)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x, y, z
      over rational field
    by ideal (x^3 - y^2*z)

julia> complement(A, Y)
Open subset
  of affine 3-space
complement to V(x^3 - y^2*z)
```
"""
complement(X::Scheme,Y::Scheme)

function complement(X::AbsAffineScheme, Z::AbsAffineScheme{<:Ring, <:MPolyRing})
  ambient_coordinate_ring(X) == ambient_coordinate_ring(Z) || error("X and Z do not compare")
  return EmptyScheme(base_ring(X))
end

function complement(X::AbsAffineScheme, Z::AbsAffineScheme{<:Ring, <:MPolyQuoRing})
  ambient_coordinate_ring(X) == ambient_coordinate_ring(Z) || error("X and Z do not compare")
  return AffineSchemeOpenSubscheme(X, modulus(OO(Z)))
end

function complement(X::AbsAffineScheme,
    Z::AbsAffineScheme{<:Ring, <:MPolyQuoLocRing};
    check::Bool=true
  )
  @check is_closed_embedding(Z, X) "not a closed embedding"
  return AffineSchemeOpenSubscheme(Y, modulus(underlying_quotient(OO(Z))))
end

########################################################################
# Conversion from AbsAffineScheme                                              #
########################################################################

AffineSchemeOpenSubscheme(X::AbsAffineScheme) = AffineSchemeOpenSubscheme(X, [one(ambient_coordinate_ring(X))], check=false)

########################################################################
# Additional constructors                                              #
########################################################################

@doc raw"""
    product(U::AffineSchemeOpenSubscheme, Z::AbsAffineScheme) -> AffineSchemeOpenSubscheme, AffineSchemeOpenSubschemeMor, AffineSchemeOpenSubschemeMor

Given a Zariski open subset `U` of an affine scheme `X`, complement to
a subscheme `Y` of `X`, and given an affine scheme `Z`, return the product
scheme $V := U \times Z$ as a Zariski open subset of the affine scheme
$X\times Z$ complement to the closed subscheme $Y\times Z$.

It is returned with the two product maps $U \to V$ and $Z \to V$.

# Examples
```jldoctest
julia> P, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> P2, (x2, y2, z2) = polynomial_ring(QQ, [:x2, :y2, :z2]);

julia> P1, (x1, y1, z1) = polynomial_ring(QQ, [:x1, :y1, :z1]);

julia> A1 = spec(P1)
Spectrum
  of multivariate polynomial ring in 3 variables x1, y1, z1
    over rational field

julia> A2 = spec(P2)
Spectrum
  of multivariate polynomial ring in 3 variables x2, y2, z2
    over rational field

julia> I1 = ideal([x1^3-y1^2*z1]);

julia> I2 = ideal([x2^2-y2^2+z2^2]);

julia> Y1 = spec(P1, I1);

julia> Y2 = spec(P2, I2);

julia> U2 = complement(A2, Y2);

julia> V, iU, iZ = product(U2, Y1)
(Complement to V(x2^2 - y2^2 + z2^2) in affine scheme with coordinates [x2, y2, z2, x1, y1, z1], Hom: complement to V(x2^2 - y2^2 + z2^2) in affine scheme with coordinates [x2, y2, z2, x1, y1, z1] -> complement to V(x2^2 - y2^2 + z2^2) in affine scheme with coordinates [x2, y2, z2], Hom: complement to V(x2^2 - y2^2 + z2^2) in affine scheme with coordinates [x2, y2, z2, x1, y1, z1] -> complement to V(1) in affine scheme with coordinates [x1, y1, z1])

julia> V
Open subset
  of affine 6-space over QQ with coordinates [x2, y2, z2, x1, y1, z1]
complement to V(x2^2 - y2^2 + z2^2)

julia> iU
Affine scheme open subscheme morphism
  from [x2, y2, z2, x1, y1, z1]  complement to V(x2^2 - y2^2 + z2^2)
  to   [x2, y2, z2]              complement to V(x2^2 - y2^2 + z2^2)
defined by the map
  affine scheme morphism
    from [x2, y2, z2, x1, y1, z1]  scheme(x1^3 - y1^2*z1) \ scheme(x2^2 - y2^2 + z2^2)
    to   [x2, y2, z2]              affine 3-space
  given by the pullback function
    x2 -> x2
    y2 -> y2
    z2 -> z2

julia> iZ
Affine scheme open subscheme morphism
  from [x2, y2, z2, x1, y1, z1]  complement to V(x2^2 - y2^2 + z2^2)
  to   [x1, y1, z1]              complement to V(1)
defined by the map
  affine scheme morphism
    from [x2, y2, z2, x1, y1, z1]  scheme(x1^3 - y1^2*z1) \ scheme(x2^2 - y2^2 + z2^2)
    to   [x1, y1, z1]              scheme(x1^3 - y1^2*z1)
  given by the pullback function
    x1 -> x1
    y1 -> y1
    z1 -> z1
```
"""
function product(U::AffineSchemeOpenSubscheme, Y::AbsAffineScheme)
  X = ambient_scheme(U)
  P, pX, pY = product(X, Y)
  V = AffineSchemeOpenSubscheme(P, lifted_numerator.(pullback(pX).(complement_equations(U))))
  res_pX = restrict(pX, V, U, check=false)
  res_pY = restrict(pY, V, AffineSchemeOpenSubscheme(Y), check=false)
  return V, res_pX, res_pY
end

function subscheme(U::AffineSchemeOpenSubscheme, I::Ideal)
  Z = subscheme(ambient_scheme(U), I) #Takes care of coercion and complains if necessary
  return AffineSchemeOpenSubscheme(Z, [g for g in complement_equations(U) if !iszero(OO(Z)(g))])
end

