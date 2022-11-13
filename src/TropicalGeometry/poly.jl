###
# Tropical polynomial computations in Oscar
# =========================================
###
export tropical_polynomial

@doc Markdown.doc"""
    tropical_polynomial(f::MPolyElem,M::Union{typeof(min),typeof(max)}=min)

Given a polynomial `f` over a field with an intrinsic valuation (i.e., a field
on which a function `valuation` is defined such as `PadicField(7,2)`),
returns the tropicalization of `f` as a polynomial over the min tropical semiring
(default) or the max tropical semiring.

# Examples
```jldoctest
julia> K = PadicField(7, 2)
Field of 7-adic numbers

julia> Kxy, (x,y) = K["x", "y"]
(Multivariate Polynomial Ring in x, y over Field of 7-adic numbers, AbstractAlgebra.Generic.MPoly{padic}[x, y])

julia> f = 7*x+y+49
(7^1 + O(7^3))*x + y + 7^2 + O(7^4)

julia> tropical_polynomial(f,min)
(1)*x + y + (2)

julia> tropical_polynomial(f,max)
(-1)*x + y + (-2)
```
"""
function tropical_polynomial(f::MPolyElem, M::Union{typeof(min),typeof(max)}=min)
  T = TropicalSemiring(M)
  if M==min
    s=1
  else
    s=-1
  end

  Tx,x = PolynomialRing(T,[repr(x) for x in gens(parent(f))])
  tropf = inf(T)

  if base_ring(parent(f)) isa NonArchLocalField
    for (c,alpha) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
      tropf = tropf + T(s*valuation(c))*monomial(Tx,alpha)
    end
  else
    for alpha in AbstractAlgebra.exponent_vectors(f)
      tropf = tropf + T(0)*monomial(Tx,alpha)
    end
  end

  return tropf
end



@doc Markdown.doc"""
    tropical_polynomial(f::MPolyElem,val::TropicalSemiringMap)

Given a polynomial `f` and a tropical semiring map `val`,
returns the tropicalization of `f` as a polynomial over the tropical semiring.

# Examples
```jldoctest
julia> R, (x,y) = PolynomialRing(QQ,["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> val = TropicalSemiringMap(QQ,7)
The 7-adic valuation on Rational Field

julia> f = 7*x+y+49
7*x + y + 49

julia> tropical_polynomial(f,val)
(1)*x + y + (2)
```
"""
function tropical_polynomial(f::MPolyElem, val::TropicalSemiringMap)
  T = TropicalSemiring(val)
  Tx,x = PolynomialRing(T,[repr(x) for x in gens(parent(f))])
  tropf = inf(T)

  for (c,alpha) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
    tropf = tropf + val(c)*monomial(Tx,alpha)
  end

  return tropf
end
