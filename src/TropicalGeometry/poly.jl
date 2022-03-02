###
# Tropical polynomial computations in Oscar
# =========================================
###

@doc Markdown.doc"""
    tropical_polynomial{M}(f::AbstractAlgebra.Generic.MPoly{Any})

Returns the tropicalization of a polynomial in form of a polynomial over the tropical numbers.
If M=min, returns a polynomial over the min-plus semiring.
If M=max, returns a polynomial over the max-plus semiring.

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
function tropical_polynomial(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, M::Union{typeof(min),typeof(max)}=min)
  T = tropical_semiring(M)
  if M==min
    s=1
  else
    s=-1
  end

  Tx,x = PolynomialRing(T,[repr(x) for x in gens(parent(f))])
  tropf = inf(T)

  if base_ring(parent(f)) isa NonArchLocalField
    for (c,alpha) in zip(coefficients(f),exponent_vectors(f))
      tropf = tropf + T(s*valuation(c))*monomial(Tx,alpha)
    end
  else
    for alpha in exponent_vectors(f)
      tropf = tropf + T(0)*monomial(Tx,alpha)
    end
  end

  return tropf
end
function tropical_polynomial(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, val::ValuationMap)
  T = val.tropical_semiring
  Tx,x = PolynomialRing(T,[repr(x) for x in gens(parent(f))])
  tropf = inf(T)

  for (c,alpha) in zip(coefficients(f),exponent_vectors(f))
    tropf = tropf + val(c)*monomial(Tx,alpha)
  end

  return tropf
end
export tropical_polynomial
