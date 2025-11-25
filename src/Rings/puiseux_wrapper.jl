########################################################################
# Wrapper functionality for `PuiseuxExpansions.lib` in Singular        #
########################################################################

# A wrapper for the raw output of `puiseux`.
# It returns a list of lists. It seems that the outer list contains 
# only one element. The inner list contains 
#   * a polynomial `h`
#   * an integer `e`
#   * more stuff that I can not make sense of at the moment; probably some multiplicities
# We interpret the result as `h` being the polynomial which evaluates 
# to the Puiseux expansion of `f` on `y^(1//e)`. 
function _puiseux(f::MPolyRingElem{T}, max_deg::Int, at_origin::Bool=true) where {T<:FieldElem}
  P = parent(f)
  SP = singular_poly_ring(P)
  Sf = SP(f)
  raw = Singular.LibPuiseuxexpansions.puiseux(Sf, max_deg, Int(at_origin))
  result = Tuple{typeof(f), Int}[(P(h), e) for (h, e, _) in raw]
end

# Method to create the default Oscar parent for the Puiseux expansion of `f`. 
function _puiseux_parent(f::MPolyRingElem, prec::Int)
  R = parent(f)
  kk = coefficient_ring(R)
  symbs = symbols(R)
  return puiseux_series_ring(kk, prec, symbs[1])[1]
end

@doc raw"""
    function puiseux_expansion(
        f::MPolyRingElem{T}, 
        precision::Int=10;
        parent = _puiseux_parent(f, precision)
      ) where {T <: FieldElem}

Compute the Puiseux expansion of `f` up to degree `precision`. 
A parent ring for the return value can be provided if needed.
```jldoctest
julia> R, (x, y) = QQ[:x, :y]
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> f = y^3 + x^2 + x^8
x^8 + x^2 + y^3

julia> h = Oscar.puiseux_expansion(f, 20)
-x^(2//3) - 1//3*x^(20//3) + O(x^(22//3))

julia> evaluate(f, [gen(parent(h)), h])
O(x^(26//3))
```
"""
function puiseux_expansion(
    f::MPolyRingElem{T}, 
    precision::Int=10;
    parent = _puiseux_parent(f, precision)
  ) where {T <: FieldElem}
  R = Oscar.parent(f)
  @assert ngens(R) == 2 "polynomial must be bivariate"
  x, y = gens(R)
  xx = gen(parent)
  prec = max_precision(parent)
  pdat = _puiseux(f, prec, true)
  @assert isone(length(pdat)) "a longer list was returned than anticipated"
  h, e = only(pdat)
  # h should have no terms involving y^k for k > 0
  hh = evaluate(h, [xx^(1//e), zero(xx)])
  return hh
end

