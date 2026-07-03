########################################################################
# Wrapper functionality for `PuiseuxExpansions.lib` in Singular        #
########################################################################

# A wrapper for Singular's `Singular.LibPuiseuxexpansions.puiseux`. 
# It returns either a list of lists with entries of the following form: 
#   * a polynomial `h`
#   * an integer `e`
#   * more stuff that is returned for internal purposes in the Singular library. 
#     We discard this information
# We interpret the result as `h` being the polynomial which evaluates 
# to the Puiseux expansion of `f` on `y^(1//e)`. 
#
# Otherwise, it returns a ring and a dictionary with additional information 
# on required extension fields for presentation of the result. Which one 
# of the two cases applies is not foreseeable from the input alone. We 
# therefore catch it with dispatch on `_process_result`. 

function _process_result(P::MPolyRing, prec::Int, SE::Singular.PolyRing, res::Dict, rest...)
  # create the necesessary field extension on the Oscar side
  kk = coefficient_ring(SE) # the extension field of QQ
  mp = Singular.modulus(kk) # the minimum polynomial
  mmp = Singular.n_transExt_to_spoly(mp) # convert into an actual polynomial
  L, t = polynomial_ring(QQ, :t; cached=false)
  kkO, alpha = extension_field(L(mmp))
  P_ext, to_P_ext = change_base_ring(kkO, P)
  Puis, xx = puiseux_series_ring(kkO, prec, symbols(P_ext)[1]; cached=false)

  # transfer the polynomials to the oscar side
  result = elem_type(Puis)[]
  for branch in res[:PE]
    h, e, _ = branch
    ctx = MPolyBuildCtx(P_ext)
    for (c, v) in zip(Singular.coefficients(h), Singular.exponent_vectors(h))
      push_term!(ctx, kkO(c), v)
    end
    hh = finish(ctx)
    push!(result, evaluate(hh, [xx^(1//e), zero(xx)]))
  end
  return result
end

function _process_result(P::MPolyRing, prec::Int, h::Singular.spoly, e::Int, rest...)
  kk = coefficient_ring(P)
  Puis, xx = puiseux_series_ring(kk, prec, symbols(P)[1]; cached=false)
  return [evaluate(P(h), [xx^(1//e), zero(xx)])]
end


@doc raw"""
    function puiseux_expansion(
        f::MPolyRingElem{T},
        max_ord::Int=10;
        precision::Int=max_ord
      ) where {T <: QQFieldElem}

Compute the Puiseux expansion of `f` up to degree `max_ord` and returns the output 
in puiseux series rings with the given `precision`. 

!!! note 
    For the moment we only support bivariate polynomials over `QQ`. Puiseux expansion may produce number fields as coefficient rings for the output. As this is not foreseeable a priori, no guarantee can be given about the parent objects for the output. 
```jldoctest
julia> R, (x, y) = QQ[:x, :y]
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> f = y^3 + x^2 + x^8
x^8 + x^2 + y^3

julia> h = Oscar.puiseux_expansion(f, 15)
1-element Vector{AbstractAlgebra.Generic.PuiseuxSeriesFieldElem{QQFieldElem}}:
 -x^(2//3) + O(x^(17//3))

julia> all(is_zero(evaluate(f, [gen(parent(h)), h])) for h in h)
true
```
"""
function puiseux_expansion(
    f::MPolyRingElem{T}, 
    max_ord::Int=10;
    precision::Int=max_ord
  ) where {T <: QQFieldElem}
  R = Oscar.parent(f)
  @assert ngens(R) == 2 "polynomial must be bivariate"

  # prepare for the Singular call
  SR = singular_poly_ring(R)
  Sf = SR(f)
  raw = Singular.LibPuiseuxexpansions.puiseux(Sf, max_ord, 1)
  return reduce(vcat, [_process_result(R, precision, data...) for data in raw])
end

