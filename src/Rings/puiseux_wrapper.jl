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
  @assert length(raw) == 1 "longer list then expected for output"
  return [_process_result(P, data...) for data in raw]
end

function _process_result(P::MPolyRing, prec::Int, SE::Singular.PolyRing, res::Dict, rest...)
  # create the necesessary field extension on the Oscar side
  kk = coefficient_ring(SE) # the extension field of QQ
  mp = Singular.modulus(kk) # the minimum polynomial
  mmp = Singular.n_transExt_to_spoly(mp) # convert into an actual polynomial
  L, t = polynomial_ring(QQ, :t; cached=false)
  kkO, alpha = extension_field(L(mmp))
  P_ext, to_P_ext = change_base_ring(kkO, P)
  Puis, xx = puiseux_series_ring(kkO, prec, symbols(P_ext)[1])

  # transfer the polynomials to the oscar side
  result = elem_type(Puis)[]
  for branch in res[:PE]
    h, e, _ = branch
    ctx = MPolyBuildCtx(P_ext)
    for (c, v) in zip(Singular.coefficients(h), Singular.exponent_vectors(h))
      push_term!(ctx, kkO(c), v)
    end
    hh = finish(ctx)
    push!(result, evaluate(hh, [xx^(1//3), zero(xx)]))
  end
  return result
end

function _process_result(P::MPolyRing, prec::Int, h::Singular.spoly, e::Int, rest...)
  kk = coefficient_ring(P)
  Puis, xx = puiseux_series_ring(kk, prec, symbols(P)[1])
  return [evaluate(P(h), [xx^(1//e), zero(xx)])]
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
        precision::Int=10
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
    precision::Int=10
  ) where {T <: FieldElem}
  R = Oscar.parent(f)
  @assert ngens(R) == 2 "polynomial must be bivariate"
  x, y = gens(R)
  xx = gen(parent)
  prec = max_precision(parent)

  # prepare for the Singular call
  SR = singular_poly_ring(R)
  Sf = SR(f)
  raw = Singular.LibPuiseuxexpansions.puiseux(Sf, prec, 1)
  return reduce(vcat, [_process_result(R, precision, data...) for data in raw])
end

