###
# Computing tropical Groebner bases in Oscar
# ==========================================
#
# For a definition of tropical Groebner basis see Section 2.4 in:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
# To see how they can be computed using standard bases see
#   T. Markwig, Y. Ren: Computing tropical varieties over fields with valuation
###


###
# temporary workarounds:
###
function symbols(Kt::AbstractAlgebra.Generic.RationalFunctionField{K} where {K})
  return Kt.S
end

function change_base_ring(Rtx::FmpzMPolyRing,I::MPolyIdeal{Ktx} where {Ktx})
  return ideal([change_base_ring(Rtx,f) for f in gens(I)])
end

# function which coerces a polynomial in QQ(t)[x_1,...,x_n] into ZZ[t,x_1,...,x_n]
# Example:
# Kt,t = RationalFunctionField(QQ,"t")
# Ktx,(x1,x2,x3) = PolynomialRing(Kt,3)
# f = x1+(1+t)*x2+(1//2+1//3*t+1//5*t^2)*x3
# Rtx,(t,y1,y2,y3) = PolynomialRing(ZZ,4)
# change_base_ring(Rtx,f)
function change_base_ring(Rtx::FmpzMPolyRing,f::AbstractAlgebra.Generic.MPoly{Kt} where {Kt})

  R = coefficient_ring(Rtx)
  fRtx = zero(Rtx)

  for i in 1:length(f)
    expvKtx = exponent_vector(f,i) # exponent vector in K(t)[x1,...,xn]
    expvRtx = vcat([0],expvKtx)    # exponent vector in R[t,x1,...,xn]
    cKt = coeff(f,i)               # coefficient in K(t)
    @assert isone(denominator(cKt)) "change_base_ring: coefficient denominators need to be 1"

    cKt = numerator(cKt)           # coefficient in K[t]
    cK = coefficients(cKt)         # vector in K
    M = lcm([denominator(c) for c in cK])
    cR = [R(M*c) for c in cK]  # vector in R

    for c in cR
      fRtx += c*monomial(Rtx,expvRtx)
      expvRtx[1] += 1
    end
  end

  return fRtx
end


function simulate_valuation(I,val_t::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K}, AbstractAlgebra.Generic.Rat{K}} where {K})

  Ktx = base_ring(I)
  Kt = base_ring(Ktx)
  K = base_ring(Kt)

  Rtx = PolynomialRing(ZZ,vcat(symbols(Kt),symbols(Ktx)))


  @warn "virtual_valuation_ring: work in progress, not fully functional yet"
  return nothing
end
export simulate_valuation


function groebner_basis(I,w,val)
  # 1: construct a valuation simulating ring

  @warn "groebner_basis: work in progress, not fully functional yet"
  return nothing
end
