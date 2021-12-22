###
# Computing tropical Groebner bases in Oscar
# ==========================================
#
# For a definition of tropical Groebner basis see Section 2.4 in:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
# To see how they can be computed using standard bases see:
#   T. Markwig, Y. Ren: Computing tropical varieties over fields with valuation
#
###


###
# temporary workarounds:
###
function symbols(Kt::AbstractAlgebra.Generic.RationalFunctionField{K} where {K})
  return Kt.S
end



#=======
function which coerces a polynomial in QQ(t)[x_1,...,x_n] into ZZ[t,x_1,...,x_n]
Example:
Kt,t = RationalFunctionField(QQ,"t")
Ktx,(x1,x2,x3) = PolynomialRing(Kt,["x1","x2","x3"])
f = x1+(1+t)*x2+(1//2+1//3*t+1//5*t^2)*x3
Rtx,(t,y1,y2,y3) = PolynomialRing(QQ,["t","x1","x2","x3"])
tropical_change_base_ring(Rtx,f)
=======#
function tropical_change_base_ring(Rtx::FmpqMPolyRing,f::AbstractAlgebra.Generic.MPoly{Kt} where {Kt})

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



#=======
function which coerces a polynomial in QQ[x_1,...,x_n] into ZZ[t,x_1,...,x_n]
Example:
Kx,(x1,x2,x3) = PolynomialRing(QQ,3)
f = x1+(1+2)*x2+(1+3*2+5*2^2)*x3
Rtx,(t,y1,y2,y3) = PolynomialRing(ZZ,4)
tropical_change_base_ring(Rtx,f)
=======#
function tropical_change_base_ring(Rtx::FmpzMPolyRing,f::fmpq_mpoly)

  # todo: rewrite to use MPolyBuildCtx

  R = coefficient_ring(Rtx)
  fRtx = zero(Rtx)

  for i in 1:length(f)
    expvKx = exponent_vector(f,i) # exponent vector in K[x1,...,xn]
    expvRtx = vcat([0],expvKx)    # exponent vector in R[t,x1,...,xn]
    cK = coeff(f,i)               # coefficient in K
    @assert isone(denominator(cK)) "change_base_ring: coefficient denominators need to be 1"

    cR = numerator(cK)            # coefficient in R
    fRtx += cR*monomial(Rtx,expvRtx)
  end

  return fRtx
end
function tropical_change_base_ring(Rtx::FmpzMPolyRing,I::MPolyIdeal{Ktx} where {Ktx})
  return ideal([tropical_change_base_ring(Rtx,f) for f in gens(I)])
end
export tropical_change_base_ring



#=======
functions which, given an ideal I in variables x1, ..., xn over a field with valuation,
returns an ideal vvI in variables t, x1, ..., xn such that tropical Groebner bases of I w.r.t. w
correspond to standard bases of I w.r.t. (-1,w)
=======#
function simulate_valuation(I,val_t::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K}, AbstractAlgebra.Generic.Rat{K}} where {K})

  Ktx = base_ring(I)
  Kt = base_ring(Ktx)
  K = base_ring(Kt)

  Rtx = PolynomialRing(K,vcat(symbols(Kt),symbols(Ktx)));
  vvI = tropical_change_base_ring(Rtx,I)

  return vvI
end
function simulate_valuation(I,val_p::ValuationMap{FlintRationalField, fmpz})

  Kx = base_ring(I)
  K = coefficient_ring(Kx)

  Rtx = PolynomialRing(ZZ,vcat([:t],symbols(Kx)))
  vvI = ideal([val_p.uniformizer-Rtx[2][1]])
  vvI = vvI+tropical_change_base_ring(Rtx[1],I)

  return vvI
end
export simulate_valuation



#=======
return true if f is homogeneous (w.r.t. total degree)
return false otherwise
=======#
function tropical_is_homogeneous(f::Union{AbstractAlgebra.Generic.MPoly{K},fmpq_mpoly,fmpz_mpoly} where {K})
  d = sum(exponent_vector(f,1))
  for i in 2:length(f)
    if d!=sum(exponent_vector(f,i))
      return false
    end
  end
  return true
end
export tropical_is_homogeneous

function tropical_is_homogeneous(I::MPolyIdeal{K} where {K})
  # todo: test whether generators are interreduced
  @warn "tropical_is_homogeneous: merely checking whether given generators are homogeneous, can result in false negative"

  for f in gens(I)
    if !tropical_is_homogeneous(f)
      return false
    end
  end
  return true
end



#=======
tropical Groebner basis
todo: proper documentation
Example:
val = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
w = [0,0,0]
I = ideal([x+2*y,y+2*z])
tropical_groebner_basis(I,val,w)
=======#
function tropical_groebner_basis(I,val,w)
  vvI = simulate_valuation(I,val)
  w = vcat([-1],w)

  Rtx = base_ring(vvI)
  # todo: replace with groebner_bases in OSCAR once more orderings are supported
  S, _ = Singular.PolynomialRing(singular_ring(base_ring(Rtx)), map(string, Nemo.symbols(Rtx)), ordering = Singular.ordering_a(w)*Singular.ordering_dp())
  SI = Singular.Ideal(S, [S(g) for g in gens(vvI)])

  vvGB = Singular.std(SI)
  return [Rtx(p) for p in Singular.gens(vvGB)]
end
export tropical_groebner_basis
