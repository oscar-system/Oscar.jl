###
# Computing (tropical) Groebner bases in Oscar
# ============================================
#
# For a definition of tropical Groebner basis see Section 2.4 in:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
# To see how they can be computed using standard bases see:
#   T. Markwig, Y. Ren: Computing tropical varieties over fields with valuation
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
pseudo_change_base_ring(Rtx,f)
=======#
function pseudo_change_base_ring(Rtx::FmpqMPolyRing,f::AbstractAlgebra.Generic.MPoly{Kt} where {Kt})

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
function pseudo_change_base_ring(Rtx::FmpqMPolyRing,I::MPolyIdeal{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Rat{K}}} where {K})
  return ideal([pseudo_change_base_ring(Rtx,f) for f in gens(I)])
end



#=======
function which coerces a polynomial in QQ[x_1,...,x_n] into ZZ[t,x_1,...,x_n]
Example:
Kx,(x1,x2,x3) = PolynomialRing(QQ,3)
f = x1+(1+2)*x2+(1+3*2+5*2^2)*x3
Rtx,(t,y1,y2,y3) = PolynomialRing(ZZ,4)
pseudo_change_base_ring(Rtx,f)
=======#
function pseudo_change_base_ring(Rtx::FmpzMPolyRing,f::fmpq_mpoly)
  fRtx = MPolyBuildCtx(Rtx)
  for (cK, expvKx) = Base.Iterators.zip(Singular.coefficients(f), Singular.exponent_vectors(f))
    cR = numerator(cK)            # coefficient in R
    expvRtx = vcat([0],expvKx)    # exponent vector in R[t,x1,...,xn]
    @assert isone(denominator(cK)) "change_base_ring: coefficient denominators need to be 1"
    push_term!(fRtx,cR,expvRtx)
  end
  return finish(fRtx)
end
function pseudo_change_base_ring(Rtx::FmpzMPolyRing,I::MPolyIdeal{Ktx} where {Ktx})
  return ideal([pseudo_change_base_ring(Rtx,f) for f in gens(I)])
end
export pseudo_change_base_ring


#=======
functions which, given an ideal I in variables x1, ..., xn over a field with valuation,
returns an ideal vvI in variables t, x1, ..., xn such that tropical Groebner bases of I w.r.t. w
correspond to standard bases of I w.r.t. (-1,w)
Example:
Kt,t = RationalFunctionField(QQ,"t")
val_t = ValuationMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
I = ideal([x+t*y,y+t*z])
simulate_valuation(I,val_t)

val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
simulate_valuation(I,val_2)
=======#
function simulate_valuation(I,val_t::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K}, AbstractAlgebra.Generic.Rat{K}} where {K})

  Ktx = base_ring(I)
  Kt = base_ring(Ktx)
  K = base_ring(Kt)

  @assert K==QQ "simulate_valuation: only function fields over QQ supported for now"

  Rtx,_ = PolynomialRing(QQ,vcat(symbols(Kt),symbols(Ktx)));
  vvI = pseudo_change_base_ring(Rtx,I)

  return vvI
end
function simulate_valuation(I,val_p::ValuationMap{FlintRationalField, fmpz})

  Kx = base_ring(I)
  K = coefficient_ring(Kx)

  Rtx = PolynomialRing(ZZ,vcat([:t],symbols(Kx)))
  vvI = ideal([val_p.uniformizer-Rtx[2][1]])
  vvI = vvI+pseudo_change_base_ring(Rtx[1],I)

  return vvI
end
export simulate_valuation


#=======
functions which, given an ideal I in variables x1, ..., xn over a field with valuation,
returns an ideal vvI in variables t, x1, ..., xn such that tropical Groebner bases of I w.r.t. w
correspond to standard bases of I w.r.t. (-1,w)
Example:
val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
vvI = simulate_valuation(I,val_2)
desimulate_valuation(vvI,val_2)
=======#
function desimulate_valuation(vvI,val_p::ValuationMap{FlintRationalField, fmpz})
  Rtx = base_ring(vvI)
  x = copy(symbols(Rtx))
  popfirst!(x)
  K = val_p.valued_field
  Kx,_ = PolynomialRing(K,x)

  vvG = [evaluate(g,[1],[val_p.uniformizer]) for g in gens(vvI)]
  G = []
  for vvg in vvG
    if !iszero(vvg) # vvI contained p-t, so one entry of vvG is 0
      g = MPolyBuildCtx(Kx)
      for (c, expvRtx) = Base.Iterators.zip(Singular.coefficients(vvg), Singular.exponent_vectors(vvg))
        expvKx = copy(expvRtx) # exponent vector in R[t,x1,...,xn]
        popfirst!(expvKx)      # exponent vector in K[x1,...,xn]
        push_term!(g,c,expvKx)
      end
      push!(G,finish(g))
    end
  end

  return ideal(Kx,G)
end
#=======
Example:
Kt,t = RationalFunctionField(QQ,"t")
val_t = ValuationMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
I = ideal([x+t*y,y+t*z])
vvI = simulate_valuation(I,val_t)
desimulate_valuation(vvI,val_t)
=======#
function desimulate_valuation(vvI,val_t::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K}, AbstractAlgebra.Generic.Rat{K}} where {K})
  Rtx = base_ring(vvI)
  x = copy(symbols(Rtx))
  popfirst!(x)
  Kt = val_t.valued_field
  t = val_t.uniformizer
  Ktx,_ = PolynomialRing(Kt,x)

  G = []
  for vvg in gens(vvI)
    g = MPolyBuildCtx(Ktx)
    for (c, expvRtx) = Base.Iterators.zip(Singular.coefficients(vvg), Singular.exponent_vectors(vvg))
      expvKtx = copy(expvRtx) # exponent vector in R[t,x1,...,xn]
      d = popfirst!(expvKtx)  # exponent vector in K(t)[x1,...,xn]
      push_term!(g,Kt(c)*t^d,expvKtx)
    end
    push!(G,finish(g))
  end

  return ideal(Ktx,G)
end
export desimulate_valuation



#=======
return true if f is homogeneous (w.r.t. total degree)
return false otherwise
=======#
function is_homogeneous(f::Union{AbstractAlgebra.Generic.MPoly{K},fmpq_mpoly,fmpz_mpoly} where {K})
  d = sum(exponent_vector(f,1))
  for i in 2:length(f)
    if d!=sum(exponent_vector(f,i))
      return false
    end
  end
  return true
end
export is_homogeneous

function is_homogeneous(I::MPolyIdeal{K} where {K})
  # todo: test whether generators are interreduced
  @warn "is_homogeneous: merely checking whether given generators are homogeneous, can result in false negative"

  for f in gens(I)
    if !is_homogeneous(f)
      return false
    end
  end
  return true
end



#=======
tropical Groebner basis
todo: proper documentation
Example:

val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
w = [0,0,0]
I = ideal([x+2*y,y+2*z])
groebner_basis(I,val_2,w)

Kt,t = RationalFunctionField(QQ,"t")
val_t = ValuationMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
w = [0,0,0]
I = ideal([x+t*y,y+t*z])
groebner_basis(I,val_t,w)
=======#
function groebner_basis(I,val::ValuationMap{valuedField,uniformizer} where{valuedField,uniformizer},w::Vector{Int}; complete_reduction::Bool=false)
  vvI = simulate_valuation(I,val)
  w = vcat([-1],w)

  Rtx = base_ring(vvI)
  # todo: replace with groebner_bases in OSCAR once more orderings are supported
  S,_ = Singular.PolynomialRing(singular_ring(base_ring(Rtx)), map(string, Nemo.symbols(Rtx)), ordering = Singular.ordering_a(w)*Singular.ordering_dp())
  SI = Singular.Ideal(S, [S(g) for g in gens(vvI)])

  vvGB = Singular.std(SI,complete_reduction=complete_reduction)
  vvIGB = ideal(Rtx,Singular.gens(vvGB))

  IGB = desimulate_valuation(vvIGB,val)
  return gens(IGB)
end
