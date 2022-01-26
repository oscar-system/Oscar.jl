###
# Valuations over exact fields for tropical geometry
# ==================================================
###

struct ValuationMap{typeofValuedField,typeofUniformizer}
  valued_field::typeofValuedField
  uniformizer::typeofUniformizer
  residue_field
  residue_map
end
export ValuationMap


###
# trivial valuation
###

# Constructor:
function ValuationMap()
  residue_map(c) = return c
  return ValuationMap{typeof(K),Nothing}(K,nothing,K,residue_map)
end

# Evaluation:
(val::ValuationMap{K,Nothing} where {K})(c) = return 0



###
# p-adic valuation on QQ
###

# Constructor:
function ValuationMap(Q::FlintRationalField,p::fmpz)
  residue_map(c) = FiniteField(p)(c)
  return ValuationMap{typeof(Q),typeof(p)}(Q,p,FiniteField(p),residue_map)
end

ValuationMap(Q::FlintRationalField,p) = ValuationMap(Q,ZZ(p)) # for other types of `p` such as `Int`

# Evaluation:
(val::ValuationMap{FlintRationalField,fmpz})(c) = valuation(c, val.uniformizer)



###
# Laurent valuation on K(t)
###

# Constructor:
function t_adic_valuation(c,t::AbstractAlgebra.Generic.Rat)
    num = numerator(c)
    nom = denominator(c)
    valnum = first(i for i in 0:degree(num) if !iszero(coeff(num, i)))
    valnom = first(i for i in 0:degree(nom) if !iszero(coeff(nom, i)))
    return valnum-valnom
end

function ValuationMap(Kt::AbstractAlgebra.Generic.RationalFunctionField,t::AbstractAlgebra.Generic.Rat)
    function residue_map(c)
        valc = t_adic_valuation(c,t)
        if (valc<0)
            error("residue_map: input has negative valuation, not in valuation ring")
        end
        return base_ring(Kt)(evaluate(c,0))
    end
    return ValuationMap{typeof(Kt),typeof(t)}(Kt,t,base_ring(Kt),residue_map)
end

# Evaluation:
(val::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K},AbstractAlgebra.Generic.Rat{K}} where {K})(c) = t_adic_valuation(val.uniformizer,c)




###
# Check whether valuation is p-adic (as in: p-adic numbers), t-adic (as in: function fields), or trivial
# ======================================================================================================
###

function is_valuation_p_adic(val::ValuationMap)
  return typeof(val.valued_field)==FlintRationalField && typeof(val.uniformizer)==fmpz
end

function is_valuation_t_adic(val::ValuationMap)
  return typeof(val.valued_field)==AbstractAlgebra.Generic.RationalFunctionField{fmpq} && typeof(val.uniformizer)==AbstractAlgebra.Generic.Rat{fmpq}
end

function is_valuation_trivial(val::ValuationMap)
  return typeof(val.uniformizer)==Nothing
end



###
# Simulating valuations for algebraic computations in tropical geometry
# =====================================================================
###



###
# temporary workarounds:
###
function symbols(Kt::AbstractAlgebra.Generic.RationalFunctionField{K} where {K})
  return Kt.S
end



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
function simulate_valuation(I::MPolyIdeal,val_t::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K}, AbstractAlgebra.Generic.Rat{K}} where {K})

  Ktx = base_ring(I)
  Kt = base_ring(Ktx)
  K = base_ring(Kt)
  @assert K==QQ "simulate_valuation: only function fields over QQ supported for now"

  Rtx,_ = PolynomialRing(QQ,vcat(symbols(Kt),symbols(Ktx)));
  R = base_ring(Rtx)
  vvG = []
  for f in gens(I)
    fRtx = MPolyBuildCtx(Rtx)
    for (cKt,expvKtx) in zip(coefficients(f),exponent_vectors(f))
      # cKt = coefficient in K(t)
      # expvKtx = exponent vector in K(t)[x1,...,xn]
      expvRtx = vcat([0],expvKtx)  # exponent vector in R[t,x1,...,xn]
      @assert isone(denominator(cKt)) "change_base_ring: coefficient denominators need to be 1"
      cKt = numerator(cKt)           # coefficient in K[t]
      cK = coefficients(cKt)         # vector in K
      M = lcm([denominator(c) for c in cK])
      cR = [R(M*c) for c in cK]      # vector in R

      for c in cR
        push_term!(fRtx,c,expvRtx)
        expvRtx[1] += 1
      end
    end
    push!(vvG,finish(fRtx))
  end

  vvI = ideal(Rtx,vvG)
  return vvI
end
function simulate_valuation(I::MPolyIdeal,val_p::ValuationMap{FlintRationalField, fmpz})

  Kx = base_ring(I)
  K = coefficient_ring(Kx)

  Rtx,tx = PolynomialRing(ZZ,vcat([:t],symbols(Kx)))
  vvG = [val_p.uniformizer-tx[1]]
  for f in gens(I)
    fRtx = MPolyBuildCtx(Rtx)
    for (cK,expvKx) = zip(coefficients(f),exponent_vectors(f))
      # cK = coefficient in K
      # expvKx = exponent vector in K[x1,...,xn]
      @assert isone(denominator(cK)) "change_base_ring: coefficient denominators need to be 1"
      cR = numerator(cK)          # coefficient in R
      expvRtx = vcat([0],expvKx)  # exponent vector in R[t,x1,...,xn]
      push_term!(fRtx,cR,expvRtx)
    end
    push!(vvG,finish(fRtx))
  end

  vvI = ideal(Rtx,vvG)
  return vvI
end
export simulate_valuation

function simulate_valuation(w::Vector,val::ValuationMap)
  # if the valuation is non-trivial, prepend -1 to the vector
  if !is_valuation_trivial(val)
    w = vcat([-1],w)
  end
  # either way, scale vector to make entries integral
  commonDenom = lcm([denominator(wi) for wi in w])
  return [numerator(commonDenom*wi) for wi in w]
end

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
