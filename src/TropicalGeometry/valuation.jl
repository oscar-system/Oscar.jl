###
# Valuations over exact fields for tropical geometry
# ==================================================
###

@doc Markdown.doc"""

    ValuationMap(K,p)

A valuation map for computations in the tropical geometry package. Currently, the only supported valuations are:
- t-adic valuation on QQ(t)
- p-adic valuations on QQ
- trivial valuation on any field

# Example for p-adic valuation on QQ
```jldoctest
julia> val_2 = ValuationMap(QQ,2);

julia> val_2(4)
2
julia> val_2(1//4)
-2
```

# Example for t-adic valuation on QQ(t)
```jldoctest
julia> Kt,t = RationalFunctionField(QQ,"t");

julia> val_t = ValuationMap(Kt,t);

julia> val_t(t^2)
2
julia> val_2(1//t^2)
-2
```

# Example for p-adic valuation on QQ
```jldoctest
julia> val = ValuationMap(QQ);

julia> val(4)
0
julia> val(1//4)
0
```
"""
struct ValuationMap{typeofValuedField,typeofUniformizer}
  valued_field::typeofValuedField
  uniformizer_field::typeofUniformizer
  valued_ring
  uniformizer_ring
  residue_field
  residue_map
  uniformizer_symbol
end
export ValuationMap


###
# trivial valuation
###

# Constructor:
function ValuationMap(K)
  residue_map(c) = return c
  return ValuationMap{typeof(K),Nothing}(K,nothing,K,nothing,K,residue_map,nothing)
end

# Evaluation:
(val::ValuationMap{K,Nothing} where {K})(c) = return 0



###
# p-adic valuation on QQ
###

# Constructor:
function ValuationMap(Q::FlintRationalField,p::fmpq)
  residue_map(c) = FiniteField(p)[1](c)
  return ValuationMap{typeof(Q),typeof(p)}(Q,p,ZZ,ZZ(p),FiniteField(ZZ(p))[1],residue_map,:p)
end

ValuationMap(Q::FlintRationalField,p) = ValuationMap(Q,QQ(p)) # for other types of `p` such as `Integer`

# Evaluation:
(val::ValuationMap{FlintRationalField,fmpq})(c) = valuation(QQ(c),val.uniformizer_ring)



###
# Laurent valuation on K(t)
###

# t-adic valuation for elements in the valued field (=rational functions):
function t_adic_valuation(c::AbstractAlgebra.Generic.Rat)
  num = numerator(c)
  nom = denominator(c)
  valnum = first(i for i in 0:degree(num) if !iszero(coeff(num, i)))
  valnom = first(i for i in 0:degree(nom) if !iszero(coeff(nom, i)))
  return valnum-valnom
end
# t-adic valuation for elements in the valued ring (=polynomials):
function t_adic_valuation(c::fmpq_poly)
  return first(i for i in 0:degree(c) if !iszero(coeff(c, i)))
end

# Constructor:
function ValuationMap(Kt::AbstractAlgebra.Generic.RationalFunctionField,t::AbstractAlgebra.Generic.Rat)
  function residue_map(c)
    valc = t_adic_valuation(c)
    if (valc<0)
      error("residue_map: input has negative valuation, not in valuation ring")
    end
    return base_ring(Kt)(evaluate(c,0))
  end
  Rt,_ = PolynomialRing(base_ring(Kt),symbols(Kt))
  return ValuationMap{typeof(Kt),typeof(t)}(Kt,t,Rt,Rt(t),base_ring(Kt),residue_map,:t)
end

# Evaluation:
(val::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K},AbstractAlgebra.Generic.Rat{K}} where {K})(c) = t_adic_valuation(c)




###
# Check whether valuation is p-adic (as in: p-adic numbers), t-adic (as in: function fields), or trivial
# ======================================================================================================
###

function is_valuation_p_adic(val::ValuationMap)
  return typeof(val.valued_field)==FlintRationalField && typeof(val.uniformizer_field)==fmpq
end

function is_valuation_t_adic(val::ValuationMap)
  return typeof(val.valued_field)==AbstractAlgebra.Generic.RationalFunctionField{fmpq} && typeof(val.uniformizer_field)==AbstractAlgebra.Generic.Rat{fmpq}
end

function is_valuation_trivial(val::ValuationMap)
  return typeof(val.uniformizer_field)==Nothing
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
K,s = RationalFunctionField(QQ,"s")
val_t = ValuationMap(K,s)
Kx,(x1,x2,x3) = PolynomialRing(K,3)
I = ideal([x1+s*x2,x2+s*x3])
simulate_valuation(I,val_t)

val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
simulate_valuation(I,val_2)
=======#
function simulate_valuation(I::MPolyIdeal, val::ValuationMap)

  # if the valuation is trivial, then nothing needs to be done
  if is_valuation_trivial(val)
    return I
  end

  R = val.valued_ring
  t = val.uniformizer_symbol

  Rtx,tx = PolynomialRing(R,vcat([t],symbols(base_ring(I))))
  vvG = [val.uniformizer_ring-tx[1]]
  for f in gens(I)
    fRtx = MPolyBuildCtx(Rtx)
    for (cK,expvKx) = zip(coefficients(f),exponent_vectors(f))
      @assert isone(denominator(cK)) "change_base_ring: coefficient denominators need to be 1"
      cR = R(numerator(cK))       # coefficient in R
      expvRtx = vcat([0],expvKx)  # exponent vector in R[t,x1,...,xn]
      push_term!(fRtx,cR,expvRtx)
    end
    push!(vvG,tighten_simulation(finish(fRtx),val))
  end

  vvI = ideal(Rtx,vvG)
  return vvI
end
export simulate_valuation

function simulate_valuation(w::Vector, val::ValuationMap)
  # if the valuation is non-trivial, prepend -1 to the vector
  if !is_valuation_trivial(val)
    w = vcat([-1],w)
  end
  # either way, scale vector to make entries integral
  commonDenom = lcm([denominator(wi) for wi in w])
  return [numerator(commonDenom*wi) for wi in w]
end

function simulate_valuation(w::Vector, u::Vector, val::ValuationMap)
  # if the valuation is non-trivial, prepend -1 to the vector
  if !is_valuation_trivial(val)
    w = vcat([-1],w)
    u = vcat([-1],u)
  end
  # either way, scale vector to make entries integral
  w_commonDenom = lcm([denominator(wi) for wi in w])
  u_commonDenom = lcm([denominator(ui) for ui in u])
  return [numerator(w_commonDenom*wi) for wi in w],[numerator(u_commonDenom*ui) for ui in u]
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

Ks,s = RationalFunctionField(QQ,"s")
val_s = ValuationMap(Ks,s)
Ksx,(x1,x2,x3) = PolynomialRing(Ks,3)
I = ideal([x1+s*x2,x2+s*x3])
vvI = simulate_valuation(I,val_s)
desimulate_valuation(vvI,val_s)
=======#
function desimulate_valuation(vvI::MPolyIdeal,val::ValuationMap)
  Rx = base_ring(vvI)
  R = coefficient_ring(Rx)
  x = copy(symbols(Rx))
  popfirst!(x)

  K = val.valued_field
  Kx,_ = PolynomialRing(K,x)

  vvG = [evaluate(g,[1],[val.uniformizer_ring]) for g in gens(vvI)]
  G = []
  for vvg in vvG
    if !iszero(vvg) # vvI contained p-t, so one entry of vvG is 0
      g = MPolyBuildCtx(Kx)
      for (c, expvRtx) = Base.Iterators.zip(coefficients(vvg), exponent_vectors(vvg))
        expvKx = copy(expvRtx) # exponent vector in R[t,x1,...,xn]
        popfirst!(expvKx)      # exponent vector in K[x1,...,xn]
        push_term!(g,K(c),expvKx)
      end
      push!(G,finish(g))
    end
  end

  return ideal(Kx,G)
end
export desimulate_valuation

function desimulate_valuation(w::Vector,val::ValuationMap)
  # if the valuation is non-trivial, scale the vector so that first entry is -1
  #   and then remove first entry
  if !is_valuation_trivial(val)
    w /= w[1]
    popfirst!(w)
  end
  return w
end


#=======
functions which reduces polynomials in variables t,x1, ..., xn simulating the valuation by p-t
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
function tighten_simulation(f::MPolyElem,val::ValuationMap)
  Rtx = parent(f)
  R = coefficient_ring(f)
  p = val.uniformizer_field
  f_tightened = MPolyBuildCtx(Rtx)
  for (c,alpha) in zip(coefficients(f),exponent_vectors(f))
    v = val(c)
    alpha[1] += v
    push_term!(f_tightened,R(c*p^-v),alpha)
  end
  return finish(f_tightened)
end
function tighten_simulation(I::MPolyIdeal,val::ValuationMap)
  return ideal([tighten_simulation(f) for f in gens(I)] + [])
end
