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

# Display:
function Base.show(io::IO, val::ValuationMap{K,Nothing} where {K})
    print(io, "The trivial valuation on $(val.valued_field)")
end


###
# p-adic valuation on QQ
###

# Constructor:
function ValuationMap(Q::FlintRationalField,p::fmpq)
  function residue_map(c)
    return FiniteField(ZZ(p))[1](ZZ(c))
  end
  return ValuationMap{typeof(Q),typeof(p)}(Q,p,ZZ,ZZ(p),FiniteField(ZZ(p))[1],residue_map,:p)
end

ValuationMap(Q::FlintRationalField,p) = ValuationMap(Q,QQ(p)) # for other types of `p` such as `Integer`

# Evaluation:
(val::ValuationMap{FlintRationalField,fmpq})(c) = valuation(QQ(c),val.uniformizer_ring)

# Display:
function Base.show(io::IO, val::ValuationMap{FlintRationalField,fmpq})
    print(io, "The $(val.uniformizer_field)-adic valuation on $(val.valued_field)")
end


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

# Display:
function Base.show(io::IO, val::ValuationMap{AbstractAlgebra.Generic.RationalFunctionField{K},AbstractAlgebra.Generic.Rat{K}} where {K})
    print(io, "The $(val.uniformizer_field)-adic valuation on $(val.valued_field)")
end



###
# Check whether valuation is p-adic (as in: p-adic numbers), t-adic (as in: function fields), or trivial
# ======================================================================================================
###

function is_valuation_p_adic(val::ValuationMap)
  return val.valued_field isa FlintRationalField && val.uniformizer_field isa fmpq
end

function is_valuation_t_adic(val::ValuationMap)
  return val.valued_field isa AbstractAlgebra.Generic.RationalFunctionField && val.uniformizer_field isa AbstractAlgebra.Generic.Rat
end

function is_valuation_trivial(val::ValuationMap)
  return typeof(val.uniformizer_field)==Nothing
end

function is_valuation_nontrivial(val::ValuationMap)
  return typeof(val.uniformizer_field)!=Nothing
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
function simulate_valuation(I::MPolyIdeal, val::ValuationMap; coefficient_field::Bool=false)
  return ideal(simulate_valuation(gens(I),val,coefficient_field=coefficient_field))
end
function simulate_valuation(G::Vector{<:MPolyElem}, val::ValuationMap; coefficient_field::Bool=false)

  # if the valuation is trivial, then nothing needs to be done
  if is_valuation_trivial(val)
    return G
  end
  if length(G)==0
    error("input vector of polynomials empty, thus ambient polynomial ring unknown")
  end

  if coefficient_field
    R = val.valued_field
  else
    R = val.valued_ring
  end
  t = val.uniformizer_symbol

  Rtx,tx = PolynomialRing(R,vcat([t],symbols(parent(G[1]))))
  vvG = [val.uniformizer_ring-tx[1]]
  for f in G
    fRtx = MPolyBuildCtx(Rtx)
    for (cK,expvKx) = zip(coefficients(f),exponent_vectors(f))
      @assert isone(denominator(cK)) "change_base_ring: coefficient denominators need to be 1"
      cR = R(numerator(cK))       # coefficient in R
      expvRtx = vcat([0],expvKx)  # exponent vector in R[t,x1,...,xn]
      push_term!(fRtx,cR,expvRtx)
    end
    push!(vvG,tighten_simulation(finish(fRtx),val))
  end

  return vvG
end
export simulate_valuation

function simulate_valuation(w::Vector, val::ValuationMap)
  # if the valuation is non-trivial, prepend -1 to the vector
  if !is_valuation_trivial(val)
    w = vcat([-1],w)
  end
  # either way, scale vector to make entries integral
  commonDenom = lcm([denominator(wi) for wi in w])
  return [Int(numerator(commonDenom*wi)) for wi in w] # casting vector entries to Int32 for Singular
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
  return [Int(numerator(w_commonDenom*wi)) for wi in w],[Int(numerator(u_commonDenom*ui)) for ui in u]
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
  return ideal(desimulate_valuation(gens(vvI),val))
end
export desimulate_valuation

function desimulate_valuation(vvG::Vector{<:MPolyElem}, val::ValuationMap)
  G = [desimulate_valuation(vvg,val) for vvg in vvG]
  return [g for g in G if !iszero(g)]
end

function desimulate_valuation(vvg::MPolyElem, val::ValuationMap)
  Rx = parent(vvg)
  R = coefficient_ring(Rx)
  x = copy(symbols(Rx))
  popfirst!(x)

  K = val.valued_field
  Kx,_ = PolynomialRing(K,x)


  vvg = evaluate(vvg,[1],[val.uniformizer_ring])
  if iszero(vvg) # vvg may be p-t
    return Kx(0)
  end

  g = MPolyBuildCtx(Kx)
  for (c, expvRtx) = Base.Iterators.zip(coefficients(vvg), exponent_vectors(vvg))
    expvKx = copy(expvRtx) # exponent vector in R[t,x1,...,xn]
    popfirst!(expvKx)      # exponent vector in K[x1,...,xn]
    push_term!(g,K(c),expvKx)
  end
  return finish(g)
end

function desimulate_valuation(w::Vector,val::ValuationMap)
  # if the valuation is non-trivial, scale the vector so that first entry is -1
  #   and then remove first entry
  if !is_valuation_trivial(val)
    w /= w[1]
    popfirst!(w)
  end
  return w
end

function desimulate_valuation(w::Vector, u::Vector, val::ValuationMap)
  # if the valuation is non-trivial, scale the vector w so that first entry is -1
  #   and then remove first entry of both w and u
  if !is_valuation_trivial(val)
    w /= w[1]
    popfirst!(w)
    popfirst!(u)
  end
  return w,u
end


#=======
Given a polynomial f in t,x_1, ..., x_n simulating the valuation:
- If f==p-t or f==t-p, returns f
- Otherwise, returns a polynomial f' in <f,p-t> such that
  * f' and f share the same monomials in x
  * all terms of f' have distinct monomials in x
  * all terms of f' have valuation 0 coefficients
Example:
val_2 = ValuationMap(QQ,2)
Rtx,(p,x1,x2,x3) = PolynomialRing(val_2.valued_ring,["p","x1","x2","x3"])
f = x1+p*x1+p^2*x1+2^2*x2+p*x2+p^2*x2+x3
tighten_simulation(f,val_2)
tighten_simulation(2^3*f,val_2)
tighten_simulation(p^3*f,val_2)

K,s = RationalFunctionField(QQ,"s")
val_s = ValuationMap(K,s)
s = val_s.uniformizer_ring
Rtx,(t,x1,x2,x3) = PolynomialRing(val_s.valued_ring,["t","x1","x2","x3"])
f = x1+t*x1+t^2*x1+s^2*x2+t*x2+t^2*x2+x3
tighten_simulation(f,val_s)
tighten_simulation(s^3*f,val_s)
tighten_simulation(t^3*f,val_s)
=======#
function tighten_simulation(f::MPolyElem,val::ValuationMap)

  # return f if f = p-t or t-p
  Rtx = parent(f)
  p = val.uniformizer_ring
  pt = p - gens(Rtx)[1]
  if f==pt || f==-pt
    return f
  end

  # subsitute first variable by uniformizer_ring so that all monomials have distinct x-monomials
  # and compute the gcd of its coefficients
  # note: if the coefficient ring is a field, gcd will always be 1
  f = evaluate(f,[1],[p]) # todo: sanity check that f is not 0
  cGcd = val.valued_field(gcd(collect(coefficients(f))))

  # next divide f by the gcd of its coefficients
  # and replace uniformizer_ring by first variable
  K = val.valued_field
  R = coefficient_ring(Rtx)
  p = val.uniformizer_field
  f_tightened = MPolyBuildCtx(Rtx)
  for (c,alpha) in zip(coefficients(f),exponent_vectors(f))
    c = K(c)//cGcd # casting c into K for the t-adic valuation case where typeof(c)=fmpq_poly
    v = val(c)
    alpha[1] += v
    push_term!(f_tightened,R(c//p^v),alpha)
  end

  return finish(f_tightened)

end
export tighten_simulation
