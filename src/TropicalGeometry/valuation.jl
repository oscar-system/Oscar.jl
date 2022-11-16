###
# Valuations over exact fields for tropical geometry
# ==================================================
###

export TropicalSemiringMap,
       simulate_valuation,
       desimulate_valuation,
       tighten_simulation

@doc Markdown.doc"""

    TropicalSemiringMap(K,p,M::Union{typeof(min),typeof(max)}=min)

Constructs a map `val` from `K` to the min tropical semiring `T` (default)
or the max tropical semiring that:
- is a semigroup homomorphism `(K,*) -> (T,+)`,
- preserves the ordering on both sides.

In other words, `val` is either a valuation on `K` with image in
`TropicalSemiring(min)` or the negative of a valuation on `K` with image in
`TropicalSemiring(max)`.

The role of `val` is to encode with respect to which valuation on `K` and
under which convention (min or max) tropical computations should take place.

Currently, the only supported valuations are:
- the $t$-adic valuation on $\mathbb{Q}(t)$
- the $p$-adic valuations on $\mathbb{Q}$
- the trivial valuation on any field

# Example ($p$-adic valuation on $\mathbb{Q}$)
```jldoctest
julia> val_2 = TropicalSemiringMap(QQ,2); # = TropicalSemiringMap(QQ,2,min)

julia> val_2(4)
(2)
julia> val_2(1//4)
(-2)
julia> val_2 = TropicalSemiringMap(QQ,2,max);

julia> val_2(4)
(-2)
julia> val_2(1//4)
(2)
```

# Example ($t$-adic valuation on $\mathbb{Q}(t)$)
```jldoctest
julia> Kt,t = RationalFunctionField(QQ,"t");

julia> val_t = TropicalSemiringMap(Kt,t);

julia> val_t(t^2)
(2)
julia> val_t(1//t^2)
(-2)
```

# Example (trivial valuation on $\mathbb{Q}$)
```jldoctest
julia> val = TropicalSemiringMap(QQ);

julia> val(4)
(0)
julia> val(1//4)
(0)
julia> val(0)
∞
```
"""
struct TropicalSemiringMap{typeofValuedField,typeofUniformizer}
  valued_field::typeofValuedField
  uniformizer_field::typeofUniformizer
  valued_ring
  uniformizer_ring
  residue_field
  residue_map
  uniformizer_symbol
  TropicalSemiring
end


################################################################################
#
#  Basic access
#
################################################################################

valued_field(val::TropicalSemiringMap) = val.valued_field
uniformizer_field(val::TropicalSemiringMap) = val.uniformizer_field
valued_ring(val::TropicalSemiringMap) = val.valued_ring
uniformizer_ring(val::TropicalSemiringMap) = val.uniformizer_ring
residue_field(val::TropicalSemiringMap) = val.residue_field
residue_map(val::TropicalSemiringMap) = val.residue_map
uniformizer_symbol(val::TropicalSemiringMap) = val.uniformizer_symbol
TropicalSemiring(val::TropicalSemiringMap) = val.TropicalSemiring
convention(val::TropicalSemiringMap) = convention(val.TropicalSemiring)

###
# trivial valuation
###

# Constructor:
function TropicalSemiringMap(K,M::Union{typeof(min),typeof(max)}=min)
  residue_map(c) = return c
  return TropicalSemiringMap{typeof(K),Nothing}(K,nothing,K,nothing,K,residue_map,nothing,TropicalSemiring(M))
end

# Evaluation:
function (val::TropicalSemiringMap{K,Nothing} where {K})(c)
  if iszero(c)
    return inf(val.TropicalSemiring)
  end
  return val.TropicalSemiring(0)
end

# Display:
function Base.show(io::IO, val::TropicalSemiringMap{K,Nothing} where {K})
    print(io, "The trivial valuation on $(val.valued_field)")
end


###
# p-adic valuation on QQ
###

# Constructor:
function TropicalSemiringMap(Q::FlintRationalField, p::fmpq, M::Union{typeof(min),typeof(max)}=min)
  function residue_map(c)
    return FiniteField(ZZ(p))[1](ZZ(c))
  end
  return TropicalSemiringMap{typeof(Q),typeof(p)}(Q,p,ZZ,ZZ(p),FiniteField(ZZ(p))[1],residue_map,:p,TropicalSemiring(M))
end
# for other types of `p` such as `Integer`
TropicalSemiringMap(Q::FlintRationalField,p::fmpz,M::Union{typeof(min),typeof(max)}=min) = TropicalSemiringMap(Q,QQ(p),M)
TropicalSemiringMap(Q::FlintRationalField,p::Int64,M::Union{typeof(min),typeof(max)}=min) = TropicalSemiringMap(Q,QQ(p),M)

# Evaluation:
function (val::TropicalSemiringMap{FlintRationalField,fmpq})(c)
  if iszero(c)
    return inf(val.TropicalSemiring)
  end
  if convention(val)==min
    return val.TropicalSemiring(valuation(QQ(c),val.uniformizer_ring))
  end
  return val.TropicalSemiring(-valuation(QQ(c),val.uniformizer_ring))
end

# Display:
function Base.show(io::IO, val::TropicalSemiringMap{FlintRationalField,fmpq})
    print(io, "The $(val.uniformizer_field)-adic valuation on $(val.valued_field)")
end


###
# Laurent valuation on K(t)
###

# t-adic valuation for elements in the valued field (=rational functions):
function t_adic_valuation(c::Generic.Rat)
  num = numerator(c)
  nom = denominator(c)
  return t_adic_valuation(num)-t_adic_valuation(nom)
end
# t-adic valuation for elements in the valued ring (=polynomials):
function t_adic_valuation(c::PolyElem)
  return first(i for i in 0:degree(c) if !iszero(coeff(c, i)))
end

# Constructor:
function TropicalSemiringMap(Kt::AbstractAlgebra.Generic.RationalFunctionField,t::AbstractAlgebra.Generic.Rat,M::Union{typeof(min),typeof(max)}=min)
  function residue_map(c)
    valc = t_adic_valuation(c)
    if (valc<0)
      error("residue_map: input has negative valuation, not in valuation ring")
    end
    return base_ring(Kt)(evaluate(c,0))
  end
  Rt,_ = PolynomialRing(base_ring(Kt),symbols(Kt))
  return TropicalSemiringMap{typeof(Kt),typeof(t)}(Kt,t,Rt,Rt(t),base_ring(Kt),residue_map,:t,TropicalSemiring(M))
end

# Evaluation:

function (val::TropicalSemiringMap{S, T})(c) where {K,
                        S <: AbstractAlgebra.Generic.RationalFunctionField{K},
                        T <: AbstractAlgebra.Generic.Rat{K}}
  if iszero(c)
    return inf(val.TropicalSemiring)
  end
  if convention(val)==min
    return val.TropicalSemiring(t_adic_valuation(c))
  end
  return val.TropicalSemiring(-t_adic_valuation(c))
end
# Display:
function Base.show(io::IO, val::TropicalSemiringMap{S, T}) where {K,
                        S <: AbstractAlgebra.Generic.RationalFunctionField{K},
                        T <: AbstractAlgebra.Generic.Rat{K}}
    print(io, "The $(val.uniformizer_field)-adic valuation on $(val.valued_field)")
end



###
# Check whether valuation is p-adic (as in: p-adic numbers), t-adic (as in: function fields), or trivial
# ======================================================================================================
###

function is_valuation_p_adic(val::TropicalSemiringMap)
  return val.valued_field isa FlintRationalField && val.uniformizer_field isa fmpq
end

function is_valuation_t_adic(val::TropicalSemiringMap)
  return val.valued_field isa AbstractAlgebra.Generic.RationalFunctionField && val.uniformizer_field isa AbstractAlgebra.Generic.Rat
end

function is_valuation_trivial(val::TropicalSemiringMap)
  return typeof(val.uniformizer_field)==Nothing
end

function is_valuation_nontrivial(val::TropicalSemiringMap)
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
val_t = TropicalSemiringMap(K,s)
Kx,(x1,x2,x3) = PolynomialRing(K,3)
I = ideal([x1+s*x2,x2+s*x3])
simulate_valuation(I,val_t)

val_2 = TropicalSemiringMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
simulate_valuation(I,val_2)
=======#
function simulate_valuation(I::MPolyIdeal, val::TropicalSemiringMap)
  return ideal(simulate_valuation(gens(I),val))
end
function simulate_valuation(G::Vector{<:MPolyElem}, val::TropicalSemiringMap)

  # if the valuation is trivial, then nothing needs to be done
  if is_valuation_trivial(val)
    return G
  end
  if length(G)==0
    error("input vector of polynomials empty, thus ambient polynomial ring unknown")
  end

  R = val.valued_ring
  t = val.uniformizer_symbol

  Rtx,tx = PolynomialRing(R,vcat([t],symbols(parent(G[1]))))
  vvG = [val.uniformizer_ring-tx[1]]
  for f in G
    fRtx = MPolyBuildCtx(Rtx)
    for (cK,expvKx) = zip(AbstractAlgebra.coefficients(f),AbstractAlgebra.exponent_vectors(f))
      @assert isone(denominator(cK)) "change_base_ring: coefficient denominators need to be 1"
      cR = R(numerator(cK))       # coefficient in R
      expvRtx = vcat([0],expvKx)  # exponent vector in R[t,x1,...,xn]
      push_term!(fRtx,cR,expvRtx)
    end
    push!(vvG,tighten_simulation(finish(fRtx),val))
  end

  return vvG
end

function simulate_valuation(w::Vector, val::TropicalSemiringMap)
  # if the valuation is non-trivial, prepend -1 to the vector
  if !is_valuation_trivial(val)
    w = vcat([-1],w)
  end
  # either way, scale vector to make entries integral
  commonDenom = lcm([denominator(wi) for wi in w])
  sw = [Int(numerator(commonDenom*wi)) for wi in w] # casting vector entries to Int32 for Singular
  if convention(val)==min
    sw *= -1
  end
  return sw
end

function simulate_valuation(w::Vector, u::Vector, val::TropicalSemiringMap)
  # if the valuation is non-trivial, prepend -1 to the vector
  if !is_valuation_trivial(val)
    w = vcat([-1],w)
    u = vcat([0],u)
  end
  # either way, scale vector to make entries integral
  w_commonDenom = lcm([denominator(wi) for wi in w])
  u_commonDenom = lcm([denominator(ui) for ui in u])
  sw = [Int(numerator(w_commonDenom*wi)) for wi in w]
  su = [Int(numerator(u_commonDenom*ui)) for ui in u]
  if convention(val)==min
    sw *= -1
    su *= -1
  end
  return sw,su
end

#=======
functions which, given an ideal I in variables x1, ..., xn over a field with valuation,
returns an ideal vvI in variables t, x1, ..., xn such that tropical Groebner bases of I w.r.t. w
correspond to standard bases of I w.r.t. (-1,w)
Example:
val_2 = TropicalSemiringMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
I = ideal([x+2*y,y+2*z])
vvI = simulate_valuation(I,val_2)
desimulate_valuation(vvI,val_2)

Ks,s = RationalFunctionField(QQ,"s")
val_s = TropicalSemiringMap(Ks,s)
Ksx,(x1,x2,x3) = PolynomialRing(Ks,3)
I = ideal([x1+s*x2,x2+s*x3])
vvI = simulate_valuation(I,val_s)
desimulate_valuation(vvI,val_s)
=======#
function desimulate_valuation(vvI::MPolyIdeal,val::TropicalSemiringMap)
  return ideal([g for g in desimulate_valuation(gens(vvI),val) if !iszero(g)])
end

function desimulate_valuation(vvG::Vector{<:MPolyElem}, val::TropicalSemiringMap)
  return [desimulate_valuation(vvg,val) for vvg in vvG]
end

function desimulate_valuation(vvg::MPolyElem, val::TropicalSemiringMap)
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
  for (c, expvRtx) = Base.Iterators.zip(AbstractAlgebra.coefficients(vvg), AbstractAlgebra.exponent_vectors(vvg))
    expvKx = copy(expvRtx) # exponent vector in R[t,x1,...,xn]
    popfirst!(expvKx)      # exponent vector in K[x1,...,xn]
    push_term!(g,K(c),expvKx)
  end
  return finish(g)
end

function desimulate_valuation(w::Vector,val::TropicalSemiringMap)
  # if the valuation is non-trivial, scale the vector so that first entry is -1
  #   and then remove first entry
  if !is_valuation_trivial(val)
    w /= w[1]
    popfirst!(w)
  end
  if convention(val)==min
    w *= -1
  end
  return w
end

function desimulate_valuation(w::Vector, u::Vector, val::TropicalSemiringMap)
  # if the valuation is non-trivial, scale the vector w so that first entry is -1
  #   and then remove first entry of both w and u
  if !is_valuation_trivial(val)
    w /= w[1]
    popfirst!(w)
    popfirst!(u)
  end
  if convention(val)==min
    w *= -1
    u *= -1
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
val_2 = TropicalSemiringMap(QQ,2)
Rtx,(p,x1,x2,x3) = PolynomialRing(val_2.valued_ring,["p","x1","x2","x3"])
f = x1+p*x1+p^2*x1+2^2*x2+p*x2+p^2*x2+x3
tighten_simulation(f,val_2)
tighten_simulation(2^3*f,val_2)
tighten_simulation(p^3*f,val_2)

K,s = RationalFunctionField(QQ,"s")
val_s = TropicalSemiringMap(K,s)
s = val_s.uniformizer_ring
Rtx,(t,x1,x2,x3) = PolynomialRing(val_s.valued_ring,["t","x1","x2","x3"])
f = x1+t*x1+t^2*x1+s^2*x2+t*x2+t^2*x2+x3
tighten_simulation(f,val_s)
tighten_simulation(s^3*f,val_s)
tighten_simulation(t^3*f,val_s)
=======#
function tighten_simulation(f::MPolyElem,val::TropicalSemiringMap)
  @assert !iszero(f)
  # return f if f = p-t or t-p
  Rtx = parent(f)
  p = val.uniformizer_ring
  pt = p - gens(Rtx)[1]
  if f==pt || f==-pt
    return f
  end

  # substitute first variable by uniformizer_ring so that all monomials have distinct x-monomials
  # and compute the gcd of its coefficients
  f = evaluate(f,[1],[p]) # todo: sanity check that f is not 0
  cGcd = val.valued_field(gcd([c for c in AbstractAlgebra.coefficients(f)]))

  # next divide f by the gcd of its coefficients
  # and replace uniformizer_ring by first variable
  K = val.valued_field
  R = val.valued_ring
  p = val.uniformizer_field
  f_tightened = MPolyBuildCtx(Rtx)
  for (c,alpha) in zip(AbstractAlgebra.coefficients(f),AbstractAlgebra.exponent_vectors(f))
    c = K(c)//cGcd # casting c into K for the t-adic valuation case where typeof(c)=fmpq_poly
    v = Int(val(c); preserve_ordering=true)
    alpha[1] += v
    push_term!(f_tightened,R(c//p^v),alpha)
  end

  return finish(f_tightened)

end


# function valuation_Int(val::TropicalSemiringMap, c)
#   assert !iszero(c)
#   vc = Int(ZZ(data(val(c))))
#   if convention(valuation)==min
#     return vc
#   else
#     return -vc
#   end
# end
