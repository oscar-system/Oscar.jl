export restrict, fiber_product, product, simplify



###########################################################
# (1) The fibre product of two morphisms of affine schemes
###########################################################

@Markdown.doc """
    fiber_product(f::SpecMor, g::SpecMor)

For morphisms ``f : Y → X`` and ``g : Z → X`` return the fiber
product ``Y ×ₓ Z`` over ``X`` together with its two canonical
projections.
"""
function fiber_product(
    f::SpecMor{SpecType, SpecType, <:Any},
    g::SpecMor{SpecType, SpecType, <:Any}
  ) where {SpecType<:StdSpec}
  Y = domain(f)
  X = codomain(f)
  X == codomain(g) || error("maps need to have the same codomain")
  Z = domain(g)
  YxZ, pY, pZ = product(Y, Z)
  RX = ambient_ring(X)
  #I = ideal(OO(YxZ), [pullback(pY)(pullback(f)(x)) - pullback(pZ)(pullback(g)(x)) for x in gens(RX)])
  W = subscheme(YxZ, [pullback(pY)(pullback(f)(x)) - pullback(pZ)(pullback(g)(x)) for x in gens(RX)])
  return W, restrict(pY, W, Y, check=false), restrict(pZ, W, Z, check=false)
end



###########################################################
# (2) The direct product of two affine schemes
###########################################################

@Markdown.doc """
    product(X::AbsSpec, Y::AbsSpec)
    
Returns a triple ``(X×Y, p₁, p₂)`` consisting of the product ``X×Y`` over 
the common base ring ``𝕜`` and the two projections ``p₁ : X×Y → X`` and 
``p₂ : X×Y → Y``.
"""
function product(X::AbsSpec, Y::AbsSpec;
    change_var_names_to::Vector{String}=["", ""]
  )
  base_ring(X) == base_ring(Y) || error("schemes are not defined over the same base ring")
  Xstd = standard_spec(X)
  Ystd = standard_spec(Y)
  XxY, prX, prY = product(Xstd, Ystd, change_var_names_to=change_var_names_to)
  return XxY, compose(prX, SpecMor(Xstd, X, gens(OO(Xstd)))), compose(prY, SpecMor(Ystd, Y, gens(OO(Ystd))))
end

function product(X::StdSpec, Y::StdSpec;
    change_var_names_to::Vector{String}=["", ""]
  )
  K = OO(X)
  L = OO(Y) 
  V = localized_ring(K)
  W = localized_ring(L)
  R = base_ring(K)
  S = base_ring(L)
  k = base_ring(R)
  k == base_ring(S) || error("varieties are not defined over the same field")

  m = length(gens(R))
  n = length(gens(S))
  new_symb = Symbol[]
  if length(change_var_names_to[1]) == 0
    new_symb = symbols(R)
  else 
    new_symb = Symbol.([change_var_names_to[1]*"$i" for i in 1:ngens(R)])
  end
  if length(change_var_names_to[2]) == 0
    new_symb = vcat(new_symb, symbols(S))
  else 
    new_symb = vcat(new_symb, Symbol.([change_var_names_to[2]*"$i" for i in 1:ngens(S)]))
  end
  RS, z = PolynomialRing(k, new_symb)
  inc1 = hom(R, RS, gens(RS)[1:m])
  inc2 = hom(S, RS, gens(RS)[m+1:m+n])
  IX = ideal(RS, inc1.(gens(modulus(quotient_ring(OO(X))))))
  IY = ideal(RS, inc2.(gens(modulus(quotient_ring(OO(Y))))))
  UX = MPolyPowersOfElement(RS, inc1.(denominators(inverted_set(OO(X)))))
  UY = MPolyPowersOfElement(RS, inc2.(denominators(inverted_set(OO(Y)))))
  XxY = Spec(RS, IX + IY, UX*UY)
  pr1 = SpecMor(XxY, X, gens(RS)[1:m], check=false)
  pr2 = SpecMor(XxY, Y, gens(RS)[m+1:m+n], check=false)
  return XxY, pr1, pr2
end



###########################################################
# (3) Simplify an affine scheme
###########################################################

@Markdown.doc """
    simplify(X::AbsSpec{<:Field})

Given an affine scheme ``X`` with coordinate ring ``R = 𝕜[x₁,…,xₙ]/I`` 
(or a localization thereof), use `Singular`'s `elimpart` to try 
to eliminate variables ``xᵢ`` to arrive at a simpler presentation 
``R ≅ R' = 𝕜[y₁,…,yₘ]/J`` for some ideal ``J``; return 
the triple ``(Y, f, g)`` where ``Y = Spec(R')`` and ``f : Y ↔ X : g``
are the identifying isomorphisms. 
"""
function simplify(X::AbsSpec{<:Field})
  L, f, g = simplify(OO(X))
  Y = Spec(L)
  YtoX = SpecMor(Y, X, f)
  XtoY = SpecMor(X, Y, g)
  set_attribute!(YtoX, :inverse, XtoY)
  set_attribute!(XtoY, :inverse, YtoX)
  return Y, YtoX, XtoY
end


########################################
# (4) Equality
########################################

function ==(f::SpecMorType, g::SpecMorType) where {SpecMorType<:AbsSpecMor}
  X = domain(f)
  X == domain(g) || return false
  codomain(f) == codomain(g) || return false
  OO(X).(pullback(f).(gens(ambient_ring(codomain(f))))) == OO(X).(pullback(f).(gens(ambient_ring(codomain(g))))) || return false
  return true
end



########################################
# (6) Display (missing?)
########################################
