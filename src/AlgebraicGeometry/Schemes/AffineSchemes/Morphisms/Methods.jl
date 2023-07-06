


###########################################################
# (1) The fibre product of two morphisms of affine schemes
###########################################################

@doc raw"""
    fiber_product(f::SpecMor{SpecType, SpecType, <:Any}, g::SpecMor{SpecType, SpecType, <:Any}) where {SpecType<:StdSpec}

For morphisms ``f : Y → X`` and ``g : Z → X`` return the fiber
product ``Y×Z`` over ``X`` together with its two canonical projections.
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
  RX = ambient_coordinate_ring(X)
  W = subscheme(YxZ, [pullback(pY)(pullback(f)(x)) - pullback(pZ)(pullback(g)(x)) for x in gens(RX)])
  return W, restrict(pY, W, Y, check=false), restrict(pZ, W, Z, check=false)
end



###########################################################
# (2) The direct product of two affine schemes
###########################################################

# First the product of the ambient spaces. Documented below.
function product(X::AbsSpec{BRT, RT}, Y::AbsSpec{BRT, RT};
    change_var_names_to::Vector{String}=["", ""]
  ) where {BRT, RT<:MPolyRing}
  K = OO(X)
  L = OO(Y)
  # V = localized_ring(K)
  # W = localized_ring(L)
  k = base_ring(K)
  k == base_ring(L) || error("varieties are not defined over the same base ring")

  m = ngens(K)
  n = ngens(L)
  new_symb = Symbol[]
  if length(change_var_names_to[1]) == 0
    new_symb = symbols(K)
  else
    new_symb = Symbol.([change_var_names_to[1]*"$i" for i in 1:ngens(L)])
  end
  if length(change_var_names_to[2]) == 0
    new_symb = vcat(new_symb, symbols(L))
  else
    new_symb = vcat(new_symb, Symbol.([change_var_names_to[2]*"$i" for i in 1:ngens(L)]))
  end
  KL, z = polynomial_ring(k, new_symb)
  XxY = Spec(KL)
  pr1 = SpecMor(XxY, X, gens(KL)[1:m], check=false)
  pr2 = SpecMor(XxY, Y, gens(KL)[m+1:m+n], check=false)
  return XxY, pr1, pr2
end

@doc raw"""
    product(X::AbsSpec, Y::AbsSpec)
    
Return a triple ``(X×Y, p₁, p₂)`` consisting of the product ``X×Y`` over
the common base ring ``𝕜`` and the two projections ``p₁ : X×Y → X`` and
``p₂ : X×Y → Y``.
"""
function product(X::AbsSpec, Y::AbsSpec;
    change_var_names_to::Vector{String}=["", ""]
  )
  # take the product of the ambient spaces and restrict
  base_ring(X) == base_ring(Y) || error("schemes are not defined over the same base ring")
  A = ambient_space(X)
  B = ambient_space(Y)
  AxB,prA, prB  = product(A, B, change_var_names_to=change_var_names_to)
  XxY = intersect(preimage(prA, X, check=false), preimage(prB, Y,check=false))
  prX = restrict(prA, XxY, X, check=false)
  prY = restrict(prB, XxY, Y, check=false)
  return XxY, prX, prY
end



#=
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

  m = ngens(R)
  n = ngens(S)
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
  RS, z = polynomial_ring(k, new_symb)
  inc1 = hom(R, RS, gens(RS)[1:m], check=false)
  inc2 = hom(S, RS, gens(RS)[m+1:m+n], check=false)
  IX = ideal(RS, inc1.(gens(modulus(underlying_quotient(OO(X))))))
  IY = ideal(RS, inc2.(gens(modulus(underlying_quotient(OO(Y))))))
  UX = MPolyPowersOfElement(RS, inc1.(denominators(inverted_set(OO(X)))))
  UY = MPolyPowersOfElement(RS, inc2.(denominators(inverted_set(OO(Y)))))
  XxY = Spec(RS, IX + IY, UX*UY)
  pr1 = SpecMor(XxY, X, gens(RS)[1:m], check=false)
  pr2 = SpecMor(XxY, Y, gens(RS)[m+1:m+n], check=false)
  return XxY, pr1, pr2
end
=#




########################################
# (4) Equality
########################################

function ==(f::SpecMorType, g::SpecMorType) where {SpecMorType<:AbsSpecMor}
  X = domain(f)
  X == domain(g) || return false
  codomain(f) == codomain(g) || return false
  OO(X).(pullback(f).(gens(ambient_coordinate_ring(codomain(f))))) == OO(X).(pullback(f).(gens(ambient_coordinate_ring(codomain(g))))) || return false
  return true
end



########################################
# (5) Display
########################################

function Base.show(io::IO, f::AbsSpecMor)
  println(io, "morphism from\n")
  println(io, "\t$(domain(f))\n")
  println(io, "to\n")
  println(io, "\t$(codomain(f))\n")
  println(io, "with coordinates\n")
  x = coordinates(codomain(f))
  print(io,"\t")
  for i in 1:length(x)-1
    print(io, "$(pullback(f)(x[i])), ")
  end
  print(io, "$(pullback(f)(last(x)))")
end

########################################################################
# (6) Base change
########################################################################

@doc raw"""
    base_change(phi::Any, f::AbsSpecMor)
        domain_map::AbsSpecMor=base_change(phi, domain(f))[2],
        codomain_map::AbsSpecMor=base_change(phi, codomain(f))[2]
      )

For a morphism ``f : X → Y`` between two schemes over a `base_ring` ``𝕜`` 
and a ring homomorphism ``φ : 𝕜 → 𝕂`` this returns a triple 
`(b₁, F, b₂)` consisting of the maps in the commutative diagram 
```
             f
    X        →    Y
    ↑ b₁          ↑ b₂
  X×ₖSpec(𝕂) → Y×ₖSpec(𝕂)
             F

The optional arguments `domain_map` and `codomain_map` can be used 
to specify the morphisms `b₁` and `b₂`, respectively. 
```
"""
function base_change(phi::Any, f::AbsSpecMor; 
    domain_map::AbsSpecMor=base_change(phi, domain(f))[2],
    codomain_map::AbsSpecMor=base_change(phi, codomain(f))[2]
  )
  X = domain(f)
  Y = codomain(f)
  XX = domain(domain_map)
  YY = domain(codomain_map)
  pbf = pullback(f)
  pb1 = pullback(domain_map)
  pb2 = pullback(codomain_map)

  R = OO(Y)
  S = OO(X)
  RR = OO(YY)
  SS = OO(XX)

  img_gens = [pb1(pbf(x)) for x in gens(R)]
  # For the pullback of F no explicit coeff_map is necessary anymore 
  # since both rings in domain and codomain have the same (extended/reduced)
  # coefficient ring by now.
  pbF = hom(RR, SS, img_gens, check=false) # TODO: Set to false after testing

  return domain_map, SpecMor(XX, YY, pbF, check=false), codomain_map # TODO: Set to false after testing
end

function _register_birationality!(f::AbsSpecMor, 
    g::AbsSpecMor, ginv::AbsSpecMor)
  set_attribute!(g, :inverse, ginv)
  set_attribute!(ginv, :inverse, g)
  return _register_birationality(f, g)
end

function _register_birationality!(f::AbsSpecMor, 
    g::AbsSpecMor
  )
  set_attribute!(f, :is_birational, true)
  set_attribute!(f, :iso_on_open_subset, g)
end

