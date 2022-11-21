export pullback, is_non_zero_divisor, find_non_zero_divisor, generic_fractions
########################################################################
# Printing                                                             #
########################################################################
function Base.show(io::IO, f::SpecOpenMor) 
  print(io, "Morphism from $(domain(f)) to $(codomain(f))")
  #given by the rational map $(generic_fractions(f))")
end

########################################################################
# Composition of maps                                                  #
########################################################################
function compose(f::SpecOpenMor, g::SpecOpenMor; check::Bool=true)
  U = domain(f)
  Cf = codomain(f)
  V = domain(g)
  if check
    issubset(Cf, V) || error("maps are not compatible")
  end
  W = codomain(g)
  X = ambient_scheme(U)
  Y = ambient_scheme(V)
  Z = ambient_scheme(W)
  pb_coords = [pullback(f)(pullback(g)(OO(W)(x))) for x in gens(ambient_coordinate_ring(Z))]
  maps_on_patches = [SpecMor(A, Z, [restrict(h, A) for h in pb_coords], check=check) for A in affine_patches(U)]
  return SpecOpenMor(U, W, maps_on_patches)
end

########################################################################
# Pullback of regular functions (deprecated)                           #
########################################################################
function pullback(f::SpecOpenMor, a::RingElem)
  U = domain(f)
  X = ambient_scheme(U)
  V = codomain(f)
  Y = ambient_scheme(V)
  R = ambient_coordinate_ring(Y)
  parent(a) == R || error("element does not belong to the correct ring")
  pb_a = [pullback(f[i])(a) for i in 1:npatches(U)]
  return SpecOpenRingElem(SpecOpenRing(X, U), pb_a)
end

########################################################################
# Equality test                                                        #
########################################################################
function ==(f::T, g::T) where {T<:SpecOpenMor} 
  (domain(f) == domain(g)) || return false
  (codomain(f) == codomain(g)) || return false
  Y = ambient_scheme(codomain(f))
  m = length(affine_patches(domain(f)))
  n = length(affine_patches(domain(g)))
  for i in 1:m
    for j in 1:n
      U = intersect(domain(f)[i], domain(g)[i])
      restrict(f[i], U, Y) == restrict(g[i], U, Y) || return false
    end
  end
  return true
end

########################################################################
# Preimages under SpecOpenMors                                         #
########################################################################
function preimage(f::SpecOpenMor, Z::AbsSpec; check::Bool=true)
  U = domain(f) 
  X = ambient_scheme(U)
  if check
    is_closed_embedding(Z, ambient_scheme(codomain(f))) || error("second argument must be closed in the codomain")
  end
  n = length(affine_patches(U))
  pbZ = [preimage(f[i], Z) for i in 1:n]
  Y = X 
  for K in pbZ
    Y = subscheme(Y, gens(modulus(underlying_quotient(OO(K)))))
  end
  return SpecOpen(Y, [g for g in gens(U) if !iszero(OO(Y)(g))])
end
function preimage(f::SpecOpenMor, W::PrincipalOpenSubset; check::Bool=true)
  V = codomain(f) 
  Y = ambient_scheme(V)
  Y === ambient_scheme(W) || error("second argument must be open in the ambient scheme of the domain of the morphism")
  h = complement_equation(W)
  pbh = pullback(f)(OO(codomain(f))(h))
  R = ambient_coordinate_ring(ambient_scheme(domain(f)))
  U = domain(f)
  X = ambient_scheme(U)
  I = ideal(R, one(R))
  for i in 1:npatches(U)
    I = intersect(I, saturated_ideal(ideal(OO(U[i]), pbh[i])))
  end
  return intersect(U, SpecOpen(X, I))
end


function preimage(f::SpecOpenMor, V::SpecOpen)
  U = domain(f)
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = ideal(R, one(R))
  for i in 1:npatches(U)
    I = intersect(I, saturated_ideal(ideal(OO(U[i]), pullback(f[i]).(gens(V)))))
  end
  return intersect(U, SpecOpen(X, I))
end

########################################################################
# Auxiliary methods                                                    #
########################################################################
function is_non_zero_divisor(f::RET, U::SpecOpen) where {RET<:RingElem}
  return all(x->(is_non_zero_divisor(f, x)), affine_patches(U))
end

function find_non_zero_divisor(U::SpecOpen)
  n = length(gens(U))
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  n == 0 && return zero(R)
  kk = base_ring(X)
  coeff = elem_type(kk)[rand(kk, 0:100) for i in 1:n]
  d = sum([coeff[i]*gens(U)[i] for i in 1:n])
  while !is_non_zero_divisor(d, U)
    d = dot([rand(kk, 0:100) for i in 1:n], gens(U))
  end
  return d
end

@Markdown.doc """
    generic_fractions(f::SpecOpenMor)

Given a morphism ``f : U → V`` of Zariski open subsets ``U ⊂ X ⊂ 𝔸ᵐ`` and ``V ⊂ Y ⊂ 𝔸ⁿ``, 
produce a tuple of fractions ``[a₁/b₁,…,aₙ/bₙ]`` such that ``f`` can be recovered 
as the maximal extension of the rational map given by 
```
   U ⊃ U' → 𝔸ⁿ,  x ↦ [a₁(x)/b₁(x),…,aₙ(x)/bₙ(x)]
```
where ``U'`` is the complement of the zero loci of the denominators ``bᵢ`` in ``U``.
In particular, this requires ``U'`` to be dense in ``U`` and this subset 
is chosen at random.
"""
function generic_fractions(f::SpecOpenMor)
  U = domain(f)
  X = ambient_scheme(U)
  V = codomain(f)
  Y = ambient_scheme(V)
  d = find_non_zero_divisor(U)
  W = hypersurface_complement(X, d)
  result = fraction.([restrict(pullback(f)(OO(V)(y)), W) for y in gens(ambient_coordinate_ring(Y))])
  return result
end

