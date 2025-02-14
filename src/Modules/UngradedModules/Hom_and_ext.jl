@doc raw"""
    iszero(f::ModuleFPHom)

Return true iff `f` is the zero map.
"""
function iszero(f::ModuleFPHom)
  return all(iszero, map(f, gens(domain(f))))
end

@doc raw"""
    hom(M::ModuleFP, N::ModuleFP; algorithm::Symbol=:maps)

Return the module `Hom(M,N)` as an object of type `SubquoModule`.

Additionally, if `H` is that object, return the map which sends an element of `H` to the corresponding homomorphism `M` $\to$ `N`.

The keyword `algorithm` can be set to `:maps` for the default algorithm or to `:matrices` for an alternative based on matrices.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = FreeMod(R, 2);

julia> V = [x*F[1], y^2*F[2]];

julia> M = quo_object(F, V)
Subquotient of submodule with 2 generators
  1: e[1]
  2: e[2]
by submodule with 2 generators
  1: x*e[1]
  2: y^2*e[2]

julia> H = hom(M, M)[1]
hom of (M, M)

julia> gens(H)
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 (e[1] -> e[1])
 (e[2] -> e[2])

julia> relations(H)
4-element Vector{FreeModElem{QQMPolyRingElem}}:
 x*(e[1] -> e[1])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])
 y^2*(e[2] -> e[2])
```
"""
function hom(M::ModuleFP, N::ModuleFP; algorithm::Symbol=:maps)
  #source: Janko's CA script: https://www.mathematik.uni-kl.de/~boehm/lehre/17_CA/ca.pdf
  if algorithm == :matrices && M isa SubquoModule && N isa SubquoModule
    if is_graded(M) && is_graded(N)
      error("This algorithm is not implemented for graded modules.")
    end
    return hom_matrices(M,N,false)
  end
  p1 = presentation(M)
  p2 = presentation(N)

  f0 = map(p1, 0)
  f1 = map(p1, 1)
  g0 = map(p2, 0)
  g1 = map(p2, 1)
  if is_graded(M) && is_graded(N)
    @assert is_graded(f0)
    @assert is_graded(f1)
    @assert is_graded(g0)
    @assert is_graded(g1)
  end

  #step 2
  H_s0_t0, mH_s0_t0 = hom(domain(f0), domain(g0))
  H_s1_t1, mH_s1_t1 = hom(domain(f1), domain(g1))
  D, pro = direct_product(H_s0_t0, H_s1_t1, task = :prod)

  H_s1_t0, mH_s1_t0 = hom(domain(f1), domain(g0))

  delta = hom(D, H_s1_t0, elem_type(H_s1_t0)[preimage(mH_s1_t0, f1*mH_s0_t0(pro[1](g))-mH_s1_t1(pro[2](g))*g1) for g = gens(D)])

  H_s0_t1, mH_s0_t1 = hom(domain(f0), domain(g1))

  rho_prime = hom(H_s0_t1, H_s0_t0, elem_type(H_s0_t0)[preimage(mH_s0_t0, mH_s0_t1(C)*g1) for C in gens(H_s0_t1)])
 
  kDelta = kernel(delta)

  projected_kernel::Vector{elem_type(H_s0_t0)} = filter(!is_zero, FreeModElem[pro[1](repres(AB)) for AB in gens(kDelta[1])])
  H = quo_object(sub_object(H_s0_t0, projected_kernel), image(rho_prime)[1])

  H_simplified, s_inj, s_proj = simplify_light(H)

  function im(x::SubquoModuleElem)
    #@assert parent(x) === H
    @assert parent(x) === H_simplified
    return hom(M, N, elem_type(N)[g0(mH_s0_t0(repres(s_inj(x)))(preimage(f0, g))) for g = gens(M)])
  end

  function pre(f::ModuleFPHom)
    @assert domain(f) === M
    @assert codomain(f) === N
    Rs0 = domain(f0)
    Rt0 = domain(g0)
    g = hom(Rs0, Rt0, elem_type(Rt0)[preimage(g0, f(f0(g))) for g = gens(Rs0)])

    return s_proj(SubquoModuleElem(repres(preimage(mH_s0_t0, g)), H))
  end
  to_hom_map = MapFromFunc(H_simplified, Hecke.MapParent(M, N, "homomorphisms"), im, pre)
  set_attribute!(H_simplified, :show => Hecke.show_hom, :hom => (M, N), :module_to_hom_map => to_hom_map)
  return H_simplified, to_hom_map
end

#=
### New and hopefully more maintainable code
# Disabled for the moment because of the book chapter. We want to keep code output stable.
# Eventually we will gradually merge this in.
function hom(F::FreeMod, G::ModuleFP)
  R = base_ring(F)
  R === base_ring(G) || error("base rings must be the same")
  
  # We construct Hom(F, G) as Hom(F, R¹) ⊗ G.
  # Then we can use coordinate independent methods to construct 
  # induced homomorphisms.
  Fdual, interp = dual(F)
  H, mult_map = tensor_product(Fdual, G; task=:with_map)

  function _elem_to_hom1(v::ModuleFPElem)
    result = hom(F, G, elem_type(G)[zero(G) for i in 1:ngens(F)]; check=false)
    for (i, c) in coordinates(v)
      # decompose the i-th generator as ψ : F → R¹ and g ∈ G
      (u, g) = preimage(mult_map, H[i])
      psi = element_to_homomorphism(u)
      # construct the homomorphism by applying it manually 
      # to the generators; psi(w)[1] is the value considered 
      # as an element of the ring R ≅ R¹.
      result = result + hom(F, G, elem_type(G)[c*psi(w)[1] * g for w in gens(F)]; check=false)
    end
    return result
  end

  function _hom_to_elem1(phi::ModuleFPHom)
    @assert domain(phi) === F
    @assert codomain(phi) === G
    return sum(mult_map((Fdual[i], phi(v))) for (i, v) in enumerate(gens(F)); init=zero(H))
  end
  
  to_hom_map1 = MapFromFunc(H, Hecke.MapParent(F, G, "homomorphisms"), _elem_to_hom1, _hom_to_elem1)
  set_attribute!(H, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map1)
  return H, to_hom_map1
end

function hom(M::SubquoModule, N::ModuleFP)
  R = base_ring(M)
  R === base_ring(N) || error("base rings must coincide")
  pres = presentation(M)
  # pres :   F1 --phi--> F0 --eps--> M --> 0
  # Construct the induced hom-sequence
  # hom(pres, N) : hom(F1, N) <--phi^T-- hom(F0, N) <--eps^T-- hom(M, N) <-- 0
  F0 = pres[0]::FreeMod
  F1 = pres[1]::FreeMod
  phi = map(pres, 1)
  @assert codomain(phi) === F0
  @assert domain(phi) === F1
  eps = map(pres, 0)

  phi_transp = hom(phi, N) # The induced map in the first argument
  hom_F0_N = domain(phi_transp)
  @assert get_attribute(hom_F0_N, :hom) === (F0, N)
  hom_F1_N = codomain(phi_transp)
  @assert get_attribute(hom_F1_N, :hom) === (F1, N)

  K, inc = kernel(phi_transp)
  
  # The output will look horribly difficult in general. Simplification is indicated.
  # H, iso = prune_with_map(K) # Still buggy
  # H, iso_inv, iso = Oscar._alt_simplify(K) # Takes too long in some cases
  H, iso, iso_inv = simplify_light(K)


  function _elem_to_hom2(v::ModuleFPElem)
    w = inc(iso(v)) 
    @assert parent(w) === hom_F0_N
    psi = element_to_homomorphism(w)
    img_gens = elem_type(N)[psi(preimage(eps, g)) for g in gens(M)]
    return hom(M, N, img_gens; check=false)
  end

  function _hom_to_elem2(phi::ModuleFPHom)
    # phi : M --> N
    @assert domain(phi) === M
    @assert codomain(phi) === N
    # The following is using that we have a 1:1-correspondence of the generators 
    # of M with those of F0 due to the implementation of `presentation`.
    phi_lift = hom(F0, N, images_of_generators(phi); check=false)::FreeModuleHom{<:FreeMod}
    w = homomorphism_to_element(hom_F0_N, phi_lift)
    return preimage(iso, preimage(inc, w))
  end

  to_hom_map2 = MapFromFunc(H, Hecke.MapParent(M, N, "homomorphisms"), _elem_to_hom2, _hom_to_elem2)
  set_attribute!(H, :show => Hecke.show_hom, :hom => (M, N), :module_to_hom_map => to_hom_map2)
  return H, to_hom_map2
end

# The induced hom on the first argument
function hom(
    phi::FreeModuleHom{<:FreeMod, <:FreeMod}, N::ModuleFP; 
    domain::ModuleFP=hom(Oscar.codomain(phi), N)[1],
    codomain::ModuleFP=hom(Oscar.domain(phi), N)[1]
  )
  R = base_ring(N)
  R === base_ring(domain) === base_ring(codomain) || error("base rings must coincide")

  img_gens = elem_type(codomain)[]
  for (i, g) in enumerate(gens(domain))
    # The map codomain(phi) --> N
    psi = element_to_homomorphism(g)
    comp = compose(phi, psi)
    push!(img_gens, homomorphism_to_element(codomain, comp))
  end
  return hom(domain, codomain, img_gens; check=false)
end

@doc raw"""
    hom(F::FreeMod, G::FreeMod)

Return a free module $S$ such that $\text{Hom}(F,G) \cong S$ along with a function 
that converts elements from $S$ into morphisms $F \to G$.

# Examples
```jldoctest
julia> R, _ = polynomial_ring(QQ, [:x, :y, :z]);

julia> F1 = free_module(R, 3)
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ

julia> F2 = free_module(R, 2)
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> V, f = hom(F1, F2)
(hom of (F1, F2), Map: V -> set of all homomorphisms from F1 to F2)

julia> f(V[1])
Map with following data
Domain:
=======
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ
Codomain:
=========
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F1 = graded_free_module(Rg, [1,2,2])
Graded free module Rg^1([-1]) + Rg^2([-2]) of rank 3 over Rg

julia> F2 = graded_free_module(Rg, [3,5])
Graded free module Rg^1([-3]) + Rg^1([-5]) of rank 2 over Rg

julia> V, f = hom(F1, F2)
(hom of (F1, F2), Map: V -> set of all homomorphisms from F1 to F2)

julia> f(V[1])
F1 -> F2
e[1] -> e[1]
e[2] -> 0
e[3] -> 0
Graded module homomorphism of degree [2]

```
"""
function hom(F::FreeMod, G::FreeMod)
  if is_graded(F) && is_graded(G)
    return _hom_graded(F, G)
  else
    return _hom_simple(F, G)
  end
end

function _hom_simple(F::FreeMod, G::FreeMod)
  R = base_ring(F)
  R === base_ring(G) || error("base rings must be the same")
  m = ngens(F)
  n = ngens(G)

  mn = m*n
  H = FreeMod(R, mn)
  # Custom printing of hom-module generators
  if rank(G) == 1
    H.S = [Symbol("($i)*") for i = F.S]
  else
    H.S = [Symbol("($i "*(is_unicode_allowed() ? "→" : "->")*" $j)") for i = F.S for j = G.S]
  end

  # We think of elements of H as matrices A ∈ Rᵐˣⁿ stored 
  # in concatenated lines: 
  # v = (a₁₁, a₁₂, …, a₁ₙ, a₂₁, …)

  function _elem_to_hom3(v::FreeModElem)
    img_gens = [sum(c*G[j] for (j, c) in coordinates(v)[(i-1)*n+1:i*n]; init=zero(G)) for i in 1:m]
    return hom(F, G, img_gens; check=false)
  end

  function _hom_to_elem3(phi::FreeModuleHom)
    @assert domain(phi) === F
    @assert codomain(phi) === G
    v = zero(H)
    for (i, g) in enumerate(images_of_generators(phi))
      v = v + sum(c*H[(i-1)*n + j] for (j, c) in coordinates(g); init=zero(H))
    end
    return v
  end
  
  to_hom_map3 = MapFromFunc(H, Hecke.MapParent(F, G, "homomorphisms"), _elem_to_hom3, _hom_to_elem3)
  set_attribute!(H, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map3)
  return H, to_hom_map3
end

function _hom_graded(F::FreeMod, G::FreeMod)
  R = base_ring(F)
  GG = grading_group(R)
  R === base_ring(G) || error("base rings must be the same")
  m = ngens(F)
  n = ngens(G)

  mn = m*n
  H = graded_free_module(R, [_degree_fast(G[j]) - _degree_fast(F[i]) for i in 1:m for j in 1:n])
  # Custom printing of hom-module generators
  if rank(G) == 1 && iszero(_degree_fast(G[1]))
    H.S = [Symbol("($i)*") for i = F.S]
  else
    H.S = [Symbol("($i "*(is_unicode_allowed() ? "→" : "->")*" $j)") for i = F.S for j = G.S]
  end

  # We think of elements of H as matrices A ∈ Rᵐˣⁿ stored 
  # in concatenated lines: 
  # v = (a₁₁, a₁₂, …, a₁ₙ, a₂₁, …)

  function _elem_to_hom4(v::FreeModElem)
    img_gens = [sum(c*G[j] for (j, c) in coordinates(v)[(i-1)*n+1:i*n]; init=zero(G)) for i in 1:m]
    return hom(F, G, img_gens; check=false)
  end

  function _hom_to_elem4(phi::FreeModuleHom)
    @assert domain(phi) === F
    @assert codomain(phi) === G
    v = zero(H)
    for (i, g) in enumerate(images_of_generators(phi))
      v = v + sum(c*H[(i-1)*n + j] for (j, c) in coordinates(g); init=zero(H))
    end
    return v
  end
  
  to_hom_map4 = MapFromFunc(H, Hecke.MapParent(F, G, "homomorphisms"), _elem_to_hom4, _hom_to_elem4)
  set_attribute!(H, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map4)
  return H, to_hom_map4
end
=#


@doc raw"""
    element_to_homomorphism(f::ModuleFPElem)

If `f` is an element of a module created via `hom(M,N)`, for some modules `M` and `N`, 
return the homomorphism `M` $\to$ `N` corresponding to `f`.
# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = FreeMod(R, 2);

julia> V = [x*F[1], y^2*F[2]];

julia> M = quo_object(F, V)
Subquotient of submodule with 2 generators
  1: e[1]
  2: e[2]
by submodule with 2 generators
  1: x*e[1]
  2: y^2*e[2]

julia> H = hom(M, M)[1];

julia> gens(H)
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 (e[1] -> e[1])
 (e[2] -> e[2])

julia> relations(H)
4-element Vector{FreeModElem{QQMPolyRingElem}}:
 x*(e[1] -> e[1])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])
 y^2*(e[2] -> e[2])

julia> a = element_to_homomorphism(H[1]+y*H[2])
Module homomorphism
  from M
  to M

julia> matrix(a)
[1   0]
[0   y]
```
"""
function element_to_homomorphism(f::ModuleFPElem)
  H = f.parent
  to_hom_map = get_attribute(H, :module_to_hom_map)
  to_hom_map === nothing && error("element doesn't live in a hom module")  
  return to_hom_map(f)
end

@doc raw"""
    homomorphism_to_element(H::ModuleFP, a::ModuleFPHom)

If the module `H` is created via `hom(M,N)`, for some modules `M` and `N`, and
`a`: `M` $\to$ `N` is a homomorphism, then return the element of `H` corresponding to `a`.
# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = FreeMod(R, 2);

julia> V = [x*F[1], y^2*F[2]];

julia> M = quo_object(F, V)
Subquotient of submodule with 2 generators
  1: e[1]
  2: e[2]
by submodule with 2 generators
  1: x*e[1]
  2: y^2*e[2]

julia> H = hom(M, M)[1];

julia> gens(H)
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 (e[1] -> e[1])
 (e[2] -> e[2])

julia> relations(H)
4-element Vector{FreeModElem{QQMPolyRingElem}}:
 x*(e[1] -> e[1])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])
 y^2*(e[2] -> e[2])

julia> W =  [M[1], y*M[2]];

julia> a = hom(M, M, W);

julia> is_welldefined(a)
true

julia> matrix(a)
[1   0]
[0   y]

julia> m = homomorphism_to_element(H, a)
(e[1] -> e[1]) + y*(e[2] -> e[2])
```
"""
function homomorphism_to_element(H::ModuleFP, a::ModuleFPHom)
  to_hom_map = get_attribute(H, :module_to_hom_map)
  to_hom_map === nothing && error("module must be a hom module")
  map_to_hom = pseudo_inv(to_hom_map)
  return map_to_hom(a)
end

@doc raw"""
    multiplication_morphism(a::RingElem, M::ModuleFP)

Return the multiplication by `a` as an endomorphism on `M`.
"""
function multiplication_morphism(a::RingElem, M::ModuleFP)
  @assert base_ring(M) === parent(a)
  return hom(M, M, [a*v for v in gens(M)]; check=false)
end

@doc raw"""
    multiplication_morphism(a::FreeModElem, M::ModuleFP)

Return the multiplication by `a` as an endomorphism on `M`. For this,
the parent of `a` must be a module of rank 1.
"""
function multiplication_morphism(a::FreeModElem, M::ModuleFP)
  @assert rank(parent(a)) == 1
  return multiplication_morphism(a[1], M)
end

@doc raw"""
    multiplication_induced_morphism(F::FreeMod, H::ModuleFP)

Let `H` be the module of endomorphisms on a module `M`. (If this is not
the case an error is thrown.) Let `F` be free of rank 1. Return the 
morphism from `F` to `H` which sends an element of `F` to its
corresponding multiplication morphism.
"""
function multiplication_induced_morphism(F::FreeMod, H::ModuleFP)
  @assert rank(F) == 1
  M_N = get_attribute(H, :hom)
  M_N === nothing && error("module must be a hom module")
  M,N = M_N
  @assert M === N
  return hom(F, H, [homomorphism_to_element(H, multiplication_morphism(F[1], M))]; check=false)
end


