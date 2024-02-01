@doc raw"""
    iszero(f::ModuleFPHom)

Return true iff `f` is the zero map.
"""
function iszero(f::ModuleFPHom)
  return all(iszero, map(f, gens(domain(f))))
end

@doc raw"""
    hom(M::ModuleFP, N::ModuleFP)

Return the module `Hom(M,N)` as an object of type `SubquoModule`.

Additionally, if `H` is that object, return the map which sends an element of `H` to the corresponding homomorphism `M` $\to$ `N`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> F = FreeMod(R, 2);

julia> V = [x*F[1], y^2*F[2]];

julia> M = quo(F, V)[1]
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y^2*e[2]

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
function hom(M::ModuleFP, N::ModuleFP, algorithm::Symbol=:maps)
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

  #step 2
  H_s0_t0, mH_s0_t0 = hom(domain(f0), domain(g0))
  H_s1_t1, mH_s1_t1 = hom(domain(f1), domain(g1))
  D, pro = direct_product(H_s0_t0, H_s1_t1, task = :prod)

  H_s1_t0, mH_s1_t0 = hom(domain(f1), domain(g0))

  delta = hom(D, H_s1_t0, Vector{elem_type(H_s1_t0)}([preimage(mH_s1_t0, f1*mH_s0_t0(pro[1](g))-mH_s1_t1(pro[2](g))*g1) for g = gens(D)]))

  H_s0_t1, mH_s0_t1 = hom(domain(f0), domain(g1))

  rho_prime = hom(H_s0_t1, H_s0_t0, Vector{elem_type(H_s0_t0)}([preimage(mH_s0_t0, mH_s0_t1(C)*g1) for C in gens(H_s0_t1)]))
 
  kDelta = kernel(delta)

  projected_kernel::Vector{elem_type(H_s0_t0)} = filter(v -> !is_zero(v), FreeModElem[pro[1](repres(AB)) for AB in gens(kDelta[1])])
  H = quo(sub(H_s0_t0, projected_kernel, :module), image(rho_prime)[1], :module)

  H_simplified, s_inj, s_proj = simplify_light(H)

  function im(x::SubquoModuleElem)
    #@assert parent(x) === H
    @assert parent(x) === H_simplified
    return hom(M, N, Vector{elem_type(N)}([g0(mH_s0_t0(repres(s_inj(x)))(preimage(f0, g))) for g = gens(M)]))
  end

  function pre(f::ModuleFPHom)
    @assert domain(f) === M
    @assert codomain(f) === N
    Rs0 = domain(f0)
    Rt0 = domain(g0)
    g = hom(Rs0, Rt0, Vector{elem_type(Rt0)}([preimage(g0, f(f0(g))) for g = gens(Rs0)]))

    return s_proj(SubquoModuleElem(repres(preimage(mH_s0_t0, g)), H))
  end
  to_hom_map = MapFromFunc(H_simplified, Hecke.MapParent(M, N, "homomorphisms"), im, pre)
  set_attribute!(H_simplified, :show => Hecke.show_hom, :hom => (M, N), :module_to_hom_map => to_hom_map)
  return H_simplified, to_hom_map
end

@doc raw"""
    element_to_homomorphism(f::ModuleFPElem)

If `f` is an element of a module created via `hom(M,N)`, for some modules `M` and `N`, 
return the homomorphism `M` $\to$ `N` corresponding to `f`.
# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> F = FreeMod(R, 2);

julia> V = [x*F[1], y^2*F[2]];

julia> M = quo(F, V)[1]
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y^2*e[2]

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
Map with following data
Domain:
=======
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y^2*e[2]
Codomain:
=========
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y^2*e[2]

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
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> F = FreeMod(R, 2);

julia> V = [x*F[1], y^2*F[2]];

julia> M = quo(F, V)[1]
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y^2*e[2]

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
  map_to_hom = to_hom_map.g
  return map_to_hom(a)
end

@doc raw"""
    multiplication_morphism(a::RingElem, M::ModuleFP)

Return the multiplication by `a` as an endomorphism on `M`.
"""
function multiplication_morphism(a::RingElem, M::ModuleFP)
  @assert base_ring(M) === parent(a)
  return hom(M, M, [a*v for v in gens(M)])
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
  return hom(F, H, [homomorphism_to_element(H, multiplication_morphism(F[1], M))])
end


