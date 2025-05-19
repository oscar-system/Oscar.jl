####################
### Chain Complexes
####################

@doc raw"""
    chain_complex(V::ModuleFPHom...; seed::Int = 0)

Given a tuple `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the chain complex defined by these homomorphisms.

    chain_complex(V::Vector{<:ModuleFPHom}; seed::Int = 0)

Given a vector `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the chain complex defined by these homomorphisms.

!!! note
    The integer `seed` indicates the lowest homological degree of a module in the complex.

!!! note
    The function checks whether successive homomorphisms indeed compose to zero.
"""
function chain_complex(V::ModuleFPHom...; seed::Int = 0)
  return ComplexOfMorphisms(ModuleFP, collect(V); typ = :chain, seed = seed)
end

function chain_complex(V::Vector{<:ModuleFPHom}; seed::Int = 0, check::Bool=true)
  return ComplexOfMorphisms(ModuleFP, V; typ = :chain, seed = seed, check=check)
end

####################

####################
### Cochain Complexes
####################

@doc raw"""
    cochain_complex(V::ModuleFPHom...; seed::Int = 0)

Given a tuple `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the cochain complex defined by these homomorphisms.

    cochain_complex(V::Vector{<:ModuleFPHom}; seed::Int = 0)

Given a vector `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the cochain complex defined by these homomorphisms.

!!! note
    The integer `seed` indicates the lowest cohomological degree of a module of the complex.

!!! note
    The function checks whether successive homomorphisms indeed compose to zero.
"""
function cochain_complex(V::ModuleFPHom...; seed::Int = 0)
  return ComplexOfMorphisms(ModuleFP, collect(V); typ = :cochain, seed = seed)
end

function cochain_complex(V::Vector{<:ModuleFPHom}; seed::Int = 0)
  return ComplexOfMorphisms(ModuleFP, V; typ = :cochain, seed = seed)
end


function hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{<:ModuleFPHom{<:ModuleFP, <:ModuleFP, Nothing}})
  tM = get_attribute(M, :tensor_product)
  tM === nothing && error("both modules must be tensor products")
  tN = get_attribute(N, :tensor_product)
  tN === nothing && error("both modules must be tensor products")
  @assert length(tM) == length(tN) == length(V)
  @assert all(i-> domain(V[i]) === tM[i] && codomain(V[i]) === tN[i], 1:length(V))
  #gens of M are M[i][j] tensor M[h][l] for i != h and all j, l
  #such a pure tensor is mapped to V[i](M[i][j]) tensor V[h](M[j][l])
  #thus need the pure map - and re-create the careful ordering of the generators as in the 
  # constructor
  #store the maps? and possibly more data, like the ordeing
  decompose_M = get_attribute(M, :tensor_generator_decompose_function)
  pure_N = get_attribute(N, :tensor_pure_function)
  function map_gen(g) # Is there something that generalizes FreeModElem and SubquoModuleElem?
    g_decomposed = decompose_M(g)
    image_as_tuple = Tuple(f(x) for (f,x) in zip(V,g_decomposed))
    res = pure_N(image_as_tuple)
    return res
  end
  return hom(M, N, Vector{elem_type(N)}(map(map_gen, gens(M))); check=false)
end

# the case of a non-trivial base change
@doc raw"""
    hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{<:ModuleFPHom})

Given modules `M`, `N` which are tensor products with the same number of factors,
say $M = M_1 \otimes \cdots \otimes M_r$, $N = N_1 \otimes \cdots \otimes N_r$,
and given a vector `V` of homomorphisms $a_i : M_i \to N_i$, return 
$a_1 \otimes \cdots \otimes a_r$.
"""
function hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{<:ModuleFPHom})
  tM = get_attribute(M, :tensor_product)
  tM === nothing && error("both modules must be tensor products")
  tN = get_attribute(N, :tensor_product)
  tN === nothing && error("both modules must be tensor products")
  @assert length(tM) == length(tN) == length(V)
  @assert all(i-> domain(V[i]) === tM[i] && codomain(V[i]) === tN[i], 1:length(V))
  #gens of M are M[i][j] tensor M[h][l] for i != h and all j, l
  #such a pure tensor is mapped to V[i](M[i][j]) tensor V[h](M[j][l])
  #thus need the pure map - and re-create the careful ordering of the generators as in the 
  # constructor
  #store the maps? and possibly more data, like the ordeing
  decompose_M = get_attribute(M, :tensor_generator_decompose_function)
  pure_N = get_attribute(N, :tensor_pure_function)
  function map_gen(g) # Is there something that generalizes FreeModElem and SubquoModuleElem?
    g_decomposed = decompose_M(g)
    image_as_tuple = Tuple(f(x) for (f,x) in zip(V,g_decomposed))
    res = pure_N(image_as_tuple)
    return res
  end
  bc_map = base_ring_map(first(V))
  bc_map !== nothing && @assert all(domain(g) === domain(bc_map) && codomain(g) === codomain(bc_map) for g in base_ring_map.(V))
  return hom(M, N, Vector{elem_type(N)}(map(map_gen, gens(M))), bc_map)
end


@doc raw"""
    hom_product(M::ModuleFP, N::ModuleFP, A::Matrix{<:ModuleFPHom{<:ModuleFP, <:ModuleFP, Nothing}})

Given modules `M` and `N` which are products with `r` respective `s` factors,  
say $M = \prod_{i=1}^r M_i$, $N = \prod_{j=1}^s N_j$, and given a $r \times s$ matrix 
`A` of homomorphisms $a_{ij} : M_i \to N_j$, return the homomorphism
$M \to N$ with $ij$-components $a_{ij}$.
"""
function hom_product(M::ModuleFP, N::ModuleFP, A::Matrix{<:ModuleFPHom{<:ModuleFP, <:ModuleFP, Nothing}})
  tM = get_attribute(M, :direct_product)
  tM === nothing && error("both modules must be direct products")
  tN = get_attribute(N, :direct_product)
  tN === nothing && error("both modules must be direct products")
  @assert length(tM) == size(A, 1) && length(tN) == size(A, 2)
  @assert all(ij -> domain(A[ij[1],ij[2]]) === tM[ij[1]] && codomain(A[ij[1],ij[2]]) === tN[ij[2]], Base.Iterators.ProductIterator((1:size(A, 1), 1:size(A, 2))))
  #need the canonical maps..., maybe store them as well?
  return hom(M,N,Vector{elem_type(N)}([sum([canonical_injection(N,j)(sum([A[i,j](canonical_projection(M,i)(g)) for i=1:length(tM)])) for j=1:length(tN)]) for g in gens(M)]))
end
# hom(prod -> X), hom(x -> prod)
# if too much time: improve the hom(A, B) in case of A and/or B are products - or maybe not...
# tensor and hom functors for chain complex
# dual: ambig: hom(M, R) or hom(M, Q(R))?

function lift_with_unit(a::FreeModElem{T}, generators::ModuleGens{T}) where {T <: MPolyRingElem}
  # TODO allow optional argument ordering
  # To do this efficiently we need better infrastructure in Singular.jl
  R = base_ring(parent(a))
  singular_assure(generators)
  if Singular.has_global_ordering(base_ring(generators.SF))
    l = lift(a, generators)
    return l, R(1)
  end
  error("Not implemented")
end


#############################
# Tor
#############################
@doc raw"""
    tensor_product(M::ModuleFP, C::ComplexOfMorphisms{ModuleFP})

Return the complex obtained by applying `M` $\otimes\;\! \bullet$ to `C`.
"""
function tensor_product(P::ModuleFP, C::Hecke.ComplexOfMorphisms{ModuleFP})
  #tensor_chain = Hecke.map_type(C)[]
  tensor_chain = valtype(C.maps)[]
  tensor_modules = [tensor_product(P, domain(map(C,first(chain_range(C)))), task=:cache_morphism)[1]]
  append!(tensor_modules, [tensor_product(P, codomain(map(C,i)), task=:cache_morphism)[1] for i in Hecke.map_range(C)])

  for i in 1:length(Hecke.map_range(C))
    A = tensor_modules[i]
    B = tensor_modules[i+1]
    success, A_fac = _is_tensor_product(A)
    @assert success
    success, B_fac = _is_tensor_product(B)
    @assert success
    @assert A_fac[1] === B_fac[1]
    
    j = Hecke.map_range(C)[i]
    @assert domain(map(C, j)) === A_fac[2]
    @assert codomain(map(C, j)) === B_fac[2]
    push!(tensor_chain, hom_tensor(A,B,[identity_map(A_fac[1]), map(C,j)]))
  end

  return Hecke.ComplexOfMorphisms(ModuleFP, tensor_chain, seed=C.seed, typ=C.typ)
end

function tensor_product(M::ModuleFP, F::FreeResolution)
  return tensor_product(M, F.C)
end


@doc raw"""
    tensor_product(C::ComplexOfMorphisms{<:ModuleFP}, M::ModuleFP)

Return the complex obtained by applying $\bullet\;\! \otimes$ `M` to `C`.
"""
function tensor_product(C::Hecke.ComplexOfMorphisms{<:ModuleFP}, P::ModuleFP)
  #tensor_chain = Hecke.map_type(C)[]
  tensor_chain = valtype(C.maps)[]
  tensor_chain = Map[]
  chain_range = Hecke.map_range(C)
  tensor_modules = [tensor_product(domain(map(C,first(chain_range))), P, task=:cache_morphism)[1]]
  append!(tensor_modules, [tensor_product(codomain(map(C,i)), P, task=:cache_morphism)[1] for i in chain_range])

  for i=1:length(chain_range)
    A = tensor_modules[i]
    B = tensor_modules[i+1]

    j = chain_range[i]
    push!(tensor_chain, hom_tensor(A,B,[map(C,j), identity_map(P)]))
  end

  return Hecke.ComplexOfMorphisms(ModuleFP, tensor_chain, seed=C.seed, typ=C.typ)
end

function tensor_product(F::FreeResolution, M::ModuleFP)
  return tensor_product(F.C, M)
end

@doc raw"""
    tor(M::ModuleFP, N::ModuleFP, i::Int)

Return $\text{Tor}_i(M,N)$.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> A = R[x; y];

julia> B = R[x^2; y^3; z^4];

julia> M = SubquoModule(A, B);

julia> F = free_module(R, 1);

julia> Q, _ = quo(F, [x*F[1]]);

julia> T0 = tor(Q, M, 0)
Subquotient of submodule with 2 generators
  1: (e[1] \otimes e[1])
  2: (e[1] \otimes e[2])
by submodule with 7 generators
  1: x*(e[1] \otimes e[1])
  2: -y*(e[1] \otimes e[1]) + x*(e[1] \otimes e[2])
  3: y^2*(e[1] \otimes e[2])
  4: y^3*(e[1] \otimes e[1])
  5: z^4*(e[1] \otimes e[1])
  6: z^4*(e[1] \otimes e[2])
  7: x*(e[1] \otimes e[2])

julia> T1 = tor(Q, M, 1)
Subquotient of submodule with 2 generators
  1: (e[1] \otimes e[1])
  2: x*(e[1] \otimes e[2])
by submodule with 6 generators
  1: x*(e[1] \otimes e[1])
  2: -y*(e[1] \otimes e[1]) + x*(e[1] \otimes e[2])
  3: y^2*(e[1] \otimes e[2])
  4: y^3*(e[1] \otimes e[1])
  5: z^4*(e[1] \otimes e[1])
  6: z^4*(e[1] \otimes e[2])

julia> T2 =  tor(Q, M, 2)
Submodule with 0 generators
represented as subquotient with no relations
```
"""
function tor(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M; length=i+2)
  lifted_resolution = tensor_product(free_res.C[first(chain_range(free_res.C)):-1:1], N) #TODO only three homs are necessary
  return simplify_light(homology(lifted_resolution,i))[1]
end

simplify_light(F::FreeMod) = (F, identity_map(F), identity_map(F))

#TODO, mF
#  (hom lift) => hom and tensor functor
#  filtrations
#  more constructors
#################################################
#
#################################################
@doc raw"""
    lift_homomorphism_contravariant(Hom_MP::ModuleFP, Hom_NP::ModuleFP, a::ModuleFPHom)

Given modules of homomorphism, say, `Hom_MP` $= \text{Hom}(M,P)$ and `Hom_NP` $= \text{Hom}(N,P)$, 
and given a homomorphism `a` $: N \to M$, return the induced homomorphism
$\text{Hom}(M,P) \to \text{Hom}(N,P)$.
"""
function lift_homomorphism_contravariant(Hom_MP::ModuleFP, Hom_NP::ModuleFP, phi::ModuleFPHom)
  # phi : N -> M
  M_P = get_attribute(Hom_MP, :hom)
  M_P === nothing && error("Both modules must be hom modules")
  N_P = get_attribute(Hom_NP, :hom)
  N_P === nothing && error("Both modules must be hom modules")
  
  @assert M_P[2] === N_P[2]
  M,P = M_P
  N,_ = N_P
  @assert domain(phi) === N
  @assert codomain(phi) === M
  
  phi_lifted = hom(Hom_MP, Hom_NP, Vector{elem_type(Hom_NP)}([homomorphism_to_element(Hom_NP, phi*element_to_homomorphism(f)) for f in gens(Hom_MP)]))
  return phi_lifted
end

@doc raw"""
    lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, a::ModuleFPHom)

Given modules of homomorphism, say, `Hom_PM` $= \text{Hom}(P,M)$ and `Hom_PN` $= \text{Hom}(P,N)$,
and given a homomorphism `a` $: M \to N$, return the induced homomorphism
$\text{Hom}(P,M) \to \text{Hom}(P,N)$.
"""
function lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, phi::ModuleFPHom)
  # phi : M -> N
  P_M = get_attribute(Hom_PM, :hom)
  P_M === nothing && error("Both modules must be hom modules")
  P_N = get_attribute(Hom_PN, :hom)
  P_N === nothing && error("Both modules must be hom modules")

  @assert P_M[1] === P_N[1]
  P,M = P_M
  _,N = P_N
  @assert domain(phi) === M
  @assert codomain(phi) === N

  if iszero(Hom_PN)
    return hom(Hom_PM, Hom_PN, Vector{elem_type(Hom_PN)}([zero(Hom_PN) for _=1:ngens(Hom_PM)]))
  end
  phi_lifted = hom(Hom_PM, Hom_PN, Vector{elem_type(Hom_PN)}([homomorphism_to_element(Hom_PN, element_to_homomorphism(f)*phi) for f in gens(Hom_PM)]))
  return phi_lifted
end

@doc raw"""
    hom(M::ModuleFP, C::ComplexOfMorphisms{ModuleFP})

Return the complex obtained by applying $\text{Hom}($`M`, $-)$ to `C`.
"""
function hom(P::ModuleFP, C::Hecke.ComplexOfMorphisms{ModuleFP})
  #hom_chain = Hecke.map_type(C)[]
  hom_chain = valtype(C.maps)[]
  chain_range = Hecke.map_range(C)
  hom_modules = [hom(P, domain(map(C,first(chain_range))))]
  append!(hom_modules, [hom(P, codomain(map(C,i))) for i in chain_range])

  for i=1:length(chain_range)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    j = chain_range[i]
    push!(hom_chain, lift_homomorphism_covariant(A,B,map(C,j)))
  end

  return Hecke.ComplexOfMorphisms(ModuleFP, hom_chain, seed=C.seed, typ=C.typ)
end

function hom(M::ModuleFP, F::FreeResolution)
  return hom(M, F.C)
end

@doc raw"""
    hom(C::ComplexOfMorphisms{ModuleFP}, M::ModuleFP)

Return the complex obtained by applying $\text{Hom}(-,$ `M`$)$ to `C`.

If `C` is a chain complex, return a cochain complex.
If `C` is a cochain complex, return a chain complex.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, [:x]);

julia> F = free_module(R, 1);

julia> A, _ = quo(F, [x^4*F[1]]);

julia> B, _ = quo(F, [x^3*F[1]]);

julia> a = hom(A, B, [x^2*B[1]]);

julia> b = hom(B, B, [x^2*B[1]]);

julia> C = chain_complex([a, b]; seed = 3);

julia> range(C)
5:-1:3

julia> D = hom(C, A);

julia> range(D)
3:5
```
"""
function hom(C::Hecke.ComplexOfMorphisms{T}, P::ModuleFP) where {T<:ModuleFP}
  #hom_chain = Hecke.map_type(C)[]
  hom_chain = valtype(C.maps)[]
  hom_chain = Map[]
  chain_range = Hecke.map_range(C)
  hom_modules = Tuple{ModuleFP, Map}[hom(domain(map(C,first(chain_range))),P)]
  append!(hom_modules, Tuple{ModuleFP, Map}[hom(codomain(map(C,i)), P) for i in chain_range])

  for i=1:length(chain_range)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    j = chain_range[i]
    push!(hom_chain, lift_homomorphism_contravariant(B,A,map(C,j)))
  end

  typ = Hecke.is_chain_complex(C) ? :cochain : :chain
  seed = C.seed
  return Hecke.ComplexOfMorphisms(ModuleFP, reverse(hom_chain), seed=seed, typ=typ)
end

function hom(F::FreeResolution, M::ModuleFP)
  return hom(F.C, M)
end

@doc raw"""
    hom_without_reversing_direction(C::ComplexOfMorphisms{ModuleFP}, M::ModuleFP)

Return the complex obtained by applying $\text{Hom}(-,$ `M`$)$ to `C`.

If `C` is a chain complex, return a chain complex.
If `C` is a cochain complex, return a cochain complex.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, [:x]);

julia> F = free_module(R, 1);

julia> A, _ = quo(F, [x^4*F[1]]);

julia> B, _ = quo(F, [x^3*F[1]]);

julia> a = hom(A, B, [x^2*B[1]]);

julia> b = hom(B, B, [x^2*B[1]]);

julia> C = chain_complex([a, b]; seed = 3);

julia> range(C)
5:-1:3

julia> D = hom_without_reversing_direction(C, A);

julia> range(D)
-3:-1:-5
```
"""
function hom_without_reversing_direction(C::Hecke.ComplexOfMorphisms{ModuleFP}, P::ModuleFP)
  #up to seed/ typ identical to the one above. Should be
  #ONE worker function with 2 interfaces.
  #hom_chain = Hecke.map_type(C)[]
  hom_chain = valtype(C.maps)[]
  m_range = Hecke.map_range(C)
  hom_modules = [hom(domain(map(C,first(m_range))),P)]
  append!(hom_modules, [hom(codomain(map(C,i)), P) for i in m_range])

  for i=1:length(m_range)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    j = m_range[i]
    push!(hom_chain, lift_homomorphism_contravariant(B,A,map(C,j)))
  end

  return Hecke.ComplexOfMorphisms(ModuleFP, reverse(hom_chain), seed=-first(chain_range(C)), typ=C.typ)
end

function hom_without_reversing_direction(F::FreeResolution, M::ModuleFP)
  return hom_without_reversing_direction(F.C, M)
end


#############################
@doc raw"""
    homology(C::ComplexOfMorphisms{<:ModuleFP})

Return the homology of `C`.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, [:x]);

julia> F = free_module(R, 1);

julia> A, _ = quo(F, [x^4*F[1]]);

julia> B, _ = quo(F, [x^3*F[1]]);

julia> a = hom(A, B, [x^2*B[1]]);

julia> b = hom(B, B, [x^2*B[1]]);

julia> C = ComplexOfMorphisms(ModuleFP, [a, b]);

julia> H = homology(C)
3-element Vector{SubquoModule{QQMPolyRingElem}}:
 Subquotient of submodule with 1 generator
  1: x*e[1]
by submodule with 1 generator
  1: x^4*e[1]
 Subquotient of submodule with 1 generator
  1: x*e[1]
by submodule with 2 generators
  1: x^3*e[1]
  2: x^2*e[1]
 Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 2 generators
  1: x^3*e[1]
  2: x^2*e[1]
```
"""
function homology(C::Hecke.ComplexOfMorphisms{<:ModuleFP})
  return [homology(C,i) for i in Hecke.range(C)]
end

function homology(C::FreeResolution)
  return homology(C.C)
end


@doc raw"""
    homology(C::ComplexOfMorphisms{<:ModuleFP}, i::Int)

Return the `i`-th homology module of `C`.

# Examples
```jldoctest
julia> R, (x,) = polynomial_ring(QQ, [:x]);

julia> F = free_module(R, 1);

julia> A, _ = quo(F, [x^4*F[1]]);

julia> B, _ = quo(F, [x^3*F[1]]);

julia> a = hom(A, B, [x^2*B[1]]);

julia> b = hom(B, B, [x^2*B[1]]);

julia> C = ComplexOfMorphisms(ModuleFP, [a, b]);

julia> H = homology(C, 1)
Subquotient of submodule with 1 generator
  1: x*e[1]
by submodule with 2 generators
  1: x^3*e[1]
  2: x^2*e[1]
```
"""
function homology(C::Hecke.ComplexOfMorphisms{<:ModuleFP}, i::Int)
  chain_range = Hecke.range(C)
  map_range = Hecke.map_range(C)
  @assert length(chain_range) > 0 #TODO we need actually only the base ring
  if i == first(chain_range)
    return kernel(map(C, first(map_range)))[1]
  elseif i == last(chain_range)
    f = map(C,last(map_range))
    return cokernel(f)    
  elseif i in chain_range
    if Hecke.is_chain_complex(C)
      return quo_object(kernel(map(C,i))[1], image(map(C,i+1))[1])
    else
      return quo_object(kernel(map(C,i))[1], image(map(C,i-1))[1])
    end
  else
    return FreeMod(base_ring(obj(C,first(chain_range))),0)
  end
end

#############################
# Ext
#############################
@doc raw"""
    ext(M::ModuleFP, N::ModuleFP, i::Int)

Return $\text{Ext}^i(M,N)$.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = FreeMod(R, 1);

julia> V = [x*F[1], y*F[1]];

julia> M = quo_object(F, V)
Subquotient of submodule with 1 generator
  1: e[1]
by submodule with 2 generators
  1: x*e[1]
  2: y*e[1]

julia> ext(M, M, 0)
Subquotient of submodule with 1 generator
  1: (e[1] -> e[1])
by submodule with 2 generators
  1: y*(e[1] -> e[1])
  2: x*(e[1] -> e[1])

julia> ext(M, M, 1)
Subquotient of submodule with 2 generators
  1: (e[1] -> e[1])
  2: (e[2] -> e[1])
by submodule with 4 generators
  1: y*(e[1] -> e[1])
  2: x*(e[1] -> e[1])
  3: y*(e[2] -> e[1])
  4: x*(e[2] -> e[1])

julia> ext(M, M, 2)
Subquotient of submodule with 1 generator
  1: (e[1] -> e[1])
by submodule with 2 generators
  1: y*(e[1] -> e[1])
  2: x*(e[1] -> e[1])

julia> ext(M, M, 3)
Submodule with 0 generators
represented as subquotient with no relations
```
"""
function ext(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M; length=i+2)
  
  lifted_resolution = hom(free_res.C[first(Hecke.map_range(free_res.C)):-1:1], N) #TODO only three homs are necessary
  return simplify_light(homology(lifted_resolution,i))[1]
end

