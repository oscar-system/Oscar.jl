###############################################################################
# FreeModuleHom constructors
###############################################################################

#=@doc raw"""
    FreeModuleHom(F::FreeMod{T}, G::S, a::Vector) where {T, S}

Construct the morphism $F \to G$ where `F[i]` is mapped to `a[i]`.
In particular, `ngens(F) == length(a)` must hold.
"""
FreeModuleHom(F::AbstractFreeMod{T}, G::S, a::Vector) where {T, S} = FreeModuleHom{T,S}(F, G, a)

@doc raw"""
    FreeModuleHom(F::FreeMod{T}, G::S, mat::MatElem{T}) where {T,S}

Construct the morphism $F \to G$ corresponding to the matrix `mat`.
"""
FreeModuleHom(F::AbstractFreeMod{T}, G::S, mat::MatElem{T}) where {T,S} = FreeModuleHom{T,S}(F, G, mat)=#

img_gens(f::FreeModuleHom) = images_of_generators(f)
images_of_generators(f::FreeModuleHom) = f.imgs_of_gens::Vector{elem_type(codomain(f))}
image_of_generator(phi::FreeModuleHom, i::Int) = phi.imgs_of_gens[i]::elem_type(codomain(phi))
base_ring_map(f::FreeModuleHom) = f.ring_map
@attr Map function base_ring_map(f::FreeModuleHom{<:SubquoModule, <:ModuleFP, Nothing})
    return identity_map(base_ring(domain(f)))
end
base_ring_map(f::SubQuoHom) = f.ring_map
@attr Map function base_ring_map(f::SubQuoHom{<:SubquoModule, <:ModuleFP, Nothing})
    return identity_map(base_ring(domain(f)))
end

@doc raw"""
    matrix(a::FreeModuleHom)

Given a homomorphism `a : F → M` of type  `FreeModuleHom`, 
return a matrix `A` over `base_ring(M)` with `rank(F)` rows and 
`ngens(M)` columns such that $a(F[i]) = \sum_j A[i,j]*M[j]$.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]];

julia> a = hom(F, G, V);

julia> matrix(a)
[y   0]
[x   y]
[0   z]
```
"""
function matrix(f::FreeModuleHom)
  if !isdefined(f, :matrix)
    D = domain(f)
    C = codomain(f)
    R = base_ring(C)
    matrix = zero_matrix(R, rank(D), ngens(C))
    for i=1:rank(D)
      image_of_gen = f(D[i])
      for j=1:ngens(C)
        matrix[i,j] = image_of_gen[j]
      end
    end
    setfield!(f, :matrix, matrix)
  end
  return f.matrix
end

(h::FreeModuleHom)(a::AbstractFreeModElem) = image(h, a)

@doc raw"""
    hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T

Given a vector `V` of `rank(F)` elements of `M`, 
return the homomorphism `F` $\to$ `M` which sends the `i`-th
basis vector of `F` to the `i`-th entry of `V`.

    hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}) where T

Given a matrix `A` with `rank(F)` rows and `ngens(M)` columns, return the
homomorphism `F` $\to$ `M` which sends the `i`-th basis vector of `F` to 
the linear combination $\sum_j A[i,j]*M[j]$ of the generators `M[j]` of `M`.

!!! note
    The module `M` may be of type `FreeMod` or `SubquoMod`. If both modules
    `F` and `M` are graded, the data must define a graded module homomorphism of some degree.
    If this degree is the zero element of the (common) grading group, we refer to
    the homomorphism under consideration as a *homogeneous module homomorphism*.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
3-element Vector{FreeModElem{QQMPolyRingElem}}:
 y*e[1]
 x*e[1] + y*e[2]
 z*e[2]

julia> a = hom(F, G, V)
Map with following data
Domain:
=======
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ
Codomain:
=========
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> a(F[2])
x*e[1] + y*e[2]

julia> B = R[y 0; x y; 0 z]
[y   0]
[x   y]
[0   z]

julia> b = hom(F, G, B)
Map with following data
Domain:
=======
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ
Codomain:
=========
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> a == b
true
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F1 = graded_free_module(Rg, 3)
Graded free module Rg^3([0]) of rank 3 over Rg

julia> G1 = graded_free_module(Rg, 2)
Graded free module Rg^2([0]) of rank 2 over Rg

julia> V1 = [y*G1[1], (x+y)*G1[1]+y*G1[2], z*G1[2]]
3-element Vector{FreeModElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y*e[1]
 (x + y)*e[1] + y*e[2]
 z*e[2]

julia> a1 = hom(F1, G1, V1)
F1 -> G1
e[1] -> y*e[1]
e[2] -> (x + y)*e[1] + y*e[2]
e[3] -> z*e[2]
Graded module homomorphism of degree [1]

julia> F2 = graded_free_module(Rg, [1,1,1])
Graded free module Rg^3([-1]) of rank 3 over Rg

julia> G2 = graded_free_module(Rg, [0,0])
Graded free module Rg^2([0]) of rank 2 over Rg

julia> V2 = [y*G2[1], (x+y)*G2[1]+y*G2[2], z*G2[2]]
3-element Vector{FreeModElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 y*e[1]
 (x + y)*e[1] + y*e[2]
 z*e[2]

julia> a2 = hom(F2, G2, V2)
F2 -> G2
e[1] -> y*e[1]
e[2] -> (x + y)*e[1] + y*e[2]
e[3] -> z*e[2]
Homogeneous module homomorphism

julia> B = Rg[y 0; x+y y; 0 z]
[    y   0]
[x + y   y]
[    0   z]

julia> b = hom(F2, G2, B)
F2 -> G2
e[1] -> y*e[1]
e[2] -> (x + y)*e[1] + y*e[2]
e[3] -> z*e[2]
Homogeneous module homomorphism

julia> a2 == b
true
```
"""
function hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}; check::Bool=true) where T
  base_ring(F) === base_ring(M) || return FreeModuleHom(F, M, V, base_ring(M))
  return FreeModuleHom(F, M, V)
end
function hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}; check::Bool=true) where T 
  base_ring(F) === base_ring(M) || return FreeModuleHom(F, M, A, base_ring(M))
  return FreeModuleHom(F, M, A)
end

@doc raw"""
    hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, h::RingMapType) where {T, RingMapType}

Given a vector `V` of `rank(F)` elements of `M` and a ring map `h`
from `base_ring(F)` to `base_ring(M)`, return the 
`base_ring(F)`-homomorphism `F` $\to$ `M` which sends the `i`-th
basis vector of `F` to the `i`-th entry of `V`, and the scalars in 
`base_ring(F)` to their images under `h`.

    hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}, h::RingMapType) where {T, RingMapType}

Given a matrix `A` over `base_ring(M)` with `rank(F)` rows and `ngens(M)` columns
and a ring map `h` from `base_ring(F)` to `base_ring(M)`, return the
`base_ring(F)`-homomorphism `F` $\to$ `M` which sends the `i`-th basis vector of `F` to 
the linear combination $\sum_j A[i,j]*M[j]$ of the generators `M[j]` of `M`, and the 
scalars in `base_ring(F)` to their images under `h`.

!!! note
    The module `M` may be of type `FreeMod` or `SubquoMod`. If both modules
    `F` and `M` are graded, the data must define a graded module homomorphism of some degree.
    If this degree is the zero element of the (common) grading group, we refer to
    the homomorphism under consideration as a *homogeneous module homomorphism*.
"""
hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, h::RingMapType; check::Bool=true) where {T, RingMapType} = FreeModuleHom(F, M, V, h)
hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}, h::RingMapType; check::Bool=true) where {T, RingMapType} = FreeModuleHom(F, M, A, h)

@doc raw"""
    identity_map(M::ModuleFP)

Return the identity map $id_M$.
"""
function identity_map(M::ModuleFP)
  phi = hom(M, M, gens(M), check=false)
  phi.generators_map_to_generators = true
  return phi
end

### type getters in accordance with the `hom`-constructors
function morphism_type(F::AbstractFreeMod, G::ModuleFP)
  base_ring(F) === base_ring(G) && return FreeModuleHom{typeof(F), typeof(G), Nothing}
  return FreeModuleHom{typeof(F), typeof(G), typeof(base_ring(G))}
end

### Careful here! Different base rings may still have the same type.
# Whenever this is the case despite a non-trivial ring map, the appropriate 
# type getter has to be called manually!
function morphism_type(::Type{T}, ::Type{U}) where {T<:AbstractFreeMod, U<:ModuleFP}
  base_ring_type(T) == base_ring_type(U) || return morphism_type(T, U, base_ring_type(U))
  return FreeModuleHom{T, U, Nothing}
end

base_ring_type(::Type{ModuleType}) where {T, ModuleType<:ModuleFP{T}} = parent_type(T)
base_ring_elem_type(::Type{ModuleType}) where {T, ModuleType<:ModuleFP{T}} = T

base_ring_type(M::ModuleType) where {ModuleType<:ModuleFP} = base_ring_type(typeof(M))
base_ring_elem_type(M::ModuleType) where {ModuleType<:ModuleFP} = base_ring_elem_type(typeof(M))

function morphism_type(F::AbstractFreeMod, G::ModuleFP, h::RingMapType) where {RingMapType}
  return FreeModuleHom{typeof(F), typeof(G), typeof(h)}
end

function morphism_type(
    ::Type{DomainType}, ::Type{CodomainType}, ::Type{RingMapType}
  ) where {DomainType<:AbstractFreeMod, CodomainType<:ModuleFP, RingMapType}
  return FreeModuleHom{DomainType, CodomainType, RingMapType}
end

function Base.show(io::IO, ::MIME"text/plain", fmh::FreeModuleHom{T1, T2, RingMapType}) where {T1 <: AbstractFreeMod, T2 <: ModuleFP, RingMapType}
  # HACK
  show(io, fmh)
end

function Base.show(io::IO, fmh::FreeModuleHom{T1, T2, RingMapType}) where {T1 <: AbstractFreeMod, T2 <: ModuleFP, RingMapType}
  compact = get(io, :compact, false)
  io_compact = IOContext(io, :compact => true)
  if is_graded(fmh)  
    print(io_compact, domain(fmh))
    print(io, " -> ")
    print(io_compact, codomain(fmh))
    if !compact
      print(io, "\n")
      for i in 1:ngens(domain(fmh))
        print(io, domain(fmh)[i], " -> ")
        print(io_compact, fmh(domain(fmh)[i]))
        print(io, "\n")
      end
      A = grading_group(fmh)
      if degree(fmh) == A[0]
        print(io, "Homogeneous module homomorphism")
      else
        print(io_compact, "Graded module homomorphism of degree ", degree(fmh))
        print(io, "\n")
      end
    end
  else
    println(io, "Map with following data")
    println(io, "Domain:")
    println(io, "=======")
    println(io, domain(fmh))
    println(io, "Codomain:")
    println(io, "=========")
    print(io, codomain(fmh))
  end
end

@doc raw"""
    hom(F::FreeMod, G::FreeMod)

Return a free module $S$ such that $\text{Hom}(F,G) \cong S$ along with a function 
that converts elements from $S$ into morphisms $F \to G$.

# Examples
```jldoctest
julia> R, _ = polynomial_ring(QQ, ["x", "y", "z"]);

julia> F1 = free_module(R, 3)
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ

julia> F2 = free_module(R, 2)
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> V, f = hom(F1, F2)
(hom of (F1, F2), Map: hom of (F1, F2) -> set of all homomorphisms from Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ to Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ)

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
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F1 = graded_free_module(Rg, [1,2,2])
Graded free module Rg^1([-1]) + Rg^2([-2]) of rank 3 over Rg

julia> F2 = graded_free_module(Rg, [3,5])
Graded free module Rg^1([-3]) + Rg^1([-5]) of rank 2 over Rg

julia> V, f = hom(F1, F2)
(hom of (F1, F2), Map: hom of (F1, F2) -> set of all homomorphisms from Graded free module Rg^1([-1]) + Rg^2([-2]) of rank 3 over Rg to Graded free module Rg^1([-3]) + Rg^1([-5]) of rank 2 over Rg)

julia> f(V[1])
F1 -> F2
e[1] -> e[1]
e[2] -> 0
e[3] -> 0
Graded module homomorphism of degree [2]

```
"""
function hom(F::FreeMod, G::FreeMod)
  @assert base_ring(F) === base_ring(G)
  ###@assert is_graded(F) == is_graded(G)
  if is_graded(F)
    d = [y - x for x in degrees(F) for y in degrees(G)]
    GH = graded_free_module(F.R, d)
  else
    GH = FreeMod(F.R, rank(F) * rank(G))
  end
  GH.S = [Symbol("($i -> $j)") for i = F.S for j = G.S]

  #list is g1 - f1, g2-f1, g3-f1, ...
  X = Hecke.MapParent(F, G, "homomorphisms")
  n = ngens(F)
  m = ngens(G)
  R = base_ring(F)
  function im(x::FreeModElem)
    return hom(F, G, Vector{elem_type(G)}([FreeModElem(x.coords[R, (i-1)*m+1:i*m], G) for i=1:n]), check=false)
  end
  function pre(h::FreeModuleHom)
    s = sparse_row(F.R)
    o = 0
    for i=1:n
      for (p,v) = h(gen(F, i)).coords
        push!(s.pos, o+p)
        push!(s.values, v)
      end
      o += m
    end
    return FreeModElem(s, GH)
  end
  to_hom_map = MapFromFunc(GH, X, im, pre)
  set_attribute!(GH, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map)
  return GH, to_hom_map
end

@doc raw"""
    kernel(a::FreeModuleHom)

Return the kernel of `a` as an object of type `SubquoModule`.

Additionally, if `K` denotes this object, return the inclusion map `K` $\to$ `domain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]];

julia> a = hom(F, G, V);

julia> kernel(a)
(Submodule with 1 generator
1 -> x*z*e[1] - y*z*e[2] + y^2*e[3]
represented as subquotient with no relations., Map with following data
Domain:
=======
Submodule with 1 generator
1 -> x*z*e[1] - y*z*e[2] + y^2*e[3]
represented as subquotient with no relations.
Codomain:
=========
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ)
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, 3);

julia> G = graded_free_module(Rg, 2);

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]];

julia> a = hom(F, G, V);

julia> kernel(a)
(Graded submodule of F
1 -> x*z*e[1] - y*z*e[2] + y^2*e[3]
represented as subquotient with no relations, Graded submodule of F
1 -> x*z*e[1] - y*z*e[2] + y^2*e[3]
represented as subquotient with no relations -> F
x*z*e[1] - y*z*e[2] + y^2*e[3] -> x*z*e[1] - y*z*e[2] + y^2*e[3]
Homogeneous module homomorphism)

```
"""
function kernel(h::FreeModuleHom)  #ONLY for free modules...
  G = domain(h)
  R = base_ring(G)
  if ngens(G) == 0
    s = sub(G, gens(G), :module)
    help = hom(s, G, gens(G), check=false)
    help.generators_map_to_generators = true
    return s, help
  end
  g = map(h, basis(G))
  if isa(codomain(h), SubquoModule)
    g = [repres(x) for x = g]
    if isdefined(codomain(h), :quo)
      append!(g, collect(codomain(h).quo.gens))
    end
  end
  #TODO allow sub-quo here as well
  ambient_free_module_codomain = ambient_free_module(codomain(h))
  b = ModuleGens(g, ambient_free_module_codomain, default_ordering(ambient_free_module_codomain))
  k = syzygy_module(b)
  if isa(codomain(h), SubquoModule)
    s = collect(k.sub.gens)
    k = sub(G, [FreeModElem(x.coords[R,1:dim(G)], G) for x = s], :module)
  else
    #the syzygie_module creates a new free module to work in
    k = sub(G, [FreeModElem(x.coords, G) for x = collect(k.sub.gens)], :module)
  end
  @assert k.F === G
  c = collect(k.sub.gens)
  return k, hom(k, parent(c[1]), c, check=false)
end

@doc raw"""
    image(a::FreeModuleHom)

Return the image of `a` as an object of type `SubquoModule`.

Additionally, if `I` denotes this object, return the inclusion map `I` $\to$ `codomain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate polynomial ring in 3 variables over QQ

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]];

julia> a = hom(F, G, V);

julia> image(a)
(Submodule with 3 generators
1 -> y*e[1]
2 -> x*e[1] + y*e[2]
3 -> z*e[2]
represented as subquotient with no relations., Map with following data
Domain:
=======
Submodule with 3 generators
1 -> y*e[1]
2 -> x*e[1] + y*e[2]
3 -> z*e[2]
represented as subquotient with no relations.
Codomain:
=========
Free module of rank 2 over Multivariate polynomial ring in 3 variables over QQ)
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> F = graded_free_module(Rg, 3);

julia> G = graded_free_module(Rg, 2);

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]];

julia> a = hom(F, G, V);

julia> image(a)
(Graded submodule of G
1 -> y*e[1]
2 -> x*e[1] + y*e[2]
3 -> z*e[2]
represented as subquotient with no relations, Graded submodule of G
1 -> y*e[1]
2 -> x*e[1] + y*e[2]
3 -> z*e[2]
represented as subquotient with no relations -> G
y*e[1] -> y*e[1]
x*e[1] + y*e[2] -> x*e[1] + y*e[2]
z*e[2] -> z*e[2]
Homogeneous module homomorphism)

```
"""
function image(h::FreeModuleHom)
  si = filter(!iszero, images_of_generators(h))
  s = sub(codomain(h), si, :module)
  phi = hom(s, codomain(h), si, check=false)
  return s, phi
end


