###############################################################################
# SubQuoHom constructors
###############################################################################

@doc raw"""
    SubQuoHom(D::SubquoModule, C::ModuleFP{T}, im::Vector{<:ModuleFPElem{T}}) where T

Return the morphism $D \to C$ for a subquotient $D$ where `D[i]` is mapped to `im[i]`.
In particular, `length(im) == ngens(D)` must hold.
"""
SubQuoHom(D::SubquoModule, C::ModuleFP{T}, im::Vector{<:ModuleFPElem{T}}; check::Bool=true) where {T} = SubQuoHom{typeof(D), typeof(C), Nothing}(D, C, im; check)
SubQuoHom(D::SubquoModule, C::ModuleFP{T}, im::Vector{<:ModuleFPElem{T}}, h::RingMapType; check::Bool=true) where {T, RingMapType} = SubQuoHom{typeof(D), typeof(C), RingMapType}(D, C, im, h; check)

@doc raw"""
    SubQuoHom(D::SubquoModule, C::ModuleFP{T}, mat::MatElem{T})

Return the morphism $D \to C$ corresponding to the given matrix, where $D$ is a subquotient.
`mat` must have `ngens(D)` many rows and `ngens(C)` many columns.
"""
function SubQuoHom(D::SubquoModule, C::ModuleFP{T}, mat::MatElem{T}; check::Bool=true) where T
  @assert nrows(mat) == ngens(D)
  @assert ncols(mat) == ngens(C)
  if C isa FreeMod
    hom = SubQuoHom(D, C, [FreeModElem(sparse_row(mat[i:i,:]), C) for i=1:ngens(D)]; check)
    return hom
  else
    hom = SubQuoHom(D, C, [SubquoModuleElem(sparse_row(mat[i:i,:]), C) for i=1:ngens(D)]; check)
    return hom
  end
end

function SubQuoHom(D::SubquoModule, C::ModuleFP{T}, mat::MatElem{T}, h::RingMapType; check::Bool=true) where {T, RingMapType}
  @assert nrows(mat) == ngens(D)
  @assert ncols(mat) == ngens(C)
  if C isa FreeMod
    hom = SubQuoHom(D, C, [FreeModElem(sparse_row(mat[i:i,:]), C) for i=1:ngens(D)], h; check)
    return hom
  else
    hom = SubQuoHom(D, C, [SubquoModuleElem(sparse_row(mat[i:i,:]), C) for i=1:ngens(D)], h; check)
    return hom
  end
end

function Base.show(io::IO, ::MIME"text/plain", fmh::SubQuoHom{T1, T2, RingMapType}) where {T1 <: AbstractSubQuo, T2 <: ModuleFP, RingMapType}
   println(terse(io), fmh)
   io = pretty(io)
   io_compact = IOContext(io, :compact => true)
   println(io_compact, Indent(), "from ", Lowercase(), domain(fmh))
   print(io_compact, "to ", Lowercase(), codomain(fmh), Dedent())

  if is_graded(fmh)
    println(io)
    print(io, "defined by", Indent())
    io = terse(io)
    domain_gens = gens(domain(fmh))
    for g in domain_gens
      println(io)
      print(io, g, " -> ", fmh(g))
    end
  end
end

function Base.show(io::IO, fmh::SubQuoHom{T1, T2, RingMapType}) where {T1 <: AbstractSubQuo, T2 <: ModuleFP, RingMapType}
  if is_terse(io)
    if is_graded(fmh)
      A = grading_group(fmh)
      if degree(fmh) == A[0]
        print(io, "Homogeneous module homomorphism")
      else
        print(io, "Graded module homomorphism of degree ", degree(fmh))
      end
    else
        print(io, "Module homomorphism")
    end
  else
    io = pretty(io)
    print(io, "Hom: ")
    io = terse(io)
    print(io, Lowercase(), domain(fmh), " -> ")
    print(io, Lowercase(), codomain(fmh))
  end
end

images_of_generators(phi::SubQuoHom) = phi.im::Vector{elem_type(codomain(phi))}
image_of_generator(phi::SubQuoHom, i::Int) = phi.im[i]::elem_type(codomain(phi))

###################################################################

@doc raw"""
    hom(M::SubquoModule{T}, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T

Given a vector `V` of `ngens(M)` elements of `N`,
return the homomorphism `M` $\to$ `N` which sends the `i`-th
generator `M[i]` of `M` to the `i`-th entry of `V`.

    hom(M::SubquoModule{T}, N::ModuleFP{T},  A::MatElem{T})) where T

Given a matrix `A` with `ngens(M)` rows and `ngens(N)` columns, return the
homomorphism `M` $\to$ `N` which sends the `i`-th generator `M[i]` of `M` to
the linear combination $\sum_j A[i,j]*N[j]$ of the generators `N[j]` of `N`.

!!! note
    The module `N` may be of type `FreeMod` or `SubquoMod`. If both modules
    `M` and `N` are graded, the data must define a graded module homomorphism of some degree.
    If this degree is the zero element of the (common) grading group, we refer to
    the homomorphism under consideration as a *homogeneous module homomorphism*.

!!! warning
    The functions do not check whether the resulting homomorphism is well-defined,
    that is, whether it sends the relations of `M` into the relations of `N`.

If you are uncertain with regard to well-definedness, use the function below.
Note, however, that the check performed by the function requires a Gröbner basis computation. This may take some time.

    is_welldefined(a::ModuleFPHom)

Return `true` if `a` is well-defined, and `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubquoModule(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]]
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 x*y^2*e[1]
 x*y*e[1]

julia> a = hom(M, N, V)
Module homomorphism
  from M
  to M

julia> is_welldefined(a)
true

julia> W = R[y^2 0; 0 x]
[y^2   0]
[  0   x]

julia> b = hom(M, N, W);

julia> a == b
true
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x; y];

julia> B = R[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> N = M;

julia> W = [y*N[1], x*N[2]]
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 x*y*e[1]
 x*y*e[1]

julia> c = hom(M, N, W);

julia> is_welldefined(c)
false
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B)
Graded subquotient of graded submodule of F with 2 generators
  1: x*e[1]
  2: y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

julia> is_welldefined(a)
true

julia> W = Rg[y^2 0; 0 x^2]
[y^2     0]
[  0   x^2]

julia> b = hom(M, N, W)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

julia> a == b
true

julia> W = [y*N[1], x*N[2]]
2-element Vector{SubquoModuleElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 x*y*e[1]
 x*y*e[1]

julia> c = hom(M, N, W)
Graded module homomorphism of degree [1]
  from M
  to M
defined by
  x*e[1] -> x*y*e[1]
  y*e[1] -> x*y*e[1]

julia> is_welldefined(c)
false
```
"""
hom(M::SubquoModule{T}, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}; check::Bool=true) where T = SubQuoHom(M, N, V; check)
hom(M::SubquoModule{T}, N::ModuleFP{T},  A::MatElem{T}; check::Bool=true) where T = SubQuoHom(M, N, A; check)


@doc raw"""
    hom(M::SubquoModule, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, h::RingMapType) where {T, RingMapType}

Given a vector `V` of `ngens(M)` elements of `N`,
return the homomorphism `M` $\to$ `N` which sends the `i`-th
generator `M[i]` of `M` to the `i`-th entry of `V`, and the
scalars in `base_ring(M)` to their images under `h`.

    hom(M::SubquoModule, N::ModuleFP{T},  A::MatElem{T}, h::RingMapType) where {T, RingMapType}

Given a matrix `A` with `ngens(M)` rows and `ngens(N)` columns, return the
homomorphism `M` $\to$ `N` which sends the `i`-th generator `M[i]` of `M` to
the linear combination $\sum_j A[i,j]*N[j]$ of the generators `N[j]` of `N`,
and the scalars in `base_ring(M)` to their images under `h`.

!!! note
    The module `N` may be of type `FreeMod` or `SubquoMod`. If both modules
    `M` and `N` are graded, the data must define a graded module homomorphism of some degree.
    If this degree is the zero element of the (common) grading group, we refer to
    the homomorphism under consideration as a *homogeneous module homomorphism*.

!!! warning
    The functions do not check whether the resulting homomorphism is well-defined,
    that is, whether it sends the relations of `M` into the relations of `N`.

If you are uncertain with regard to well-definedness, use the function below.
Note, however, that the check performed by the function requires a Gröbner basis computation. This may take some time.

    is_welldefined(a::ModuleFPHom)

Return `true` if `a` is well-defined, and `false` otherwise.

"""
hom(M::SubquoModule, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, h::RingMapType; check::Bool=true) where {T, RingMapType} = SubQuoHom(M, N, V, h; check)
hom(M::SubquoModule, N::ModuleFP{T}, A::MatElem{T}, h::RingMapType; check::Bool=true) where {T, RingMapType} = SubQuoHom(M, N, A, h; check)

function is_welldefined(H::Union{FreeModuleHom,FreeModuleHom_dec})
  return true
end

function is_welldefined(H::SubQuoHom)
  M = domain(H)
  pres = presentation(M)
  # is a short exact sequence with maps
  # M <--eps-- F0 <--g-- F1
  # and H : M -> N
  eps = map(pres, 0)
  g = map(pres, 1)
  F0 = pres[0]
  N = codomain(H)
  # the induced map phi : F0 --> N
  phi = hom(F0, N, elem_type(N)[H(eps(v)) for v in gens(F0)]; check=false)
  # now phi ∘ g : F1 --> N has to be zero.
  return iszero(compose(g, phi))

  C = present_as_cokernel(M).quo
  n = ngens(C)
  m = rank(C.F)
  ImH = images_of_generators(H)
  for i=1:n
    if !iszero(sum([C[i][j]*ImH[j] for j=1:m]; init=zero(codomain(H))))
      return false
    end
  end
  return true
end

# No ring map
function (==)(f::ModuleFPHom{D, C, Nothing}, g::ModuleFPHom{D, C, Nothing}) where {D, C}
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return all(f(v) == g(v) for v in gens(domain(f)))
end

# With ring map
function (==)(f::ModuleFPHom, g::ModuleFPHom)
  if isnothing(ring_map(f))
    isnothing(ring_map(g)) || is_trivial(ring_map(g)) || return false
  end
  if isnothing(ring_map(g))
    isnothing(ring_map(f)) || is_trivial(ring_map(f)) || return false
  end
  # TODO: Catch the case where == defaults to === and this returns false, 
  # but the ring_maps are identical in the mathematical sense nevertheless.
  ring_map(f) == ring_map(g) || return false # Will throw if ring maps do not compare
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  return all(f(v) == g(v) for v in gens(domain(f)))
end

# TODO: Move to Hecke?
function (==)(f::Map, g::Map)
  f === g && return true
  error("comparison of maps of type $(typeof(f)) and $(typeof(g)) not implemented")
end

function Base.hash(f::ModuleFPHom{T}, h::UInt) where {U<:FieldElem, S<:MPolyRingElem{U}, T<:ModuleFP{S}}
  b = 0x535bbdbb2bc54b46 % UInt
  h = hash(typeof(f), h)
  h = hash(domain(f), h)
  h = hash(codomain(f), h)
  for g in images_of_generators(f)
    h = hash(g, h)
  end
  return xor(h, b)
end

function Base.hash(f::ModuleFPHom, h::UInt)
  b = 0x535bbdbb2bc54b46 % UInt
  h = hash(typeof(f), h)
  h = hash(domain(f), h)
  h = hash(codomain(f), h)
  # We can not assume that the images of generators
  # have a hash in general
  return xor(h, b)
end
###################################################################

@doc raw"""
    matrix(a::SubQuoHom)

Given a homomorphism `a` of type  `SubQuoHom` with domain `M`
and codomain `N`, return a matrix `A` with `ngens(M)` rows and
`ngens(N)` columns such that `a == hom(M, N, A)`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over R

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubquoModule(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]];

julia> a = hom(M, N, V);

julia> A = matrix(a)
[y^2   0]
[  0   x]

julia> a(M[1])
x*y^2*e[1]
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B)
Graded subquotient of graded submodule of F with 2 generators
  1: x*e[1]
  2: y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

julia> matrix(a)
[y^2     0]
[  0   x^2]

```
"""
function matrix(f::SubQuoHom)
  if !isdefined(f, :matrix)
    D = domain(f)
    C = codomain(f)
    R = base_ring(C)
    matrix = zero_matrix(R, ngens(D), ngens(C))
    for i=1:ngens(D), j=1:ngens(C)
      matrix[i,j] = f.im[i][j]
    end
    f.matrix = matrix
  end
  return f.matrix
end

function show_morphism(f::ModuleFPHom)
  show(stdout, "text/plain", matrix(f))
end


@doc raw"""
    image(a::SubQuoHom, m::SubquoModuleElem)

Return the image $a(m)$.
"""
function image(f::SubQuoHom, a::SubquoModuleElem)
  @assert a.parent === domain(f)
  iszero(a) && return zero(codomain(f))
  # The code in the comment below was an attempt to make
  # evaluation of maps faster. However, it turned out that
  # for the average use case the comparison was more expensive
  # than the gain for mappings. The flag should be set by constructors
  # nevertheless when applicable.
  #if f.generators_map_to_generators === nothing
  #  f.generators_map_to_generators = images_of_generators(f) == gens(codomain(f))
  #end
  f.generators_map_to_generators === true && return codomain(f)(map_entries(base_ring_map(f), coordinates(a)))
  h = base_ring_map(f)
  return sum(h(b)*image_of_generator(f, i) for (i, b) in coordinates(a); init=zero(codomain(f)))
end

function image(f::SubQuoHom{<:SubquoModule, <:ModuleFP, Nothing}, a::SubquoModuleElem)
  # TODO matrix vector multiplication
  @assert a.parent === domain(f)
 #if f.generators_map_to_generators === nothing
 #  f.generators_map_to_generators = images_of_generators(f) == gens(codomain(f))
 #end
  f.generators_map_to_generators === true && return codomain(f)(coordinates(a))
  return sum(c*image_of_generator(f, i) for (i, c) in coordinates(a); init=zero(codomain(f)))
end

@doc raw"""
    image(f::SubQuoHom, a::FreeModElem)

Return $f(a)$. `a` must represent an element in the domain of `f`.
"""
function image(f::SubQuoHom, a::FreeModElem)
  return image(f, SubquoModuleElem(a, domain(f)))
end

function image(f::SubQuoHom{<:SubquoModule, <:ModuleFP, Nothing}, a::FreeModElem)
  return image(f, SubquoModuleElem(a, domain(f)))
end

@doc raw"""
    preimage(f::SubQuoHom, a::Union{SubquoModuleElem,FreeModElem})

Compute a preimage of `a` under `f`.
"""
function preimage(f::SubQuoHom{<:SubquoModule, <:ModuleFP}, a::Union{SubquoModuleElem,FreeModElem})
  @assert parent(a) === codomain(f)
  phi = base_ring_map(f)
  D = domain(f)
  i = zero(D)
  b = coordinates(a isa FreeModElem ? a : repres(a), image(f)[1])
  bb = map_entries(x->(preimage(phi, x)), b)
  for (p,v) = bb
    i += v*gen(D, p)
  end
  return i
end

function preimage(f::SubQuoHom{<:SubquoModule, <:ModuleFP, Nothing},
        a::Union{SubquoModuleElem,FreeModElem})
  @assert parent(a) === codomain(f)
  D = domain(f)
  i = zero(D)
  b = coordinates(a isa FreeModElem ? a : repres(a), image(f)[1])
  for (p,v) = b
    i += v*gen(D, p)
  end
  return i
end

(f::SubQuoHom)(a::FreeModElem) = image(f, SubquoModuleElem(a, domain(f)))
(f::SubQuoHom)(a::SubquoModuleElem) = image(f, a)

@doc raw"""
    image(a::SubQuoHom)

Return the image of `a` as an object of type `SubquoModule`.

Additionally, if `I` denotes this object, return the inclusion map `I` $\to$ `codomain(a)`.
"""
@attr Tuple{<:SubquoModule, <:ModuleFPHom} function image(h::SubQuoHom)
  s = sub_object(codomain(h), images_of_generators(h))
  inc = hom(s, codomain(h), images_of_generators(h), check=false)
  return s, inc
end

@doc raw"""
    image(a::ModuleFPHom)

Return the image of `a` as an object of type `SubquoModule`.

Additionally, if `I` denotes this object, return the inclusion map `I` $\to$ `codomain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> F = free_module(R, 3);

julia> G = free_module(R, 2);

julia> W = R[y 0; x y; 0 z]
[y   0]
[x   y]
[0   z]

julia> a = hom(F, G, W);

julia> I, incl = image(a);

julia> I
Submodule with 3 generators
  1: y*e[1]
  2: x*e[1] + y*e[2]
  3: z*e[2]
represented as subquotient with no relations

julia> incl
Module homomorphism
  from I
  to G
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> F = free_module(R, 1);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubquoModule(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]]
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 x*y^2*e[1]
 x*y*e[1]

julia> a = hom(M, N, V);

julia> I, incl = image(a);

julia> I
Subquotient of submodule with 2 generators
  1: x*y^2*e[1]
  2: x*y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> incl
Module homomorphism
  from I
  to M
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B)
Graded subquotient of graded submodule of F with 2 generators
  1: x*e[1]
  2: y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

julia> image(a)
(Graded subquotient of graded submodule of F with 2 generators
  1: x*y^2*e[1]
  2: x^2*y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1], Hom: graded subquotient of graded submodule of F with 2 generators
  1: x*y^2*e[1]
  2: x^2*y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1] -> M)
```
"""
function image(a::ModuleFPHom)
 error("image is not implemented for the given types.")
end

@doc raw"""
    kernel(a::SubQuoHom)

Return the kernel of `a` as an object of type `SubquoModule`.

Additionally, if `K` denotes this object, return the inclusion map `K` $\to$ `domain(a)`.
"""
function kernel(h::SubQuoHom)
  D = domain(h)
  R = base_ring(D)
  is_graded(h) ? F = graded_free_module(R, degrees_of_generators(D)) : F = FreeMod(R, ngens(D))
  hh = hom(F, codomain(h), images_of_generators(h), check=false)
  K, inc_K = kernel(hh)
  @assert domain(inc_K) === K
  @assert codomain(inc_K) === F
  v = gens(D)
  imgs = Vector{elem_type(D)}(filter(!iszero, [sum(a*v[i] for (i, a) in coordinates(g); init=zero(D)) for g in images_of_generators(inc_K)]))
  k = sub_object(D, imgs)
  return k, hom(k, D, imgs, check=false)
end

@doc raw"""
    kernel(a::ModuleFPHom)

Return the kernel of `a` as an object of type `SubquoModule`.

Additionally, if `K` denotes this object, return the inclusion map `K` $\to$ `domain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> F = free_module(R, 3);

julia> G = free_module(R, 2);

julia> W = R[y 0; x y; 0 z]
[y   0]
[x   y]
[0   z]

julia> a = hom(F, G, W);

julia> K, incl = kernel(a);

julia> K
Submodule with 1 generator
  1: x*z*e[1] - y*z*e[2] + y^2*e[3]
represented as subquotient with no relations

julia> incl
Module homomorphism
  from K
  to F
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> F = free_module(R, 1);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubquoModule(F, A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]]
2-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 x*y^2*e[1]
 x*y*e[1]

julia> a = hom(M, N, V);

julia> K, incl = kernel(a);

julia> K
Subquotient of submodule with 3 generators
  1: (-x + y^2)*e[1]
  2: x*y*e[1]
  3: -x*y*e[1]
by submodule with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> incl
Module homomorphism
  from K
  to M
```

```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B)
Graded subquotient of graded submodule of F with 2 generators
  1: x*e[1]
  2: y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x^2*N[2]];

julia> a = hom(M, N, V)
Graded module homomorphism of degree [2]
  from M
  to M
defined by
  x*e[1] -> x*y^2*e[1]
  y*e[1] -> x^2*y*e[1]

julia> kernel(a)
(Graded subquotient of graded submodule of F with 2 generators
  1: y*e[1]
  2: -x*y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1], Hom: graded subquotient of graded submodule of F with 2 generators
  1: y*e[1]
  2: -x*y*e[1]
by graded submodule of F with 3 generators
  1: x^2*e[1]
  2: y^3*e[1]
  3: z^4*e[1] -> M)

```
"""
function kernel(a::ModuleFPHom)
 error("kernel is not implemented for the given types.")
end

#TODO
#  replace the +/- for the homs by proper constructors for homs and direct sums
#  relshp to store the maps elsewhere

@doc raw"""
    *(a::ModuleFPHom, b::ModuleFPHom)

Return the composition `b` $\circ$ `a`.
"""
function *(h::ModuleFPHom{T1, T2, Nothing}, g::ModuleFPHom{T2, T3, Nothing}) where {T1, T2, T3}
  @assert codomain(h) === domain(g)
  return hom(domain(h), codomain(g), Vector{elem_type(codomain(g))}([g(h(x)) for x = gens(domain(h))]), check=false)
end

function *(h::ModuleFPHom{T1, T2, <:Map}, g::ModuleFPHom{T2, T3, <:Map}) where {T1, T2, T3}
  @assert codomain(h) === domain(g)
  return hom(domain(h), codomain(g), Vector{elem_type(codomain(g))}([g(h(x)) for x = gens(domain(h))]), compose(base_ring_map(h), base_ring_map(g)), check=false)
end

function *(h::ModuleFPHom{T1, T2, <:Any}, g::ModuleFPHom{T2, T3, <:Any}) where {T1, T2, T3}
  @assert codomain(h) === domain(g)
  return hom(domain(h), codomain(g),
             Vector{elem_type(codomain(g))}([g(h(x)) for x = gens(domain(h))]),
             MapFromFunc(base_ring(domain(h)),
                         base_ring(codomain(g)),
                         x->(base_ring_map(g)(base_ring_map(h)(x)))),
             check=false
            )

end

compose(h::ModuleFPHom, g::ModuleFPHom) = h*g

-(h::ModuleFPHom{D, C, Nothing}) where {D, C} = hom(domain(h), codomain(h), elem_type(codomain(h))[-h(x) for x in gens(domain(h))], check=false)
-(h::ModuleFPHom{D, C, T}) where {D, C, T} = hom(domain(h), codomain(h), elem_type(codomain(h))[-h(x) for x in gens(domain(h))], base_ring_map(h), check=false)

function -(h::ModuleFPHom{D, C, T}, g::ModuleFPHom{D, C, T}) where {D, C, T}
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  @assert base_ring_map(h) === base_ring_map(g)
  return hom(domain(h), codomain(h), elem_type(codomain(h))[h(x) - g(x) for x in gens(domain(h))], base_ring_map(h), check=false)
end

function -(h::ModuleFPHom{D, C, Nothing}, g::ModuleFPHom{D, C, Nothing}) where {D, C}
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  return hom(domain(h), codomain(h), elem_type(codomain(h))[h(x) - g(x) for x in gens(domain(h))], check=false)
end

function +(h::ModuleFPHom{D, C, T}, g::ModuleFPHom{D, C, T}) where {D, C, T}
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  @assert base_ring_map(h) === base_ring_map(g)
  return hom(domain(h), codomain(h), elem_type(codomain(h))[h(x) + g(x) for x in gens(domain(h))], base_ring_map(h), check=false)
end

function +(h::ModuleFPHom{D, C, Nothing}, g::ModuleFPHom{D, C, Nothing}) where {D, C}
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  return hom(domain(h), codomain(h), elem_type(codomain(h))[h(x) + g(x) for x in gens(domain(h))], check=false)
end

function *(a::AdmissibleModuleFPRingElem, g::ModuleFPHom{D, C, Nothing}) where {D, C}
  @assert base_ring(codomain(g)) === parent(a)
  return hom(domain(g), codomain(g), elem_type(codomain(g))[a*g(x) for x in gens(domain(g))], check=false)
end

function *(a::AdmissibleModuleFPRingElem, g::ModuleFPHom{D, C, T}) where {D, C, T}
  @assert base_ring(codomain(g)) === parent(a)
  return hom(domain(g), codomain(g), elem_type(codomain(g))[a*g(x) for x in gens(domain(g))], base_ring_map(g), check=false)
end


@doc raw"""
    restrict_codomain(H::ModuleFPHom, M::SubquoModule)

Return, if possible, a homomorphism, which is mathematically identical to `H`,
but has codomain `M`. `M` has to be a submodule of the codomain of `H`.
"""
function restrict_codomain(H::ModuleFPHom, M::SubquoModule)
  D = domain(H)
  return hom(D, M, map(v -> SubquoModuleElem(v, M), map(x -> repres(H(x)), gens(D))), check=false)
end

@doc raw"""
    restrict_domain(H::SubQuoHom, M::SubquoModule)

Restrict the morphism `H` to `M`. For this `M` has to be a submodule
of the domain of `H`. The relations of `M` must be the relations of
the domain of `H`.
"""
function restrict_domain(H::SubQuoHom, M::SubquoModule)
  for (cod, t) in M.outgoing
    if cod === domain(H)
      return _recreate_morphism(M, cod, t)*H
    end
  end
  # else there is no cached map
  if ngens(M) > 0
    @assert M.quo == domain(H).quo
  end
  _, i = sub(domain(H), map(m -> SubquoModuleElem(repres(m), domain(H)), gens(M)), cache_morphism=true)
  return i*H
end

@doc raw"""
    induced_map(f::FreeModuleHom, M::SubquoModule, check::Bool = true)

Return the map which sends an element `v` of `M` to `f(repres(v))`.
If `check` is set to true the well-definedness of the map is checked.
"""
function induced_map(f::FreeModuleHom, M::SubquoModule, check::Bool = true)
  @assert ambient_free_module(M) === domain(f)
  ind_f = hom(M, codomain(f), [f(repres(v)) for v in gens(M)], check=false)
  if check
    @assert is_welldefined(ind_f)
  end
  return ind_f
end

@doc raw"""
    inv(a::ModuleFPHom)

If `a` is bijective, return its inverse.
"""
function inv(H::ModuleFPHom)
  if isdefined(H, :inverse_isomorphism)
    return H.inverse_isomorphism
  end
  @assert is_bijective(H)
  N = domain(H)
  M = codomain(H)

  Hinv = hom(M,N, Vector{elem_type(N)}([preimage(H,m) for m in gens(M)]), check=false)
  Hinv.inverse_isomorphism = H
  H.inverse_isomorphism = Hinv

  return Hinv
end

######################################
# Migrating test
######################################
@doc raw"""
    projection(F::FreeMod, indices::AbstractArray)

Return the canonical projection from $F = R^I$ to $R^(\texttt{indices})$ where $\texttt{indices} \subset I$.
"""
function projection(F::FreeMod, indices::AbstractArray)
  @assert all(<=(ngens(F)), indices)
  @assert length(Set(indices)) == length(indices) # unique indices
  R = base_ring(F)
  G = FreeMod(R, length(indices))
  return hom(F, G, Vector{elem_type(G)}([i in indices ? G[findfirst(==(i),indices)] : zero(G) for i=1:ngens(F)]), check=false)
end

@doc raw"""
    preimage(H::SubQuoHom,N::SubquoModule{T}, task::Symbol = :none) where {T}

Return the preimage of the submodule `N` under the morphism `H`
as a subquotient, as well as the injection homomorphism into the domain of $H$.
"""
function preimage(H::SubQuoHom,N::SubquoModule{T}, task::Symbol = :none) where {T}
  inclusion = get_attribute(N, :canonical_inclusion)
  if inclusion !== nothing && codomain(inclusion) === codomain(H)
    elems = [inclusion(v) for v in gens(N)]
  else
    elems = [SubquoModuleElem(repres(v),codomain(H)) for v in gens(N)]
  end
  return preimage(H,elems,task)
end

@doc raw"""
    preimage(H::SubQuoHom,elems::Vector{SubquoModuleElem{T}}, task::Symbol = :none) where {T}

Return the preimage of the submodule generated by the Elements `elems` under $H$
as a subquotient, as well as the injection homomorphism into the domain of $H$.
"""
function preimage(H::SubQuoHom,elems::Vector{SubquoModuleElem{T}}, task::Symbol = :none) where {T}
  if length(elems)==0
      k,emb = kernel(H)
      if task == :none
        return k
      else
        return k,emb
      end
  end
  @assert all(x->parent(x)===codomain(H),elems)
  cod_coker,i_cod_coker_inv = present_as_cokernel(codomain(H), :with_morphism)
  i_cod_coker = inv(i_cod_coker_inv) # this is cheap
  elems_in_coker = map(i_cod_coker, elems)
  cokernel_modulo_elmes,projection = quo(cod_coker,elems_in_coker)
  preimage, emb = kernel(H*i_cod_coker*projection)

  if task != :none
    return preimage, emb
  else
    return preimage
  end
end

@doc raw"""
    matrix_kernel(A::MatElem)

Compute the kernel of `A` where `A` is considered as the corresponding morphism
between free modules.
"""
function matrix_kernel(A::MatElem)
  R = base_ring(A)
  F_domain = FreeMod(R, nrows(A))
  F_codomain = FreeMod(R, ncols(A))

  phi = FreeModuleHom(F_domain, F_codomain, A)
  _, inclusion = kernel(phi)
  return matrix(inclusion)
end

@doc raw"""
    simplify_light(M::SubquoModule)

Simplify the given subquotient `M` and return the simplified subquotient `N` along
with the injection map $N \to M$ and the projection map $M \to N$. These maps are
isomorphisms.
The only simplifications which are done are the following:
- Remove all generators which are represented by the zero element in the ambient
  free module.
- Remove all generators which are in the generating set of the relations.
- Remove all duplicates in the generators and relations sets.
"""
function simplify_light(M::SubquoModule)
  M_gens = ambient_representatives_generators(M)
  M_rels = relations(M)

  N_rels = unique!(filter(!is_zero, M_rels))
  N_gens = unique!(setdiff!(filter(!is_zero, M_gens), N_rels))

  N = length(N_rels) == 0 ? SubquoModule(ambient_free_module(M), N_gens) : SubquoModule(ambient_free_module(M), N_gens, N_rels)

  index_of_N_in_M = indexin(N_gens, M_gens)
  inj = hom(N, M, Vector{elem_type(M)}([M[index_of_N_in_M[i]] for i in 1:ngens(N)]), check=false)

  index_of_M_in_N = indexin(M_gens, N_gens)
  proj = hom(M, N, Vector{elem_type(N)}([index_of_M_in_N[i] === nothing ? zero(N) : N[index_of_M_in_N[i]] for i in 1:ngens(M)]), check=false)

  return N, inj, proj
end

@doc raw"""
    simplify_with_same_ambient_free_module(M::SubquoModule)

Simplify the given subquotient `M` and return the simplified subquotient `N` along
with the injection map $N \to M$ and the projection map $M \to N$. These maps are
isomorphisms. The ambient free module of `N` is the same as that of `M`.
"""
function simplify_with_same_ambient_free_module(M::SubquoModule)
  _, to_M, from_M = _old_simplify(M)
  N, N_to_M = image(to_M)
  return N, N_to_M, hom(M, N, [N(coordinates(from_M(g))) for g in gens(M)], check=false)
  #return N, N_to_M, hom(M, N, [N(repres(g)) for g in gens(M)])
end

### Old version of simplify which is trying to do things directly on the 
# subquotient without a presentation. Does not respect gradings and is, 
# hence, not fully functional; see issue #3108.
function _old_simplify(M::SubquoModule)
  respect_grading = is_graded(M)
  function standard_unit_vector_in_relations(i::Int, M::SubquoModule)
    F = ambient_free_module(M)
    !isdefined(M, :quo) && return iszero(F[i])
    return in(F[i], M.quo)
  end

  function delete_rows(A::MatElem, to_delete::Vector{Int})
    Mat = A[setdiff(1:nrows(A),to_delete),:]
    return Mat
  end
  function delete_columns(A::MatElem, to_delete::Vector{Int})
    return transpose(delete_rows(transpose(A), to_delete))
  end

  function assign_row!(A::MatElem, v::Vector, row_index::Int)
    if length(v) != size(A)[2]
      throw(DimensionMismatch("Different row lengths"))
    end
    for i=1:length(v)
      A[row_index,i] = v[i]
    end
    return A
  end

  function assign_row!(A::MatElem, v::MatElem, row_index::Int)
    if size(v)[1] > 1
      throw(DimensionMismatch("Expected row vector"))
    end
    if length(v) != size(A)[2]
      throw(DimensionMismatch("Different row lengths"))
    end
    for i=1:length(v)
      A[row_index,i] = v[1,i]
    end
    return A
  end

  function rows_to_delete(A::MatElem, max_index::Int, M::SubquoModule, respect_grading::Bool=false)
    to_delete_indices::Vector{Int} = []
    corresponding_row_index::Vector{Int} = []
    if max_index < nrows(A)
      A = vcat(A[(max_index+1):nrows(A),:],A[1:max_index,:])
    end
    K = matrix_kernel(A)
    if max_index < nrows(A)
      K = hcat(K[:,(ncols(K)-max_index+1):ncols(K)],K[:,1:(ncols(K)-max_index)])
    end
    for i=1:size(K)[1], j=1:max_index
      #if is_unit(K[i,j]) && (!respect_grading || degree(M.sub.O[j]) == degree(M.quo.O[i]))
      if is_unit(K[i,j])
        deletion_possible = true
        for k in to_delete_indices
          if !iszero(K[i,k])
            deletion_possible = false
            break
          end
        end
        if deletion_possible
          push!(to_delete_indices, j)
          push!(corresponding_row_index, i)
        end
      end
    end
    return to_delete_indices, corresponding_row_index, K
  end

  R = base_ring(M)
  #remove columns

  M_generators = generator_matrix(M.sub)
  M_relations = isdefined(M, :quo) ? generator_matrix(M.quo) : zero_matrix(R, 1,ncols(M_generators))

  to_delete::Vector{Int} = []
  for i=1:size(M_relations)[2]
    if standard_unit_vector_in_relations(i, M)
      push!(to_delete, i)
    end
  end

  new_generators = delete_columns(M_generators, to_delete)
  new_relations = delete_columns(M_relations, to_delete)

  to_delete,_,_ = rows_to_delete(transpose(vcat(new_generators, new_relations)), size(new_relations)[2], M, respect_grading)

  new_generators = delete_columns(new_generators, to_delete)
  new_relations = delete_columns(new_relations, to_delete)

  #remove rows
  #simplify relations
  to_delete,_,_ = rows_to_delete(new_relations, size(new_relations)[1], M, respect_grading)

  new_relations = delete_rows(new_relations, to_delete)

  #simplify generators
  to_delete, corresponding_row, K_gen = rows_to_delete(vcat(new_generators, new_relations), size(new_generators)[1], M, respect_grading)

  injection_matrix = delete_rows(identity_matrix(R, size(M_generators)[1]), to_delete)
  projection_matrix = zero_matrix(R, size(M_generators)[1], size(K_gen)[2]-length(to_delete))
  for i=1:size(M_generators)[1]
    if i in to_delete
      index = findfirst(==(i), to_delete)
      assign_row!(projection_matrix, R(-1)*R(inv(coeff(K_gen[corresponding_row[index],i], 1)))*delete_columns(K_gen[corresponding_row[index]:(corresponding_row[index]),:], to_delete), i)
    else
      standard_unit_vector_index = i-length(filter(<(i), to_delete))
      standard_unit_vector = [j == standard_unit_vector_index ? R(1) : R(0) for j=1:size(projection_matrix)[2]]
      assign_row!(projection_matrix, standard_unit_vector, i)
    end
  end

  new_generators = delete_rows(new_generators, to_delete)

  if length(new_generators)==0
    zero_module = FreeMod(R,0)
    injection = FreeModuleHom(zero_module, M, Vector{elem_type(M)}())
    projection = SubQuoHom(M, zero_module, [zero(zero_module) for i=1:ngens(M)])
    # TODO early return or register morphisms?
    return zero_module,injection,projection
  else
    SQ = iszero(new_relations) ? SubquoModule(SubModuleOfFreeModule(new_generators)) : SubquoModule(new_generators, new_relations)
    injection = SubQuoHom(SQ, M, injection_matrix)
    projection = SubQuoHom(M, SQ, projection_matrix[:,1:size(projection_matrix)[2]-size(new_relations)[1]])
  end
  register_morphism!(injection)
  register_morphism!(projection)
  injection.inverse_isomorphism = projection
  projection.inverse_isomorphism = injection

  return SQ,injection,projection
end

######################################
# Matrix to morphism
######################################
@doc raw"""
    map(F::FreeMod{T}, A::MatrixElem{T}) where T

Converts a given $n \times m$-matrix into the corresponding morphism $A : R^n \to F$,
with `rank(F) == m`.
"""
function map(F::FreeMod{T}, A::MatrixElem{T}) where {T <: RingElement}
  if is_graded(F)
    return graded_map(F,A)
  end
  R = base_ring(F)
  F_domain = FreeMod(R, nrows(A))

  phi = FreeModuleHom(F_domain, F, A)
  return phi
end

@doc raw"""
    map(A::MatElem)

Converts a given $n \times m$-matrix into the corresponding morphism $A : R^n \to R^m$.
"""
function map(A::MatElem)
  R = base_ring(A)
  F_codomain = FreeMod(R, ncols(A))
  return map(F_codomain,A)
end

@doc raw"""
    is_injective(f::ModuleFPHom)

Test if `f` is injective.
"""
function is_injective(f::ModuleFPHom)
  return iszero(kernel(f)[1])
end

@doc raw"""
    is_surjective(f::ModuleFPHom)

Test if `f` is surjective.
"""
function is_surjective(f::ModuleFPHom)
  return image(f)[1] == codomain(f)
end

@doc raw"""
    is_bijective(f::ModuleFPHom)

Test if `f` is bijective.
"""
function is_bijective(f::ModuleFPHom)
  return is_injective(f) && is_surjective(f)
end

