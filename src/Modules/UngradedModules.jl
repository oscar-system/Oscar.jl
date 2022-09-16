export presentation, coords, coeffs, repres, cokernel, index_of_gen, sub,
      quo, presentation, present_as_cokernel, is_equal_with_morphism, 
      std_basis, groebner_basis, reduced_groebner_basis, leading_module, 
      reduce, hom_tensor, hom_product, coordinates, 
      represents_element, free_resolution, free_resolution_via_kernels,
      element_to_homomorphism, homomorphism_to_element, generator_matrix, 
      restrict_codomain, restrict_domain, direct_product, tensor_product, 
      free_module, tor, ext, lift_homomorphism_contravariant, 
      lift_homomorphism_covariant, hom_without_reversing_direction, transport, 
      find_morphism, find_morphisms, is_canonically_isomorphic, 
      is_canonically_isomorphic_with_map, register_morphism!, dense_row,
      show_subquo, show_morphism, show_morphism_as_map,
      simplify_light, simplify_with_same_ambient_free_module,
      matrix_kernel, simplify, map, is_injective,
      is_surjective, is_bijective, is_welldefined, subquotient,
      multiplication_morphism, multiplication_induced_morphism,
      ambient_free_module, ambient_module, ambient_representative,
      ambient_representatives_generators, relations, img_gens, is_complete,
      chain_complex, cochain_complex

# TODO replace asserts by error messages?

#TODO: "fix" to allow QuoElem s as well...
# this requires
#  re-typeing of FreeModule
#  typing of BiModArray
# ... and all the rest.
# parametrization has to be by elem_type(coeff_ring) and not, like currently, the bottom coeff ring
# Also: qring is a Singular native. So it needs to be added to the ring creation

###############################################################################
# FreeMod constructors
###############################################################################

@doc Markdown.doc"""
    FreeMod(R::Ring, n::Int, name::String = "e"; cached::Bool = false)

Construct a free module over the ring `R` with rank `n`.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for 
the standard unit vectors.
"""
function FreeMod(R::Ring, n::Int, name::String = "e"; cached::Bool = false) # TODO cached?
  return FreeMod{elem_type(R)}(n, R, [Symbol("$name[$i]") for i=1:n])
end

@doc Markdown.doc"""
    free_module(R::MPolyRing, p::Int, name::String = "e"; cached::Bool = false)

Return the free $R$-module $R^p$, created with its basis of standard unit vectors.

The string `name` specifies how the basis vectors are printed. 

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y over Rational Field

julia> F[2]
e[2]

julia> typeof(F)
FreeMod{fmpq_mpoly}
```
"""
free_module(R::MPolyRing, p::Int, name::String = "e"; cached::Bool = false) = FreeMod(R, p, name, cached = cached)

#=XXX this cannot be as it is inherently ambigous
  - FreeModule(R, n)
  - direct sum of rings, ie. a ring
  - set of n-th powers of R
thus the "category" needs to be set explicitly

=#


function (F::FreeMod)()
  return FreeModElem(sparse_row(base_ring(F)), F)
end

function show(io::IO, F::FreeMod)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Free module of rank $(F.n) over ")
  print(IOContext(io, :compact =>true), F.R)
end

@doc Markdown.doc"""
    rank(F::FreeMod)
    ngens(F::AbstractFreeMod)
    dim(F::AbstractFreeMod)

Return the rank of `F`.
"""
dim(F::AbstractFreeMod) = rank(F)
rank(F::FreeMod) = F.n
ngens(F::AbstractFreeMod) = rank(F)

@doc Markdown.doc"""
    ==(F::FreeMod, G::FreeMod)

Return  `true` if `F` and `G` are equal, `false` otherwise.

Here, `F` and `G` are equal iff their base rings, ranks, and names for printing the basis elements are equal.
"""
function (==)(F::FreeMod, G::FreeMod)
  # two free modules are equal if the rank and the ring are
  # TODO it this enough or e.g. stored morphisms also be considered?
  return F.R == G.R && rank(F) == rank(G) && F.S == G.S
end

function hash(F::FreeMod, h::UInt)
  return hash((base_ring(F), rank(F), F.S), h)
end

@doc Markdown.doc"""
    is_isomorphic(F::FreeMod, G::FreeMod)

Return  `true` if `F` and `G` are isomorphic, `false` otherwise.

Here, `F` and `G` are isomorphic iff their base rings and ranks are equal.
"""
function is_isomorphic(F::FreeMod, G::FreeMod)
  return F.R == G.R && rank(F) == rank(G)
end

@doc Markdown.doc"""
    iszero(F::AbstractFreeMod)

Return `true` if `F` is the zero module, `false` otherwise.
"""
function iszero(F::AbstractFreeMod)
  return rank(F) == 0
end

@doc Markdown.doc"""
    canonical_isomorphism(F::FreeMod{T}, G::FreeMod{T})

For `F` and `G` which have equal rank (otherwise an error is thrown)
return the canonical isomorphism, that is map the i-th basis vector of the
canonical basis of `F` to the i-th basis vector of the canonical basis of `G`.
"""
function canonical_isomorphism(F::FreeMod{T}, G::FreeMod{T}) where T
  @assert rank(F) == rank(G)
  return hom(F, G, gens(G))
end

###############################################################################
# FreeModElem constructors
###############################################################################

@doc Markdown.doc"""
    FreeModElem(c::SRow{T}, parent::FreeMod{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
FreeModElem(c::SRow{T}, parent::FreeMod{T}) where T = FreeModElem{T}(c, parent)

@doc Markdown.doc"""
    FreeModElem(c::Vector{T}, parent::FreeMod{T}) where T

Return the element of `F` whose coefficients with respect to the basis of
standard unit vectors of `F` are given by the entries of `c`.
"""
function FreeModElem(c::Vector{T}, parent::FreeMod{T}) where T
  @assert length(c) == rank(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:rank(parent)), c)
  return FreeModElem{T}(sparse_coords,parent)
end

#@doc Markdown.doc"""
#    (F::FreeMod{T})(c::SRow{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#"""
function (F::FreeMod{T})(c::SRow{T}) where T
  return FreeModElem(c, F)
end

#@doc Markdown.doc"""
#    (F::FreeMod{T})(c::Vector{T}) where T
#
#Return the element of `F` whose coefficients with respect to the basis of
#standard unit vectors of `F` are given by the entries of `c`.
#
## Examples
#```jldoctest
#julia> R, (x,y) = PolynomialRing(QQ, ["x", "y"])
#(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])
#
#julia> F = FreeMod(R,3)
#Free module of rank 3 over Multivariate Polynomial Ring in x, y over Rational Field
#
#julia> V = [x, zero(R), y]
#3-element Vector{fmpq_mpoly}:
# x
# 0
# y
#
#julia> f = F(V)
#x*e[1] + y*e[3]
#```
#"""
function (F::FreeMod{T})(c::Vector{T}) where T 
 return FreeModElem(c, F)
end

function (F::AbstractFreeMod{T})(v::AbstractFreeModElem{T}) where T
  @assert parent(v) === F
  return v
end

function in(v::AbstractFreeModElem, F::AbstractFreeMod)
  return parent(v) === F
end

function in(v::AbstractFreeModElem, M::SubQuo)
  return represents_element(v, M)
end

@doc Markdown.doc"""
    coords(v::AbstractFreeModElem)

Return the entries (with respect to the standard basis) of `v` as a sparse row.
"""
function coords(v::AbstractFreeModElem)
  return v.coords
end

@doc Markdown.doc"""
    coeffs(v::AbstractFreeModElem)

Return the entries (with respect to the standard basis) of `v` as a sparse row.
"""
function coeffs(v::AbstractFreeModElem)
  return coords(v)
end

#######################################################
@doc Markdown.doc"""
    coefficients(f::AbstractFreeModElem)

Return the coefficients of `f` with respect to the basis of standard unit vectors. 

The result is returned as a sparse row.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = FreeMod(R,3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y over Rational Field

julia> f = x*gen(F,1)+y*gen(F,3)
x*e[1] + y*e[3]

julia> coefficients(f)
Sparse row with positions [1, 3] and values fmpq_mpoly[x, y]
```
"""
coefficients(f::AbstractFreeModElem) = coeffs(f)
#########################################################

@doc Markdown.doc"""
    repres(v::AbstractFreeModElem)

Return just `v`. This function exists for compatibility (with subquotient elements) reasons.
"""
function repres(v::AbstractFreeModElem)
  return v
end

function getindex(v::AbstractFreeModElem, i::Int)
  if isempty(coords(v))
    return zero(base_ring(parent(v)))
  end
  return coords(v)[i]
end

elem_type(::Type{FreeMod{T}}) where {T} = FreeModElem{T}
parent_type(::Type{FreeModElem{T}}) where {T} = FreeMod{T}
elem_type(::FreeMod{T}) where {T} = FreeModElem{T}
parent_type(::FreeModElem{T}) where {T} = FreeMod{T}

function generator_symbols(F::FreeMod)
  return F.S
end

function expressify(e::AbstractFreeModElem; context = nothing)
  sum = Expr(:call, :+)
  for (pos, val) in e.coords
     # assuming generator_symbols(parent(e)) is an array of strings/symbols
     push!(sum.args, Expr(:call, :*, expressify(val, context = context), generator_symbols(parent(e))[pos]))
  end
  return sum
end
@enable_all_show_via_expressify FreeModElem

@doc Markdown.doc"""
    Vector(e::FreeModElem)

Return the coefficients of `e` as a Vector.
"""
function Vector(e::FreeModElem)
   return [e[i] for i in 1:rank(parent(e))]
end

@doc Markdown.doc"""
    basis(F::AbstractFreeMod)

Return the standard basis of `F`.
"""
function basis(F::AbstractFreeMod)
  bas = elem_type(F)[]
  for i=1:dim(F)
    s = Hecke.sparse_row(base_ring(F), [(i, base_ring(F)(1))])
    push!(bas, FreeModElem(s, F))
  end
  return bas
end


@doc Markdown.doc"""
    gens(F::AbstractFreeMod)

Return the (canonical) generators of the free module `F`.
"""
gens(F::AbstractFreeMod) = basis(F)

@doc Markdown.doc"""
    basis(F::AbstractFreeMod, i::Int)

    gen(F::AbstractFreeMod, i::Int)

Return the `i`th basis vector of `F`, that is, return the `i`th standard unit vector.
"""
function basis(F::AbstractFreeMod, i::Int)
  @assert 0 < i <= ngens(F)
  s = Hecke.sparse_row(base_ring(F), [(i, base_ring(F)(1))])
  return FreeModElem(s, F)
end
gen(F::AbstractFreeMod, i::Int) = basis(F,i)

function getindex(F::AbstractFreeMod, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

@doc Markdown.doc"""
    base_ring(F::AbstractFreeMod)

Return the underlying ring of `F`.
"""
base_ring(F::FreeMod) = F.R

#TODO: Parent - checks everywhere!!!

# the negative of a free module element
-(a::AbstractFreeModElem) = FreeModElem(-coords(a), parent(a))

# Addition of free module elements
function +(a::AbstractFreeModElem, b::AbstractFreeModElem)
   check_parent(a, b)
   return FreeModElem(coords(a)+coords(b), parent(a))
end

# Subtraction of free module elements
function -(a::AbstractFreeModElem, b::AbstractFreeModElem)
    check_parent(a,b)
    return FreeModElem(coords(a)-coords(b), parent(a))
end

# Equality of free module elements
function (==)(a::AbstractFreeModElem, b::AbstractFreeModElem) 
  if parent(a) !== parent(b)
    return false
  end
  return a.coords == b.coords
end

function hash(a::AbstractFreeModElem, h::UInt)
  return hash(tuple(parent(a), coords(a)), h)
end

# scalar multiplication with polynomials, integers
function *(a::MPolyElem_dec, b::AbstractFreeModElem)
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
  end
  return FreeModElem(a*coords(b), parent(b))
end
function *(a::MPolyElem, b::AbstractFreeModElem) 
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  return FreeModElem(a*coords(b), parent(b))
end
function *(a::RingElem, b::AbstractFreeModElem) 
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  return FreeModElem(a*coords(b), parent(b))
end
*(a::Int, b::AbstractFreeModElem) = FreeModElem(a*coords(b), parent(b))
*(a::Integer, b::AbstractFreeModElem) = FreeModElem(base_ring(parent(b))(a)*coords(b), parent(b))
*(a::fmpq, b::AbstractFreeModElem) = FreeModElem(base_ring(parent(b))(a)*coords(b), parent(b))

@doc Markdown.doc"""
    zero(F::AbstractFreeMod)

Return the zero element of  `F`.
"""
zero(F::AbstractFreeMod) = FreeModElem(sparse_row(base_ring(F), Tuple{Int, elem_type(base_ring(F))}[]), F)

@doc Markdown.doc"""
    parent(a::AbstractFreeModElem)

Return the free module where `a` lives in.
"""
parent(a::AbstractFreeModElem) = a.parent

@doc Markdown.doc"""
    iszero(f::AbstractFreeModElem)

Return `true` if `f` is zero, `false` otherwise.
"""
iszero(f::AbstractFreeModElem) = iszero(coords(f))

###############################################################################
# ModuleGens constructors
###############################################################################

@doc Markdown.doc"""
    ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, SF::Singular.FreeMod) where T

Construct `ModuleGens` from an array of free module elements, specifying the free module 
and Singular free module. 
This function is only useful indirectly.
"""
ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, SF::Singular.FreeMod) where T = ModuleGens{T}(O, F, SF)

@doc Markdown.doc"""
    ModuleGens(F::FreeMod{S}, s::Singular.smodule) where {S}

Construct `ModuleGens` from a given Singular submodule.
"""
ModuleGens(F::FreeMod{S}, s::Singular.smodule) where {S} = ModuleGens{S}(F, s)

@doc Markdown.doc"""
    ModuleGens(O::Vector{<:FreeModElem})

Construct `ModuleGens` from an array of free module elements.

!!! note 
    
    The array must not be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem})
  # TODO Empty generating set
  @assert length(O) > 0
  #SF = singular_module(parent(O[1]))
  #return ModuleGens(O, SF)
  return ModuleGens(O, parent(O[1]))
end

@doc Markdown.doc"""
    ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}) where {T}

Construct `ModuleGens` from an array of free module elements, specifying the free module.

!!! note

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}) where {T}
  #SF = singular_module(F)
  return ModuleGens{T}(O, F)
end

@doc Markdown.doc"""
    ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, ordering::ModuleOrdering) where {T}

Construct `ModuleGens` from an array of free module elements, specifying the free module.
Moreover, the ordering is defined by `ordering`. 

!!! note

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, ordering::ModuleOrdering) where {T}
  SF = singular_module(F, ordering)
  M = ModuleGens{T}(O, F, SF)
  M.ordering = ordering
  return M
end

@doc Markdown.doc"""
    ModuleGens(O::Vector{<:FreeModElem}, SF::Singular.FreeMod)

Construct `ModuleGens` from an array of free module elements, specifying the Singular free module.

!!! note 

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, SF::Singular.FreeMod)
  return ModuleGens{elem_type(base_ring(parent(O[1])))}(O, parent(O[1]), SF)
end

@doc Markdown.doc"""
    base_ring(M::ModuleGens)

Return the base ring of `M` (that is, if `M` is an `R`-module, return `R`).
"""
function base_ring(M::ModuleGens)
  return base_ring(M.F)
end

@doc Markdown.doc"""
    singular_generators(M::ModuleGens)

Return the generators of `M` from Singular side.
"""
function singular_generators(M::ModuleGens)
  singular_assure(M)
  return M.S
end

@doc Markdown.doc"""
    oscar_generators(M::ModuleGens)  

Return the generators of `M` from the Oscar side.
"""
function oscar_generators(M::ModuleGens)
  oscar_assure(M)
  return M.O
end

@doc Markdown.doc"""
    iszero(M::ModuleGens)

Check if `M` is zero.
"""
function iszero(M::ModuleGens)
  return iszero(singular_generators(M))
end

function show(io::IO, F::ModuleGens)
  if F.isGB
    show_gb(io, F)
  else
    print(io, "Module generating system of length ", length(F))
    for i=1:length(F)
      if isassigned(F.O, i)
        print(io, "\n", i, " -> ", F.O[i])
      end
    end
  end
end

@doc Markdown.doc"""
    length(F::ModuleGens)

Return the number of elements of the module generating set.

!!! note

    This is not the length in the mathematical sense!
"""
length(F::ModuleGens) = length(oscar_generators(F))

@doc Markdown.doc"""
    ngens(F::ModuleGens)

Return the number of elements of the module generating set.
"""
ngens(F::ModuleGens) = length(oscar_generators(F))

# i-th entry of module generating set on Oscar side
# Todo: clean up, convert or assure
function getindex(F::ModuleGens, ::Val{:O}, i::Int)
  if !isassigned(F.O, i)
    F.O[i] = F.F(singular_generators(F)[i])
  end
  return oscar_generators(F)[i]
end

# i-th entry of module generating set on Singular side
# Todo: clean up, convert or assure
function getindex(F::ModuleGens, ::Val{:S}, i::Int)
  singular_assure(F)
  if !isdefined(F, :S)
    F.S = Singular.Module(base_ring(F.SF), [F.SF(x) for x = oscar_generators(F)]...)
  end
  return F.S[i]
end

@doc Markdown.doc"""
    oscar_assure(F::ModuleGens)

If fields of `F` from the Oscar side are not defined, they
are computed, given the Singular side.
"""
function oscar_assure(F::ModuleGens)
  if !isdefined(F, :O)
    F.O = [F.F(singular_generators(F)[i]) for i=1:Singular.ngens(singular_generators(F))]
  end
end

@doc Markdown.doc"""
    singular_assure(F::ModuleGens)

If fields of `F` from the Singular side are not defined, they
are computed, given the Oscar side.
"""
function singular_assure(F::ModuleGens)
  if !isdefined(F, :S) || !isdefined(F, :SF)
    SF = singular_module(F.F)
    sr = base_ring(SF)
    F.SF = SF
    if length(F) == 0
      F.S = Singular.Module(sr, Singular.vector(sr, sr(0)))
      return 
    end
    F.S = Singular.Module(base_ring(F.SF), [F.SF(x) for x = oscar_generators(F)]...)
    return
  end
  #F[Val(:S), 1]
end

# i-th entry of module generating set (taken from Oscar side)
getindex(F::ModuleGens, i::Int) = getindex(F, Val(:O), i)

@doc Markdown.doc"""
    union(M::ModuleGens, N::ModuleGens)

Compute the union of `M` and `N`.
"""
function union(M::ModuleGens, N::ModuleGens)
  @assert M.F === N.F
  O = vcat(M.O, N.O)
  return ModuleGens(M.F, O)
end

@doc Markdown.doc"""
    singular_module(F::FreeMod)

Create a Singular module from a given free module.
"""
function singular_module(F::FreeMod)
  Sx = singular_poly_ring(base_ring(F), keep_ordering=false)
  return Singular.FreeModule(Sx, dim(F))
end

@doc Markdown.doc"""
    singular_module(F::FreeMod, ordering::ModuleOrdering)

Create a Singular module from a given free module over the given Singular polynomial ring.
"""
function singular_module(F::FreeMod, ordering::ModuleOrdering)
  Sx = singular_ring(base_ring(F), singular(ordering))
  return Singular.FreeModule(Sx, dim(F))
end

@doc Markdown.doc"""
    (SF::Singular.FreeMod)(m::FreeModElem)

Convert a free module element to the Singular side.
"""
function (SF::Singular.FreeMod)(m::FreeModElem)
  g = Singular.gens(SF)
  e = SF()
  Sx = base_ring(SF)
  for (p,v) = m.coords
    e += Sx(v)*g[p]
  end
  return e
end

@doc Markdown.doc"""
    (F::FreeMod)(s::Singular.svector)

Convert a Singular vector to a free module element.
"""
function (F::FreeMod)(s::Singular.svector)
  pos = Int[]
  values = []
  Rx = base_ring(F)
  R = base_ring(Rx)
  for (i, e, c) = s
    f = Base.findfirst(x->x==i, pos)
    if f === nothing
      push!(values, MPolyBuildCtx(base_ring(F)))
      f = length(values)
      push!(pos, i)
    end
    push_term!(values[f], R(c), e)
  end
  pv = Tuple{Int, elem_type(Rx)}[(pos[i], base_ring(F)(finish(values[i]))) for i=1:length(pos)]
  return FreeModElem(sparse_row(base_ring(F), pv), F)
end

###############################################################################
# FreeModuleHom constructors
###############################################################################

#=@doc Markdown.doc"""
    FreeModuleHom(F::FreeMod{T}, G::S, a::Vector) where {T, S}

Construct the morphism $F \to G$ where `F[i]` is mapped to `a[i]`.
In particular, `ngens(F) == length(a)` must hold.
"""
FreeModuleHom(F::AbstractFreeMod{T}, G::S, a::Vector) where {T, S} = FreeModuleHom{T,S}(F, G, a)

@doc Markdown.doc"""
    FreeModuleHom(F::FreeMod{T}, G::S, mat::MatElem{T}) where {T,S}

Construct the morphism $F \to G$ corresponding to the matrix `mat`.
"""
FreeModuleHom(F::AbstractFreeMod{T}, G::S, mat::MatElem{T}) where {T,S} = FreeModuleHom{T,S}(F, G, mat)=#

img_gens(f::FreeModuleHom) = gens(image(f)[1])
base_ring_map(f::FreeModuleHom) = f.ring_map

@doc Markdown.doc"""
    matrix(a::FreeModuleHom)

Given a homomorphism `a` of type  `FreeModuleHom` with domain `F`
and codomain `M`, return a matrix `A` with `rank(F)` rows and 
`ngens(M)` columns such that `a == hom(F, M, A)`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]];

julia> a = hom(F, G, V);

julia> matrix(a)
[y   0]
[x   y]
[0   z]
```
"""
function matrix(a::FreeModuleHom)
  if !isdefined(a, :matrix)
    D = domain(a)
    C = codomain(a)
    R = base_ring(D)
    matrix = zero_matrix(R, rank(D), ngens(C))
    for i=1:rank(D)
      image_of_gen = a(D[i])
      for j=1:ngens(C)
        matrix[i,j] = image_of_gen[j]
      end
    end
    setfield!(a, :matrix, matrix)
  end
  return a.matrix
end

(h::FreeModuleHom)(a::AbstractFreeModElem) = image(h, a)

@doc Markdown.doc"""
    hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T

Given a vector `V` of `rank(F)` elements of `M`, 
return the homomorphism `F` $\to$ `M` which sends the `i`-th
basis vector of `F` to the `i`-th entry of `V`.

    hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}) where T

Given a matrix `A` with `rank(F)` rows and `ngens(M)` columns, return the
homomorphism `F` $\to$ `M` which sends the `i`-th basis vector of `F` to 
the linear combination $\sum_j A[i,j]*M[j]$ of the generators `M[j]` of `M`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> V = [y*G[1], x*G[1]+y*G[2], z*G[2]]
3-element Vector{FreeModElem{fmpq_mpoly}}:
 y*e[1]
 x*e[1] + y*e[2]
 z*e[2]

julia> a = hom(F, G, V)
Map with following data
Domain:
=======
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field
Codomain:
=========
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> a(F[2])
x*e[1] + y*e[2]

julia> W = matrix(R, [y 0; x y; 0 z])
[y   0]
[x   y]
[0   z]

julia> b = hom(F, G, W)
Map with following data
Domain:
=======
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field
Codomain:
=========
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> a == b
true
```
"""
function hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T
  base_ring(F) === base_ring(M) || return FreeModuleHom(F, M, V, base_ring(M))
  return FreeModuleHom(F, M, V)
end
function hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}) where T 
  base_ring(F) === base_ring(M) || return FreeModuleHom(F, M, A, base_ring(M))
  return FreeModuleHom(F, M, A)
end

@doc Markdown.doc"""
    hom(F::FreeMod, M::ModuleFP, V::Vector{<:ModuleFPElem}, h::RingMapType) where {RingMapType}

Given a vector `V` of `rank(F)` elements of `M` and a ring map `h`
from `base_ring(F)` to `base_ring(M)`, return the 
`base_ring(F)`-homomorphism `F` $\to$ `M` which sends the `i`-th
basis vector of `F` to the `i`-th entry of `V`.

    hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}, h::RingMapType) where {T, RingMapType}

Given a matrix `A` over `base_ring(M)` with `rank(F)` rows and `ngens(M)` columns
and a ring map `h` from `base_ring(F)` to `base_ring(M)`, return the
`base_ring(F)`-homomorphism `F` $\to$ `M` which sends the `i`-th basis vector of `F` to 
the linear combination $\sum_j A[i,j]*M[j]$ of the generators `M[j]` of `M`.
"""
hom(F::FreeMod, M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, h::RingMapType) where {T, RingMapType} = FreeModuleHom(F, M, V, h)
hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}, h::RingMapType) where {T, RingMapType} = FreeModuleHom(F, M, A, h)

@doc Markdown.doc"""
    identity_map(M::ModuleFP)

Return the identity map $id_M$.
"""
function identity_map(M::ModuleFP)
  return hom(M, M, gens(M))
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


###############################################################################
# SubModuleOfFreeModule
###############################################################################
@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModElem}) where {R}

Construct the submodule of `F` generated by the elements of `gens` (the elements of 
`gens` must live in `F`, if not, an execption is thrown).
"""
function SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModElem}) where {R}
  @assert all(x -> parent(x) === F, gens)
  return SubModuleOfFreeModule(F, ModuleGens(gens, F))
end

@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModElem}, default_ordering::ModuleOrdering) where {R}

Construct the submodule of `F` generated by the elements of `gens` (the elements of 
`gens` must live in `F`, if not, an execption is thrown). 
Moreover, set the default ordering to `default_ordering`.
"""
function SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModElem}, default_ordering::ModuleOrdering) where {R}
  @assert all(x -> parent(x) === F, gens)
  return SubModuleOfFreeModule(F, ModuleGens(gens, F, default_ordering))
end

function SubModuleOfFreeModule(F::FreeMod{R}, singular_module::Singular.smodule) where {R} 
  return SubModuleOfFreeModule(F, ModuleGens(F, singular_module))
end

@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{R}, gens::ModuleGens{R}) where {R}

Construct the submodule of `F` generated by `gens`.
"""
function SubModuleOfFreeModule(F::FreeMod{R}, gens::ModuleGens{R}) where {R} 
  subModule = SubModuleOfFreeModule{R}(F)
  subModule.gens = gens
  return subModule
end

@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}) where {L}

Construct the submodule generated by the rows of `A`. The embedding free
module is `F`. In particular, `rank(F) == ncols(A)` must hold.
"""
function SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}) where {L} 
  subModule = SubModuleOfFreeModule{L}(F)
  O = [FreeModElem(sparse_row(A[i,:]), F) for i in 1:nrows(A)]
  subModule.gens = ModuleGens(O, F)
  subModule.matrix = A
  return subModule
end

@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}, default_ordering::ModuleOrdering) where {L}

Construct the submodule generated by the rows of `A`. The embedding free
module is `F`. In particular, `rank(F) == ncols(A)` must hold.
Morevoer, set the default ordering to `default_ordering`.
"""
function SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}, default_ordering::ModuleOrdering) where {L} 
  subModule = SubModuleOfFreeModule{L}(F)
  O = [FreeModElem(sparse_row(A[i,:]), F) for i in 1:nrows(A)]
  subModule.gens = ModuleGens(O, F, default_ordering)
  subModule.matrix = A
  return subModule
end

@doc Markdown.doc"""
    SubModuleOfFreeModule(A::MatElem{L}) where {L} 

Construct the submodule generated by the rows of `A`.

!!! note

    The ambient free module of the submodule is constructed by the function and therefore
    not compatible with free modules that are defined by the user or by other functions. 
    For compatibility, use `SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}) where {L}`.
"""
function SubModuleOfFreeModule(A::MatElem{L}) where {L} 
  R = base_ring(A)
  F = FreeMod(R, ncols(A))
  return SubModuleOfFreeModule(F, A)
end


@doc Markdown.doc"""
    SubModuleOfFreeModule(A::MatElem{L}, default_ordering::ModuleOrdering) where {L} 

Construct the submodule generated by the rows of `A`.

!!! note

    The ambient free module of the submodule is constructed by the function and therefore
    not compatible with free modules that are defined by the user or by other functions. 
    For compatibility, use `SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}) where {L}`.
"""
function SubModuleOfFreeModule(A::MatElem{L}, default_ordering::ModuleOrdering) where {L} 
  R = base_ring(A)
  F = FreeMod(R, ncols(A))
  return SubModuleOfFreeModule(F, A, default_ordering)
end

function getindex(M::SubModuleOfFreeModule, i::Int)
  return oscar_generators(M.gens)[i]
end

@doc Markdown.doc"""
    iszero(M::SubModuleOfFreeModule)

Check if `M` is zero.
"""
function iszero(M::SubModuleOfFreeModule)
  return iszero(M.gens)
end

@doc Markdown.doc"""
    base_ring(M::SubModuleOfFreeModule)

Return the base ring of `M`.
"""
function base_ring(M::SubModuleOfFreeModule)
  return base_ring(M.F)
end

@doc Markdown.doc"""
    ambient_free_module(M::SubModuleOfFreeModule)

Return the ambient free module of `M`.
"""
function ambient_free_module(M::SubModuleOfFreeModule)
  return M.F
end

@doc Markdown.doc"""
    default_ordering(M::SubModuleOfFreeModule)

Get the default ordering of `M`.
"""
function default_ordering(M::SubModuleOfFreeModule)
  if !isdefined(M, :default_ordering)
    ord = default_ordering(ambient_free_module(M))
    set_default_ordering!(M, ord)
    return M.default_ordering
  end
  return M.default_ordering
end

@doc Markdown.doc"""
    set_default_ordering!(M::SubModuleOfFreeModule, ord::ModuleOrdering)

Set the default ordering in `M` to `ord`.
"""
function set_default_ordering!(M::SubModuleOfFreeModule, ord::ModuleOrdering)
  M.default_ordering = ord
end

@doc Markdown.doc"""
    std_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod))

Compute a standard basis of `submod` with respect to the given `odering``.
The return type is `ModuleGens`.
"""
function std_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod))
  gb = get!(submod.groebner_basis, ordering) do
    return compute_std_basis(submod, ordering)
  end::ModuleGens
  return gb
end

@doc Markdown.doc"""
    groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod))

Compute a Gröbner of `submod` with respect to the given `ordering`.
The ordering must be global. The return type is `ModuleGens`.
"""
function groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod))
  @assert is_global(ordering)
  return std_basis(submod, ordering)
end

@doc Markdown.doc"""
    reduced_groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod))

Compute a reduced Gröbner basis with respect to the given `ordering`. The return type is `ModuleGens`.
"""
function reduced_groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod))
  @assert is_global(ordering)

  gb = get!(submod.groebner_basis, ordering) do
    return compute_std_basis(submod, ordering, true)
  end::ModuleGens
  gb.is_reduced && return gb
  return get_attribute!(gb, :reduced_groebner_basis) do
    return compute_std_basis(submod, ordering, true)
  end::ModuleGens
end

function leading_module(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod))
  gb = std_basis(submod, ordering)
  return SubModuleOfFreeModule(submod.F, leading_monomials(gb))
end

@doc Markdown.doc"""
    compute_std_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod), reduced::Bool=false)

Compute a standard basis of `submod` with respect to the given ordering.
Allowed orderings are those that are allowed as orderings for Singular polynomial rings.
In case `reduced` is `true` and the ordering is global, a reduced Gröbner basis is computed.
"""
function compute_std_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering = default_ordering(submod), reduced::Bool=false)
  if reduced
    @assert is_global(ordering)
  end
  mg = ModuleGens(oscar_generators(submod.gens), submod.F , ordering)
  gb = std_basis(mg, reduced)
  oscar_assure(gb)
  gb.isGB = true
  gb.S.isGB = true
  gb.ordering = ordering
  gb.is_reduced = reduced
  return gb
end

function show_gb(io::IO, GB::ModuleGens)
  if isdefined(GB, :quo_GB)
    show_relative_groebner_basis(io, GB, GB.quo_GB, GB.is_reduced)
  else
    show_non_relative_groebner_basis(io, GB, GB.is_reduced)
  end
end

function show_groebner_basis_helper(io::IO, sub::ModuleGens, init::String)
  F = sub.F

  print(io, init)

  # TODO Remove this / use Oscar orderings when available
  singular_assure(sub)
  vectors = Vector{Vector{elem_type(F)}}(undef, length(sub.O))
  for i in 1:length(sub.O)
    monomials = []
    v = sub.S[i]
    l = Singular.lead(v)
    while !iszero(l)
      push!(monomials, l)
      v = v - l
      l = Singular.lead(v)
    end
    monomials = Vector{elem_type(F)}(map(x -> F(x), monomials))
    vectors[i] = monomials
  end

  for v in vectors
    print(io, "\n")
    print(io, v[1])
    for i in 2:length(v)
      print(io, " + ", v[i])
    end
  end
end

function show_non_relative_groebner_basis(io::IO, sub::ModuleGens, reduced::Bool = false)
  init = "Gröbner basis with elements"
  if reduced
    init = "Reduced " * init
  end
  show_groebner_basis_helper(io, sub, init)
end

function show_relative_groebner_basis(io::IO, sub::ModuleGens, quo::ModuleGens, reduced::Bool = false)
  init = "Gröbner basis with elements"
  if reduced
    init = "Reduced relative " * init
  else
    init = "Relative " * init
  end
  show_groebner_basis_helper(io, sub, init)
  print(io, "\n")
  show_groebner_basis_helper(io, quo, "\nwith quotient Gröbner basis defined by the elements")
end

@doc Markdown.doc"""
    generator_matrix(submod::SubModuleOfFreeModule)

Return the generators of `submod` in matrix-form, that is the rows of the 
matrix generate `submod`.
"""
function generator_matrix(submod::SubModuleOfFreeModule)
  if !isdefined(submod, :matrix)
    R = base_ring(submod)
    matrix = zero_matrix(R, length(submod.gens), rank(submod.F))
    for i = 1:nrows(matrix), j = 1:ncols(matrix)
      matrix[i,j] = submod.gens[i][j] 
    end
    submod.matrix = matrix
  end
  return submod.matrix
end

@doc Markdown.doc"""
    is_generated_by_standard_unit_vectors(M::SubModuleOfFreeModule)

Check if `M` is generated by the standard unit vectors.
"""
function is_generated_by_standard_unit_vectors(M::SubModuleOfFreeModule)
  return issubset(gens(M.F), gens(M))
end

function show(io::IO, M::SubModuleOfFreeModule)
  if ngens(M) == 1
    print(io, "Submodule with ", ngens(M), " generator")
  else
    print(io, "Submodule with ", ngens(M), " generators")
  end
  for i=1:ngens(M)
    if isassigned(M.gens.O, i)
      print(io, "\n", i, " -> ", M[i])
    end
  end
end

function length(M::SubModuleOfFreeModule)
  error("use ngens() instead of length()")
  return length(M.gens)
end

@doc Markdown.doc"""
    ngens(M::SubModuleOfFreeModule)

Return the number of generators of `M`.
"""
function ngens(M::SubModuleOfFreeModule)
  return ngens(M.gens)
end

@doc Markdown.doc"""
    gens(M::SubModuleOfFreeModule)

Return the generators of `M` as an array of `FreeModElem`s.
"""
function gens(M::SubModuleOfFreeModule)
  return oscar_generators(M.gens)
end

@doc Markdown.doc"""
    gen(M::SubModuleOfFreeModule, i::Int)

Return the `i`th generator of `M`.
"""
function gen(M::SubModuleOfFreeModule, i::Int)
  return oscar_generators(M.gens)[i]
end

@doc Markdown.doc"""
    sum(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Compute $M+N$.
"""
function sum(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  @assert M.F === N.F
  return SubModuleOfFreeModule(M.F, vcat(collect(M.gens), collect(N.gens)))
end

@doc Markdown.doc"""
    issubset(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Check if `M` is a subset of `N`. For this their embedding free modules must be 
identical (`===`).
"""
function issubset(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  if M.F !== N.F
    return false
  end
  M_mod_N = reduce(M, N)
  return iszero(M_mod_N)
end

@doc Markdown.doc"""
    ==(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Check for equality. For two submodules of free modules to be equal their embedding 
free modules must be identical (`===`) and the generators must generate equal submodules.
"""
function (==)(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  if M === N
    return true
  end
  if M.F !== N.F
    return false
  end
  #TODO should there be a check for === up to permutation in order to avoid std-computation?
  M_mod_N = reduce(M, N)
  N_mod_M = reduce(N, M)
  return iszero(M_mod_N) && iszero(N_mod_M)
end

@doc Markdown.doc"""
    is_canonically_isomorphic(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Check if `M` are canonically isomorphic. This means that if `F = ambient_free_module(M)` is 
isomorphic to `G = ambient_free_module(N)` and the image of `M` under `canonical_isomorphism(F, G)`
is equal to `N`.
"""
function is_canonically_isomorphic(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  F = ambient_free_module(M)
  G = ambient_free_module(N)
  if !is_isomorphic(F, G)
    return false
  end
  f = canonical_isomorphism(F,G)
  return SubModuleOfFreeModule(G, [f(v) for v in gens(M)]) == N
end

#+(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule) = sum(M, N)

###############################################################################
# SubQuo constructors
###############################################################################

@doc Markdown.doc"""
    SubQuo(sub::SubModuleOfFreeModule{R}) where {R}

Construct the module `sub` as a subquotient.
"""
function SubQuo(sub::SubModuleOfFreeModule{R}) where {R}
  subquo = SubQuo{R}(ambient_free_module(sub))
  subquo.sub = sub
  subquo.sum = sub
  return subquo
end

@doc Markdown.doc"""
    SubQuo(sub::SubModuleOfFreeModule{R}, quo::SubModuleOfFreeModule{R}) where {R}

Construct the subquotient module $\texttt{sub} + texttt{quo} / \texttt{quo}$. 
`sub` and `quo` must be submodules of the same free module.
"""
function SubQuo(sub::SubModuleOfFreeModule{R}, quo::SubModuleOfFreeModule{R}) where {R}
  F = ambient_free_module(sub)
  @assert F === ambient_free_module(quo)
  subquo = SubQuo{R}(F)
  subquo.sub = sub
  subquo.quo = quo
  subquo.sum = sum(subquo.sub, subquo.quo)
  return subquo
end

@doc Markdown.doc"""
    SubQuo(F::FreeMod{R}, O::Vector{<:FreeModElem}) where {R}

Construct the module generated by the elements of `O` as a subquotient.
The elements of `O` must live in `F`.

# Examples
```jldoctest
julia> R, (x,y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = FreeMod(R,2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y over Rational Field

julia> O = [x*F[1]+F[2],y*F[2]]
2-element Vector{FreeModElem{fmpq_mpoly}}:
 x*e[1] + e[2]
 y*e[2]

julia> M = SubQuo(F, O)
Submodule with 2 generators
1 -> x*e[1] + e[2]
2 -> y*e[2]
represented as subquotient with no relations.

```
"""
function SubQuo(F::FreeMod{R}, O::Vector{<:FreeModElem}) where {R} 
  sub = SubModuleOfFreeModule(F, O)
  return SubQuo(sub)
end

@doc Markdown.doc"""
    SubQuo(S::SubQuo{L}, O::Vector{<:FreeModElem}) where {L}

Construct a subquotient where the generators are those of `S` and the relations are 
the union of `O` and the relations of `S`.
"""
function SubQuo(S::SubQuo{L}, O::Vector{<:FreeModElem}) where {L}
  F = ambient_free_module(S)
  subquo = SubQuo{L}(F)
  subquo.sub = S.sub
  O_as_submodule = SubModuleOfFreeModule(F, O)
  subquo.quo = isdefined(S,:quo) ? sum(S.quo,O_as_submodule) : O_as_submodule
  subquo.sum = sum(subquo.sub, subquo.quo)
  return subquo
end

function SubQuo(F::FreeMod{R}, s::Singular.smodule) where {R}
  subquo = SubQuo{R}(F)
  subquo.sub = SubModuleOfFreeModule(F, s)
  subquo.sum = subquo.sub
  return subquo
end

function SubQuo(F::FreeMod{R}, s::Singular.smodule, t::Singular.smodule) where {R}
  subquo = SubQuo{R}(F)
  subquo.sub = SubModuleOfFreeModule(F, s)
  subquo.quo = SubModuleOfFreeModule(F, t)
  subquo.sum = sum(subquo.sub, subquo.quo)
  return subquo
end

@doc Markdown.doc"""
    SubQuo(F::FreeMod{R}, A::MatElem{R}, B::MatElem{R}) where {R}

Given matrices `A` and `B` with entries in a ring `R` representing maps 
of free $R$-modules with the same codomain `F`, return the subquotient 
$(\text{im } A + \text{im }  B)/\text{im }  B$.
"""
function SubQuo(F::FreeMod{R}, A::MatElem{R}, B::MatElem{R}) where {R}
  @assert ncols(A) == ncols(B) == rank(F)
  return SubQuo(SubModuleOfFreeModule(F, A), SubModuleOfFreeModule(F, B))
end

@doc Markdown.doc"""
    SubQuo(A::MatElem{R}, B::MatElem{R}) where {R}

Given matrices `A` and `B` with entries in a ring `R` 
representing maps of free $R$-modules with the same codomain,
return the subquotient $(\text{im } A + \text{im }  B)/\text{im }  B.$ 

!!! note

    The ambient free module of the subquotient is constructed by the function and therefore
    not compatible with free modules defined by the user or by other functions. 
    For compatibility, use `SubQuo(F::FreeMod{R}, A::MatElem{R}, B::MatElem{R})`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubQuo(A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 4 generators
1 -> x^2*e[1]
2 -> x*y*e[1]
3 -> y^2*e[1]
4 -> z^4*e[1]
```
"""
function SubQuo(A::MatElem{R}, B::MatElem{R}) where {R}
  @assert ncols(A) == ncols(B)
  S = base_ring(A)
  F = FreeMod(S, ncols(A))
  return SubQuo(SubModuleOfFreeModule(F, A), SubModuleOfFreeModule(F, B))
end

@doc Markdown.doc"""
    SubQuo(F::FreeMod{T}, g::Vector{FreeModElem{T}}, q::Vector{FreeModElem{T}}) where {T<:RingElem} 

Construct the subquotient with ambient free module `F`, generators `g`
and relations `q`.
"""
function SubQuo(F::FreeMod{T}, g::Vector{FreeModElem{T}}, q::Vector{FreeModElem{T}}) where {T<:RingElem} 
  return SubQuo(SubModuleOfFreeModule(F, g), SubModuleOfFreeModule(F, q))
end

#######################################################
@doc Markdown.doc"""
    subquotient(a::FreeModuleHom, b::FreeModuleHom)

Given homomorphisms `a` and `b` between free modules such that 
`codomain(a) === codomain(b)`, 
return $(\text{im } a + \text{im } b)/\text{im } b$.

    subquotient(F::FreeMod{T}, A::MatElem{T}, B::MatElem{T}) where T

Given matrices `A` and `B` with rank `F` columns, return  
$(\text{im } a + \text{im } b)/\text{im } b$,
where `a` and `b` are free module homomorphisms with codomain `F` represented by `A` and `B`.

    subquotient(A::MatElem{T}, B::MatElem{T}) where T

Given matrices `A` and `B` with the same number of columns, create a free module `F` whose rank 
is that number, and return $(\text{im } a + \text{im } b)/\text{im } b$, where `a` and `b` are 
free module homomorphisms with codomain `F` represented by `A` and `B`.

# Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
```
"""
function subquotient(a::FreeModuleHom, b::FreeModuleHom)
  F = codomain(a)
  @assert F === codomain(b)
  A = matrix(a)
  B = matrix(b)
  return SubQuo(F, A, B)
end
subquotient(F::FreeMod{T}, A::MatElem{T}, B::MatElem{T}) where {T} = SubQuo(F, A, B)
subquotient(A::MatElem{T}, B::MatElem{T}) where {T} = SubQuo(A, B)
#######################################################

function show(io::IO, SQ::SubQuo)
  @show_name(io, SQ)
  @show_special(io, SQ)

  if isdefined(SQ, :quo)
    print(io, "Subquotient of ", SQ.sub, "\nby ", SQ.quo)
  else
    #println(io, "Subquotient by ", SQ.sub)
    print(io, SQ.sub, "\n")
    print(io, "represented as subquotient with no relations.")
  end
end

@doc Markdown.doc"""
    show_subquo(SQ::SubQuo)

Show `SQ` as a subquotient of *matrices* `A` and `B`.
"""
function show_subquo(SQ::SubQuo)
  #@show_name(io, SQ)
  #@show_special(io, SQ)

  if isdefined(SQ, :quo)
    if is_generated_by_standard_unit_vectors(SQ.sub)
      println("Cokernel of")
      display(generator_matrix(SQ.quo))
    else
      println("Subquotient with of image of")
      display(generator_matrix(SQ.sub))
      println("by image of")
      display(generator_matrix(SQ.quo))
    end
  else
    println("Image of")
    display(generator_matrix(SQ.sub))
  end
end

function show_morphism_as_map(f::ModuleFPHom, print_non_zero_only = false)
  D = domain(f)
  for i in 1:ngens(D)
    generator = gen(D, i)
    im = f(generator)
    if !(print_non_zero_only && iszero(im))
      print(generator, " -> ", f(generator))
      if i < ngens(D)
        print("\n")
      end
    end
  end
  return
end

@doc Markdown.doc"""
    cokernel(a::ModuleFPHom)

Return the cokernel of `a` as an object of type `SubQuo`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> W = matrix(R, [y 0; x y; 0 z])
[y   0]
[x   y]
[0   z]

julia> a = hom(F, G, W);

julia> cokernel(a)
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 3 generators
1 -> y*e[1]
2 -> x*e[1] + y*e[2]
3 -> z*e[2]
```

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F1 = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F1, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]]
2-element Vector{SubQuoElem{fmpq_mpoly}}:
 x*y^2*e[1]
 x*y*e[1]

julia> b = hom(M, N, V);

julia> cokernel(b)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 5 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
4 -> x*y^2*e[1]
5 -> x*y*e[1]
```
"""
function cokernel(f::ModuleFPHom)
  return quo(codomain(f), image(f)[1], :module)
end

function cokernel(F::FreeMod{R}, A::MatElem{R}) where R
  return cokernel(map(F,A))
end

function cokernel(A::MatElem)
  return cokernel(map(A))
end

@doc Markdown.doc"""
    default_ordering(M::SubQuo)    

Return the default ordering of `M`.
"""
function default_ordering(M::SubQuo)
  if !isdefined(M.sub, :default_ordering)
    ord = default_ordering(ambient_free_module(M))
    set_default_ordering!(M, ord)
    return default_ordering(M.sub)
  end
  return default_ordering(M.sub)
end

@doc Markdown.doc"""
    set_default_ordering!(M::SubQuo, ord::ModuleOrdering)

Set the default ordering in `M` to `ord`.
"""
function set_default_ordering!(M::SubQuo, ord::ModuleOrdering)
  set_default_ordering!(M.sub, ord)
  if isdefined(M, :quo)
    set_default_ordering!(M.quo, ord)
  end
  set_default_ordering!(M.sum, ord)
end

function std_basis(M::SubQuo, ord::ModuleOrdering = default_ordering(M))
  if !haskey(M.groebner_basis, ord)
    if isdefined(M, :quo)
      quo_gb = std_basis(M.quo, ord)
      sub_union_gb_of_quo = SubModuleOfFreeModule(M.F, ModuleGens(vcat(M.sub.gens.O, quo_gb.O), M.F))
      gb = compute_std_basis(sub_union_gb_of_quo, ord)
      rel_gb_list = Vector{elem_type(ambient_free_module(M))}()

      for i in 1:length(gb.O)
        v = gb.S[i]
        if !iszero(_reduce(v, quo_gb.S))
          push!(rel_gb_list, gb.O[i])
        end
      end

      rel_gb_ModuleGens = ModuleGens(rel_gb_list, M.F)
      rel_gb_ModuleGens.isGB = true
      rel_gb_ModuleGens.ordering = ord
      rel_gb_ModuleGens.quo_GB = quo_gb
      M.groebner_basis[ord] = rel_gb_ModuleGens
    else
      M.groebner_basis[ord] = std_basis(M.sub, ord)
    end
  end
  return M.groebner_basis[ord]
end

function groebner_basis(M::SubQuo, ord::ModuleOrdering = default_ordering(M))
  @assert is_global(ord)
  return std_basis(M, ord)
end

function reduced_groebner_basis(M::SubQuo, ord::ModuleOrdering = default_ordering(M))
  if !haskey(M.groebner_basis, ord) || !(M.groebner_basis[ord].is_reduced)
    if !isdefined(M, :quo)
      M.groebner_basis[ord] = reduced_groebner_basis(M.sub, ord)
    else
      quo_gb = reduced_groebner_basis(M.quo, ord)
      sub_union_gb_of_quo = SubModuleOfFreeModule(M.F, ModuleGens(vcat(M.sub.gens.O, quo_gb.O), M.F))
      gb = reduced_groebner_basis(sub_union_gb_of_quo, ord)
      rel_gb_list = Vector{elem_type(ambient_free_module(M))}()

      for i in 1:length(gb.O)
        v = gb.S[i]
        if !iszero(_reduce(v, quo_gb.S))
          push!(rel_gb_list, gb.O[i])
        end
      end

      rel_gb_ModuleGens = ModuleGens(rel_gb_list, M.F)
      rel_gb_ModuleGens.isGB = true
      rel_gb_ModuleGens.is_reduced = true
      rel_gb_ModuleGens.quo_GB = quo_gb
      if !haskey(M.groebner_basis, ord)
        M.groebner_basis[ord] = rel_geb_ModuleGens
      else
        set_attribute!(M.groebner_basis[ord], :reduced_groebner_basis => rel_gb_ModuleGens)
        return rel_geb_ModuleGens
      end
    end
  end
  return M.groebner_basis[ord]
end

function leading_module(M::SubQuo, ord::ModuleOrdering = default_ordering(M))
  @assert !isdefined(M, :quo)
  return SubQuo(leading_module(M.sub, ord))
end

@doc Markdown.doc"""
    issubset(M::SubQuo{T}, N::SubQuo{T}) where T

Given subquotients `M` and `N` such that `ambient_module(M) == ambient_module(N)`,
return `true` if `M` is contained in `N`, where `M` and `N` are regarded as submodules 
of the common ambient module.

# Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> AM = R[x;]
[x]

julia> BM = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, AM, BM)
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> AN = R[x; y]
[x]
[y]

julia> BN = R[x^2+y^4; y^3; z^4]
[x^2 + y^4]
[      y^3]
[      z^4]

julia> N = SubQuo(F, AN, BN)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> (x^2 + y^4)*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> issubset(M, N)
true
```
"""
function issubset(M::SubQuo{T}, N::SubQuo{T}) where T
  if !isdefined(M, :quo) 
    if !isdefined(N, :quo)
      return issubset(M.sub, N.sub)
    else
      return iszero(N.quo) && issubset(M.sub, N.sub)
    end
  else
    if !isdefined(N, :quo)
      return iszero(M.quo) && issubset(M.sub, N.sub)
    else
      return M.quo == N.quo && issubset(M.sum, N.sum)
    end
  end
end

function compare_helper(M::SubQuo{T}, N::SubQuo{T}, comparer::Function) where {T}
  if !isdefined(M, :quo) 
    if !isdefined(N, :quo)
      return comparer(M.sub, N.sub)
    else
      return iszero(N.quo) && comparer(M.sub, N.sub)
    end
  else
    if !isdefined(N, :quo)
      return iszero(M.quo) && comparer(M.sub, N.sub)
    else
      return comparer(M.quo, N.quo) && comparer(M.sum, N.sum)
    end
  end
end

@doc Markdown.doc"""
    ==(M::SubQuo{T}, N::SubQuo{T}) where {T}

Given subquotients `M` and `N` such that `ambient_module(M) == ambient_module(N)`,
return `true` if `M` equals `N`, where `M` and `N` are regarded as submodules 
of the common ambient module.

Here, `ambient_module(M) == ambient_module(N)` if

- `ambient_free_module(M) === ambient_free_module(N)`, and
- the submodules of the common ambient free module generated by the relations of `M` and `N`, respectively, are equal.

# Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> AM = R[x;]
[x]

julia> BM = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, AM, BM)
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> AN = R[x; y]
[x]
[y]

julia> BN = R[x^2+y^4; y^3; z^4]
[x^2 + y^4]
[      y^3]
[      z^4]

julia> N = SubQuo(F, AN, BN)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> (x^2 + y^4)*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> M == N
false
```
"""
function (==)(M::SubQuo{T}, N::SubQuo{T}) where {T} # TODO replace implementation by two inclusion checks?
  return compare_helper(M, N, (==))
end

@doc Markdown.doc"""
    is_canonically_isomorphic(M::SubQuo{T}, N::SubQuo{T}) where {T}

Check if `M` and `N` are isomorphic under `canonical_isomorphism(F, G)` where
`F` and `G` are the ambient free modules of `M` and `N` respectively.
Return `false` if the ambient free modules are not isomorphic.
"""
function is_canonically_isomorphic(M::SubQuo{T}, N::SubQuo{T}) where {T}
  F = ambient_free_module(M)
  G = ambient_free_module(N)
  if !is_isomorphic(F, G)
    return false
  end
  return compare_helper(M, N, is_canonically_isomorphic)
end

@doc Markdown.doc"""
    is_canonically_isomorphic_with_map(M::SubQuo{T}, N::SubQuo{T}) where {T}

Check if `M` and `N` are isomorphic under `canonical_isomorphism(F, G)` where
`F` and `G` are the ambient free modules of `M` and `N` respectively.
Moreover, if `M` and `N` are canonically isomorphic then return also the isomorphism, 
otherwise return the zero map.
"""
function is_canonically_isomorphic_with_map(M::SubQuo{T}, N::SubQuo{T}) where {T}
  is_canonically_iso = is_canonically_isomorphic(M, N)
  if is_canonically_iso
    F = ambient_free_module(M)
    G = ambient_free_module(N)
    f = canonical_isomorphism(F, G)
    f = induced_map(f, M, false) # is certainly well-defined
    f = restrict_codomain(f, N)
    return true, f
  else
    return false, hom(M, N, [zero(N) for _ in 1:gens(M)])
  end
end

@doc Markdown.doc"""
    sum(M::SubQuo{T},N::SubQuo{T}) where T

Given subquotients `M` and `N` such that `ambient_module(M) == ambient_module(N)`,
return the sum of `M` and `N` regarded as submodules of the common ambient module.

Additionally, return the inclusion maps `M` $\to$ `M + N` and `N` $\to$ `M + N`.

# Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> AM = R[x;]
[x]

julia> BM = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, AM, BM)
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> AN = R[y;]
[y]

julia> BN = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> N = SubQuo(F, AN, BN)
Subquotient of Submodule with 1 generator
1 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> O = sum(M, N);

julia> O[1]
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> O[2]
Map with following data
Domain:
=======
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
Codomain:
=========
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> O[3]
Map with following data
Domain:
=======
Subquotient of Submodule with 1 generator
1 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
Codomain:
=========
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
```
"""
function sum(M::SubQuo{T},N::SubQuo{T}) where T
  @assert ambient_free_module(M) === ambient_free_module(N)

  if !isdefined(M, :quo) && !isdefined(N, :quo)
    SQ = SubQuo(sum(M.sub,N.sub))
  else

    M_quo = isdefined(M, :quo) ? M.quo : SubModuleOfFreeModule(ambient_free_module(M), Vector{elem_type(ambient_free_module(M))}())
    N_quo = isdefined(N, :quo) ? N.quo : SubModuleOfFreeModule(ambient_free_module(N), Vector{elem_type(ambient_free_module(N))}())
  
    if M_quo == N_quo
      SQ = SubQuo(sum(M.sub,N.sub),M_quo)
    else
      error("Relations are not equal")
    end
  end

  iM = SubQuoHom(M,SQ,[SQ[i] for i=1:ngens(M)])
  iN = SubQuoHom(N,SQ,[SQ[i] for i=ngens(M)+1:ngens(SQ)])

  register_morphism!(iM)
  register_morphism!(iN)

  SQ_simplified, _, s_proj = simplify_light(SQ)
  return SQ_simplified, iM*s_proj, iN*s_proj
end

@doc Markdown.doc"""
    +(M::SubQuo{T},N::SubQuo{T}) where T

Given subquotients `M` and `N` such that `ambient_module(M) == ambient_module(N)`,
return the sum of `M` and `N` regarded as submodules of the common ambient module. 

# Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> AM = R[x;]
[x]

julia> BM = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, AM, BM)
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> AN = R[y;]
[y]

julia> BN = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> N = SubQuo(F, AN, BN)
Subquotient of Submodule with 1 generator
1 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> O = M + N
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
```
"""
function +(M::SubQuo{T},N::SubQuo{T}) where T
  return sum(M,N)[1]
end

@doc Markdown.doc"""
    intersect(M::SubQuo{T}, N::SubQuo{T}) where T

Given subquotients `M` and `N` such that `ambient_module(M) == ambient_module(N)`,
return the intersection of `M` and `N` regarded as submodules of the common ambient module.

Additionally, return the inclusion maps `M` $\cap$ `N` $\to$ `M` and `M` $\cap$ `N` $\to$ `N`.

# Examples

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> AM = R[x;]
[x]

julia> BM = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, AM, BM)
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> AN = R[y;]
[y]

julia> BN = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> N = SubQuo(F, AN, BN)
Subquotient of Submodule with 1 generator
1 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> intersect(M, N)
(Subquotient of Submodule with 3 generators
1 -> x^2*e[1]
2 -> -x*y*e[1]
3 -> x*z^4*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1], Map with following data
Domain:
=======
Subquotient of Submodule with 3 generators
1 -> x^2*e[1]
2 -> -x*y*e[1]
3 -> x*z^4*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
Codomain:
=========
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1], Map with following data
Domain:
=======
Subquotient of Submodule with 3 generators
1 -> x^2*e[1]
2 -> -x*y*e[1]
3 -> x*z^4*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
Codomain:
=========
Subquotient of Submodule with 1 generator
1 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1])
```
"""
function intersect(M::SubQuo{T}, N::SubQuo{T}) where T
  #TODO allow task as argument?
  @assert ambient_free_module(M) === ambient_free_module(N)
  M_quo = isdefined(M, :quo) ? M.quo : SubModuleOfFreeModule(ambient_free_module(M), Vector{elem_type(ambient_free_module(M))}())
  N_quo = isdefined(N, :quo) ? N.quo : SubModuleOfFreeModule(ambient_free_module(N), Vector{elem_type(ambient_free_module(N))}())
  R = base_ring(M)

  if M_quo == N_quo

    F1 = FreeMod(R, ngens(M.sub) + ngens(N.sub) + ngens(M_quo))
    F2 = ambient_free_module(M)
    phi = FreeModuleHom(F1,F2,vcat(gens(M.sub),gens(N.sub),gens(M_quo)))
    K,i = kernel(phi)

    intersection_gens_array_with_zeros = [sum([repres(k)[i]*M.sub[i] for i=1:ngens(M.sub)]) for k in gens(K)]
    iszero_array = map(!iszero, intersection_gens_array_with_zeros)

    intersection_gens = SubModuleOfFreeModule(ambient_free_module(M), intersection_gens_array_with_zeros[iszero_array] )
    SQ = SubQuo(intersection_gens,M_quo)

    m = ngens(M)
    M_hom = SubQuoHom(SQ,M,[sum([repres(k)[i]*M[i] for i=1:m]) for k in gens(K)][iszero_array])
    N_hom = SubQuoHom(SQ,N,[sum([repres(k)[i]*N[i-m] for i=m+1:m+ngens(N)]) for k in gens(K)][iszero_array])

    register_morphism!(M_hom)
    register_morphism!(N_hom)

    SQ_simplified, s_inj, _ = simplify_light(SQ)
    return SQ_simplified, s_inj*M_hom, s_inj*N_hom
  end
  throw(ArgumentError("Relations of M and N are not equal."))
end

###############################################################################
# SubQuoElem constructors
###############################################################################

@doc Markdown.doc"""
    SubQuoElem(v::SRow{R}, SQ::SubQuo) where {R}

Return the element $\sum_i v[i] \cdot SQ[i]$.
"""
SubQuoElem(v::SRow{R}, SQ::SubQuo) where {R} = SubQuoElem{R}(v, SQ)

@doc Markdown.doc"""
    SubQuoElem(a::FreeModElem{R}, SQ::SubQuo) where {R}

Construct an element $v \in SQ$ that is represented by $a$.
"""
SubQuoElem(a::FreeModElem{R}, SQ::SubQuo) where {R} = SubQuoElem{R}(a, SQ) 

elem_type(::SubQuo{T}) where {T} = SubQuoElem{T}
parent_type(::SubQuoElem{T}) where {T} = SubQuo{T}
elem_type(::Type{SubQuo{T}}) where {T} = SubQuoElem{T}
parent_type(::Type{SubQuoElem{T}}) where {T} = SubQuo{T}

function in(v::SubQuoElem, M::SubQuo)
  ambient_free_module(parent(v)) === ambient_free_module(M) || return false
  return represents_element(repres(v), M)
end

@doc Markdown.doc"""
    getindex(v::SubQuoElem, i::Int)

Let $v \in M$ with $v = \sum_i a[i] \cdot M[i]$. Return $a[i]$
"""
function getindex(v::SubQuoElem, i::Int)
  if isempty(coeffs(v))
    return zero(base_ring(v.parent))
  end
  return coeffs(v)[i]
end

@doc Markdown.doc"""
    coeffs(v::SubQuoElem)

Let $v \in M$. Return an `SRow` `a` such that $\sum_i a[i] \cdot M[i] = v$.
"""
function coeffs(v::SubQuoElem)
  return v.coeffs
end

#######################################################
@doc Markdown.doc"""
    coefficients(m::SubQuoElem)

Given an element `m` of a subquotient $M$ over a ring $R$, say,
return the coefficients of an $R$-linear combination of the generators of $M$
which gives $m$.

Return the coefficients of `m` with respect to the basis of standard unit vectors. 

The result is returned as a sparse row.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubQuo(A, B);

julia> m = z*M[1] + M[2]
(x*z + y)*e[1]

julia> coefficients(m)
Sparse row with positions [1, 2] and values fmpq_mpoly[z, 1]
```
"""
coefficients(m::SubQuoElem) = coeffs(m)
#########################################################

@doc Markdown.doc"""
    repres(v::SubQuoElem)

Return a free module element that is a representative of `v`.
"""
function repres(v::SubQuoElem)
  return v.repres
end

#######################################################
@doc Markdown.doc"""
    ambient_representative(m::SubQuoElem)

Given an element `m` of a subquotient $M$, say, return 

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubQuo(A, B);

julia> m = z*M[1] + M[2]
(x*z + y)*e[1]

julia> typeof(m)
SubQuoElem{fmpq_mpoly}

julia> fm = ambient_representative(m)
(x*z + y)*e[1]

julia> typeof(fm)
FreeModElem{fmpq_mpoly}

julia> parent(fm) == ambient_free_module(M)
true
```
"""
ambient_representative(m::SubQuoElem) = repres(m)
#######################################################

@doc Markdown.doc"""
    Vector(v::SubQuoElem)

Return the coefficients of a representative of `v` as a Vector.
"""
function Vector(v::SubQuoElem)
  return Vector(repres(v))
end

@doc Markdown.doc"""
    std_basis(F::ModuleGens{T}, reduced::Bool=false) where {T <: MPolyElem}

Return a standard basis of `F` as an object of type `ModuleGens`.
If `reduced` is set to `true` and the ordering of the underlying ring is global, 
a reduced Gröbner basis is computed.
"""
function std_basis(F::ModuleGens{T}, reduced::Bool=false) where {T <: MPolyElem}
  singular_assure(F)
  if reduced
    @assert Singular.has_global_ordering(base_ring(F.SF))
  end
  if singular_generators(F).isGB && !reduced
    return F
  end
  return ModuleGens(F.F, Singular.std(singular_generators(F), complete_reduction=reduced))
end

@doc Markdown.doc"""
    lift_std(M::ModuleGens{T}) where {T <: MPolyElem}

Return a standard basis `G` of `F` as an object of type `ModuleGens` along with 
a transformation matrix `T` such that `T*matrix(M) = matrix(G)`.
"""
function lift_std(M::ModuleGens{T}) where {T <: MPolyElem}
  singular_assure(M)
  R = base_ring(M)
  G,Trans_mat = Singular.lift_std(singular_generators(M)) # When Singular supports reduction add it also here
  mg = ModuleGens(M.F, G)
  mg.isGB = true
  mg.S.isGB = true
  mg.ordering = default_ordering(M.F)
  mat = map_entries(R, transpose(Trans_mat))
  set_attribute!(mg, :transformation_matrix => mat)
  return mg, mat
end

@doc Markdown.doc"""
    lift_std(M::ModuleGens{T}, ordering::ModuleOrdering) where {T <: MPolyElem}

Return a standard basis `G` of `F` with respect to the given `ordering`
as an object of type `ModuleGens` along with a transformation 
matrix `T` such that `T*matrix(M) = matrix(G)`.
"""
function lift_std(M::ModuleGens{T}, ordering::ModuleOrdering) where {T <: MPolyElem}
  M = ModuleGens(M.O, M.F, ordering)
  mg, mat = lift_std(M)
  mg.ordering = ordering
  return mg, mat
end

@doc Markdown.doc"""
    leading_monomials(F::ModuleGens)

Return the leading monomials of `F` as an object of type `ModuleGens`.
The leading module is with respect to the ordering defined on the Singular side.
"""
function leading_monomials(F::ModuleGens)
  # TODO
  # The following doesn't work yet. When comparison / lead for module elements
  # is implemented this should be uncommented.
  #if !isdefined(F, :S)
    #return ModuleGens(F.F, [lead(g) for g in F.O])
  #end
  singular_assure(F)
  singular_gens = singular_generators(F)
  return ModuleGens(F.F, Singular.lead(singular_gens))
end

function show(io::IO, b::SubQuoElem)
  print(io, b.repres)
end

@doc Markdown.doc"""
    parent(b::SubQuoElem)

Let $b \in M$. Return $M$.
"""
parent(b::SubQuoElem) = b.parent

@doc Markdown.doc"""
    (M::SubQuo{T})(f::FreeModElem{T}) where T

Given an element `f` of the ambient free module of `M` which represents an element of `M`, 
return the represented element.
"""
function (M::SubQuo{T})(f::FreeModElem{T}) where T
  coords = coordinates(f, M)
  if coords === nothing
    error("not in the module")
  end
  return SubQuoElem(coords, M)
end

@doc Markdown.doc"""
    (M::SubQuo{T})(c::SRow{T}) where T   

Return the subquotient element $\sum_i a[i] \cdot M[i]\in M$.
"""
function (R::SubQuo)(a::SRow)
  return SubQuoElem(a, R)
end

@doc Markdown.doc"""
    SubQuoElem(c::Vector{T}, parent::SubQuo{T}) where T
    
Return the element of  `parent`  defined as a linear combination
of the generators of $parent$ with coefficients given by the entries of `c`.
"""
function SubQuoElem(c::Vector{T}, parent::SubQuo{T}) where T
  @assert length(c) == ngens(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:ngens(parent)), c)
  return SubQuoElem{T}(sparse_coords,parent)
end


@doc Markdown.doc"""
    (M::SubQuo{T})(c::Vector{T}) where T

Return the element of `M` defined as a linear combination
of the generators of $M$ with coefficients given by the entries of `c`.
"""
function (M::SubQuo{T})(c::Vector{T}) where T
 return SubQuoElem(c, M)
end

@doc Markdown.doc"""
    (R::SubQuo)(a::SubQuoElem)

Return `a` if it lives in `R`.
"""
function (R::SubQuo)(a::SubQuoElem)
  if parent(a) == R
    return a
  end
  error("illegal coercion")
end

@doc Markdown.doc"""
    index_of_gen(v::SubQuoElem)

Let $v \in G$ with $v$ the `i`th generator of $G$. Return `i`.
"""
function index_of_gen(v::SubQuoElem)
  @assert length(coeffs(v).pos) == 1
  @assert isone(coeffs(v).values[1])
  return coeffs(v).pos[1]
end

# function to check whether a free module element is in a particular free module
function check_parent(a::Union{AbstractFreeModElem,SubQuoElem}, b::Union{AbstractFreeModElem,SubQuoElem})
  if parent(a) !== parent(b)
    error("elements not compatible")
  end  
end

function +(a::SubQuoElem, b::SubQuoElem)
  check_parent(a,b)
  return SubQuoElem(coeffs(a)+coeffs(b), a.parent)
end
function -(a::SubQuoElem, b::SubQuoElem) 
  check_parent(a,b)
  return SubQuoElem(coeffs(a)-coeffs(b), a.parent)
end
-(a::SubQuoElem) = SubQuoElem(-coeffs(a), a.parent)
function *(a::MPolyElem_dec, b::SubQuoElem) 
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  return SubQuoElem(a*coeffs(b), b.parent)
end
function *(a::MPolyElem, b::SubQuoElem) 
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  return SubQuoElem(a*coeffs(b), b.parent)
end
function *(a::RingElem, b::SubQuoElem) 
  if parent(a) !== base_ring(parent(b))
    return base_ring(parent(b))(a)*b # this will throw if conversion is not possible
  end
  return SubQuoElem(a*coeffs(b), b.parent)
end
*(a::Int, b::SubQuoElem) = SubQuoElem(a*coeffs(b), b.parent)
*(a::Integer, b::SubQuoElem) = SubQuoElem(a*coeffs(b), b.parent)
*(a::fmpq, b::SubQuoElem) = SubQuoElem(a*coeffs(b), b.parent)
function (==)(a::SubQuoElem, b::SubQuoElem) 
  if parent(a) !== parent(b)
    return false
  end
  return iszero(a-b)
end

@doc Markdown.doc"""
    sub(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of elements of `F`, return the submodule of `F` generated by these elements.

Put more precisely, return the submodule as an object of type `SubQuo`. 

Additionally, if `N` denotes this object, 

- return the inclusion map `N` $\to$ `F` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `F` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T
  s = SubQuo(F, V)
  emb = hom(s, F, V)
  set_attribute!(s, :canonical_inclusion => emb)
  return return_sub_wrt_task(s, emb, task)
end

@doc Markdown.doc"""
    sub(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T} 

Given a matrix `A`, return the submodule of `F` generated by the rows of `A`.

Put more precisely, return this submodule as an object of type `SubQuo`. 

Additionally, if `N` denotes this submodule, 

- return the inclusion map `N` $\to$ `F` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `F` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T}
  M = SubQuo(SubModuleOfFreeModule(F, A)) 
  #M = SubQuo(F, A, zero(MatrixSpace(base_ring(F), 1, rank(F))))
  emb = hom(M, F, ambient_representatives_generators(M))
  emb.matrix = A
  set_attribute!(M, :canonical_inclusion => emb)
  return return_sub_wrt_task(M, emb, task)
end

@doc Markdown.doc"""
    sub(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T

Return `S` as a submodule of `F`, where `S` is generated by `O`.
The embedding module of the parent of the elements of `O` must be `F`.
If `task` is set to `:none` or to `:module` return only `S`.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the canonical injection morphism
$S \to F$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function sub(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T
  s = SubQuo(F, [x.repres for x = O])
  return sub(F, s, task)
end

@doc Markdown.doc"""
    sub(F::FreeMod{T}, s::SubQuo{T}, task::Symbol = :with_morphism) where T

Return `s` as a submodule of `F`, that is the embedding free module of `s` must 
be `F` and `s` has no relations.
If `task` is set to `:none` or to `:module` return only `s`.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the canonical injection morphism
$s \to F$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function sub(F::FreeMod{T}, s::SubQuo{T}, task::Symbol = :with_morphism) where T
  @assert !isdefined(s, :quo)
  @assert s.F === F
  emb = hom(s, F, Vector{elem_type(F)}([repres(x) for x in gens(s)]))
  #emb = hom(s, F, [FreeModElem(x.repres.coords, F) for x in gens(s)])
  set_attribute!(s, :canonical_inclusion => emb)
  return return_sub_wrt_task(s, emb, task)
end

@doc Markdown.doc"""
    sub(M::SubQuo{T}, V::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of elements of `M`, return the submodule of `M` generated by these elements.

Put more precisely, return this submodule as an object of type `SubQuo`. 

Additionally, if `N` denotes this object,

- return the inclusion map `N` $\to$ `M` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `M` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(M::SubQuo{T}, V::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T
  @assert all(x -> x.parent === M, V)
  t = SubQuo(M.F, FreeModElem[repres(x) for x in V])
  if isdefined(M, :quo)
    t.quo = M.quo
    t.sum = sum(t.sub, t.quo)
  end
  emb = hom(t, M, V)
  set_attribute!(t, :canonical_inclusion => emb)
  return return_sub_wrt_task(t, emb, task)
end

@doc Markdown.doc"""
    sub(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of elements of `M`, return the submodule of `M` generated by these elements.

Put more precisely, return this submodule as an object of type `SubQuo`. 

Additionally, if `N` denotes this object,

- return the inclusion map `N` $\to$ `M` if `task = :with_morphism` (default),
- return and cache the inclusion map `N` $\to$ `M` if `task = :cache_morphism`,
- do none of the above if `task = :none`.

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, task::Symbol = :with_morphism) where T
 error("sub is not implemented for the given types.")
end

@doc Markdown.doc"""
    return_sub_wrt_task(M::SubQuo, emb::SubQuoHom, task::Symbol)

This helper function returns `M`, `emb` or both according
to `task`.
"""
function return_sub_wrt_task(M::SubQuo, emb::SubQuoHom, task::Symbol)
  (task == :none || task == :module) && return M
  task == :cache_morphism && register_morphism!(emb)
  task == :only_morphism && return emb 
  (task == :cache_morphism || task == :both || task == :with_morphism) && return M, emb
  error("No valid option for task.")
end

@doc Markdown.doc"""
    quo(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of elements of `F`, return the quotient of `F` by the submodule of `F` which is generated by these elements.

Put more precisely, return this quotient as an object of type `SubQuo`. 

Additionally, if `N` denotes this object,

- return the projection map `F` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `F` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(F::FreeMod{T}, V::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T
  S = SubQuo(F, basis(F))
  Q = SubQuo(S, V)

  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
    quo(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T}

Given a vector `A`, return the quotient of `F` by the submodule of `F` which is generated by 
the rows of `A`.

Put more precisely, return this quotient as an object of type `SubQuo`. 

Additionally, if `N` denotes this object,

- return the projection map `F` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `F` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(F::FreeMod{T}, A::MatElem{T}, task::Symbol = :with_morphism) where {T}
  E = one(MatrixSpace(base_ring(F), rank(F), rank(F)))
  Q = SubQuo(F, E, A)

  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
    quo(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T

Compute $F / T$, where $T$ is generated by $O$.
The embedding free module of the parent of the elements of `O` must be `F`.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the 
canonical projection morphism $F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
If `task` is set to `:none` or `:module` return only the module.
"""
function quo(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T
  S = SubQuo(F, basis(F))
  Q = SubQuo(S, [x.repres for x = O])
  
  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
    quo(F::SubQuo{T}, O::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T

Compute $F / T$, where $T$ is generated by $O$.
The elements of `O` must be elements of the embedding free module of `S`.
If `task` is set to `:with_morphism` (default) or to `:both` return also the 
canonical projection morphism $F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
If `task` is set to `:none` or `:module` return only the module.
"""
function quo(F::SubQuo{T}, O::Vector{<:FreeModElem{T}}, task::Symbol = :with_morphism) where T
  if length(O) > 0
    @assert parent(O[1]) === F.F
  end
  if isdefined(F, :quo)
    oscar_assure(F.quo.gens)
    singular_assure(F.quo.gens)
    s = Singular.Module(base_ring(F.quo.gens.SF), [F.quo.gens.SF(x) for x = [O; oscar_generators(F.quo.gens)]]...)
    Q = SubQuo(F.F, singular_generators(F.sub.gens), s)
    return return_quo_wrt_task(F, Q, task)
  end
  Q = SubQuo(F, O)
  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
    quo(M::SubQuo{T}, V::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of elements of `M`, return the quotient of `M` by the submodule of `M` which is generated by these elements.

Put more precisely, return the quotient as an object of type `SubQuo`. 

Additionally, if `N` denotes this object,

- return the projection map `M` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(M::SubQuo{T}, V::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T
  return quo(M, [x.repres for x = V], task)
end

@doc Markdown.doc"""
    quo(M::SubQuo{T}, V::Vector{<:SubQuoElem{T}}, task::Symbol = :with_morphism) where T

Given a vector `V` of elements of `M`, return the quotient of `M` by the submodule of `M` which is generated by these elements.

Put more precisely, return the quotient as an object of type `SubQuo`. 

Additionally, if `N` denotes this object,

- return the projection map `M` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(M::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, task::Symbol = :with_morphism) where T
 error("quo is not implemented for the given types.")
end


@doc Markdown.doc"""
    quo(M::SubQuo{T}, U::SubQuo{T}, task::Symbol = :with_morphism) where T

Return the quotient $M / U$.

Put more precisely, if `N` denotes this quotient, return `N` as an object of type `SubQuo`. Additionally,

- return the projection map `M` $\to$ `N` if `task = :with_morphism` (default),
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` or `task = :module`.

If `task = :only_morphism`, return only the projection map.
"""
function quo(M::SubQuo{T}, U::SubQuo{T}, task::Symbol = :with_morphism) where T
  if isdefined(M, :quo) && isdefined(U, :quo)
    @assert M.quo == U.quo
  else
    @assert !isdefined(M, :quo) && !isdefined(U, :quo)
  end
  Q = SubQuo(M, oscar_generators(U.sub.gens))
  return return_quo_wrt_task(M, Q, task)
end

@doc Markdown.doc"""
    quo(F::FreeMod{R}, T::SubQuo{R}, task::Symbol = :with_morphism) where R

Compute $F / T$.
If `task` is set to `:with_morphism` (default option) or to `:both` return also the 
canonical projection morphism $F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
If `task` is set to `:none` or `:module` return only the module.
"""
function quo(F::FreeMod{R}, T::SubQuo{R}, task::Symbol = :with_morphism) where R
  @assert !isdefined(T, :quo)
  return quo(F, gens(T), task)
end

@doc Markdown.doc"""
    return_quo_wrt_task(M::ModuleFP, Q::ModuleFP, task)

This helper function returns the module `Q = M / N` for some `N` 
along with the canonical projection morphism $M \to Q$ according to the 
given `task`.
"""
function return_quo_wrt_task(M::ModuleFP, Q::ModuleFP, task)
  if task == :none || task == :module
    return Q
  else
    pro = hom(M, Q, gens(Q))
    task == :cache_morphism && register_morphism!(pro)
    task == :only_morphism && return pro
    return Q, pro
  end
end

@doc Markdown.doc"""
    syzygy_module(F::ModuleGens; sub = FreeMod(base_ring(F.F), length(oscar_generators(F))))
"""
function syzygy_module(F::ModuleGens{T}; sub = FreeMod(base_ring(F.F), length(oscar_generators(F)))) where {T <: MPolyElem}
  singular_assure(F)
  # TODO Obtain the Gröbner basis and cache it
  s = Singular.syz(singular_generators(F))
  return SubQuo(sub, s)
end

@doc Markdown.doc"""
    gens(M::SubQuo{T}) where T

Return the generators of `M`.
"""
function gens(M::SubQuo{T}) where T
  return SubQuoElem{T}[gen(M,i) for i=1:ngens(M)]
end

@doc Markdown.doc"""
    gen(M::SubQuo{T}, i::Int) where T

Return the `i`th generator of `M`.
"""
function gen(M::SubQuo{T}, i::Int) where T
  R = base_ring(M)
  v::SRow{T} = sparse_row(R)
  v.pos = [i]
  v.values = [R(1)]
  return SubQuoElem{T}(v, M)
end

@doc Markdown.doc"""
    ngens(M::SubQuo)

Return the number of generators of `M`.
"""
ngens(M::SubQuo) = ngens(M.sub)

@doc Markdown.doc"""
    base_ring(M::SubQuo)

Given an `R`-module `M`, return `R`.
"""
base_ring(M::SubQuo) = base_ring(M.F)

@doc Markdown.doc"""
    zero(M::SubQuo)

Return the zero element of `M`.
"""
zero(M::SubQuo) = SubQuoElem(SRow(base_ring(M)), M)

@doc Markdown.doc"""
    iszero(M::SubQuo)

Return `true` if `M` is the zero module, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x^2+y^2;]
[x^2 + y^2]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, A, B)
Subquotient of Submodule with 1 generator
1 -> (x^2 + y^2)*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> iszero(M)
false
```
"""
function iszero(M::SubQuo)
  return all(iszero, gens(M))
end

@doc Markdown.doc"""
    getindex(F::SubQuo, i::Int)

Return the `i`th generator of `F`.
"""
function getindex(F::SubQuo, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

function iterate(F::ModuleGens, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end
eltype(::ModuleGens{T}) where {T} = FreeModElem{T} 

#??? A scalar product....
function *(a::FreeModElem, b::Vector{FreeModElem})
  @assert dim(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) = a.coords
    s += v*b[p]
  end
  return s
end

####################
### Chain Complexes
####################

@doc Markdown.doc"""
    chain_complex(V::ModuleFPHom...; start::Int = 0)

Given a tuple `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the chain complex defined by these homomorphisms.

    chain_complex(V::Vector{<:ModuleFPHom}; start::Int = 0)

Given a vector `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the chain complex defined by these homomorphisms.

!!! note
    The integer `start` indicates the lowest homological degree of a module in the complex.

!!! note
    The function checks whether succesive homomorphisms indeed compose to zero.
"""
function chain_complex(V::ModuleFPHom...; start::Int = 0)
  return ChainComplex(collect(V); direction = :left, start = start)
end

function chain_complex(V::Vector{<:ModuleFPHom}; start::Int = 0)
  return ChainComplex(V; direction = :left, start = start)
end

####################

####################
### Cochain Complexes
####################

@doc Markdown.doc"""
    cochain_complex(V::ModuleFPHom...; start::Int = 0)

Given a tuple `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the cochain complex defined by these homomorphisms.

    cochain_complex(V::Vector{<:ModuleFPHom}; start::Int = 0)

Given a vector `V` of module homorphisms between successive modules over a multivariate polynomial ring, 
return the cochain complex defined by these homomorphisms.

!!! note
    The integer `start` indicates the lowest homological degree of a module of the complex.

!!! note
    The function checks whether succesive homomorphisms indeed compose to zero.
"""
function cochain_complex(V::ModuleFPHom...; start::Int = 0)
  return ChainComplex(collect(V); direction = :right, start = start)
end

function cochain_complex(V::Vector{<:ModuleFPHom}; start::Int = 0)
  return ChainComplex(V; direction = :right, start = start)
end

####################

@doc Markdown.doc"""
    presentation(M::SubQuo)
ive
Return a free presentation of `M`. 
"""
function presentation(SQ::SubQuo)
  #A+B/B is generated by A and B
  #the relations are A meet B? written wrt to A
  R = base_ring(SQ)
  F = FreeMod(R, ngens(SQ.sub))
  set_attribute!(F,  :name => "br^$(ngens(SQ.sub))")
  q = elem_type(F)[]
  if is_generated_by_standard_unit_vectors(SQ.sub)
    if isdefined(SQ, :quo)
      q = [FreeModElem(coords(g), F) for g in gens(SQ.quo)]
    end
  else
    s, _ = kernel(hom(FreeMod(R,ngens(SQ.sum)), ambient_free_module(SQ), gens(SQ.sum)))
    #s = syzygy_module(SQ.sum.gens)
    #TODO: wait for Hans to release Modulo(A, B) that does exactly this
    c = collect(s.sub.gens)
    #q = elem_type(F)[]

    for x = c
      b = sparse_row(R)
      e = zero(SQ.F)
      for (i,v) = x.coords
        if i>ngens(SQ)
          break
        end
        e += v*gen(SQ, i).repres
        push!(b.pos, i)
        push!(b.values, v)
      end
      if length(b) == 0
        continue
      end
      push!(q, FreeModElem(b, F))
    end
  end
  #want R^a -> R^b -> SQ -> 0
  #TODO sort decoration and fix maps, same decoration should be bundled (to match pretty printing)
  G = FreeMod(R, length(q))
  set_attribute!(G, :name => "br^$(length(q))")
  h_G_F = hom(G, F, q)
  h_F_SQ = hom(F, SQ, gens(SQ)) # DO NOT CHANGE THIS LINE, see present_as_cokernel and preimage

  Z = FreeMod(F.R, 0)
  set_attribute!(Z, :name => "0")
  h_SQ_Z = hom(SQ, Z, Vector{elem_type(Z)}([zero(Z) for i=1:ngens(SQ)]))
  return Hecke.ChainComplex(ModuleFP, ModuleFPHom[h_G_F, h_F_SQ, h_SQ_Z], check = false, start=-2)
end

@doc Markdown.doc"""
    presentation(F::FreeMod)

Return a free presentation of $F$.
"""
function presentation(F::FreeMod)
  Z = FreeMod(F.R, 0)
  set_attribute!(Z, :name => "0")
  return Hecke.ChainComplex(ModuleFP, ModuleFPHom[hom(Z, F, Vector{elem_type(F)}()), hom(F, F, gens(F)), hom(F, Z, Vector{elem_type(Z)}([zero(Z) for i=1:ngens(F)]))], check = false, start=-2)
end

@doc Markdown.doc"""
    presentation(M::ModuleFP)

Return a free presentation of $M$.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> A = R[x; y];

julia> B = R[x^2; y^3; z^4];

julia> M = SubQuo(A, B);

julia> P = presentation(M)
br^5 ----> br^2 ----> C_-1 ----> 0
where:
	C_-1 = Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
```
"""
function presentation(M::ModuleFP)
 error("presentation is not implemented for the given types.")
end

@doc Markdown.doc"""
    present_as_cokernel(SQ::SubQuo, task::Symbol = :none)

Return a subquotient $M = R^n / im(f) $, i.e. $M = \text{coker}(f)$, such that
$M \cong SQ$.
If `task` is set to `:with_morphism` or to `:both` then return also an isomorphism $M \to SQ$. Calling `inv()`
on this isomorphism is cheap.
If `task` is set to `:cache_morphism` then the isomorphism is cached.
If `task` is set to `:only_morphism` then return only the isomorphism.

# Examples
"""
function present_as_cokernel(SQ::SubQuo, task::Symbol = :none)
  chainComplex = presentation(SQ)
  R_b = obj(chainComplex, 0)
  f = map(chainComplex, 1)
  g = map(chainComplex, 0)
  presentation_module = quo(R_b, image(f)[1], :module)

  if task == :none
    return presentation_module
  end
  
  # The isomorphism is just the identity matrix
  isomorphism = hom(presentation_module, SQ, Vector{elem_type(SQ)}([g(x) for x in gens(R_b)]))
  inverse_isomorphism = hom(SQ, presentation_module, Vector{elem_type(presentation_module)}([presentation_module[i] for i=1:ngens(SQ)]))
  isomorphism.inverse_isomorphism = inverse_isomorphism

  if task == :cache_morphism
    register_morphism!(isomorphism)
    register_morphism!(inverse_isomorphism)
  end
  task == :only_morphism && return isomorphism
  
  return presentation_module, isomorphism
end

@doc Markdown.doc"""
    is_equal_with_morphism(M::SubQuo{T}, N::SubQuo{T}, task::Symbol = :none) where {T}

If $M = N$ (mathematically, but with (possibly) different generating systems), return $\phi : M \to N$ 
which is mathematically the identity. 
If `task == :inverse` also the inverse map is computed and cached (in the morphism).
If `task == :cache_morphism` the inverse map is also cached in `M` and `N`.
"""
function is_equal_with_morphism(M::SubQuo{T}, N::SubQuo{T}, task::Symbol = :none) where {T}
  @assert M == N

  M_to_N = hom(M, N, Vector{elem_type(N)}([SubQuoElem(coordinates(m.repres, N), N) for m in gens(M)]))

  if task == :cache_morphism || task == :inverse
    N_to_M = hom(N, M, Vector{elem_type(M)}([SubQuoElem(coordinates(n.repres, M), M) for n in gens(N)]))
    M_to_N.inverse_isomorphism = N_to_M
    N_to_M.inverse_isomorphism = M_to_N

    if task == :cache_morphism
      register_morphism!(M_to_N) 
      register_morphism!(N_to_M)
    end
  end
  
  return M_to_N
end

###############################################################################
# SubQuoHom constructors
###############################################################################

@doc Markdown.doc"""
    SubQuoHom(D::SubQuo, C::ModuleFP{T}, im::Vector{<:ModuleFPElem{T}}) where T

Return the morphism $D \to C$ for a subquotient $D$ where `D[i]` is mapped to `im[i]`.
In particular, `length(im) == ngens(D)` must hold.
"""
SubQuoHom(D::SubQuo, C::ModuleFP{T}, im::Vector{<:ModuleFPElem{T}}) where {T} = SubQuoHom{typeof(D), typeof(C), Nothing}(D, C, im)
SubQuoHom(D::SubQuo, C::ModuleFP{T}, im::Vector{<:ModuleFPElem{T}}, h::RingMapType) where {T, RingMapType} = SubQuoHom{typeof(D), typeof(C), RingMapType}(D, C, im, h)

@doc Markdown.doc"""
    SubQuoHom(D::SubQuo, C::ModuleFP{T}, mat::MatElem{T})

Return the morphism $D \to C$ corresponding to the given matrix, where $D$ is a subquotient.
`mat` must have `ngens(D)` many rows and `ngens(C)` many columns.
"""
function SubQuoHom(D::SubQuo, C::ModuleFP{T}, mat::MatElem{T}) where T
  @assert nrows(mat) == ngens(D)
  @assert ncols(mat) == ngens(C)
  if C isa FreeMod
    hom = SubQuoHom(D, C, [FreeModElem(sparse_row(mat[i,:]), C) for i=1:ngens(D)])
    return hom
  else
    hom = SubQuoHom(D, C, [SubQuoElem(sparse_row(mat[i,:]), C) for i=1:ngens(D)])
    return hom
  end
end

###################################################################

@doc Markdown.doc"""
    hom(M::SubQuo{T}, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T

Given a vector `V` of `ngens(M)` elements of `N`, 
return the homomorphism `M` $\to$ `N` which sends the `i`-th
generator `M[i]` of `M` to the `i`-th entry of `V`.

    hom(M::SubQuo{T}, N::ModuleFP{T},  A::MatElem{T})) where T

Given a matrix `A` with `ngens(M)` rows and `ngens(N)` columns, return the
homomorphism `M` $\to$ `N` which sends the `i`-th generator `M[i]` of `M` to 
the linear combination $\sum_j A[i,j]*N[j]$ of the generators `N[j]` of `N`.

!!! warning
    The functions do not check whether the resulting homomorphism is well-defined,
    that is, whether it sends the relations of `M` into the relations of `N`. 

If you are uncertain with regard to well-definedness, use the function below.
Note, however, that the check performed by the function requires a Gröbner basis computation. This may take some time.

    is_welldefined(a::ModuleFPHom)

Return `true` if `a` is well-defined, and `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]]
2-element Vector{SubQuoElem{fmpq_mpoly}}:
 x*y^2*e[1]
 x*y*e[1]

julia> a = hom(M, N, V)
Map with following data
Domain:
=======
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
Codomain:
=========
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> is_welldefined(a)
true

julia> W = matrix(R,  [y^2 0; 0 x])
[y^2   0]
[  0   x]

julia> b = hom(M, N, W);

julia> a == b
true
```

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y];

julia> B = R[x^2; y^3; z^4];

julia> M = SubQuo(F, A, B);

julia> N = M;

julia> W = [y*N[1], x*N[2]]
2-element Vector{SubQuoElem{fmpq_mpoly}}:
 x*y*e[1]
 x*y*e[1]

julia> c = hom(M, N, W);

julia> is_welldefined(c)
false
```
"""
hom(M::SubQuo, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}) where T = SubQuoHom(M, N, V) 
hom(M::SubQuo, N::ModuleFP{T}, V::Vector{<:ModuleFPElem{T}}, h::RingMapType) where {RingMapType, T} = SubQuoHom(M, N, V, h) 
hom(M::SubQuo, N::ModuleFP{T},  A::MatElem{T}) where T = SubQuoHom(M, N, A)
function is_welldefined(H::ModuleFPHom)
  if H isa Union{FreeModuleHom,FreeModuleHom_dec}
    return true
  end
  M = domain(H)
  C = present_as_cokernel(M).quo
  n = ngens(C)
  m = rank(C.F)
  ImH = map(x -> H(x), gens(M))
  for i=1:n
    if !iszero(sum([C[i][j]*ImH[j] for j=1:m]))
      return false
    end
  end
  return true
end

function (==)(f::ModuleFPHom, g::ModuleFPHom)
  domain(f) === domain(g) || return false
  codomain(f) === codomain(g) || return false
  M = domain(f)
  for v in gens(M)
    f(v) == g(v) || return false
  end
  return true
end
###################################################################

@doc Markdown.doc"""
    matrix(a::SubQuoHom)

Given a homomorphism `a` of type  `SubQuoHom` with domain `M`
and codomain `N`, return a matrix `A` with `ngens(M)` rows and 
`ngens(N)` columns such that `a == hom(M, N, A)`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]];

julia> a = hom(M, N, V);

julia> A = matrix(a)
[y^2   0]
[  0   x]

julia> a(M[1])
x*y^2*e[1]
```
"""
function matrix(f::SubQuoHom)
  if !isdefined(f, :matrix)
    D = domain(f)
    C = codomain(f)
    R = base_ring(D)
    matrix = zero_matrix(R, ngens(D), ngens(C))
    for i=1:ngens(D), j=1:ngens(C)
      matrix[i,j] = f.im[i][j]
    end
    f.matrix = matrix
  end
  return f.matrix
end

function show_morphism(f::ModuleFPHom)
  display(matrix(f))
end

@doc Markdown.doc"""
    hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{ <: ModuleFPHom})

Given modules `M`, `N` which are tensor products with the same number of factors,
say $M = M_1 \otimes \cdots \otimes M_r$, $N = N_1 \otimes \cdots \otimes N_r$,
and given a vector `V` of homomorphisms $a_i : M_i \to N_i$, return 
$a_1 \otimes \cdots \otimes a_r$.
"""
function hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{ <: ModuleFPHom})
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
  function map_gen(g) # Is there something that generalizes FreeModElem and SubQuoElem?
    g_decomposed = decompose_M(g)
    image_as_tuple = Tuple(f(x) for (f,x) in zip(V,g_decomposed))
    res = pure_N(image_as_tuple)
    return res
  end
  return hom(M, N, Vector{elem_type(N)}(map(map_gen, gens(M))))
end

@doc Markdown.doc"""
    hom_product(M::ModuleFP, N::ModuleFP, A::Matrix{<:ModuleFPHom})

Given modules `M`, `N` which are products with the same number of factors,  
say $M = \prod_{i=1}^r M_i$, $N = \prod_{j=1}^r N_j$, and given a matrix 
`A` of homomorphisms $a_{ij} : M_i \to N_j$, return the homomorphism
$M \rightarrow N$ with $ij$-components $a_{ij}$.
"""
function hom_product(M::ModuleFP, N::ModuleFP, A::Matrix{<:ModuleFPHom})
  tM = get_attribute(M, :direct_product)
  tM === nothing && error("both modules must be direct products")
  tN = get_attribute(N, :direct_product)
  tN === nothing && error("both modules must be direct products")
  @assert length(tM) == size(A, 1) && length(tN) == size(A, 2)
  @assert all(ij -> domain(A[ij[1],ij[2]]) === tM[ij[1]] && codomain(A[ij[1],ij[2]]) === tN[ij[2]], Base.Iterators.ProductIterator((1:size(A, 1), 1:size(A, 2))))
  #need the canonical maps..., maybe store them as well?
  return hom(M,N,Vector{elem_type(N)}([sum([Hecke.canonical_injection(N,j)(sum([A[i,j](Hecke.canonical_projection(M,i)(g)) for i=1:length(tM)])) for j=1:length(tN)]) for g in gens(M)]))
end
# hom(prod -> X), hom(x -> prod)
# if too much time: improve the hom(A, B) in case of A and/or B are products - or maybe not...
# tensor and hom functors for chain complex
# dual: ambig: hom(M, R) or hom(M, Q(R))?

function lift_with_unit(a::FreeModElem{T}, generators::ModuleGens{T}) where {T <: MPolyElem}
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

function lift(a::FreeModElem{T}, generators::ModuleGens{T}) where {T <: MPolyElem}
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  singular_assure(generators)
  S = singular_generators(generators)
  b = ModuleGens([a], generators.SF)
  singular_assure(b)
  s, r = Singular.lift(S, singular_generators(b))
  if Singular.ngens(s) == 0 || iszero(s[1])
    error("The free module element is not liftable to the given generating system.")
  end
  Rx = base_ring(generators)
  return sparse_row(Rx, s[1], 1:ngens(generators))
end

@doc Markdown.doc"""
    coordinates(a::FreeModElem{T}, generators::ModuleGens{T})

Compute a sparse row `r` such that `a = sum([r[i]*gen(generators,i) for i in 1:ngens(generators)])`.
If no such `r` exists, an execption is thrown.
"""
function coordinates(a::FreeModElem{T}, generators::ModuleGens{T}) where {T <: MPolyElem}
  singular_assure(generators)
  if !Singular.has_global_ordering(base_ring(generators.SF))
    error("Ordering must be global")
  end
  return lift(a, generators)
end

@doc Markdown.doc"""
    coordinates_via_transform(a::FreeModElem{T}, generators::ModuleGens{T}) where T

Let `generators` be a Gröbner basis and let `A*M = generators`.
Compute a sparse row `r` such that `a = sum([r[i]*gen(M,i) for i in 1:ngens(M)])`.
If no such `r` exists, an exception is thrown.
"""
function coordinates_via_transform(a::FreeModElem{T}, generators::ModuleGens{T}) where T
  A = get_attribute(generators, :transformation_matrix)
  A === nothing && error("No transformation matrix in the Gröbner basis.")
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  @assert generators.isGB
  if base_ring(generators) isa Union{MPolyQuo,MPolyRing}
    if !is_global(generators.ordering)
      error("Ordering is not global")
    end
  end

  singular_assure(generators)
  S = singular_generators(generators)
  S.isGB = generators.isGB
  b = ModuleGens([a], generators.SF)
  singular_assure(b)
  s, r = Singular.lift(S, singular_generators(b)) # Possibly use division with remainder
  if Singular.ngens(s) == 0 || iszero(s[1])
    error("The free module element is not liftable to the given generating system.")
  end
  Rx = base_ring(generators)
  coords_wrt_groebner_basis = sparse_row(Rx, s[1], 1:ngens(generators))

  return mul(coords_wrt_groebner_basis,sparse_matrix(A))
end

@doc Markdown.doc"""
    coordinates(a::FreeModElem, M::SubModuleOfFreeModule, task::Symbol = :auto)

Compute a sparse row `r` such that `a = sum([r[i]*gen(M,i) for i in 1:ngens(M)])`.
If no such `r` exists, an exception is thrown.
For `task` there are the following options:
- `:auto` (default option): Use `:via_transform` if coefficient ring of base ring is a field, 
else use `:via_lift` 
- `:via_transform`: Compute first a Gröbner basis with a transformation matrix representing
the Gröbner basis in terms of the generators and cache the data
- `:via_lift`: Compute the lift and do not cache the auxiliary data
Note: `:via_lift` is typically faster than `:via_transform` for a single vector while the latter
is faster if many vectors are lifted
"""
function coordinates(a::FreeModElem, M::SubModuleOfFreeModule, task::Symbol = :auto)
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  if task == :auto
    if base_ring(base_ring(parent(a))) isa Field
      task = :via_transform
    else
      task = :via_lift
    end
  end
  for i in 1:ngens(M)
    g = gen(M,i)
    if a == g
      R = base_ring(M)
      return sparse_row(R, [(i,R(1))])
    end
  end
  if task == :via_transform
    std, _ = lift_std(M)
    return coordinates_via_transform(a, std)
  elseif task == :via_lift
    return coordinates(a, M.gens)
  else
    error("Invalid task given.")
  end
end

function lift_std(M::SubModuleOfFreeModule)
  if haskey(M.groebner_basis, default_ordering(M))
    gb = M.groebner_basis[default_ordering(M)]
    transform = get_attribute(gb, :transformation_matrix)
    if transform !== nothing
      return gb, transform
    end
  end
  for gb in values(M.groebner_basis)
    transform = get_attribute(gb, :transformation_matrix)
    if transform !== nothing
      return gb, transform
    end
  end
  gb, transform = lift_std(M.gens, default_ordering(M))
  M.groebner_basis[default_ordering(M)] = gb
  return gb, transform
end

@doc Markdown.doc"""
    coordinates(a::FreeModElem, SQ::SubQuo, task::Symbol = :auto)

Compute a sparse row `r` such that `a` is a representative of `SubQuoElem(r, SQ)`.
If no such `r` exists, an exception is thrown.
For `task` there are the following options:
- `:auto` (default option): Use `:via_transform` if coefficient ring of base ring is a field, 
else use `:via_lift` 
- `:via_transform`: Compute first a Gröbner basis with a transformation matrix representing
the Gröbner basis in terms of the generators and cache the data
- `:via_lift`: Compute the lift and do not cache the auxiliary data
Note: `:via_lift` is typically faster than `:via_transform` for a single vector while the latter
is faster if many vectors are lifted
"""
function coordinates(a::FreeModElem, SQ::SubQuo, task::Symbol = :auto)
  coords = coordinates(a, SQ.sum, task)
  return coords[1:ngens(SQ)]
end

@doc Markdown.doc"""
    in(a::FreeModElem, M::SubModuleOfFreeModule)

Check if `a` is an element of `M`.
"""
function in(a::FreeModElem, M::SubModuleOfFreeModule)
  F = ambient_free_module(M)
  return iszero(reduce(a, std_basis(M, default_ordering(F))))
end

@doc Markdown.doc"""
    represents_element(a::FreeModElem, SQ::SubQuo)

Check if `a` represents an element `SQ`.
"""
function represents_element(a::FreeModElem, SQ::SubQuo)
  return in(a, SQ.sum)
end

@doc Markdown.doc"""
    normal_form(M::ModuleGens, GB::ModuleGens)

Compute a normal form of `M` (that is of each element of `M`) with respect to the Gröbner basis `GB`.
"""
function normal_form(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyElem}
  @assert M.F === GB.F
  @assert GB.isGB # TODO When Singular.jl can handle reduce with non-GB remove this

  P = isdefined(GB, :quo_GB) ? union(GB, GB.quo_GB) : GB

  singular_assure(P)
  singular_assure(M)

  red = _reduce(M.S, P.S)
  res = ModuleGens(M.F, red)
  oscar_assure(res)
  return res
end

@doc Markdown.doc"""
    normal_form_with_unit(M::ModuleGens, GB::ModuleGens)

Compute a normal form of `M` (that is of each element of `M`) with respect to the Gröbner basis `GB`.
Moreover, return a vector `U` of unit elements such that 
`U[i]*M[i]` is the `i`th element of the normal form `ModuleGens`.
"""
function normal_form_with_unit(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyElem}
  @assert M.F === GB.F
  @assert GB.isGB # TODO When Singular.jl can handle reduce/nf with non-GB remove this
  if !is_global(GB.ordering)
    error("normal_form_with_unit not yet implemented for non-global orderings") # This function doesn't exist yet in Singular.jl
  end
  R = base_ring(M)

  P = isdefined(GB, :quo_GB) ? union(GB, GB.quo_GB) : GB

  singular_assure(P)
  singular_assure(M)

  red = _reduce(M.S, P.S)
  res = ModuleGens(M.F, red)
  oscar_assure(res)
  return res, [R(1) for _ in 1:ngens(M)]
end

function normal_form_with_unit_and_coefficients(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyElem}
  # TODO requires additional functionality in Singular.jl
  error("Not yet implemented")
end

@doc Markdown.doc"""
    normal_form(v::AbstractFreeModElem, GB::ModuleGens)

Compute a normal_form of `v` with respect to the Gröbner basis `GB`.
"""
function normal_form(v::AbstractFreeModElem, GB::ModuleGens)
  @assert GB.isGB
  return normal_form(ModuleGens([v], parent(v)), GB).O[1]
end

@doc Markdown.doc"""
    normal_form_with_unit(v::AbstractFreeModElem, GB::ModuleGens)

Compute a normal form of `v` with respect to the Gröbner basis `GB`.
Moreover, return a unit `u` such that `u*v` is the normal form.
"""
function normal_form_with_unit(v::AbstractFreeModElem, GB::ModuleGens)
  @assert GB.isGB
  if !is_global(GB.ordering)
    error("normal_form_with_unit not yet implemented for non-global orderings") # This function doesn't exist yet in Singular.jl
  end
  red, unit = normal_form_with_unit(ModuleGens([v], parent(v)), GB)
  return red[1], unit[1]
end

@doc Markdown.doc"""
    reduce(M::ModuleGens, GB::ModuleGens)

Reduce `M` with respect to the Gröbner basis `GB`.
"""
function reduce(M::ModuleGens{T}, GB::ModuleGens{T}) where {T <: MPolyElem}
  @assert GB.isGB
  @assert is_global(GB.ordering)
  return normal_form(M, GB)
end

function normal_form(M::SubModuleOfFreeModule{T}, N::SubModuleOfFreeModule{T}) where {T <: MPolyElem}
  @assert is_global(default_ordering(N))
  # TODO reduced flag to be implemented in Singular.jl
  #return SubModuleOfFreeModule(M.F, normal_form(M.gens, std_basis(N), reduced = true))
  error("Not yet implemented")
end

function normal_form(v::AbstractFreeModElem{T}, N::SubModuleOfFreeModule{T}) where {T <: MPolyElem}
  @assert is_global(default_ordering(N))
  # TODO reduced flag to be implemented in Singular.jl
  #return normal_form(v, N.gens, reduced = true)
  error("Not yet implemented")
end

@doc Markdown.doc"""
    reduce(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Reduce `M` with respect to `N`, that is with respect to a Gröbner basis of `N` (the Gröbner basis is computed for this).
"""
function reduce(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  return SubModuleOfFreeModule(M.F, reduce(M.gens, groebner_basis(N, default_ordering(M))))
end

@doc Markdown.doc"""
    reduce(v::FreeModElem, GB::ModuleGens)

Reduce the element `v` with respect to the Gröbner basis `GB`.
"""
function reduce(v::AbstractFreeModElem, GB::ModuleGens)
  @assert GB.isGB
  @assert is_global(GB.ordering)
  return normal_form(v, GB)
end

@doc Markdown.doc"""
    reduce(v::FreeModElem, N::SubModuleOfFreeModule)

Reduce the element `v` with respect to a Gröbner basis of `N`.
"""
function reduce(v::AbstractFreeModElem, N::SubModuleOfFreeModule)
  return reduce(v, groebner_basis(N))
end

@doc Markdown.doc"""
    reduce(v::FreeModElem, N::SubQuo)

Let `N` be a submodule of a free module (that is, `N` has no quotient).
Reduce the element `v` with respect to a Gröbner basis of `N`.
"""
function reduce(v::AbstractFreeModElem, N::SubQuo)
  @assert !(isdefined(N, :quo))
  return reduce(v, groebner_basis(N))
end

@doc Markdown.doc"""
    image(a::SubQuoHom, m::SubQuoElem)

Return the image $a(m)$.
"""
function image(f::SubQuoHom, a::SubQuoElem)
  # TODO matrix vector multiplication
  @assert a.parent === domain(f)
  i = zero(codomain(f))
  b = coeffs(a)
  for (p,v) = b
    i += v*f.im[p]
  end
  return i
end

@doc Markdown.doc"""
    image(f::SubQuoHom, a::FreeModElem)

Return $f(a)$. `a` must represent an element in the domain of `f`.
"""
function image(f::SubQuoHom, a::FreeModElem)
  return image(f, SubQuoElem(a, domain(f)))
end

@doc Markdown.doc"""
    preimage(f::SubQuoHom, a::Union{SubQuoElem,FreeModElem})

Compute a preimage of `a` under `f`.
"""
function preimage(f::SubQuoHom, a::Union{SubQuoElem,FreeModElem})
  @assert parent(a) === codomain(f)
  D = domain(f)
  i = zero(D)
  b = coordinates(a isa FreeModElem ? a : a.repres, image(f)[1])
  for (p,v) = b
    i += v*gen(D, p)
  end
  return i
end

(f::SubQuoHom)(a::FreeModElem) = image(f, a)
(f::SubQuoHom)(a::SubQuoElem) = image(f, a)

@doc Markdown.doc"""
    iszero(m::SubQuoElem)

Return `true` if `m` is zero, `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> iszero(M[1])
false

julia> iszero(x*M[1])
true
```
"""
function iszero(m::SubQuoElem)
  C = parent(m)
  if !isdefined(C, :quo)
    return iszero(repres(m))
  end
  x = reduce(repres(m), C.quo)
  return iszero(x)
end

@doc Markdown.doc"""
    hom(F::FreeMod, G::FreeMod)

Return a subquotient $S$ such that $\text{Hom}(F,G) \cong S$ along with a function 
that converts elements from $S$ into morphisms $F \to G$.
"""
function hom(F::FreeMod, G::FreeMod)
  @assert base_ring(F) === base_ring(G)
  GH = FreeMod(F.R, rank(F) * rank(G))
  GH.S = [Symbol("($i -> $j)") for i = F.S for j = G.S]

  #list is g1 - f1, g2-f1, g3-f1, ...
  X = Hecke.MapParent(F, G, "homomorphisms")
  n = ngens(F)
  m = ngens(G)
  R = base_ring(F)
  function im(x::FreeModElem)
    return hom(F, G, Vector{elem_type(G)}([FreeModElem(x.coords[R, (i-1)*m+1:i*m], G) for i=1:n]))
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
  to_hom_map = Hecke.MapFromFunc(im, pre, GH, X)
  set_attribute!(GH, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map)
  return GH, to_hom_map
end

@doc Markdown.doc"""
    kernel(a::FreeModuleHom)

Return the kernel of `a` as an object of type `SubQuo`.

Additionally, if `K` denotes this object, return the inclusion map `K` $\rightarrow$ `domain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

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
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field)
```
"""
function kernel(h::FreeModuleHom)  #ONLY for free modules...
  G = domain(h)
  R = base_ring(G)
  if ngens(G) == 0
    s = sub(G, gens(G), :module)
    return s, hom(s, G, gens(G))
  end
  g = map(h, basis(G))
  if isa(codomain(h), SubQuo)
    g = [x.repres for x = g]
    if isdefined(codomain(h), :quo)
      append!(g, collect(codomain(h).quo.gens))
    end
  end
  #TODO allow sub-quo here as well
  ambient_free_module_codomain = ambient_free_module(codomain(h))
  b = ModuleGens(g, ambient_free_module_codomain, default_ordering(ambient_free_module_codomain))
  k = syzygy_module(b)
  if isa(codomain(h), SubQuo)
    s = collect(k.sub.gens)
    k = sub(G, [FreeModElem(x.coords[R,1:dim(G)], G) for x = s], :module)
  else
    #the syzygie_module creates a new free module to work in
    k = sub(G, [FreeModElem(x.coords, G) for x = collect(k.sub.gens)], :module)
  end
  @assert k.F === G
  c = collect(k.sub.gens)
  return k, hom(k, parent(c[1]), c)
end

@doc Markdown.doc"""
    image(a::FreeModuleHom)

Return the image of `a` as an object of type `SubQuo`.

Additionally, if `I` denotes this object, return the inclusion map `I` $\rightarrow$ `codomain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

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
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field)
```
"""
function image(h::FreeModuleHom)
  si = [x for x = map(h, basis(domain(h))) if !iszero(x)]
  s = sub(codomain(h), si, :module)
  return s, hom(s, codomain(h), si)
end

@doc Markdown.doc"""
    image(a::SubQuoHom)

Return the image of `a` as an object of type `SubQuo`.

Additionally, if `I` denotes this object, return the inclusion map `I` $\rightarrow$ `codomain(a)`.
"""
function image(h::SubQuoHom)
  h_image_vector::Vector{elem_type(codomain(h))} = h.im
  s = sub(codomain(h), h_image_vector, :module)
  return s, hom(s, codomain(h), h_image_vector)
end

@doc Markdown.doc"""
    image(a::ModuleFPHom)

Return the image of `a` as an object of type `SubQuo`.

Additionally, if `I` denotes this object, return the inclusion map `I` $\rightarrow$ `codomain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> W = matrix(R, [y 0; x y; 0 z])
[y   0]
[x   y]
[0   z]

julia> a = hom(F, G, W);

julia> I, iota = image(a)
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
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field)
```

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F1 = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F1, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]]
2-element Vector{SubQuoElem{fmpq_mpoly}}:
 x*y^2*e[1]
 x*y*e[1]

julia> b = hom(M, N, V);

julia> image(b)
(Subquotient of Submodule with 2 generators
1 -> x*y^2*e[1]
2 -> x*y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1], Map with following data
Domain:
=======
Subquotient of Submodule with 2 generators
1 -> x*y^2*e[1]
2 -> x*y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
Codomain:
=========
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1])
```
"""
function image(a::ModuleFPHom)
 error("kernel is not implemented for the given types.")
end

@doc Markdown.doc"""
    kernel(a::SubQuoHom)

Return the kernel of `a` as an object of type `SubQuo`.

Additionally, if `K` denotes this object, return the inclusion map `K` $\rightarrow$ `domain(a)`.
"""
function kernel(h::SubQuoHom)
  D = domain(h)
  R = base_ring(D)
  F = FreeMod(R, ngens(D))
  hh = hom(F, codomain(h), Vector{elem_type(codomain(h))}(map(h, gens(D))))
  k = kernel(hh)
  @assert domain(k[2]) === k[1]
  @assert codomain(k[2]) === F
  hh = hom(F, D, gens(D))
  im::Vector{elem_type(D)} = filter(x -> !iszero(x), map(x->hh(k[2](x)), gens(k[1])))
  k = sub(D, im, :module)
  return k, hom(k, D, im)
end

@doc Markdown.doc"""
    kernel(a::ModuleFPHom)

Return the kernel of `a` as an object of type `SubQuo`.

Additionally, if `K` denotes this object, return the inclusion map `K` $\rightarrow$ `domain(a)`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> G = free_module(R, 2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> W = matrix(R, [y 0; x y; 0 z])
[y   0]
[x   y]
[0   z]

julia> a = hom(F, G, W);

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
Free module of rank 3 over Multivariate Polynomial Ring in x, y, z over Rational Field)
```

```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> F1 = free_module(R, 1)
Free module of rank 1 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; y^3; z^4]
[x^2]
[y^3]
[z^4]

julia> M = SubQuo(F1, A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]

julia> N = M;

julia> V = [y^2*N[1], x*N[2]]
2-element Vector{SubQuoElem{fmpq_mpoly}}:
 x*y^2*e[1]
 x*y*e[1]

julia> b = hom(M, N, V);

julia> K, iota = kernel(b)
(Subquotient of Submodule with 3 generators
1 -> (-x + y^2)*e[1]
2 -> x*y*e[1]
3 -> -x*y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1], Map with following data
Domain:
=======
Subquotient of Submodule with 3 generators
1 -> (-x + y^2)*e[1]
2 -> x*y*e[1]
3 -> -x*y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1]
Codomain:
=========
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 3 generators
1 -> x^2*e[1]
2 -> y^3*e[1]
3 -> z^4*e[1])
```
"""
function kernel(a::ModuleFPHom)
 error("kernel is not implemented for the given types.")
end

@doc Markdown.doc"""
    free_resolution(F::FreeMod)

Return a free resolution of `F`.

# Examples
"""
function free_resolution(F::FreeMod)
  res = presentation(F)
  #set_attribute!(res, :show => Hecke.free_show, :free_res => F)
  return res
end

@doc Markdown.doc"""
    is_complete(FR::FreeResolution)

Return `true` if the free resolution `fr` is complete, otherwise return `false`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubQuo(A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 4 generators
1 -> x^2*e[1]
2 -> x*y*e[1]
3 -> y^2*e[1]
4 -> z^4*e[1]

julia> fr = free_resolution(M, length=1)

rank   | 6  2
-------|------
degree | 1  0

julia> is_complete(fr)
false

julia> fr = free_resolution(M)

rank   | 0  2  6  6  2
-------|---------------
degree | 4  3  2  1  0

julia> is_complete(fr)
true

```
"""
is_complete(FR::FreeResolution) = FR.C.complete

#= Fill functions (and helpers) for Hecke ChainComplexes in terms of free resolutions =#
function _get_last_map_key(cc::Hecke.ChainComplex)
  ctr = cc.start
  while haskey(cc.maps, ctr)
    ctr -= 1
  end
  return ctr+1
end

function _extend_free_resolution(cc::Hecke.ChainComplex, idx::Int; algorithm::Symbol=:fres)
  key = _get_last_map_key(cc)
  if cc.complete == true
    return cc.maps[key]
  end

  kernel_entry          = image(cc.maps[key])[1]
  br                    = base_ring(kernel_entry)
  singular_free_module  = singular_module(ambient_free_module(kernel_entry))
  singular_kernel_entry = Singular.Module(base_ring(singular_free_module),
                              [singular_free_module(repres(g)) for g in gens(kernel_entry)]...)
  singular_kernel_entry.isGB = true

  len = idx + key - cc.start + 1
  if algorithm == :fres
    res = Singular.fres(singular_kernel_entry, len, "complete")
  elseif algorithm == :sres
    res = Singular.fres(singular_kernel_entry, len)
  elseif algorithm == :lres
    error("LaScala's method is not yet available in Oscar.")
  else
    error("Unsupported algorithm $algorithm")
  end

  if Singular.length(res) < len
    cc.complete = true
  end
  dom = domain(cc.maps[key])
  j   = 2
  while j <= Singular.length(res)
    key -= 1
    codom = dom
    dom   = free_module(br, Singular.ngens(res[j]))
    SM    = SubModuleOfFreeModule(codom, res[j])
    generator_matrix(SM)
    map = hom(dom, codom, SM.matrix)
    cc.maps[key] = map
    j += 1
  end
  # Finalize maps.
  if cc.complete == true
    key -= 1
    Z = FreeMod(br, 0)
    set_attribute!(Z, :name => "0")
    cc.maps[key] = hom(Z, domain(cc.maps[key+1]), Vector{elem_type(domain(cc.maps[key+1]))}())
  end
  return cc.maps[key]
end

@doc Markdown.doc"""
    free_resolution(M::SubQuo; ordering::ModuleOrdering = default_ordering(M),
        length::Int=0, algorithm::Symbol=:fres)

Return a free resolution of `M`.

If `length != 0`, the free resolution is only computed up to the `length`-th free module.
`algorithm` can be set to `:sres` or `:fres`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"])
(Multivariate Polynomial Ring in x, y, z over Rational Field, fmpq_mpoly[x, y, z])

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubQuo(A, B)
Subquotient of Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]
by Submodule with 4 generators
1 -> x^2*e[1]
2 -> x*y*e[1]
3 -> y^2*e[1]
4 -> z^4*e[1]

julia> fr = free_resolution(M, length=1)

rank   | 6  2
-------|------
degree | 1  0

julia> is_complete(fr)
false

julia> fr[4]
Free module of rank 0 over Multivariate Polynomial Ring in x, y, z over Rational Field

julia> fr

rank   | 0  2  6  6  2
-------|---------------
degree | 4  3  2  1  0

julia> is_complete(fr)
true

julia> fr = free_resolution(M, algorithm=:sres)

rank   | 0  2  6  6  2
-------|---------------
degree | 4  3  2  1  0
```
"""
function free_resolution(M::SubQuo; ordering::ModuleOrdering = default_ordering(M),
        length::Int=0, algorithm::Symbol=:fres)

  coefficient_ring(base_ring(M)) isa AbstractAlgebra.Field ||
      error("Must be defined over a field.")

  cc_complete = false

  #= Start with presentation =#
  pm = presentation(M)
  maps = [map(pm, j) for j in range(pm)]

  br = base_ring(M)
  kernel_entry          = image(pm.maps[1])[1]
  singular_free_module  = singular_module(ambient_free_module(kernel_entry))
  singular_kernel_entry = Singular.Module(base_ring(singular_free_module),
                              [singular_free_module(repres(g)) for g in gens(kernel_entry)]...)

  singular_kernel_entry.isGB = true

  #= This is the single computational hard part of this function =#
  if algorithm == :fres
    res = Singular.fres(singular_kernel_entry, length, "complete")
  elseif algorithm == :sres
    res = Singular.fres(singular_kernel_entry, length)
  elseif algorithm == :lres
    error("LaScala's method is not yet available in Oscar.")
  else
    error("Unsupported algorithm $algorithm")
  end

  if length == 0 || Singular.length(res) < length
    cc_complete = true
  end

  #= Add maps from free resolution computation, start with second entry
   = due to inclusion of presentation(M) at the beginning. =#
  dom = domain(pm.maps[1])
  j   = 2
  while j <= Singular.length(res)
    codom = dom
    rk    = Singular.ngens(res[j])
    dom   = free_module(br, rk)
    SM    = SubModuleOfFreeModule(codom, res[j])
    generator_matrix(SM)
    set_attribute!(dom, :name => "br^$rk")
    insert!(maps, 1, hom(dom, codom, SM.matrix))
    j += 1
  end
  if cc_complete == true
    # Finalize maps.
    Z = FreeMod(br, 0)
    set_attribute!(Z, :name => "0")
    insert!(maps, 1, hom(Z, domain(maps[1]), Vector{elem_type(domain(maps[1]))}()))
  end

  cc = Hecke.ChainComplex(Oscar.ModuleFP, maps, check = false, start = -2)
  cc.fill     = _extend_free_resolution
  cc.complete = cc_complete

  return FreeResolution(cc)
end


@doc Markdown.doc"""
    free_resolution_via_kernels(M::SubQuo, limit::Int = -1)

Return a free resolution of `M`.

If `limit != -1`, the free resolution
is only computed up to the `limit`-th free module.

# Examples
"""
function free_resolution_via_kernels(M::SubQuo, limit::Int = -1)
  p = presentation(M)
  mp = [map(p, j) for j in range(p)]
  while true
    k, mk = kernel(mp[1])
    nz = findall(x->!iszero(x), gens(k))
    if length(nz) == 0 
      Z = FreeMod(base_ring(M), 0)
      set_attribute!(Z, :name => "0")
      h = hom(Z, domain(mp[1]), Vector{elem_type(domain(mp[1]))}())
      insert!(mp, 1, h)
      break
    elseif limit != -1 && length(mp) > limit
      break
    end
    F = FreeMod(base_ring(M), length(nz))
    g = hom(F, codomain(mk), collect(k.sub.gens)[nz])
    insert!(mp, 1, g)
  end
  C = Hecke.ChainComplex(ModuleFP, mp, check = false)
  #set_attribute!(C, :show => Hecke.free_show, :free_res => M) # doesn't work
  return C
end

function Hecke.ring(I::MPolyIdeal)
  return parent(gen(I, 1))
end

function *(I::MPolyIdeal{T}, M::ModuleFP{T}) where {T<:RingElem}
  return sub(M, [g*e for g in gens(I) for e in gens(M)])
end

@doc Markdown.doc"""
    ideal_to_module(I::MPolyIdeal{T}, F::FreeMod{T})

Convert the ideal `I` into a submodule of `F`. If the rank of `F`
is not equal to 1 an error is thrown.
"""
function ideal_to_module(I::MPolyIdeal{T}, F::FreeMod{T}) where T
  @assert rank(F) == 1
  return sub(F, [x*gen(F,1) for x in gens(I)], :module)
end

@doc Markdown.doc"""
    ideal_to_module(I::MPolyIdeal)

Convert the ideal `I` into a submodule.
"""
function ideal_to_module(I::MPolyIdeal)
  F = FreeMod(base_ring(I), 1)
  return ideal_to_module(I, F)
end

@doc Markdown.doc"""
    free_resolution(I::MPolyIdeal)

Compute a free resolution of `I`.

# Examples
"""
function free_resolution(I::MPolyIdeal)
  F = free_module(Hecke.ring(I), 1)
  S = sub(F, [x * gen(F, 1) for x = gens(I)], :module)
  n = Hecke.find_name(I)
  if n !== nothing
    AbstractAlgebra.set_name!(S, string(n))
  end
  return free_resolution(S)
end

@doc Markdown.doc"""
    free_resolution(Q::MPolyQuo)

Compute a free resolution of `Q`.

# Examples
"""
function free_resolution(Q::MPolyQuo)
  F = free_module(Q, 1)
  q = quo(F, [x * gen(F, 1) for x = gens(Q.I)], :module)
  n = Hecke.find_name(Q)
  if n !== nothing
    AbstractAlgebra.set_name!(q, String(n))
  end
  return free_resolution(q)
end

@doc Markdown.doc"""
    iszero(f::ModuleFPHom)

Return true iff `f` is the zero map.
"""
function iszero(f::ModuleFPHom)
  return all(iszero, map(f, gens(domain(f))))
end

@doc Markdown.doc"""
    hom(M::ModuleFP, N::ModuleFP)

Return the module `Hom(M,N)` as an object of type `SubQuo`.

Additionally, if `H` is that object, return the map which sends an element of `H` to the corresponding homomorphism `M` $\to$ `N`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

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
hom of (Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y^2*e[2], Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y^2*e[2])

julia> gens(H)
4-element Vector{SubQuoElem{fmpq_mpoly}}:
 (e[1] -> e[1])
 (e[2] -> e[2])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])

julia> relations(H)
4-element Vector{FreeModElem{fmpq_mpoly}}:
 x*(e[1] -> e[1])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])
 y^2*(e[2] -> e[2])
```
"""
function hom(M::ModuleFP, N::ModuleFP, alg::Symbol=:maps)
  #source: Janko's CA script: https://www.mathematik.uni-kl.de/~boehm/lehre/17_CA/ca.pdf
  if alg == :matrices && M isa SubQuo && N isa SubQuo
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

  function im(x::SubQuoElem)
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

    return s_proj(SubQuoElem(repres(preimage(mH_s0_t0, g)), H))
  end
  to_hom_map = MapFromFunc(im, pre, H_simplified, Hecke.MapParent(M, N, "homomorphisms"))
  set_attribute!(H_simplified, :show => Hecke.show_hom, :hom => (M, N), :module_to_hom_map => to_hom_map)
  return H_simplified, to_hom_map
end

@doc Markdown.doc"""
    element_to_homomorphism(f::ModuleFPElem)

If `f` is an element of a module created via `hom(M,N)`, for some modules `M` and `N`, 
return the homomorphism `M` $\to$ `N` corresponding to `f`.
# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

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
4-element Vector{SubQuoElem{fmpq_mpoly}}:
 (e[1] -> e[1])
 (e[2] -> e[2])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])

julia> relations(H)
4-element Vector{FreeModElem{fmpq_mpoly}}:
 x*(e[1] -> e[1])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])
 y^2*(e[2] -> e[2])

julia> a = homomorphism(H[1]+y*H[2])
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

@doc Markdown.doc"""
    homomorphism_to_element(H::ModuleFP, a::ModuleMap)

If the module `H` is created via `hom(M,N)`, for some modules `M` and `N`, and
`a`: `M` $\to$ `N` is a homomorphism, then return the element of `H` corresponding to `a`.
# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

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
4-element Vector{SubQuoElem{fmpq_mpoly}}:
 (e[1] -> e[1])
 (e[2] -> e[2])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])

julia> relations(H)
4-element Vector{FreeModElem{fmpq_mpoly}}:
 x*(e[1] -> e[1])
 y^2*(e[1] -> e[2])
 x*(e[2] -> e[1])
 y^2*(e[2] -> e[2])

julia> W =  [M[1], y*M[2]];

julia> a = hom(M, M, W);

julia> iswelldefined(a)
true

julia> matrix(a)
[1   0]
[0   y]

julia> m = module_elem(H, a)
(e[1] -> e[1]) + y*(e[2] -> e[2])
```
"""
function homomorphism_to_element(H::ModuleFP, phi::ModuleMap)
  to_hom_map = get_attribute(H, :module_to_hom_map)
  to_hom_map === nothing && error("module must be a hom module")
  map_to_hom = to_hom_map.g
  return map_to_hom(phi)
end

@doc Markdown.doc"""
    multiplication_morphism(a::RingElem, M::ModuleFP)

Return the multiplication by `a` as an endomorphism on `M`.
"""
function multiplication_morphism(a::RingElem, M::ModuleFP)
  @assert base_ring(M) === parent(a)
  return hom(M, M, [a*v for v in gens(M)])
end

@doc Markdown.doc"""
    multiplication_morphism(a::FreeModElem, M::ModuleFP)

Return the multiplication by `a` as an endomorphism on `M`. For this,
the parent of `a` must be a module of rank 1.
"""
function multiplication_morphism(a::FreeModElem, M::ModuleFP)
  @assert rank(parent(a)) == 1
  return multiplication_morphism(a[1], M)
end

@doc Markdown.doc"""
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

#TODO
#  replace the +/- for the homs by proper constructors for homs and direct sums
#  relshp to store the maps elsewhere

@doc Markdown.doc"""
    *(h::ModuleFPHom, g::ModuleFPHom)

Return the composition $g \circ h$.
"""
function *(h::ModuleFPHom, g::ModuleFPHom)
  @assert codomain(h) === domain(g)
  return hom(domain(h), codomain(g), Vector{elem_type(codomain(g))}([g(h(x)) for x = gens(domain(h))]))
end

compose(h::ModuleFPHom, g::ModuleFPHom) = h*g

-(h::ModuleFPHom) = hom(domain(h), codomain(h), [-h(x) for x in gens(domain(h))])
function -(h::ModuleFPHom, g::ModuleFPHom)
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  return hom(domain(h), codomain(h), Vector{elem_type(codomain(h))}([h(x) - g(x) for x in gens(domain(h))]))
end
function +(h::ModuleFPHom, g::ModuleFPHom)
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  return hom(domain(h), codomain(h), Vector{elem_type(codomain(h))}([h(x) + g(x) for x in gens(domain(h))]))
end
function *(a::RingElem, g::ModuleMap)
  @assert base_ring(codomain(g)) === parent(a)
  return hom(domain(g), codomain(g), Vector{elem_type(codomain(g))}([a*g(x) for x in gens(domain(g))]))
end


@doc Markdown.doc"""
    restrict_codomain(H::ModuleFPHom, M::SubQuo)

Return, if possible, a homomorphism, which is mathematically identical to `H`,
but has codomain `M`. `M` has to be a submodule of the codomain of `H`.
"""
function restrict_codomain(H::ModuleFPHom, M::SubQuo)
  D = domain(H)
  return hom(D, M, map(v -> SubQuoElem(v, M), map(x -> repres(H(x)), gens(D))))
end

@doc Markdown.doc"""
    restrict_domain(H::SubQuoHom, M::SubQuo)

Restrict the morphism `H` to `M`. For this `M` has to be a submodule
of the domain of `H`. The relations of `M` must be the relations of 
the domain of `H`.
"""
function restrict_domain(H::SubQuoHom, M::SubQuo)
  for i in M.outgoing_morphisms
    if codomain(i) === domain(H)
      return i*H
    end
  end
  # else there is no cached map
  if ngens(M) > 0
    @assert M.quo == domain(H).quo
  end
  i = sub(domain(H), map(m -> SubQuoElem(repres(m), domain(H)), gens(M)), :cache_morphism)[2]
  return i*H
end

@doc Markdown.doc"""
    induced_map(f::FreeModuleHom, M::SubQuo, check::Bool = true)

Return the map which sends an element `v` of `M` to `f(repres(v))`.
If `check` is set to true the well-definedness of the map is checked.
"""
function induced_map(f::FreeModuleHom, M::SubQuo, check::Bool = true)
  @assert ambient_free_module(M) === domain(f)
  ind_f = hom(M, codomain(f), [f(repres(v)) for v in gens(M)])
  if check
    @assert is_welldefined(ind_f)
  end
  return ind_f
end

@doc Markdown.doc"""
    inv(H::ModuleFPHom)

Compute $H^{-1}$. `H` must be bijective.
"""
function inv(H::ModuleFPHom)
  if isdefined(H, :inverse_isomorphism)
    return H.inverse_isomorphism
  end
  @assert is_bijective(H)
  N = domain(H)
  M = codomain(H)

  Hinv = hom(M,N, Vector{elem_type(N)}([preimage(H,m) for m in gens(M)]))
  Hinv.inverse_isomorphism = H
  H.inverse_isomorphism = Hinv

  return Hinv
end

##################################################
# direct product
##################################################
@doc Markdown.doc"""
    direct_product(F::FreeMod{T}...; task::Symbol = :prod) where T

Given free modules $F_1\dots F_n$, say, return the direct product $\prod_{i=1}^n F_i$.

Additionally, return
- a vector containing the canonical projections  $\prod_{i=1}^n F_i\rightarrow F_i$ if `task = :prod` (default),
- a vector containing the canonical injections  $F_i\rightarrow\prod_{i=1}^n F_i$ if `task = :sum`,
- two vectors containing the canonical projections and injections, respectively, if `task = :both`,
- none of the above maps if `task = :none`.
"""
function direct_product(F::FreeMod{T}...; task::Symbol = :prod) where {T}
  R = base_ring(F[1])
  G = FreeMod(R, sum([rank(f) for f = F]))
  G.S = []
  for i = 1:length(F)
    s = "("
    for j=1:i-1
      s *= "0, "
    end
    e = ""
    if i<length(F)
      e*=", "
    end
    for j=i+1:length(F)-1
      e *= "0, "
    end
    if i<length(F)
      e *= "0"
    end
    e*=")"

    for t = F[i].S
      push!(G.S, Symbol(s*string(t)*e))
    end
  end
  set_attribute!(G, :show => Hecke.show_direct_product, :direct_product => F)
  emb = []
  pro = []
  projection_dictionary = IdDict{Int,ModuleFPHom}()
  injection_dictionary = IdDict{Int,ModuleFPHom}()
  set_attribute!(G, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  i = 0
  for f = F
    if task in [:sum, :both]
      push!(emb, hom(f, G, Vector{elem_type(G)}([gen(G, j+i) for j=1:ngens(f)])))
      injection_dictionary[length(emb)] = emb[length(emb)]
    end
    if task in [:prod, :both]
      push!(pro, hom(G, f, vcat(elem_type(f)[zero(f) for j=1:i], gens(f), elem_type(f)[zero(f) for j=i+ngens(f)+1:ngens(G)])))
      projection_dictionary[length(pro)] = pro[length(pro)]
    end
    i += ngens(f)
  end
  if task == :none
    return G
  elseif task == :sum
    return G, emb
  elseif task == :prod
    return G, pro
  elseif task == :both
    return G, pro, emb
  end
end

@doc Markdown.doc"""
    direct_product(M::ModuleFP{T}...; task::Symbol = :prod) where T

Given modules $M_1\dots M_n$, say, return the direct product $\prod_{i=1}^n M_i$.

Additionally, return
- a vector containing the canonical projections  $\prod_{i=1}^n M_i\rightarrow M_i$ if `task = :prod` (default),
- a vector containing the canonical injections  $M_i\rightarrow\prod_{i=1}^n M_i$ if `task = :sum`,
- two vectors containing the canonical projections and injections, respectively, if `task = :both`,
- none of the above maps if `task = :none`.
"""
function direct_product(M::ModuleFP{T}...; task::Symbol = :prod) where T
  F, pro, mF = direct_product([ambient_free_module(x) for x = M]..., task = :both)
  s, emb_sF = sub(F, vcat([[mF[i](y) for y = ambient_representatives_generators(M[i])] for i=1:length(M)]...), :both)
  q::Vector{elem_type(F)} = vcat([[mF[i](y) for y = rels(M[i])] for i=1:length(M)]...)
  pro_quo = nothing
  if length(q) != 0
    s, pro_quo = quo(s, q, :both)
  end
  set_attribute!(s, :show => Hecke.show_direct_product, :direct_product => M)
  projection_dictionary = IdDict{Int,ModuleFPHom}()
  injection_dictionary = IdDict{Int,ModuleFPHom}()
  set_attribute!(s, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  if task == :none
    return s
  end
  if task == :prod || task != :sum
    if pro_quo === nothing
      for i=1:length(pro)
        pro[i] = hom(s, M[i], Vector{elem_type(M[i])}([M[i](pro[i](emb_sF(gen))) for gen in gens(s)])) # TODO distinction between pro on the left and pro on the right side!
        projection_dictionary[i] = pro[i]
      end
    else
      for i=1:length(pro)
        pro[i] = hom(s, M[i], Vector{elem_type(M[i])}([M[i](pro[i](emb_sF(preimage(pro_quo,gen)))) for gen in gens(s)]))
        projection_dictionary[i] = pro[i]
      end
    end
    if task == :prod
      return s, pro
    end
  end
  if task == :sum || task != :prod
    if pro_quo === nothing
      for i=1:length(mF)
        mF[i] = hom(M[i], s, Vector{elem_type(s)}([preimage(emb_sF, mF[i](repres(g))) for g in gens(M[i])]))
        injection_dictionary[i] = mF[i]
      end
    else
      for i=1:length(mF)
        mF[i] = hom(M[i], s, Vector{elem_type(s)}([pro_quo(preimage(emb_sF, mF[i](repres(g)))) for g in gens(M[i])]))
        injection_dictionary[i] = mF[i]
      end
    end
    if task == :sum
      return s, mF
    else
      return s, pro, mF
    end
  end
end
##################################################
# direct sum
##################################################
@doc Markdown.doc"""
    direct_sum(M::ModuleFP{T}...; task::Symbol = :sum) where T

Given modules $M_1\dots M_n$, say, return the direct sum $\bigoplus_{i=1}^n M_i$.  
 
Additionally, return 
- a vector containing the canonical injections  $M_i\rightarrow\bigoplus_{i=1}^n M_i$ if `task = :sum` (default),
- a vector containing the canonical projections  $\bigoplus_{i=1}^n M_i\rightarrow M_i$ if `task = :prod`,
- two vectors containing the canonical injections and projections, respectively, if `task = :both`,
- none of the above maps if `task = :none`.
"""
function direct_sum(M::ModuleFP{T}...; task::Symbol = :sum) where {T}
  res = direct_product(M...; task)
  if task == :sum || task == :prod
    ds, f = res
    set_attribute!(ds, :show => Hecke.show_direct_sum, :direct_sum => M)
    return ds, f
  elseif task == :both
    ds, p, i = res
    set_attribute!(ds, :show => Hecke.show_direct_sum, :direct_sum => M)
    return ds, i, p
  else
    set_attribute!(res, :show => Hecke.show_direct_sum, :direct_sum => M)
    return res
  end
end

⊕(M::ModuleFP...) = direct_sum(M..., task = :none)


@doc Markdown.doc"""
    Hecke.canonical_injection(G::ModuleFP, i::Int)

Return the canonical injection $G_i \to G$ where $G = G_1 \oplus \cdot \oplus G_n$.
"""
function Hecke.canonical_injection(G::ModuleFP, i::Int)
  H = get_attribute(G, :direct_product)
  if H === nothing
    error("module not a direct product")
  end
  injection_dictionary = get_attribute(G, :injection_morphisms)
  if haskey(injection_dictionary, i)
    return injection_dictionary[i]
  end
  0<i<= length(H) || error("index out of bound")
  j = i == 1 ? 0 : sum(ngens(H[l]) for l=1:i-1)
  emb = hom(H[i], G, Vector{elem_type(G)}([G[l+j] for l = 1:ngens(H[i])]))
  injection_dictionary[i] = emb
  return emb
end

@doc Markdown.doc"""
    Hecke.canonical_projection(G::ModuleFP, i::Int)

Return the canonical projection $G \to G_i$ where $G = G_1 \oplus \cdot \oplus G_n$.
"""
function Hecke.canonical_projection(G::ModuleFP, i::Int)
  H = get_attribute(G, :direct_product)
  if H === nothing
    error("module not a direct product")
  end
  projection_dictionary = get_attribute(G, :projection_morphisms)
  if haskey(projection_dictionary, i)
    return projection_dictionary[i]
  end
  0<i<= length(H) || error("index out of bound")
  j = i == 1 ? 0 : sum(ngens(H[l]) for l=1:i-1) 
  pro = hom(G, H[i], Vector{elem_type(H[i])}(vcat([zero(H[i]) for l=1:j], gens(H[i]), [zero(H[i]) for l=1+j+ngens(H[i]):ngens(G)])))
  projection_dictionary[i] = pro
  return pro
end
    
##################################################
# Tensor
##################################################

@doc Markdown.doc"""
    tensor_product(F::FreeMod...; task::Symbol = :none)

Given a collection of free modules, say, $F_1, \dots, F_n$ over a ring $R$, return $F_1\otimes_R \cdots \otimes_R F_n$.


If `task = :map`, additionally return the map which sends a tuple $(f_1,\dots, f_n)$ of elements $f_i\in F_i$ to the pure tensor $f_1\otimes\dots\otimes f_n$.
"""
function tensor_product(G::FreeMod...; task::Symbol = :none)
  s = G[1].S
  t = [[x] for x = 1:ngens(G[1])]
  for H = G[2:end]
    s = [Symbol("$x \\otimes $y") for x = s  for y = H.S]
    t = [push!(deepcopy(x), y) for x = t  for y = 1:ngens(H)]
  end

  F = FreeMod(G[1].R, prod([rank(g) for g in G]))
  F.S = s
  set_attribute!(F, :show => Hecke.show_tensor_product, :tensor_product => G)

  function pure(g::FreeModElem...)
    @assert length(g) == length(G)
    @assert all(i -> parent(g[i]) === G[i], 1:length(G))
    z = [[x] for x = g[1].coords.pos]
    zz = g[1].coords.values
    for h = g[2:end]
      zzz = Vector{Int}[]
      zzzz = elem_type(F.R)[]
      for i = 1:length(z)
        for (p, v) = h.coords
          push!(zzz, push!(deepcopy(z[i]), p))
          push!(zzzz, zz[i]*v)
        end
      end
      z = zzz
      zz = zzzz
    end
    indices = Vector{Int}([findfirst(x->x == y, t) for y = z])
    return FreeModElem(sparse_row(F.R, indices, zz), F)
  end
  function pure(T::Tuple)
    return pure(T...)
  end
  function inv_pure(e::FreeModElem)
    if length(e.coords.pos) == 0
      return Tuple(zero(g) for g = G)
    end
    @assert length(e.coords.pos) == 1
    @assert isone(e.coords.values[1])
    return Tuple(gen(G[i], t[e.coords.pos[1]][i]) for i = 1:length(G))
  end

  set_attribute!(F, :tensor_pure_function => pure, :tensor_generator_decompose_function => inv_pure)

  if task == :none
    return F
  end

  return F, MapFromFunc(pure, inv_pure, Hecke.TupleParent(Tuple([g[0] for g = G])), F)
end

⊗(G::ModuleFP...) = tensor_product(G..., task = :none)

@doc Markdown.doc"""
    ambient_free_module(F::FreeMod)

Just return `F`. This function exists only for compatibility reasons.
"""
function ambient_free_module(F::FreeMod)
  return F
end

@doc Markdown.doc"""
    ambient_module(F::FreeMod, task = :none)

Just return `F`. This function exists only for compatibility reasons.

If task is set to :with_morphism, return also the identity map.
"""
function ambient_module(F::FreeMod, task = :none)
  if task == :none
    return F
  else
    return F, identity_map(F)
  end
end

@doc Markdown.doc"""
    ambient_free_module(M::SubQuo)

Return the ambient free module of `M`.
"""
function ambient_free_module(M::SubQuo)
  return M.F
end

@doc Markdown.doc"""
    ambient_module(M::SubQuo, task = :none)

If $M = (P + Q) / Q$ then return $F / Q$ where `F == ambient_free_module(M)`.
If $M$ is a submodule of a free module, then the (ambient) free module is returned.

If `task == :with_morphism`, return also the canonical inclusion.
"""
function ambient_module(M::SubQuo, task = :none)
  if !isdefined(M, :quo)
    if task == :none
      return F
    else
      return F, hom(M,F,[repres(v) for v in gens(M)])
    end
  end
  F = ambient_free_module(M)
  g = SubModuleOfFreeModule(F,basis(F))
  SQ = SubQuo(g, M.quo)
  if task == :none
    return SQ
  else
    return SQ, hom(M,SQ,[SubQuoElem(coords(repres(v)), SQ) for v in gens(M)])
  end
end

@doc Markdown.doc"""
    ambient_representatives_generators(F::FreeMod)

Return the generators of `F`. This function exists only for 
compatibility reasons.
"""
function ambient_representatives_generators(F::FreeMod)
  return gens(F)
end

rels(F::FreeMod) = elem_type(F)[]

@doc Markdown.doc"""
    ambient_representatives_generators(M::SubQuo)

Return elements of the ambient free module of `M` which represent the generators of `M`.
"""
function ambient_representatives_generators(M::SubQuo)
  G = ambient_free_module(M)
  return [FreeModElem(coords(repres(x)), G) for x in gens(M)]
end


@doc Markdown.doc"""
    rels(M::SubQuo)

Return the relations of `M`.
"""
rels(M::SubQuo) = isdefined(M, :quo) ? collect(M.quo.gens) : elem_type(M.F)[]

@doc Markdown.doc"""
    relations(M::SubQuo)

Return the relations of `M`.
"""
relations(M::SubQuo) = rels(M)

@doc Markdown.doc"""
    tensor_product(M::ModuleFP...; task::Symbol = :none)

Given a collection of modules, say, $M_1, \dots, M_n$ over a ring $R$, return $M_1\otimes_R \cdots \otimes_R M_n$.

If `task = :map`, additionally return the map which sends a tuple $(m_1,\dots, m_n)$ of elements $m_i\in M_i$ to the pure tensor $m_1\otimes\dots\otimes m_n$.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> F = free_module(R, 1);

julia> A = R[x; y];

julia> B = R[x^2; y^3; z^4];

julia> M = SubQuo(F, A, B);

julia> gens(M)
2-element Vector{SubQuoElem{fmpq_mpoly}}:
 x*e[1]
 y*e[1]

julia> T, t = tensor_product(M, M; task = :map);

julia> gens(T)
4-element Vector{SubQuoElem{fmpq_mpoly}}:
 x^2*e[1] \otimes e[1]
 x*y*e[1] \otimes e[1]
 x*y*e[1] \otimes e[1]
 y^2*e[1] \otimes e[1]

julia> domain(t)
parent of tuples of type Tuple{SubQuoElem{fmpq_mpoly}, SubQuoElem{fmpq_mpoly}}

julia> t((M[1], M[2]))
x*y*e[1] \otimes e[1]
```
"""
function tensor_product(G::ModuleFP...; task::Symbol = :none)
  F, mF = tensor_product([ambient_free_module(x) for x = G]..., task = :map)
  # We want to store a dict where the keys are tuples of indices and the values
  # are the corresponding pure vectors (i.e. a tuple (2,1,5) represents the 
  # 2nd, 1st and 5th generator of the 1st, 2nd and 3rd module, which we are 
  # tensoring. The corresponding value is then G[1][2] ⊗ G[2][1] ⊗ G[2][5]). 
  corresponding_tuples_as_indices = vec([x for x = Base.Iterators.ProductIterator(Tuple(1:ngens(x) for x = G))])
  # In corresponding_tuples we store tuples of the actual generators, so in 
  # the example above we would store (G[1][2], G[2][1], G[2][5]).
  corresponding_tuples = map(index_tuple -> Tuple(map(index -> G[index][index_tuple[index]],1:length(index_tuple))), corresponding_tuples_as_indices)

  generating_tensors::Vector{elem_type(F)} = map(mF, map(tuple -> map(x -> parent(x) isa FreeMod ? x : x.repres, tuple), corresponding_tuples))
  s, emb = sub(F, generating_tensors, :with_morphism)
  #s, emb = sub(F, vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(gens(x, ambient_free_module(x)) for x = G))]), :with_morphism)
  q::Vector{elem_type(F)} = vcat([vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(i == j ? rels(G[i]) : gens(ambient_free_module(G[i])) for i=1:length(G)))]) for j=1:length(G)]...) 
  local projection_map
  if length(q) != 0
    s, projection_map = quo(s, q, :with_morphism)
  end

  tuples_pure_tensors_dict = IdDict(zip(corresponding_tuples_as_indices, gens(s)))
  set_attribute!(s, :show => Hecke.show_tensor_product, :tensor_product => G)

  
  function pure(tuple_elems::Union{SubQuoElem,FreeModElem}...)
    coeffs_tuples = vec([x for x = Base.Iterators.ProductIterator(Tuple(coeffs(x) for x = tuple_elems))])
    res = zero(s)
    for coeffs_tuple in coeffs_tuples
      indices = map(x -> x[1], coeffs_tuple)
      coeff_for_pure = prod(map(x -> x[2], coeffs_tuple))
      res += coeff_for_pure*tuples_pure_tensors_dict[indices]
    end
    return res
  end
  function pure(T::Tuple)
    return pure(T...)
  end

  decompose_generator = function(v::SubQuoElem)
    i = index_of_gen(v)
    return corresponding_tuples[i]
  end

  set_attribute!(s, :tensor_pure_function => pure, :tensor_generator_decompose_function => decompose_generator)

  if task == :none
    return s
  end

  return s, MapFromFunc(pure, Hecke.TupleParent(Tuple([g[0] for g = G])), s)
end

#############################
# Tor
#############################
@doc Markdown.doc"""
    tensor_product(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})

Apply $P \otimes \bullet$ to `C`.
"""
function tensor_product(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})
  #tensor_chain = Hecke.map_type(C)[]
  tensor_chain = valtype(C.maps)[]
  tensor_modules = [tensor_product(P, domain(map(C,first(range(C)))), task=:cache_morphism)[1]]
  append!(tensor_modules, [tensor_product(P, codomain(map(C,i)), task=:cache_morphism)[1] for i in range(C)])

  for i in 1:length(range(C))
    A = tensor_modules[i]
    B = tensor_modules[i+1]

    j = range(C)[i]
    push!(tensor_chain, hom_tensor(A,B,[identity_map(P), map(C,j)]))
  end

  start = Hecke.is_chain_complex(C) ? C.start - length(tensor_chain) - 1 : C.start
  return Hecke.ChainComplex(ModuleFP, tensor_chain, start=start, direction=C.direction)
end

@doc Markdown.doc"""
    tensor_product(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)

Apply $\bullet \otimes P$ to `C`.
"""
function tensor_product(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
  #tensor_chain = Hecke.map_type(C)[]
  tensor_chain = valtype(C.maps)[]
  chain_range = range(C)
  tensor_modules = [tensor_product(domain(map(C,first(chain_range))), P, task=:cache_morphism)[1]]
  append!(tensor_modules, [tensor_product(codomain(map(C,i)), P, task=:cache_morphism)[1] for i in chain_range])

  for i=1:length(chain_range)
    A = tensor_modules[i]
    B = tensor_modules[i+1]

    j = chain_range[i]
    push!(tensor_chain, hom_tensor(A,B,[map(C,j), identity_map(P)]))
  end

  start = Hecke.is_chain_complex(C) ? C.start - length(tensor_chain) - 1 : C.start
  return Hecke.ChainComplex(ModuleFP, tensor_chain, start=start, direction=C.direction)
end

@doc Markdown.doc"""
    tor(M::ModuleFP, N::ModuleFP, i::Int)

Return $\text{Tor}_i(M,N)$.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> A = R[x; y];

julia> B = R[x^2; y^3; z^4];

julia> M = SubQuo(A, B);

julia> F = free_module(R, 1);

julia> Q, _ = quo(F, [x*F[1]]);

julia> T0 = tor(Q, M, 0)
Subquotient of Submodule with 2 generators
1 -> x*e[1] \otimes e[1]
2 -> y*e[1] \otimes e[1]
by Submodule with 5 generators
1 -> x^2*e[1] \otimes e[1]
2 -> y^3*e[1] \otimes e[1]
3 -> z^4*e[1] \otimes e[1]
4 -> x^2*e[1] \otimes e[1]
5 -> x*y*e[1] \otimes e[1]

julia> T1 = tor(Q, M, 1)
Subquotient of Submodule with 1 generator
1 -> -x*e[1] \otimes e[1]
by Submodule with 3 generators
1 -> x^2*e[1] \otimes e[1]
2 -> y^3*e[1] \otimes e[1]
3 -> z^4*e[1] \otimes e[1]

julia> T2 =  tor(Q, M, 2)
Submodule with 0 generators
represented as subquotient with no relations.
```
"""
function tor(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M; length=i+2)
  lifted_resolution = tensor_product(free_res.C[first(range(free_res.C)):-1:1], N) #TODO only three homs are necessary
  return simplify_light(homology(lifted_resolution,i))[1]
end

#TODO, mF
#  (hom lift) => hom and tensor functor
#  filtrations
#  more constructors
#################################################
#
#################################################
@doc Markdown.doc"""
    lift_homomorphism_contravariant(Hom_MP::ModuleFP, Hom_NP::ModuleFP, phi::ModuleFPHom)

Let `Hom_MP` $= \text{Hom}(M,P)$, `Hom_NP` $= \text{Hom}(N,P)$ and `phi` $= \phi : N \to M$ a morphism.
Compute $\phi^{\ast} : \text{Hom}(M,P) \to \text{Hom}(N,P)$.
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

@doc Markdown.doc"""
    lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, phi::ModuleFPHom)

Let `Hom_PM` $= \text{Hom}(P,M)$, `Hom_PN` $= \text{Hom}(P,N)$ and `phi` $= \phi : M \to N$ a morphism.
Compute $\phi_{\ast} : \text{Hom}(P,M) \to \text{Hom}(P,N)$.
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

@doc Markdown.doc"""
    hom(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})

Apply $\text{Hom}(P,-)$ to `C`. Return the lifted chain complex.
"""
function hom(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})
  #hom_chain = Hecke.map_type(C)[]
  hom_chain = valtype(C.maps)[]
  chain_range = range(C)
  hom_modules = [hom(P, domain(map(C,first(chain_range))))]
  append!(hom_modules, [hom(P, codomain(map(C,i))) for i in chain_range])

  for i=1:length(chain_range)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    j = chain_range[i]
    push!(hom_chain, lift_homomorphism_covariant(A,B,map(C,j)))
  end

  start = Hecke.is_chain_complex(C) ? C.start - length(hom_chain) - 1 : C.start
  return Hecke.ChainComplex(ModuleFP, hom_chain, start=start, direction=C.direction)
end

@doc Markdown.doc"""
    hom(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)

Apply $\text{Hom}(-,P)$ to `C`. Return the lifted chain complex.
"""
function hom(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
  #hom_chain = Hecke.map_type(C)[]
  hom_chain = valtype(C.maps)[]
  chain_range = range(C)
  hom_modules = [hom(domain(map(C,first(chain_range))),P)]
  append!(hom_modules, [hom(codomain(map(C,i)), P) for i in chain_range])

  for i=1:length(chain_range)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    j = chain_range[i]
    push!(hom_chain, lift_homomorphism_contravariant(B,A,map(C,j)))
  end

  direction = Hecke.is_chain_complex(C) ? :right : :left
  start = Hecke.is_chain_complex(C) ? last(range(C))-1 : first(range(C)) - 1
  return Hecke.ChainComplex(ModuleFP, reverse(hom_chain), start=start, direction=direction)
end

@doc Markdown.doc"""
    hom_without_reversing_direction(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
Apply $\text{Hom}(-,P)$ to `C`. If `C` is a chain complex, return a chain complex
and accordingly if `C` is a cochain complex, return a cochain complex.
"""
function hom_without_reversing_direction(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
  #hom_chain = Hecke.map_type(C)[]
  hom_chain = valtype(C.maps)[]
  chain_range = range(C)
  hom_modules = [hom(domain(map(C,first(chain_range))),P)]
  append!(hom_modules, [hom(codomain(map(C,i)), P) for i in chain_range])

  for i=1:length(chain_range)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    j = chain_range[i]
    push!(hom_chain, lift_homomorphism_contravariant(B,A,map(C,j)))
  end

  start = Hecke.is_chain_complex(C) ? -first(chain_range)-1 : -last(chain_range)-1  
  return Hecke.ChainComplex(ModuleFP, reverse(hom_chain), start=start, direction=C.direction)
end

#############################
@doc Markdown.doc"""
    homology(C::ChainComplex{<:ModuleFP})

Return the homology of `C`.

# Examples
```jldoctest
julia> R, (x,) = PolynomialRing(QQ, ["x"]);

julia> F = free_module(R, 1);

julia> A, _ = quo(F, [x^4*F[1]]);

julia> B, _ = quo(F, [x^3*F[1]]);

julia> a = hom(A, B, [x^2*B[1]]);

julia> b = hom(B, B, [x^2*B[1]]);

julia> C = ChainComplex(ModuleFP, [a, b])
C_2 ----> C_1 ----> C_0
where:
	C_2 = Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> x^4*e[1]
	C_1 = Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> x^3*e[1]
	C_0 = Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> x^3*e[1]

julia> H = homology(C)
3-element Vector{SubQuo{fmpq_mpoly}}:
 Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 1 generator
1 -> x^4*e[1]
 Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 2 generators
1 -> x^3*e[1]
2 -> x^2*e[1]
 Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 2 generators
1 -> x^3*e[1]
2 -> x^2*e[1]
```
"""
function homology(C::Hecke.ChainComplex{<:ModuleFP})
  return [homology(C,i) for i in Hecke.object_range(C)]
end

@doc Markdown.doc"""
    homology(C::ChainComplex{<:ModuleFP}, i::Int)

Return the `i`-th homology module of `C`.

# Examples
```jldoctest
julia> R, (x,) = PolynomialRing(QQ, ["x"]);

julia> F = free_module(R, 1);

julia> A, _ = quo(F, [x^4*F[1]]);

julia> B, _ = quo(F, [x^3*F[1]]);

julia> a = hom(A, B, [x^2*B[1]]);

julia> b = hom(B, B, [x^2*B[1]]);

julia> C = ChainComplex(ModuleFP, [a, b])
C_2 ----> C_1 ----> C_0
where:
	C_2 = Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> x^4*e[1]
	C_1 = Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> x^3*e[1]
	C_0 = Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 1 generator
1 -> x^3*e[1]

julia> H = homology(C, 1)
Subquotient of Submodule with 1 generator
1 -> x*e[1]
by Submodule with 2 generators
1 -> x^3*e[1]
2 -> x^2*e[1]
```
"""
function homology(C::Hecke.ChainComplex{<:ModuleFP}, i::Int)
  chain_range = Hecke.object_range(C)
  @assert length(chain_range) > 0 #TODO we need actually only the base ring
  if i == first(chain_range)
    i = Hecke.is_chain_complex(C) ? i : i+1
    return kernel(map(C,i))[1]
  elseif i == last(chain_range)
    i = Hecke.is_chain_complex(C) ? i+1 : i
    f = map(C,i)
    return cokernel(f)    
  elseif i in chain_range
    next_index = Hecke.is_chain_complex(C) ? i : i+1
    i = Hecke.is_chain_complex(C) ? i+1 : i
    return quo(kernel(map(C,next_index))[1], image(map(C,i))[1], :module)
  else
    return FreeMod(base_ring(obj(C,first(chain_range))),0)
  end
end

#############################
# Ext
#############################
@doc Markdown.doc"""
    ext(M::ModuleFP, N::ModuleFP, i::Int)

Return $\text{Ext}^i(M,N)$.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> F = FreeMod(R, 1);

julia> V = [x*F[1], y*F[1]];

julia> M = quo(F, V)[1]
Subquotient of Submodule with 1 generator
1 -> e[1]
by Submodule with 2 generators
1 -> x*e[1]
2 -> y*e[1]

julia> ext(M, M, 0)
Subquotient of Submodule with 1 generator
1 -> (e[1] -> e[1])
by Submodule with 2 generators
1 -> x*(e[1] -> e[1])
2 -> y*(e[1] -> e[1])

julia> ext(M, M, 1)
Subquotient of Submodule with 2 generators
1 -> (e[2] -> e[1])
2 -> (e[1] -> e[1])
by Submodule with 5 generators
1 -> x*(e[1] -> e[1])
2 -> y*(e[1] -> e[1])
3 -> x*(e[2] -> e[1])
4 -> y*(e[2] -> e[1])
5 -> 0

julia> ext(M, M, 2)
Subquotient of Submodule with 1 generator
1 -> (e[1] -> e[1])
by Submodule with 4 generators
1 -> x*(e[1] -> e[1])
2 -> y*(e[1] -> e[1])
3 -> -y*(e[1] -> e[1])
4 -> 0

julia> ext(M, M, 3)
Subquotient of Submodule with 0 generators
by Submodule with 1 generator
1 -> 0
```
"""
function ext(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M; length=i+2)
  lifted_resolution = hom(free_res.C[first(range(free_res.C)):-1:1], N) #TODO only three homs are necessary
  return simplify_light(homology(lifted_resolution,i))[1]
end

#############################
# TODO ?
#############################
@doc Markdown.doc"""
    find_sequence_of_morphisms(N::SubQuo, M::SubQuo)

Compute a path from `N` to `M` in the graph of cached (canonical) morphisms.
Return it as Vector of maps where the first element has domain `N` and the last 
has codomain `M`.
If there exists no path an error is thrown.
"""
function find_sequence_of_morphisms(N::SubQuo, M::SubQuo)
  if M===N
    return [identity_map(M)]
  end
  parent_hom = IdDict{SubQuo,ModuleFPHom}()
  modules = [M]
  found_N = false
  for A in modules
    for H in A.incoming_morphisms
      B = domain(H)
      if B!==A # on trees "B!==A" is enough!
        if findfirst(x->x===B,modules) == nothing #if !(B in modules) doesn't work since it uses == instead of ===
          parent_hom[B] = H
          push!(modules,B)
        end
      end
      if B===N
        found_N = true
        break
      end
    end
    if found_N
      break
    end
  end
  if !found_N
    throw(DomainError("There is no path of canonical homomorphisms between the modules!"))
  end
  morphisms = Vector{ModuleFPHom}()
  A = N
  while A !== M
    f = parent_hom[A]
    push!(morphisms, f)
    A = codomain(f)
  end
  return morphisms
end

@doc Markdown.doc"""
    transport(M::SubQuo, v::SubQuoElem)

Map the element `v` to an element of the module `M` using cached 
canonical homomorphisms between the parent Module of `v` and `M`.
If this is not possible an error is thrown.
"""
function transport(M::SubQuo, v::SubQuoElem)
  N = parent(v)
  morphisms = find_sequence_of_morphisms(N, M)

  return foldl((x,f) -> f(x), morphisms; init=v)
end

@doc Markdown.doc"""
    find_morphism(M::SubQuo, N::SubQuo)

Find a morphism from `M` to `N` in the graph of cached (canonical) morphisms
and return it.
If there exists no such morphism an error is thrown.
"""
function find_morphism(M::SubQuo, N::SubQuo)
  morphisms = find_sequence_of_morphisms(M, N)
  return reduce(*, morphisms)
end

@doc Markdown.doc"""
    find_morphisms(N::SubQuo, M::SubQuo)

Traverse the graph of all cached morphisms originating in `N` and 
return the compositions of all loop-free paths from `N` to `M`.
"""
function find_morphisms(N::SubQuo, M::SubQuo)
  # from N to M

  all_paths = []

  function helper_dfs!(U::SubQuo, D::SubQuo, visited::Vector{<:ModuleFPHom}, path::Vector)
    if U === D
      push!(all_paths, path)
      return
    end
    for neighbor_morphism in U.outgoing_morphisms
      if findfirst(x->x===neighbor_morphism, visited) === nothing #if !(neighbor_morphism in visited) doesn't work since it uses == instead of ===
        helper_dfs!(codomain(neighbor_morphism), D, vcat(visited, [neighbor_morphism]), union(path, [neighbor_morphism]))
      end
    end
  end

  helper_dfs!(N, M, Vector{ModuleFPHom}(), [])

  morphisms = Vector{ModuleFPHom}()
  for path in all_paths
    phi = identity_map(N)
    for h in path
      phi = phi*h
    end
    push!(morphisms, phi)
  end

  return morphisms
end

#############################
# Useful functions
#############################

@doc Markdown.doc"""
    register_morphism!(f::ModuleFPHom)

Cache the morphism `f` in the corresponding caches of the domain and codomain of `f`.
"""
function register_morphism!(f::ModuleFPHom)
  push!(domain(f).outgoing_morphisms, f)
  push!(codomain(f).incoming_morphisms, f)
end

#############################
#TODO move to Hecke
#  re-evaluate and use or not

function getindex(r::Hecke.SRow, u::UnitRange)
  R = base_ring(r)
  s = sparse_row(R)
  shift = 1-first(u)
  for (p,v) = r
    if p in u
      push!(s.pos, p+shift)
      push!(s.values, v)
    end
  end
  return s
end

function getindex(r::Hecke.SRow, R::AbstractAlgebra.Ring, u::UnitRange)
  s = sparse_row(R)
  shift = 1-first(u)
  for (p,v) = r
    if p in u
      push!(s.pos, p+shift)
      push!(s.values, v)
    end
  end
  return s
end

function getindex(a::Hecke.SRow, b::AbstractVector{Int})
  if length(a.pos) == 0
    return a
  end
  m = minimum(b)
  b = sparse_row(parent(a.values[1]))
  for (k,v) = a
    if k in b
      push!(b.pos, k-b+1)
      push!(b.values, v)
    end
  end
  return b
end

@doc Markdown.doc"""
    sparse_row(A::MatElem)

Convert `A` to a sparse row. 
`nrows(A) == 1` must hold.
"""
function sparse_row(A::MatElem)
  @assert nrows(A) == 1
  if ncols(A) == 0
    return sparse_row(base_ring(A))
  end
  return Hecke.sparse_matrix(A)[1]
end

@doc Markdown.doc"""
    dense_row(r::Hecke.SRow, n::Int)

Convert `r[1:n]` to a dense row, that is an AbstractAlgebra matrix.
"""
function dense_row(r::Hecke.SRow, n::Int)
  R = base_ring(r)
  A = zero_matrix(R, 1, n)
  for i in intersect(r.pos, 1:n)
    A[1,i] = r[i]
  end
  return A
end

function default_ordering(F::FreeMod)
  if iszero(F)
    return default_ordering(base_ring(F))*ModuleOrdering(F, Orderings.ModOrdering(Vector{Int}(), :lex))
  end
  return default_ordering(base_ring(F))*lex(gens(F))
end

##############################
#should be in Singular.jl
function Singular.intersection(a::Singular.smodule, b::Singular.smodule)
  c = base_ring(a)
  return Singular.Module(c, Singular.libSingular.id_Intersection(a.ptr, b.ptr, c.ptr))
end

function _reduce(a::Singular.smodule, b::Singular.smodule)
  @assert b.isGB
  p = Singular.libSingular.p_Reduce(a.ptr, b.ptr, base_ring(b).ptr)
  return Singular.Module(base_ring(b), p)
end

function _reduce(a::Singular.svector, b::Singular.smodule)
  @assert b.isGB
  p = _reduce(Singular.Module(base_ring(b), a), b)[1] # TODO When available in Singular.jl use reduce(::svector, ::smodule)
  return Singular.Module(base_ring(b), p)[1]
end

#TODO: tensor_product from Raul's H is broken


######################################
# Migrating test
######################################
@doc Markdown.doc"""
    projection(F::FreeMod, indices::AbstractArray)

Return the canonical projection from $F = R^I$ to $R^(\texttt{indices})$ where $\texttt{indices} \subset I$.
"""
function projection(F::FreeMod, indices::AbstractArray)
  @assert all(x -> x <= ngens(F), indices)
  @assert length(Set(indices)) == length(indices) # unique indices
  R = base_ring(F)
  G = FreeMod(R, length(indices))
  return hom(F, G, Vector{elem_type(G)}([i in indices ? G[findfirst(x->x==i,indices)] : zero(G) for i=1:ngens(F)]))
end

@doc Markdown.doc"""
    preimage(H::SubQuoHom,N::SubQuo{T}, task::Symbol = :none) where {T}

Return the preimage of the submodule `N` under the morphism `H` 
as a subquotient, as well as the injection homomorphism into the domain of $H$.
"""
function preimage(H::SubQuoHom,N::SubQuo{T}, task::Symbol = :none) where {T}
  inclusion = get_attribute(N, :canonical_inclusion)
  if inclusion != nothing && codomain(inclusion) === codomain(H)
    elems = [inclusion(v) for v in gens(N)]
  else
    elems = [SubQuoElem(repres(v),codomain(H)) for v in gens(N)]
  end
  return preimage(H,elems,task)
end

@doc Markdown.doc"""
    preimage(H::SubQuoHom,elems::Vector{SubQuoElem{T}}, task::Symbol = :none) where {T}

Return the preimage of the submodule generated by the Elements `elems` under $H$
as a subquotient, as well as the injection homomorphism into the domain of $H$.
"""
function preimage(H::SubQuoHom,elems::Vector{SubQuoElem{T}}, task::Symbol = :none) where {T}
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
  elems_in_coker = map(x->i_cod_coker(x),elems)
  cokernel_modulo_elmes,projection = quo(cod_coker,elems_in_coker,:with_morphism)
  preimage, emb = kernel(H*i_cod_coker*projection)
  
  if task != :none
    return preimage, emb
  else
    return preimage
  end
end

@doc Markdown.doc"""
    matrix_kernel(A::MatElem)

Compute the kernel of `A` where `A` is considered as the correponding morphism
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

@doc Markdown.doc"""
    simplify_light(M::SubQuo)

Simplify the given subquotient `M` and return the simplified subquotient `N` along
with the injection map $N \to M$ and the projection map $M \to N$. These maps are 
isomorphisms.
The only simplifications which are done are the following: 
- Remove all generators which are represented by the zero element in the ambient 
free module.
- Remove all generators which are in the generating set of the relations.
- Remove all duplicates in the generators and relations sets.
"""
function simplify_light(M::SubQuo)
  M_gens = ambient_representatives_generators(M)
  M_rels = relations(M)

  N_rels = unique(filter(x -> !iszero(x), M_rels))
  N_gens = unique(setdiff(filter(x -> !iszero(x), M_gens), N_rels))

  N = length(N_rels) == 0 ? SubQuo(ambient_free_module(M), N_gens) : SubQuo(ambient_free_module(M), N_gens, N_rels)

  index_of_N_in_M = indexin(N_gens, M_gens)
  inj = hom(N, M, Vector{elem_type(M)}([M[index_of_N_in_M[i]] for i in 1:ngens(N)]))

  index_of_M_in_N = indexin(M_gens, N_gens)
  proj = hom(M, N, Vector{elem_type(N)}([index_of_M_in_N[i] === nothing ? zero(N) : N[index_of_M_in_N[i]] for i in 1:ngens(M)]))

  return N, inj, proj
end

@doc Markdown.doc"""
    simplify_with_same_ambient_free_module(M::SubQuo)

Simplify the given subquotient `M` and return the simplified subquotient `N` along
with the injection map $N \to M$ and the projection map $M \to N$. These maps are 
isomorphisms. The ambient free module of `N` is the same as that of `M`.
"""
function simplify_with_same_ambient_free_module(M::SubQuo)
  _, to_M, from_M = simplify(M)
  N, N_to_M = image(to_M)
  return N, N_to_M, hom(M, N, [N(coeffs(from_M(g))) for g in gens(M)])
  #return N, N_to_M, hom(M, N, [N(repres(g)) for g in gens(M)])
end

@doc Markdown.doc"""
    simplify(M::SubQuo)

Simplify the given subquotient `M` and return the simplified subquotient `N` along
with the injection map $N \to M$ and the projection map $M \to N$. These maps are 
isomorphisms. 
The simplifcation is heuristical and includes steps like for example removing
zero-generators or removing the i-th component of all vectors if those are 
reduced by a relation.
"""
function simplify(M::SubQuo)
  function standard_unit_vector_in_relations(i::Int, M::SubQuo)
    F = ambient_free_module(M)
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

  function rows_to_delete(A::MatElem, max_index::Int)
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

  to_delete,_,_ = rows_to_delete(transpose(vcat(new_generators, new_relations)),size(new_relations)[2])

  new_generators = delete_columns(new_generators, to_delete)
  new_relations = delete_columns(new_relations, to_delete)

  #remove rows
  #simplify relations
  to_delete,_,_ = rows_to_delete(new_relations, size(new_relations)[1])

  new_relations = delete_rows(new_relations, to_delete)

  #simplify generators
  to_delete, corresponding_row, K_gen = rows_to_delete(vcat(new_generators, new_relations), size(new_generators)[1])


  injection_matrix = delete_rows(identity_matrix(R, size(M_generators)[1]), to_delete)
  projection_matrix = zero_matrix(R, size(M_generators)[1], size(K_gen)[2]-length(to_delete))
  for i=1:size(M_generators)[1]
    if i in to_delete
      index = findfirst(x -> x==i, to_delete)
      assign_row!(projection_matrix, R(-1)*R(inv(coeff(K_gen[corresponding_row[index],i], 1)))*delete_columns(K_gen[corresponding_row[index],:], to_delete), i)
    else
      standard_unit_vector_index = i-length(filter(x -> x < i, to_delete))
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
    SQ = iszero(new_relations) ? SubQuo(SubModuleOfFreeModule(new_generators)) : SubQuo(new_generators, new_relations)
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
# Not only for testing
######################################
@doc Markdown.doc"""
    map(F::FreeMod{T}, A::MatrixElem{T}) where T

Converts a given $n \times m$-matrix into the corresponding morphism $A : R^n \to F$, 
with `rank(F) == m`.
"""
function map(F::FreeMod{T}, A::MatrixElem{T}) where {T <: RingElement}
  R = base_ring(F)
  F_domain = FreeMod(R, nrows(A))

  phi = FreeModuleHom(F_domain, F, A)
  return phi
end

@doc Markdown.doc"""
    map(A::MatElem)

Converts a given $n \times m$-matrix into the corresponding morphism $A : R^n \to R^m$.
"""
function map(A::MatElem)
  R = base_ring(A)
  F_codomain = FreeMod(R, ncols(A))
  return map(F_codomain,A)
end

@doc Markdown.doc"""
    is_injective(f::ModuleFPHom)

Test if `f` is injective.
"""
function is_injective(f::ModuleFPHom)
  return iszero(kernel(f)[1])
end

@doc Markdown.doc"""
    is_surjective(f::ModuleFPHom)

Test if `f` is surjective.
"""
function is_surjective(f::ModuleFPHom)
  return image(f)[1] == codomain(f)
end

@doc Markdown.doc"""
    is_bijective(f::ModuleFPHom)

Test if `f` is bijective.
"""
function is_bijective(f::ModuleFPHom)
  return is_injective(f) && is_surjective(f)
end


#############################################
# Test hom
#############################################
function to_julia_matrix(A::Union{MatElem})
  return eltype(A)[A[i, j] for i = 1:nrows(A), j = 1:ncols(A)]
end

function copy_and_reshape(M::MatElem, n, m)
  julia_matrix = to_julia_matrix(M)
  julia_matrix = reshape(julia_matrix, n, m)
  R = base_ring(M)
  mat_space = MatrixSpace(R, n, m)
  return mat_space(julia_matrix)
end

function hom_matrices_helper(f1::MatElem{T}, g1::MatElem{T}) where T
  R = base_ring(f1)
  s1, s0 = size(f1)
  t1, t0 = size(g1)

  g2 = matrix_kernel(g1)
  n = s1*t0
  m = s0*t0 + s1*t1
  delta::MatrixElem{T} = zero_matrix(R, m,n)
  for j=1:s0*t0
    b_vector::MatrixElem{T} = zero_matrix(R, 1,s0*t0)
    b_vector[1,j] = R(1)
    A = copy_and_reshape(b_vector, s0, t0)
    res = f1*A
    delta[j,:] = copy_and_reshape(res, 1, n)
  end
  for j=s0*t0+1:m
    b_vector::MatrixElem{T} = zero_matrix(R, 1,m-s0*t0)
    b_vector[1,j-s0*t0] = R(1)
    B = copy_and_reshape(b_vector, s1, t1)
    res = -B*g1
    delta[j,:] = copy_and_reshape(res, 1,n)
  end

  gamma = matrix_kernel(delta)
  gamma = gamma[:,1:s0*t0]

  rho::MatrixElem{T} = zero_matrix(R, s0*t1, s0*t0)

  for j=1:s0*t1
    b_vector = zero_matrix(R, 1,s0*t1)
    b_vector[1,j] = R(1)
    C = copy_and_reshape(b_vector, s0, t1)
    res1 = C*g1
    rho[j,1:length(res1)] = copy_and_reshape(res1, 1,length(res1))
  end

  M = SubQuo(gamma, rho)

  function convert_to_matrix(v::SubQuoElem{T})
    if parent(v) !== M
      throw(DomainError("v does not represent a homomorphism"))
    end
    R = base_ring(M)
    A = copy_and_reshape(dense_row(v.repres.coords[R, 1:s0*t0], s0*t0), s0, t0)
    return A
  end

  return M, convert_to_matrix
end


# Hom and hom_matrices implement the same mathematical algorithm, but the implementations 
# differ a lot e.g. in the data structures. As a result performance differs depending 
# on the example favoring the one or the other. So it makes sense to offer both. 
# With option :matrices in hom() hom_matrices is used.
@doc Markdown.doc"""
    hom_matrices(M::SubQuo{T},N::SubQuo{T}) where T

Return a subquotient $S$ such that $\text{Hom}(M,N) \cong S$
"""
function hom_matrices(M::SubQuo{T},N::SubQuo{T},simplify_task=true) where T
  f1 = generator_matrix(present_as_cokernel(M).quo)
  g1 = generator_matrix(present_as_cokernel(N).quo)
  R = base_ring(M)
  SQ, convert_to_matrix = hom_matrices_helper(f1,g1)
  if simplify_task
    SQ2, i, p = simplify(SQ)
    to_homomorphism = function(elem::SubQuoElem{T})
      elem2 = i(elem)
      A = convert_to_matrix(elem2)
      return SubQuoHom(M,N,A)
    end
    to_subquotient_elem = function(H::ModuleFPHom)
      m = length(matrix(H))
      v = copy_and_reshape(matrix(H),1,m)
      v = FreeModElem(sparse_row(v), FreeMod(R, length(v)))
      return p(SQ(v))
    end

    to_hom_map = MapFromFunc(to_homomorphism, to_subquotient_elem, SQ2, Hecke.MapParent(M, N, "homomorphisms"))
    set_attribute!(SQ2, :hom => (M, N), :module_to_hom_map => to_hom_map)

    return SQ2, to_hom_map
  else
    to_subquotient_elem = function(H::ModuleFPHom)
      m = length(matrix(H))
      v = copy_and_reshape(matrix(H),1,m)
      v = FreeModElem(sparse_row(v), FreeMod(R, length(v)))
      return SQ(v)
    end
    to_homomorphism = function(elem::SubQuoElem{T})
      A = convert_to_matrix(elem)
      return SubQuoHom(M,N,A)
    end

    to_hom_map = MapFromFunc(to_homomorphism, to_subquotient_elem, SQ, Hecke.MapParent(M, N, "homomorphisms"))
    set_attribute!(SQ, :hom => (M, N), :module_to_hom_map => to_hom_map)

    return SQ, to_hom_map
  end
end

function change_base_ring(S::Ring, F::FreeMod)
  R = base_ring(F)
  r = ngens(F)
  FS = FreeMod(S, r)
  map = hom(F, FS, gens(FS), MapFromFunc(x->S(x), R, S))
  return FS, map
end

function change_base_ring(f::Hecke.Map{DomType, CodType}, F::FreeMod) where {DomType<:Ring, CodType<:Ring}
  domain(f) == base_ring(F) || error("ring map not compatible with the module")
  S = codomain(f)
  r = ngens(F)
  FS = FreeMod(S, r)
  map = hom(F, FS, gens(FS), f)
  return FS, map
end

function change_base_ring(S::Ring, M::SubQuo)
  F = ambient_free_module(M)
  R = base_ring(M)
  FS, mapF = change_base_ring(S, F)
  g = ambient_representatives_generators(M)
  rels = relations(M)
  MS = SubQuo(FS, mapF.(g), mapF.(rels))
  map = SubQuoHom(M, MS, gens(MS), MapFromFunc(x->S(x), R, S))
  return MS, map
end

function change_base_ring(f::Hecke.Map{DomType, CodType}, M::SubQuo) where {DomType<:Ring, CodType<:Ring}
  domain(f) == base_ring(M) || error("ring map not compatible with the module")
  S = codomain(f)
  F = ambient_free_module(M)
  R = base_ring(M)
  FS, mapF = change_base_ring(S, F)
  g = ambient_representatives_generators(M)
  rels = relations(M)
  MS = SubQuo(FS, mapF.(g), mapF.(rels))
  map = SubQuoHom(M, MS, gens(MS), f)
  return MS, map
end
