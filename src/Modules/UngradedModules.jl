export presentation, coords, coeffs, repres, cokernel, index_of_gen, sub,
      quo, presentation, present_as_cokernel, is_equal_with_morphism, 
      groebner_basis, reduced_groebner_basis, leading_module, 
      show_morphism, hom_tensor, hom_prod_prod, coordinates, 
      represents_element, free_resolution, homomorphism, module_elem, generator_matrix,
      restrict_codomain, restrict_domain, direct_product, tensor_product, 
      free_module, tor, lift_homomorphism_contravariant, lift_homomorphism_covariant, 
      ext, map_canonically, all_canonical_maps, register_morphism!, dense_row, 
      matrix_kernel, simplify, map, isinjective, issurjective, isbijective, iswelldefined,
      subquotient, ambient_free_module, ambient_module, ambient_representative, 
      ambient_representatives_generators, relations

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

Return the free module $R^p$, created with its basis of standard unit vectors.

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

@doc Markdown.doc"""
    iszero(F::AbstractFreeMod)

Return `true` if `F` is the zero module, `false` otherwise.
"""
function iszero(F::AbstractFreeMod)
  return rank(F) == 0
end

###############################################################################
# FreeModElem constructors
###############################################################################

@doc Markdown.doc"""
    FreeModElem(c::SRow{T}, parent::FreeMod{T}) where T

Return the element of  `F`  whose coefficients with respect to the basis of 
standard unit vectors of `F` are given by the entries of `c`.
"""
FreeModElem(c::SRow{T}, parent::FreeMod{T}) where T = FreeModElem{T}(c, parent)

@doc Markdown.doc"""
    FreeModElem(c::Vector{T}, parent::FreeMod{T}) where T
    
Return the element of  `F`  whose coefficients with respect to the basis of 
standard unit vectors of `F` are given by the entries of `c`.
"""
function FreeModElem(c::Vector{T}, parent::FreeMod{T}) where T
  @assert length(c) == rank(parent)
  sparse_coords = sparse_row(base_ring(parent), collect(1:rank(parent)), c)
  return FreeModElem{T}(sparse_coords,parent)
end

@doc Markdown.doc"""
    (F::FreeMod{T})(c::SRow{T}) where T
    
Return the element of  `F`  whose coefficients with respect to the basis of 
standard unit vectors of `F` are given by the entries of `c`.
"""
function (F::FreeMod{T})(c::SRow{T}) where T
  return FreeModElem(c, F)
end

@doc Markdown.doc"""
    (F::FreeMod{T})(c::Vector{T}) where T

Return the element of  `F`  whose coefficients with respect to the basis of 
standard unit vectors of `F` are given by the entries of `c`.

# Examples
```jldoctest
julia> R, (x,y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = FreeMod(R,3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y over Rational Field

julia> V = [x, zero(R), y]
3-element Vector{fmpq_mpoly}:
 x
 0
 y

julia> f = F(V)
x*e[1] + y*e[3]
```
"""
function (F::FreeMod{T})(c::Vector{T}) where T 
 return FreeModElem(c, F)
end

function (F::AbstractFreeMod{T})(v::AbstractFreeModElem{T}) where T
  @assert parent(v) === F
  return v
end

function in(v::AbstractFreeModElem, M::ModuleFP)
  return parent(v) === M
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

# scalar multiplication with polynomials, integers
function *(a::MPolyElem_dec, b::AbstractFreeModElem)
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
  end
  return FreeModElem(a*coords(b), parent(b))
end
function *(a::MPolyElem, b::AbstractFreeModElem) 
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
  end
  return FreeModElem(a*coords(b), parent(b))
end
function *(a::RingElem, b::AbstractFreeModElem) 
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
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

Construct `ModuleGens` from an array of Oscar free module elements, specifying the Oscar free module 
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

Construct `ModuleGens` from an array of Oscar free module elements.

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

Construct `ModuleGens` from an array of Oscar free module elements, specifying the Oscar free module.

!!! note

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}) where {T}
  #SF = singular_module(F)
  #return ModuleGens{T}(O, F, SF)
  return ModuleGens{T}(O, F)
end

@doc Markdown.doc"""
    ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, ordering::ModuleOrdering) where {T}

Construct `ModuleGens` from an array of Oscar free module elements, specifying the Oscar free module.
Moreover, the ordering is defined by `ordering1` and `ordering2`. Allowed symbols are those that are 
allowed for Singular polynomial rings.

!!! note

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, F::FreeMod{T}, ordering::ModuleOrdering) where {T}
  #SF = singular_module(F, ordering)
  #return ModuleGens{T}(O, F, SF)
  return ModuleGens{T}(O, F)
end

@doc Markdown.doc"""
    ModuleGens(O::Vector{<:FreeModElem}, SF::Singular.FreeMod)

Construct `ModuleGens` from an array of Oscar free module elements, specifying the Singular free module.

!!! note 

    The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModElem}, SF::Singular.FreeMod)
  return ModuleGens{elem_type(base_ring(parent(O[1])))}(O, parent(O[1]), SF)
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
    singular_generators(M::ModuleGB)

Return the elements of `M` from Singular side.
"""
function singular_generators(M::ModuleGB)
  return singular_generators(M.groebner_basis)
end

@doc Markdown.doc"""
    oscar_generators(M::ModuleGens)  

Return the generators of `M` from the OSCAR side.
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
  print(io, "Array of length ", length(F))
  for i=1:length(F)
    if isassigned(F.O, i)
      print(io, "\n", i, " -> ", F.O[i])
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

If fields of `F` from the OSCAR side are not defined, they
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
are computed, given the OSCAR side.
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

Create a Singular module from an OSCAR free module.
"""
function singular_module(F::FreeMod)
  Sx = singular_poly_ring(base_ring(F), keep_ordering=false)
  return Singular.FreeModule(Sx, dim(F))
end

@doc Markdown.doc"""
    singular_module(F::FreeMod, ordering::ModuleOrdering)

Create a Singular module from an OSCAR free module over the given Singular polynomial ring.
"""
function singular_module(F::FreeMod, ordering::ModuleOrdering)
  # @assert exactly_one_module_ordering(ordering) # Remove
  Sx = singular_ring(base_ring(F), singular(ordering))
  return Singular.FreeModule(Sx, dim(F))
end

@doc Markdown.doc"""
    (SF::Singular.FreeMod)(m::FreeModElem)

Convert an OSCAR free module element to the Singular side.
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

Convert a Singular vector to a free module element on the OSCAR side.
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

im_gens(f::FreeModuleHom) = gens(image(f)[1])
base_ring_map(f::FreeModuleHom) = f.ring_map

@doc Markdown.doc"""
    matrix(a::FreeModuleHom)

Given a homomorphism `a` of type  `FreeModuleHom` with domain `F`
and codomain `M`, return a matrix `A` with `rank(F)` rows and 
`ngens(M)` columns such that `a == hom(F, M, A)`.
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
    hom(F::FreeMod, M::ModuleFP, V::Vector{<:ModuleFPElem})

Given a vector `V` of `rank(F)` elements of `M`, 
return the homomorphism `F` $\to$ `M` which sends the `i`-th
basis vector of `F` to the `i`-th entry of `V`.

    hom(F::FreeMod, M::ModuleFP{T}, A::MatElem{T}) where T

Given a matrix `A` with `rank(F)` rows and `ngens(M)` columns, return the
homomorphism `F` $\to$ `M` which sends the `i`-th basis vector of `F` to 
the linear combination $\sum_j A[i,j]*M[j]$ of the generators `M[j]` of `M`.
"""
function hom(F::FreeMod, M::ModuleFP, V::Vector{<:ModuleFPElem}) 
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
hom(F::FreeMod, M::ModuleFP, V::Vector{<:ModuleFPElem}, h::RingMapType) where {RingMapType} = FreeModuleHom(F, M, V, h)
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
# Groebner basis constructors
###############################################################################
@doc Markdown.doc"""
    ModuleGB(gb::ModuleGens{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering)
  
Create a Groebner basis, where `gb` contains the Groebner basis, 
`leading_monomial` is the corresponding leading monomials. All is with respect to
the ordering `ordering`.
"""
ModuleGB(gb::ModuleGens{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering) where {T} =
          ModuleGB{T}(gb, leading_monomials, ordering)

@doc Markdown.doc"""
    ModuleGB(gb::ModuleGens{T}, quotient_gb::ModuleGB{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering)
  
Create a relative Groebner basis (relative to `quotient_gb`), where `gb` contains the Groebner basis, 
`leading_monomial` is the (to `gb`) corresponding leading monomials. All is with respect to
the ordering `ordering`.
"""
ModuleGB(gb::ModuleGens{T}, quotient_gb::ModuleGB{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering) where {T} =
          ModuleGB{T}(gb, quotient_gb, leading_monomials, ordering)


@doc Markdown.doc"""
    ReducedModuleGB(gb::ModuleGens{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering)
      
Create a reduced Groebner basis, where `gb` contains the reduced Groebner basis, 
`leading_monomial` is the corresponding leading monomials. All is with respect to
the ordering `ordering`.
"""
function ReducedModuleGB(gb::ModuleGens{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering) where {T}
  gb = ModuleGB{T}(gb, leading_monomials, ordering)
  gb.reduced_groebner_basis = gb.groebner_basis
  if isdefined(gb, :quo_groebner_basis)
    gb.quo_groebner_basis.reduced_groebner_basis = gb.quo_groebner_basis.groebner_basis
  end
  return gb
end

@doc Markdown.doc"""
    ReducedModuleGB(gb::ModuleGens{T}, quotient_gb::ModuleGB{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering)
      
Create a reduced relative Groebner basis (relative to `quotient_gb`), where `gb` contains the reduced Groebner basis, 
`leading_monomial` is the corresponding leading monomials. All is with respect to
the ordering `ordering`.
"""
function ReducedModuleGB(gb::ModuleGens{T}, quotient_gb::ModuleGB{T}, leading_monomials::ModuleGens{T}, ordering::ModuleOrdering) where {T}
  gb = ModuleGB{T}(gb, quotient_gb, leading_monomials, ordering)
  gb.reduced_groebner_basis = gb.groebner_basis
  if isdefined(gb, :quo_groebner_basis)
    gb.quo_groebner_basis.reduced_groebner_basis = gb.quo_groebner_basis.groebner_basis
  end
  return gb
end

function exactly_one_module_ordering(ordering::Singular.sordering)
  count_module_ordering = 0
  for ord in ordering.data
    count_module_ordering += Int(ord.order == Singular.ringorder_c || ord.order == Singular.ringorder_C)
  end
  return count_module_ordering == 1
end

function ordering_compatible_with_ring(R::MPolyRing, ordering::Singular.sordering)
  return sum([ord.size for ord in ordering.data]) == nvars(R)
end

function ordering_is_global(R::MPolyRing, ordering::Singular.sordering)
  return Singular.has_global_ordering(singular_ring(R, ordering))
end


###############################################################################
# SubModuleOfFreeModule
###############################################################################
@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModElem}, default_ordering::ModuleOrdering = default_ordering(F)) where {R}

Construct the submodule of `F` generated by the elements of `gens` (the elements of 
`gens` must live in `F`).
"""
SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModElem}, default_ordering::ModuleOrdering = default_ordering(F)) where {R} = 
                        SubModuleOfFreeModule{R}(F, gens, default_ordering)    

SubModuleOfFreeModule(F::FreeMod{R}, singular_module::Singular.smodule, default_ordering::ModuleOrdering = default_ordering(F)) where {R} = 
                        SubModuleOfFreeModule{R}(F, singular_module, default_ordering)

@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{R}, gens::ModuleGens{R}, default_ordering::ModuleOrdering = default_ordering(F)) where {R}

Construct the submodule of `F` generated by `gens`.
"""
SubModuleOfFreeModule(F::FreeMod{R}, gens::ModuleGens{R}, default_ordering::ModuleOrdering = default_ordering(F)) where {R} = 
                        SubModuleOfFreeModule{R}(F, gens, default_ordering)

@doc Markdown.doc"""
    SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}, default_ordering::ModuleOrdering = default_ordering(F)) where {L}

Construct the submodule generated by the rows of `A`. The embedding free
module is `F`. In particular, `rank(F) == ncols(A)` must hold.
"""
SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}, default_ordering::ModuleOrdering = default_ordering(F)) where {L} = 
                        SubModuleOfFreeModule{L}(F, A, default_ordering)

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
  ord = default_ordering(F)
  return SubModuleOfFreeModule{L}(F, A, ord)
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
  return SubModuleOfFreeModule{L}(F, A, default_ordering)
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
    get_default_ordering(M::SubModuleOfFreeModule)

Get the default ordering of `M`.
"""
function default_ordering(M::SubModuleOfFreeModule)
  return M.default_ordering
end

@doc Markdown.doc"""
    set_default_ordering!(M::SubModuleOfFreeModule, ord::ModuleOrdering)

Set the default ordering in `M` to `ord`.
"""
function set_default_ordering!(M::SubModuleOfFreeModule, ord::ModuleOrdering)
  # @assert exactly_one_module_ordering(ord) # Remove
  # @assert ordering_compatible_with_ring(base_ring(M), ord) # Remove
  M.default_ordering = ord
end

@doc Markdown.doc"""
    std_basis(submod::SubModuleOfFreeModule)

Compute a standard basis of `submod` with respect to its default odering.
"""
function std_basis(submod::SubModuleOfFreeModule)
  return std_basis(submod, default_ordering(submod))
end

@doc Markdown.doc"""
    groebner_basis(submod::SubModuleOfFreeModule)

Compute a Gröbner basis of `submod` with respect to its default odering.
For this, the default ordering must be global.
"""
function groebner_basis(submod::SubModuleOfFreeModule)
  return groebner_basis(submod, default_ordering(submod))
end

@doc Markdown.doc"""
    std_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering)

Compute a standard basis of `submod` with respect to the given `ordering`.
"""
function std_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering)
  if !haskey(submod.groebner_basis, ordering)
    gb = std_basis_as_ModuleGens(submod, ordering)
    lm = leading_module_as_ModuleGens(submod, ordering)
    submod.groebner_basis[ordering] = ModuleGB(gb, lm, ordering)
  end
  return submod.groebner_basis[ordering]
end

@doc Markdown.doc"""
    groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering)

Compute a Gröbner of `submod` with respect to the given `ordering`.
The ordering must be global.
"""
function groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering)
  @assert isglobal(ordering)
  return std_basis(submod, ordering)
end

@doc Markdown.doc"""
    reduced_groebner_basis(submod::SubModuleOfFreeModule)

Compute a reduced Gröbner basis. The return type is `ModuleGens`.
"""
function reduced_groebner_basis(submod::SubModuleOfFreeModule)
  return reduced_groebner_basis(submod, default_ordering(submod))
end

@doc Markdown.doc"""
    reduced_groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering)

Compute a reduced Gröbner basis with respect to the given `ordering`. The return type is `ModuleGens`.
"""
function reduced_groebner_basis(submod::SubModuleOfFreeModule, ordering::ModuleOrdering)
  @assert isglobal(ordering)
  
  if !haskey(submod.groebner_basis, ordering)
    gb = std_basis_as_ModuleGens(submod, ordering, true)
    lm = leading_module_as_ModuleGens(submod, ordering)
    submod.groebner_basis[ordering] = ReducedModuleGB(gb, lm, ordering)
  else
    if !isdefined(submod.groebner_basis[ordering], :reduced_groebner_basis)
      # TODO compute reduced GB from given GB
      red_gb = std_basis_as_ModuleGens(submod, ordering, true)
      submod.groebner_basis[ordering].reduced_groebner_basis = red_gb
    end
  end
  return submod.groebner_basis[ordering]
end

function leading_module(submod::SubModuleOfFreeModule)
  return leading_module(submod, default_ordering(submod))
end

function leading_module(submod::SubModuleOfFreeModule, ordering::ModuleOrdering)
  gb = std_basis(submod, ordering)
  lm = gb.leading_monomials
  return SubModuleOfFreeModule(submod.F, lm, ordering)
end

@doc Markdown.doc"""
    std_basis_as_ModuleGens(submod::SubModuleOfFreeModule)

Compute a standard basis of `submod` and return it as a `ModuleGens`.
"""
function std_basis_as_ModuleGens(submod::SubModuleOfFreeModule)
  return std_basis_as_ModuleGens(submod, default_ordering(submod))
end

@doc Markdown.doc"""
    std_basis_as_ModuleGens(submod::SubModuleOfFreeModule, ordering::ModuleOrdering, reduced::Bool=false)

Compute a stanard basis of `submod` with respect to the given ordering.
Allowed orderings are those that are allowed as orderings for Singular polynomial rings.
In case `reduced` is `true` and the ordering is global, a reduced Gröbner basis is computed.
"""
function std_basis_as_ModuleGens(submod::SubModuleOfFreeModule, ordering::ModuleOrdering, reduced::Bool=false)
  if reduced
    @assert isglobal(ordering)
  end
  mg = ModuleGens(oscar_generators(submod.gens), submod.F , ordering)
  gb = std_basis(mg, reduced)
  oscar_assure(gb)
  return gb
end

function leading_module_as_ModuleGens(M::SubModuleOfFreeModule)
  ord = default_ordering(M)
  return leading_module_as_ModuleGens(M, ord)
end

function leading_module_as_ModuleGens(M::SubModuleOfFreeModule, ordering::ModuleOrdering)
  mg = ModuleGens(oscar_generators(M.gens), M.F , ordering)
  lm = leading_module(mg)
  oscar_assure(lm)
  return lm
end

function show(io::IO, GB::ModuleGB)
  if isdefined(GB, :quo_groebner_basis)
    show_relative_groebner_basis(io, GB.groebner_basis, GB.quo_groebner_basis.groebner_basis)
    if isdefined(GB, :reduced_groebner_basis)
      print(io, "\n")
      show_relative_groebner_basis(io, GB.reduced_groebner_basis, GB.quo_groebner_basis.reduced_groebner_basis)
    end
  else
    show_groebner_basis(io, GB.groebner_basis)
    if isdefined(GB, :reduced_groebner_basis)
      print(io, "\n\n")
      show_groebner_basis(io, GB.reduced_groebner_basis, true)
    end
  end
end

function show_groebner_basis_helper(io::IO, sub::ModuleGens, init::String)
  F = sub.F

  print(io, init)

  # This shall be removed when Oscar orderings are used
  vectors = [[] for _ in 1:length(sub.O)]
  for i in 1:length(sub.O)
    monomials = []
    v = sub.S[i]
    l = Singular.lead(v)
    while !iszero(l)
      push!(monomials, l)
      v = v - l
      l = Singular.lead(v)
    end
    monomials = Vector{FreeModElem}(map(x -> F(x), monomials))
    vectors[i] = monomials
  end

  for v in vectors
    print(io, "\n")
    print(io, v[1])
    for i in 2:length(v)
      print(io, " + ", v[i])
    end
  end

  #=for v in sub.O    
    print(io, "\n", v)
  end=#
end

function show_groebner_basis(io::IO, sub::ModuleGens, reduced::Bool = false)
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
    init = "Relative "
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
    isgenerated_by_standard_unit_vectors(M::SubModuleOfFreeModule)

Check if `M` is generated by the standard unit vectors.
"""
function isgenerated_by_standard_unit_vectors(M::SubModuleOfFreeModule)
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
  if M.F !== N.F
    return false
  end
  #TODO should there be a check for === up to permutation in order to avoid std-computation?
  M_mod_N = reduce(M, N)
  N_mod_M = reduce(N, M)
  return iszero(M_mod_N) && iszero(N_mod_M)
end

#+(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule) = sum(M, N)

###############################################################################
# SubQuo constructors
###############################################################################

@doc Markdown.doc"""
    SubQuo(sub::SubModuleOfFreeModule{R}) where {R}

Construct the module `sub` as a subquotient.
"""
SubQuo(sub::SubModuleOfFreeModule{R}) where {R} = SubQuo{R}(sub)

@doc Markdown.doc"""
    SubQuo(sub::SubModuleOfFreeModule{R}, quo::SubModuleOfFreeModule{R}) where {R}

Construct the subquotient module $\texttt{sub} + texttt{quo} / \texttt{quo}$. 
`sub` and `quo` must be submodules of the same free module.
"""
SubQuo(sub::SubModuleOfFreeModule{R}, quo::SubModuleOfFreeModule{R}) where {R} = SubQuo{R}(sub, quo)

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
SubQuo(F::FreeMod{R}, O::Vector{<:FreeModElem}) where {R} = SubQuo{R}(F, O)

@doc Markdown.doc"""
    SubQuo(S::SubQuo{L}, O::Vector{<:FreeModElem}) where {L}

Construct a subquotient where the generators are those of `S` and the relations are 
the union of `O` and the relations of `S`.
"""
SubQuo(S::SubQuo{L}, O::Vector{<:FreeModElem}) where {L} = SubQuo{L}(S, O)

SubQuo(F::FreeMod{R}, s::Singular.smodule) where {R} = SubQuo{R}(F, s)

SubQuo(F::FreeMod{R}, s::Singular.smodule, t::Singular.smodule) where {R} = SubQuo{R}(F, s, t)

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

#######################################################
@doc Markdown.doc"""
    subquotient(a::FreeModuleHom, b::FreeModuleHom)

Given homomorphisms `a` and `b` between free modules such that 
`codomain(a) === codomain(b)`, 
return $(\text{im } a + \text{im } b)/\text{im } b$.

    subquotient(F::FreeMod{T}, A::MatElem{T}, B::MatElem{T}) where T

Given matrices `A` and `B` with rank `F` columns, return  
$(\text{im } a + \text{im } b)/\text{im } b$,
where `a` and `b` are the free module homomorphisms with codomain `F` represented by `A` and `B`.

    subquotient(A::MatElem{T}, B::MatElem{T}) where T

Given matrices `A` and `B` with the same number of columns, create a free module `F`, say, whose rank 
is that number, and return $(\text{im } a + \text{im } b)/\text{im } b$, where `a` and `b` are as above.

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
function subquotient(
    a::FreeModuleHom{<:Any, <:Any, Nothing}, 
    b::FreeModuleHom{<:Any, <:Any, Nothing}
  )
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
    if isgenerated_by_standard_unit_vectors(SQ.sub)
      println("Cokernel of ", generator_matrix(SQ.quo))
    else
      println("Subquotient with of image of")
      display(generator_matrix(SQ.sub))
      println("by image of")
      display(generator_matrix(SQ.quo))
    end
  else
    println("Image of ", generator_matrix(SQ.sub))
  end
end

@doc Markdown.doc"""
    cokernel(a::FreeModuleHom)

Given a homomorphism `a` between free modules, return the cokernel of `a` as an object of type `SubQuo`.

    cokernel(F::FreeMod{T}, A::MatElem{T}) where T

Given a matrix `A` with with rank `F` columns, return `cokernel(a)`, 
where `a` is the free module homomorphism with codomain `F` represented by `A`.

    cokernel(A::MatElem)

Create a free module `F`, say, with rank `F` columns, and return `cokernel(F, A)`.

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

julia> cokernel(a)
Subquotient of Submodule with 2 generators
1 -> e[1]
2 -> e[2]
by Submodule with 3 generators
1 -> y*e[1]
2 -> x*e[1] + y*e[2]
3 -> z*e[2]
```
"""
function cokernel(f::FreeModuleHom)
  @assert codomain(f) isa FreeMod
  return quo(codomain(f), image(f)[1])
end

@doc Markdown.doc"""
    cokernel(F::FreeMod{R}, A::MatElem{R}) where R

Let $F = R^m$ and $A$ an $n \times m$-matrix. Return the subquotient $F / \im(A)$.
"""
function cokernel(F::FreeMod{R}, A::MatElem{R}) where R
  return cokernel(map(F,A))
end

@doc Markdown.doc"""
    cokernel(A::MatElem)

Let $A$ be an $n \times m$-matrix over the ring $R$. Then return the subquotient $R^m / \im(A)$. 

!!! note

    The free module $R^m$ is constructed by the function and therefore not compatible with 
    free modules $R^m$ that are defined by the user or by other functions. For compatibility,
    use `cokernel(F::FreeMod{R}, A::MatElem{R})`.
"""
function cokernel(A::MatElem)
  return cokernel(map(A))
end

@doc Markdown.doc"""
    default_ordering(M::SubQuo)    

Return the default ordering of `M`.
"""
function default_ordering(M::SubQuo)
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

function std_basis(M::SubQuo, ord::ModuleOrdering)
  if !haskey(M.groebner_basis, ord)
    if isdefined(M, :quo)
      quo_gb = std_basis(M.quo, ord)
      sub_union_gb_of_quo = SubModuleOfFreeModule(M.F, ModuleGens(vcat(M.sub.gens.O, quo_gb.groebner_basis.O), M.F))
      gb = std_basis(sub_union_gb_of_quo, ord)
      rel_gb_list = Vector{FreeModElem}()

      for i in 1:length(gb.groebner_basis.O)
        v = gb.groebner_basis.S[i]
        if !iszero(_reduce(v, quo_gb.groebner_basis.S))
          push!(rel_gb_list, gb.groebner_basis.O[i])
        end
      end

      rel_gb_ModuleGens = ModuleGens(rel_gb_list, M.F)
      M.groebner_basis[ord] = ModuleGB(rel_gb_ModuleGens, quo_gb, gb.leading_monomials, ord)
    else
      M.groebner_basis[ord] = std_basis(M.sub, ord)
    end
  end
  return M.groebner_basis[ord]
end

function groebner_basis(M::SubQuo, ord::ModuleOrdering)
  @assert isglobal(ord)
  return std_basis(M, ord)
end

function std_basis(M::SubQuo)
  return std_basis(M, default_ordering(M))
end

function groebner_basis(M::SubQuo)
  return groebner_basis(M, default_ordering(M))
end

function reduced_groebner_basis(M::SubQuo)
  return reduced_groebner_basis(M, default_ordering(M))
end

function reduced_groebner_basis(M::SubQuo, ord::ModuleOrdering)
  if !haskey(M.groebner_basis, ord) || !isdefined(M.groebner_basis[ord], :reduced_groebner_basis)
    if !isdefined(M, :quo)
      M.groebner_basis[ord] = reduced_groebner_basis(M.sub, ord)
    else
      quo_gb = reduced_groebner_basis(M.quo, ord)
      sub_union_gb_of_quo = SubModuleOfFreeModule(M.F, ModuleGens(vcat(M.sub.gens.O, quo_gb.reduced_groebner_basis.O), M.F))
      gb = reduced_groebner_basis(sub_union_gb_of_quo, ord)
      rel_gb_list = Vector{FreeModElem}()

      for i in 1:length(gb.groebner_basis.O)
        v = gb.groebner_basis.S[i]
        if !iszero(_reduce(v, quo_gb.groebner_basis.S))
          push!(rel_gb_list, gb.groebner_basis.O[i])
        end
      end

      rel_gb_ModuleGens = ModuleGens(rel_gb_list, M.F)
      if !haskey(M.groebner_basis, ord)
        M.groebner_basis[ord] = ReducedModuleGB(rel_gb_ModuleGens, quo_gb, gb.leading_monomials, ord)
      else
        M.groebner_basis[ord].reduced_groebner_basis = rel_gb_ModuleGens
        M.groebner_basis[ord].quo_groebner_basis.reduced_groebner_basis = quo_gb.groebner_basis
      end
    end
  end
  return M.groebner_basis[ord]
end

function leading_module(M::SubQuo, ord::ModuleOrdering)
  @assert !isdefined(M, :quo)
  return SubQuo(leading_module(M.sub, ord))
end

function leading_module(M::SubQuo)
  return leading_module(M, default_ordering(M))
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

@doc Markdown.doc"""
    ==(M::SubQuo{T}, N::SubQuo{T}) where {T}

Given subquotients `M` and `N` such that `ambient_module(M) == ambient_module(N)`,
return `true` if `M` equals `N`, where `M` and `N` are regarded as submodules 
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

julia> M == N
false
```
"""
function (==)(M::SubQuo{T}, N::SubQuo{T}) where {T} # TODO replace implementation by two inclusion checks?
  if !isdefined(M, :quo) 
    if !isdefined(N, :quo)
      return M.sub == N.sub
    else
      return iszero(N.quo) && M.sub == N.sub
    end
  else
    if !isdefined(N, :quo)
      return iszero(M.quo) && M.sub == N.sub
    else
      return M.quo == N.quo && M.sum == N.sum
    end
  end
end

@doc Markdown.doc"""
    sum(M::SubQuo{T},N::SubQuo{T}) where T

Given subquotients `M` and `N` such that `ambient_module(M) == ambient_module(N)`,
return the sum of `M` and `N` regarded as submodules of the common ambient module,
together with the inclusion maps `M` $\to$ `M + N` and `N` $\to$ `M + N`.

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
  #TODO use SubModuleOfFreeModule instead of matrices
  M_quo = isdefined(M, :quo) ? M.quo : SubModuleOfFreeModule(ambient_free_module(M), Vector{FreeModElem}())
  N_quo = isdefined(N, :quo) ? N.quo : SubModuleOfFreeModule(ambient_free_module(N), Vector{FreeModElem}())
  
  if M_quo == N_quo
    SQ = SubQuo(sum(M.sub,N.sub),M_quo)
    iM = SubQuoHom(M,SQ,[SQ[i] for i=1:ngens(M)])
    iN = SubQuoHom(N,SQ,[SQ[i] for i=ngens(M)+1:ngens(SQ)])

    register_morphism!(iM)
    register_morphism!(iN)

    return SQ, iM, iN
  end
  throw(ArgumentError("img(M.relations) != img(N.relations)"))
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
return the intersection of `M` and `N` regarded as submodules of the common ambient module,
together with the inclusion maps `M` $\cap$ `N` $\to$ `M` and `M` $\cap$ `N` $\to$ `N`.

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
  M_quo = isdefined(M, :quo) ? M.quo : Oscar.SubModuleOfFreeModule(ambient_free_module(M), Vector{FreeModElem}())
  N_quo = isdefined(N, :quo) ? N.quo : Oscar.SubModuleOfFreeModule(ambient_free_module(N), Vector{FreeModElem}())
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

    return SQ,M_hom,N_hom
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

function in(v::SubQuoElem, M::ModuleFP)
  return parent(v) === M
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
    std_basis(F::ModuleGens, reduced::Bool=false)

Return a standard basis of `F` as an object of type `ModuleGens.
If `reduced` is set to `true` and the ordering of the underlying ring is global, 
a reduced Gröbner basis is computed.
"""
function std_basis(F::ModuleGens, reduced::Bool=false)
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
    leading_module(F::ModuleGens)

Return the leading module basis of `F` as an object of type `ModuleGens.
The leading module is with respect to the ordering defined on the Singular side.
In particular, something on the Singular side must be defined.
"""
function leading_module(F::ModuleGens)
  #if !isdefined(F, :SF) && !isdefined(F, :S)
  #  error("No singular ring and therefore no ordering defined.")
  #end
  singular_assure(F)
  singular_gens = singular_generators(F)
  if singular_gens.isGB
    return ModuleGens(F.F, Singular.lead(singular_gens))
  end
  return ModuleGens(F.F, Singular.lead(Singular.std(singular_gens)))
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
    (M::SubQuo{T})(f::FreeModElem{T}; check::Bool = true) where T

Given an element `f` of the ambient free module of `M` which represents an element of `M`, 
return the represented element. If `check == true` (default) it is checked that `f` represents
indeed an element of `M`.
"""
function (M::SubQuo{T})(f::FreeModElem{T}; check::Bool = true) where T
  if check
    singular_assure(M.sum.gens)
    b = M.sum.gens.SF(f)
    c = _reduce(b, singular_generators(std_basis(M.sum)))
    iszero(c) || error("not in the module")
  end
  return SubQuoElem(f, M)
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
    error("elements not compatible")
  end
  return SubQuoElem(a*coeffs(b), b.parent)
end
function *(a::MPolyElem, b::SubQuoElem) 
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
  end
  return SubQuoElem(a*coeffs(b), b.parent)
end
function *(a::RingElem, b::SubQuoElem) 
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
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
    sub(F::FreeMod, V::Vector{<:FreeModElem}, task::Symbol = :none)

Given a vector `V` of elements of `F`, return the submodule of `F` generated by these elements.

Put more precisely, if `N` denotes this submodule, return `N` as an object of type `SubQuo`. Additionally,

- return the inclusion map `N` $\to$ `F` if `task = :with_morphism`,
- return and cache the inclusion map `N` $\to$ `F` if `task = :cache_morphism`,
- do none of the above if `task = :none` (default).

If `task = :only_morphism`, return only the inclusion map.
"""
function sub(F::FreeMod, V::Vector{<:FreeModElem}, task::Symbol = :none)
  s = SubQuo(F, V)
  emb = hom(s, F, V)
  set_attribute!(s, :canonical_inclusion => emb)
  (task == :none || task == :module) && return s
  task == :cache_morphism && register_morphism!(emb)
  task == :only_morphism && return emb
  (task == :cache_morphism || task == :both || task == :with_morphism) && return s, emb
  error("No valid option for task.")
end

@doc Markdown.doc"""
    sub(F::FreeMod, O::Vector{<:SubQuoElem}, task::Symbol = :none)

Return `S` as a submodule of `F`, where `S` is generated by `O`.
The embedding module of the parent of the elements of `O` must be `F`.
If `task` is set to `:none` (default option) or to `:module` return only `S`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical injection morphism
$S \to F$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function sub(F::FreeMod, O::Vector{<:SubQuoElem}, task::Symbol = :none)
  s = SubQuo(F, [x.repres for x = O])
  return sub(F, s, task)
end

@doc Markdown.doc"""
    sub(F::FreeMod, s::SubQuo, task::Symbol = :none)

Return `s` as a submodule of `F`, that is the embedding free module of `s` must 
be `F` and `s` has no relations.
If `task` is set to `:none` (default option) or to `:module` return only `s`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical injection morphism
$s \to F$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function sub(F::FreeMod, s::SubQuo, task::Symbol = :none)
  @assert !isdefined(s, :quo)
  @assert s.F === F
  emb = hom(s, F, Vector{ModuleFPElem}([repres(x) for x in gens(s)]))
  #emb = hom(s, F, [FreeModElem(x.repres.coords, F) for x in gens(s)])
  set_attribute!(s, :canonical_inclusion => emb)
  (task == :none || task == :module) && return s
  task == :cache_morphism && register_morphism!(emb)
  task == :only_morphism && return emb 
  (task == :cache_morphism || task == :both || task == :with_morphism) && return s, emb
  error("No valid option for task.")
end

@doc Markdown.doc"""
    sub(M::SubQuo, V::Vector{<:SubQuoElem}, task::Symbol = :none, check = true)

Given a vector `V` of elements of `M`, return the submodule of `M` generated by these elements.

Put more precisely, if `N` denotes this submodule, return `N` as an object of type `SubQuo`. Additionally,

- return the inclusion map `N` $\to$ `M` if `task = :with_morphism`,
- return and cache the inclusion map `N` $\to$ `M` if `task = :cache_morphism`,
- do none of the above if `task = :none` (default).

If `task = :only_morphism`, return only the inclusion map.

If `check = false`, the condition on `V` is not checked by the function.
"""
function sub(M::SubQuo, V::Vector{<:SubQuoElem}, task::Symbol = :none, check = true)
  if check
    @assert all(x -> x.parent === M, V)
  end
  t = SubQuo(M.F, [x.repres for x in V])
  if isdefined(M, :quo)
    t.quo = M.quo
    t.sum = sum(t.sub, t.quo)
  end
  emb = hom(t, M, V)
  set_attribute!(t, :canonical_inclusion => emb)
  (task == :none || task == :module) && return t
  task == :cache_morphism && register_morphism!(emb)
  task == :only_morphism && return emb 
  (task == :cache_morphism || task == :both || task == :with_morphism) && return t, emb
  error("No valid option for task.")
end

@doc Markdown.doc"""
    quo(F::FreeMod, V::Vector{<:FreeModElem}, task::Symbol = :none)

Given a vector `V` of elements of `F`, return the quotient of `F` by the submodule of `F` which is generated by these elements.

Put more precisely, if `N` denotes this quotient, return `N` as an object of type `SubQuo`. Additionally,

- return the projection map `F` $\to$ `N` if `task = :with_morphism`,
- return and cache the projection map `F` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` (default).

If `task = :only_morphism`, return only the projection map.
"""
function quo(F::FreeMod, V::Vector{<:FreeModElem}, task::Symbol = :none)
  S = SubQuo(F, basis(F))
  Q = SubQuo(S, V)

  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
    quo(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :none) where T

Compute $F / T$, where $T$ is generated by $O$.
The embedding free module of the parent of the elements of `O` must be `F`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function quo(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :none) where T
  S = SubQuo(F, basis(F))
  Q = SubQuo{T}(S, [x.repres for x = O])
  
  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
    quo(F::SubQuo, O::Vector{<:FreeModElem}, task::Symbol = :none)

Compute $F / T$, where $T$ is generated by $O$.
The elements of `O` must be elements of the embedding free module of `S`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function quo(F::SubQuo, O::Vector{<:FreeModElem}, task::Symbol = :none)
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
    quo(M::SubQuo, V::Vector{<:SubQuoElem}, task::Symbol = :none)

Given a vector `V` of elements of `M`, return the quotient of `M` by the submodule of `M` which is generated by these elements.

Put more precisely, if `N` denotes this quotient, return `N` as an object of type `SubQuo`. Additionally,

- return the projection map `M` $\to$ `N` if `task = :with_morphism`,
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` (default).

If `task = :only_morphism`, return only the projection map.
"""
function quo(M::SubQuo, V::Vector{<:SubQuoElem}, task::Symbol = :none)
  return quo(M, [x.repres for x = V], task)
end

@doc Markdown.doc"""
    quo(M::SubQuo, U::SubQuo, task::Symbol = :none)

Put more precisely, if `N` denotes this quotient, return `N` as an object of type `SubQuo`. Additionally,

- return the projection map `M` $\to$ `N` if `task = :with_morphism`,
- return and cache the projection map `M` $\to$ `N` if `task = :cache_morphism`,
- do none of the above if `task = :none` (default).

If `task = :only_morphism`, return only the projection map.
"""
function quo(M::SubQuo, U::SubQuo, task::Symbol = :none)
#  @assert !isdefined(U, :quo)
  # TODO @assert M.quo == U.quo or too expensive?
  Q = SubQuo(M, oscar_generators(U.sum.gens))
  return return_quo_wrt_task(M, Q, task)
end

@doc Markdown.doc"""
    quo(F::FreeMod{R}, T::SubQuo{R}, task::Symbol = :none) where R

Compute $F / T$.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$F \to F/T$.
If `task` is set to `:cache_morphism` the morphism is also cached.
If `task` is set to `:only_morphism` return only the morphism.
"""
function quo(F::FreeMod{R}, T::SubQuo{R}, task::Symbol = :none) where R
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
function syzygy_module(F::ModuleGens; sub = FreeMod(base_ring(F.F), length(oscar_generators(F))))
  F[Val(:S), 1] #to force the existence of F.S
  s = Singular.syz(singular_generators(F)) # TODO syz is sometimes too slow, example [8*x^2*y^2*z^2 + 13*x*y*z^2 12*x^2 + 7*y^2*z; 13*x*y^2 + 12*y*z^2 4*x^2*y^2*z + 8*x*y*z; 9*x*y^2 + 4*z 12*x^2*y*z^2 + 9*x*y^2*z]
  return SubQuo(sub, s)
end

@doc Markdown.doc"""
    gens(M::SubQuo{T}) where T

Return the generators of `M`.
"""
function gens(M::SubQuo{T}) where T
  return [gen(M,i) for i=1:ngens(M)]
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
zero(M::SubQuo) = SubQuoElem(zero(M.F), M)

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

@doc Markdown.doc"""
    presentation(M::SubQuo)

Return a free presentation of `M`. 

# Examples
"""
function presentation(SQ::SubQuo)
  #A+B/B is generated by A and B
  #the relations are A meet B? written wrt to A
  R = base_ring(SQ)
  F = FreeMod(R, ngens(SQ.sub))
  q = elem_type(F)[]
  if isgenerated_by_standard_unit_vectors(SQ.sub)
    if isdefined(SQ, :quo)
      q = [FreeModElem(coords(g), F) for g in gens(SQ.quo)]
    end
  else
    s = syzygy_module(SQ.sum.gens)
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
  h_G_F = hom(G, F, q)
  h_F_SQ = hom(F, SQ, gens(SQ)) # DO NOT CHANGE THIS LINE, see present_as_cokernel and preimage

  Z = FreeMod(F.R, 0)
  set_attribute!(Z, :name => "Zero")
  h_SQ_Z = hom(SQ, Z, Vector{ModuleFPElem}([zero(Z) for i=1:ngens(SQ)]))
  return Hecke.ChainComplex(Oscar.ModuleFP, Oscar.ModuleMap[h_G_F, h_F_SQ, h_SQ_Z], check = false)
end

@doc Markdown.doc"""
    presentation(F::FreeMod)

Return a free presentation of $F$.

# Examples
"""
function presentation(F::FreeMod)
  Z = FreeMod(F.R, 0)
  set_attribute!(Z, :name => "Zero")
  return Hecke.ChainComplex(ModuleFP, ModuleMap[hom(Z, F, FreeModElem[]), hom(F, F, gens(F)), hom(F, Z, Vector{ModuleFPElem}([zero(Z) for i=1:ngens(F)]))], check = false)
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
  R_b = obj(chainComplex, 1)
  f = map(chainComplex, 2)
  g = map(chainComplex, 1)
  presentation_module = quo(R_b, image(f)[1])

  if task == :none
    return presentation_module
  end
  
  # The isomorphism is just the identity matrix
  isomorphism = hom(presentation_module, SQ, Vector{ModuleFPElem}([g(x) for x in gens(R_b)]))
  inverse_isomorphism = hom(SQ, presentation_module, Vector{ModuleFPElem}([presentation_module[i] for i=1:ngens(SQ)]))
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

  M_to_N = hom(M, N, Vector{SubQuoElem}([SubQuoElem(coordinates(m.repres, N), N) for m in gens(M)]))

  if task == :cache_morphism || task == :inverse
    N_to_M = hom(N, M, Vector{SubQuoElem}([SubQuoElem(coordinates(n.repres, M), M) for n in gens(N)]))
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
    SubQuoHom(D::SubQuo, C::ModuleFP, im::Vector)

Return the morphism $D \to C$ for a subquotient $D$ where `D[i]` is mapped to `im[i]`.
In particular, `length(im) == ngens(D)` must hold.
"""
SubQuoHom(D::SubQuo, C::ModuleFP, im::Vector) = SubQuoHom{typeof(D), typeof(C)}(D, C, im)

@doc Markdown.doc"""
    SubQuoHom(D::SubQuo, C::ModuleFP, mat::MatElem)

Return the morphism $D \to C$ corresponding to the given matrix, where $D$ is a subquotient.
`mat` must have `ngens(D)` many rows and `ngens(C)` many columns.
"""
function SubQuoHom(D::SubQuo, C::ModuleFP, mat::MatElem)
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

    iswelldefined(a::ModuleMap)

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

julia> iswelldefined(a)
true

julia> W = [y*N[1], x*N[2]]
2-element Vector{SubQuoElem{fmpq_mpoly}}:
 x*y*e[1]
 x*y*e[1]

julia> b = hom(M, N, W)
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

julia> iswelldefined(b)
false
```
"""
hom(M::SubQuo{T}, N::ModuleFP{T}, V::Vector{<:ModuleFPElem}) where T = SubQuoHom(M, N, V) 
hom(M::SubQuo{T}, N::ModuleFP{T},  A::MatElem{T}) where T = SubQuoHom(M, N, A)
function iswelldefined(H::ModuleMap)
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

function show_morphism(f::ModuleMap)
  display(matrix(f))
end

@doc Markdown.doc"""
    hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{ <: ModuleMap})

Given modules `M`, `N` which are tensor products with the same number of factors,
say $M = M_1 \otimes \cdots \otimes M_r$, $N = N_1 \otimes \cdots \otimes N_r$,
and given a vector `V` of homomorphisms $a_i : M_i \to N_i$, return 
$a_1 \otimes \cdots \otimes a_r$.
"""
function hom_tensor(M::ModuleFP, N::ModuleFP, V::Vector{ <: ModuleMap})
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
  return hom(M, N, map(map_gen, gens(M)))
end

@doc Markdown.doc"""
    hom_prod_prod(M::ModuleFP, N::ModuleFP, A::Matrix{<:ModuleMap})

Given modules `M`, `N` which are products with the same number of factors,  
say $M = \prod_{i=1}^r M_i$, $N = \prod_{j=1}^r N_j$, and given a matrix 
`A` of homomorphisms $a_{ij} : M_i \to N_j$, return the homomorphism
$M \rightarrow N$ with $ij$-components $a_{ij}$.
"""
function hom_prod_prod(M::ModuleFP, N::ModuleFP, A::Matrix{<:ModuleMap})
  tM = get_attribute(M, :direct_product)
  tM === nothing && error("both modules must be direct products")
  tN = get_attribute(N, :direct_product)
  tN === nothing && error("both modules must be direct products")
  @assert length(tM) == size(A, 1) && length(tN) == size(A, 2)
  @assert all(ij -> domain(A[ij[1],ij[2]]) === tM[ij[1]] && codomain(A[ij[1],ij[2]]) === tN[ij[2]], Base.Iterators.ProductIterator((1:size(A, 1), 1:size(A, 2))))
  #need the canonical maps..., maybe store them as well?
  return hom(M,N,Vector{ModuleFPElem}([sum([Hecke.canonical_injection(N,j)(sum([A[i,j](Hecke.canonical_projection(M,i)(g)) for i=1:length(tM)])) for j=1:length(tN)]) for g in gens(M)]))
end
# hom(prod -> X), hom(x -> prod)
# if too much time: improve the hom(A, B) in case of A and/or B are products - or maybe not...
# tensor and hom functors for chain complex
# dual: ambig: hom(M, R) or hom(M, Q(R))?

@doc Markdown.doc"""
    coordinates(a::FreeModElem, SQ::SubQuo)::Union{Nothing,Oscar.SRow}

Compute a sparse row `r` such that `a` is a representative of `SubQuoElem(r, SQ)`.
If no such `r` exists, `nothing` is returned.
"""
function coordinates(a::FreeModElem, SQ::SubQuo)::Union{Nothing,Oscar.SRow}
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  oscar_assure(SQ.sub.gens)
  if isdefined(SQ, :quo)
    oscar_assure(SQ.quo.gens)
    generators = sum(SQ.sub, SQ.quo)
  else
    generators = SQ.sub
  end
  S = singular_generators(generators.gens)
  singular_assure(SQ.sum.gens)
  b = ModuleGens([a], SQ.sum.gens.SF)
  singular_assure(b)
  s, r = Singular.lift(S, singular_generators(b))
  if Singular.ngens(s) == 0 || iszero(s[1])
    return nothing
  end
  Rx = base_ring(SQ)
  return sparse_row(Rx, s[1], 1:ngens(SQ))
end

@doc Markdown.doc"""
    represents_element(a::FreeModElem, SQ::SubQuo)

Check if `a` represents an element `SQ`.
"""
function represents_element(a::FreeModElem, SQ::SubQuo)
  return !isnothing(coordinates(a,SQ))
end

@doc Markdown.doc"""
    reduce(M::ModuleGens, GB::ModuleGB)

Reduce `M` with respect to the Gröbner basis `GB`.
"""
function reduce(M::ModuleGens, GB::ModuleGB)
  @assert M.F === GB.groebner_basis.F

  P = isdefined(GB, :quo_groebner_basis) ? union(GB.groebner_basis, GB.quo_groebner_basis.groebner_basis) : GB.groebner_basis

  singular_assure(P)
  singular_assure(M)

  red = _reduce(M.S, P.S)
  res = ModuleGens(M.F, red)
  oscar_assure(res)
  return res
end

@doc Markdown.doc"""
    reduce(M::SubModuleOfFreeModule, GB::ModuleGB)

Reduce `M` with respect to the Gröbner basis `GB`.
"""
function reduce(M::SubModuleOfFreeModule, GB::ModuleGB)
  return SubModuleOfFreeModule(M.F, reduce(M.gens, GB), default_ordering(M))
end

@doc Markdown.doc"""
    reduce(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Reduce `M` with respect to `N`, that is with respect to a Gröbner basis of `N` (the Gröbner basis is computed for this).
"""
function reduce(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  return reduce(M, groebner_basis(N, default_ordering(M)))
end

@doc Markdown.doc"""
    reduce(v::FreeModElem, GB::ModuleGB)

Reduce the element `v` with respect to the Gröbner basis `GB` 
"""
function reduce(v::AbstractFreeModElem, GB::ModuleGB)
  return reduce(ModuleGens(parent(v), [v]), GB).O[1]
end

@doc Markdown.doc"""
    reduce(v::FreeModElem, N::SubModuleOfFreeModule)

Reduce the element `v` with respect to a Gröbner basis of `N`.
"""
function reduce(v::AbstractFreeModElem, N::SubModuleOfFreeModule)
  return reduce(v, groebner_basis(N))
end

hom(D::SubQuo, C::ModuleFP, A::Vector) = SubQuoHom(D, C, A)

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
    return iszero(m.repres)
  end
  singular_assure(C.quo.gens)
  x = _reduce(C.quo.gens.SF(m.repres), singular_generators(std_basis(C.quo)))
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
    return hom(F, G, Vector{FreeModElem}([FreeModElem(x.coords[R, (i-1)*m+1:i*m], G) for i=1:n]))
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

Additionally, if `K` denotes this kernel, return the inclusion map `K` $\rightarrow$ `domain(a)`.

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
    s = sub(G, gens(G))
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
  b = ModuleGens(g)
  k = syzygy_module(b)
  if isa(codomain(h), SubQuo)
    s = collect(k.sub.gens)
    k = sub(G, [FreeModElem(x.coords[R,1:dim(G)], G) for x = s])
  else
    #the syzygie_module creates a new free module to work in
    k = sub(G, [FreeModElem(x.coords, G) for x = collect(k.sub.gens)])
  end
  @assert k.F === G
  c = collect(k.sub.gens)
  return k, hom(k, parent(c[1]), c)
end

@doc Markdown.doc"""
    image(a::FreeModuleHom)

Return the image of `a` as an object of type `SubQuo`.

Additionally, if `I` denotes this image, return the inclusion map `I` $\rightarrow$ `codomain(a)`.

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
  s = sub(codomain(h), si)
  return s, hom(s, codomain(h), si)
end

@doc Markdown.doc"""
    image(a::SubQuoHom)

Return the image of `a` as an object of type `SubQuo`.

Additionally, if `I` denotes this image, return the inclusion map `I` $\rightarrow$ `codomain(a)`.
"""
function image(h::SubQuoHom)
  s = sub(codomain(h), h.im)
  return s, hom(s, codomain(h), h.im)
end

@doc Markdown.doc"""
    kernel(a::SubQuoHom)

Return the kernel of `a` as an object of type `SubQuo`.

Additionally, if `K` denotes this kernel, return the inclusion map `K` $\rightarrow$ `domain(a)`.
"""
function kernel(h::SubQuoHom)
  D = domain(h)
  R = base_ring(D)
  F = FreeMod(R, ngens(D))
  hh = hom(F, codomain(h), map(h, gens(D)))
  k = kernel(hh)
  @assert domain(k[2]) === k[1]
  @assert codomain(k[2]) === F
  hh = hom(F, D, gens(D))
  im::Vector{SubQuoElem} = map(x->hh(k[2](x)), gens(k[1]))
  k = sub(D, im)
  return k, hom(k, D, im)
end

@doc Markdown.doc"""
    free_resolution(F::FreeMod)

Return a free resolution of `F`.

# Examples
"""
function free_resolution(F::FreeMod)
  return presentation(F)
end

@doc Markdown.doc"""
    free_resolution(M::SubQuo, limit::Int = -1)

Return a free resolution of `M`. 

If `limit != -1`, the free resolution
is only computed up to the `limit`-th free module.

# Examples
"""
function free_resolution(M::SubQuo, limit::Int = -1)
  p = presentation(M)
  mp = [map(p, j) for j=1:length(p)]
  while true
    k, mk = kernel(mp[1])
    nz = findall(x->!iszero(x), gens(k))
    if length(nz) == 0 
      Z = FreeMod(base_ring(M), 0)
      set_attribute!(Z, :name => "Zero")
      h = hom(Z, domain(mp[1]), FreeModElem[])
      insert!(mp, 1, h)
      break
    elseif limit != -1 && length(mp) > limit
      break
    end
    F = FreeMod(base_ring(M), length(nz))
    g = hom(F, codomain(mk), collect(k.sub.gens)[nz])
    insert!(mp, 1, g)
  end
  return Hecke.ChainComplex(ModuleFP, mp, check = false, direction = :right)
end

function Hecke.ring(I::MPolyIdeal)
  return parent(gen(I, 1))
end

@doc Markdown.doc"""
    free_resolution(I::MPolyIdeal)

Compute a free resolution of `I`.

# Examples
"""
function free_resolution(I::MPolyIdeal)
  F = free_module(Hecke.ring(I), 1)
  S = sub(F, [x * gen(F, 1) for x = gens(I)])
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
  q = quo(F, [x * gen(F, 1) for x = gens(Q.I)])
  n = Hecke.find_name(Q)
  if n !== nothing
    AbstractAlgebra.set_name!(q, String(n))
  end
  return free_resolution(q)
end

@doc Markdown.doc"""
    iszero(f::ModuleMap)

Return true iff `f` is the zero map.
"""
function iszero(f::ModuleMap)
  return all(iszero, map(f, gens(domain(f))))
end

@doc Markdown.doc"""
    hom(M::ModuleFP, N::ModuleFP, alg::Symbol=:maps)

Return a subquotient $S$ such that $\text{Hom}(M,N) \cong S$ along with a function 
that converts elements from $S$ into morphisms $M \to N$.
If `alg` is `:matrices` a different implementation that is using matrices instead of maps is used.
"""
function hom(M::ModuleFP, N::ModuleFP, alg::Symbol=:maps)  
  d1 = 1
  d2 = 2
  if alg == :matrices && M isa SubQuo && N isa SubQuo
    return hom_matrices(M,N,false)
  end
  p1 = presentation(M)
  p2 = presentation(N)
  k, mk = kernel(map(p2, d2))
  #Janko: have R^t1 -- g1 = map(p2, 0) -> R^t0 -> G
  #kernel g1: k -> R^t1
  #source: Janko's CA script: https://www.mathematik.uni-kl.de/~boehm/lehre/17_CA/ca.pdf
  F = FreeMod(base_ring(M), ngens(k))
  g2 = hom(F, codomain(mk), collect(k.sub.gens)) #not clean - but maps not (yet) working
  #step 2
  H_s0_t0, mH_s0_t0 = hom(domain(map(p1, d1)), domain(map(p2, d1)))
  H_s1_t1, mH_s1_t1 = hom(domain(map(p1, d2)), domain(map(p2, d2)))
  D, pro, emb = direct_product(H_s0_t0, H_s1_t1, task = :both)

  H_s1_t0, mH_s1_t0 = hom(domain(map(p1, d2)), domain(map(p2, d1)))

  delta = hom(D, H_s1_t0, Vector{ModuleFPElem}([preimage(mH_s1_t0, map(p1, d2)*mH_s0_t0(pro[1](g))-mH_s1_t1(pro[2](g))*map(p2, d2)) for g = gens(D)]))

  H_s0_t1, mH_s0_t1 = hom(domain(map(p1, d1)), domain(map(p2, d2)))
  H_s1_t2, mH_s1_t2 = hom(domain(map(p1, d2)), F)

  E, pr = direct_product(H_s0_t1, H_s1_t2, task = :prod)

  rho = hom(E, D, Vector{ModuleFPElem}([emb[1](preimage(mH_s0_t0, mH_s0_t1(pr[1](g))*map(p2, d2))) + 
                  emb[2](preimage(mH_s1_t1, map(p1, d2)*mH_s0_t1(pr[1](g)) - mH_s1_t2(pr[2](g))*g2)) for g = gens(E)]))
  #need quo(kern(delta), image(rho))
 
  kDelta = kernel(delta)

  psi = kDelta[2]*pro[1]
  #psi = hom(kDelta[1], H_s0_t0, [psi(g) for g = gens(kDelta[1])])

  H = quo(sub(D, kDelta[1]), image(rho)[1])

  #x in ker delta: mH_s0_t0(pro[1](x)) should be a hom from M to N
  function im(x::SubQuoElem)
    @assert parent(x) === H
    #return SubQuoHom(M, N, mH_s0_t0(pro[1](x.repres)).matrix)
    return hom(M, N, Vector{ModuleFPElem}([map(p2, d1)(mH_s0_t0(pro[1](x.repres))(preimage(map(p1, d1), g))) for g = gens(M)]))
  end

  function pre(f::ModuleMap)
    @assert domain(f) === M
    @assert codomain(f) === N
    Rs0 = domain(map(p1, d1))
    Rt0 = domain(map(p2, d1))
    g = hom(Rs0, Rt0, Vector{ModuleFPElem}([preimage(map(p2, d1), f(map(p1, d1)(g))) for g = gens(Rs0)]))

    #return H(preimage(psi, (preimage(mH_s0_t0, g))).repres)
    return SubQuoElem(repres(preimage(psi, (preimage(mH_s0_t0, g)))), H)
    #return SubQuoElem(preimage(kDelta[2], emb[1](preimage(mH_s0_t0, g))).repres, H)
    #return SubQuoElem(emb[1](preimage(mH_s0_t0, g)), H) #???
  end
  to_hom_map = MapFromFunc(im, pre, H, Hecke.MapParent(M, N, "homomorphisms"))
  set_attribute!(H, :show => Hecke.show_hom, :hom => (M, N), :module_to_hom_map => to_hom_map)
  return H, to_hom_map
end

@doc Markdown.doc"""
    homomorphism(f::Union{SubQuoElem,FreeModElem})

If `f` is an element in a module created via `hom(M,N)` for some `M` and `N`, 
return the morphism $\phi : M \to N$ that corresponds to `f`.
"""
function homomorphism(f::Union{AbstractSubQuoElem,AbstractFreeModElem})
  H = f.parent
  to_hom_map = get_attribute(H, :module_to_hom_map)
  to_hom_map === nothing && error("element doesn't live in a hom module")  
  return to_hom_map(f)
end

@doc Markdown.doc"""
    module_elem(H::ModuleFP, phi::ModuleMap)

Let `H` be created via `hom(M,N)` for some `M` and `N`. Return 
the element in `H` corresponding to `phi`.
"""
function module_elem(H::ModuleFP, phi::ModuleMap)
  to_hom_map = get_attribute(H, :module_to_hom_map)
  to_hom_map === nothing && error("module must be a hom module")
  map_to_hom = to_hom_map.g
  return map_to_hom(phi)
end

#TODO
#  replace the +/- for the homs by proper constructors for homs and direct sums
#  relshp to store the maps elsewhere

@doc Markdown.doc"""
    *(h::ModuleMap, g::ModuleMap)

Return the composition $g \circ h$.
"""
function *(h::ModuleMap, g::ModuleMap)
  @assert codomain(h) === domain(g)
  return hom(domain(h), codomain(g), Vector{ModuleFPElem}([g(h(x)) for x = gens(domain(h))]))
end

-(h::ModuleMap) = hom(domain(h), codomain(h), [-h(x) for x in gens(domain(h))])
function -(h::ModuleMap, g::ModuleMap)
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  return hom(domain(h), codomain(h), Vector{ModuleFPElem}([h(x) - g(x) for x in gens(domain(h))]))
end
function +(h::ModuleMap, g::ModuleMap)
  @assert domain(h) === domain(g)
  @assert codomain(h) === codomain(g)
  return hom(domain(h), codomain(h), Vector{ModuleFPElem}([h(x) + g(x) for x in gens(domain(h))]))
end


@doc Markdown.doc"""
    restrict_codomain(H::ModuleMap, M::SubQuo)

Return, if possible, a homomorphism, which is mathematically identical to `H`,
but has codomain `M`. `M` has to be a submodule of the codomain of `H`.
"""
function restrict_codomain(H::ModuleMap, M::SubQuo)
  @assert codomain(H) isa SubQuo
  return hom(domain(H), M, map(v -> SubQuoElem(v, M), map(x -> H(x).repres, gens(domain(H)))))
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
  i = sub(domain(H), map(m -> SubQuoElem(repres(m), domain(H)), gens(M)), :cache_morphism, false)[2]
  return i*H
end

@doc Markdown.doc"""
    inv(H::ModuleMap)

Compute $H^{-1}$. `H` must be bijective.
"""
function inv(H::ModuleMap)
  if isdefined(H, :inverse_isomorphism)
    return H.inverse_isomorphism
  end
  @assert isbijective(H)
  N = domain(H)
  M = codomain(H)

  Hinv = hom(M,N, Vector{ModuleFPElem}([preimage(H,m) for m in gens(M)]))
  Hinv.inverse_isomorphism = H
  H.inverse_isomorphism = Hinv

  return Hinv
end

##################################################
# direct sum
##################################################
@doc Markdown.doc"""
    direct_sum(F::FreeMod{T}...; task::Symbol = :sum) where T

Given free modules $F_1\dots F_n$, say, return the direct sum $\bigoplus_{i=1}^n F_i$.

Additionally, return 
- a vector containing the canonical injections  $F_i\rightarrow\bigoplus_{i=1}^n F_i$ if `task = :sum` (default),
- a vector containing the canonical projections  $\bigoplus_{i=1}^n F_i\rightarrow F_i$ if `task = :prod`,
- two vectors containing the canonical injections and projections, respectively, if `task = :both`,
- none of the above if `task = :none`.
"""
function direct_sum(F::FreeMod{T}...; task::Symbol = :sum) where {T}
  return direct_product(F...; task = task)
end

##################################################
# direct product
##################################################
@doc Markdown.doc"""
    direct_product(F::FreeMod{T}...; task::Symbol = :prod) where T

Given free modules $F_1\dots F_n$, say, return the direct product $\prod_{i=1}^n F_i$.

Additionally, return a vector containing
- a vector containing the canonical projections  $\prod_{i=1}^n F_i\rightarrow F_i$ if `task = :prod` (default),
- a vector containing the canonical injections  $F_i\rightarrow\prod_{i=1}^n F_i$ if `task = :sum`,
- two vectors containing the canonical projections and injections, respectively, if `task = :both`,
- none of the above if `task = :none`.
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
  projection_dictionary = IdDict{Int,ModuleMap}()
  injection_dictionary = IdDict{Int,ModuleMap}()
  set_attribute!(G, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  i = 0
  for f = F
    if task in [:sum, :both]
      push!(emb, hom(f, G, Vector{ModuleFPElem}([gen(G, j+i) for j=1:ngens(f)])))
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
    direct_product(M::ModuleFP{T}...; task::Symbol = :none) where T

Given modules $M_1\dots M_n$, say, return the direct product $\prod_{i=1}^n M_i$.

Additionally, return a vector containing
- a vector containing the canonical projections  $\prod_{i=1}^n M_i\rightarrow M_i$ if `task = :prod` (default),
- a vector containing the canonical injections  $M_i\rightarrow\prod_{i=1}^n M_i$ if `task = :sum`,
- two vectors containing the canonical projections and injections, respectively, if `task = :both`,
- none of the above if `task = :none`.
"""
function direct_product(M::ModuleFP{T}...; task::Symbol = :none) where T
  F, pro, mF = direct_product([ambient_free_module(x) for x = M]..., task = :both)
  s, emb_sF = sub(F, vcat([[mF[i](y) for y = ambient_representatives_generators(M[i])] for i=1:length(M)]...), :both)
  q::Vector{FreeModElem} = vcat([[mF[i](y) for y = rels(M[i])] for i=1:length(M)]...)
  pro_quo = nothing
  if length(q) != 0
    s, pro_quo = quo(s, q, :both)
  end
  set_attribute!(s, :show => Hecke.show_direct_product, :direct_product => M)
  projection_dictionary = IdDict{Int,ModuleMap}()
  injection_dictionary = IdDict{Int,ModuleMap}()
  set_attribute!(s, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  if task == :none
    return s
  end
  if task == :prod || task != :sum
    if pro_quo === nothing
      for i=1:length(pro)
        pro[i] = hom(s, M[i], Vector{ModuleFPElem}([M[i](pro[i](emb_sF(gen))) for gen in gens(s)])) # TODO distinction between pro on the left and pro on the right side!
        projection_dictionary[i] = pro[i]
      end
    else
      for i=1:length(pro)
        pro[i] = hom(s, M[i], Vector{ModuleFPElem}([M[i](pro[i](emb_sF(preimage(pro_quo,gen)))) for gen in gens(s)]))
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
        mF[i] = hom(M[i], s, Vector{ModuleFPElem}([preimage(emb_sF, mF[i](repres(g))) for g in gens(M[i])]))
        injection_dictionary[i] = mF[i]
      end
    else
      for i=1:length(mF)
        mF[i] = hom(M[i], s, Vector{ModuleFPElem}([pro_quo(preimage(emb_sF, mF[i](repres(g)))) for g in gens(M[i])]))
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
    direct_sum(M::ModuleFP{T}...; task::Symbol = :none) where T

Given modules $M_1\dots M_n$, say, return the direct sum $\bigoplus_{i=1}^n M_i$.   
Additionally, return 
- a vector containing the canonical injections  $M_i\rightarrow\bigoplus_{i=1}^n M_i$ if `task = :sum` (default),
- a vector containing the canonical projections  $\bigoplus_{i=1}^n M_i\rightarrow M_i$ if `task = :prod`,
- two vectors containing the canonical injections and projections, respectively, if `task = :both`,
- none of the above if `task = :none`.
"""
function direct_sum(M::ModuleFP{T}...; task::Symbol = :none) where {T}
  return direct_product(M...; task)
end

⊕(M::ModuleFP...) = direct_product(M..., task = :none)


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
  emb = hom(H[i], G, Vector{ModuleFPElem}([G[l+j] for l = 1:ngens(H[i])]))
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
  pro = hom(G, H[i], Vector{ModuleFPElem}(vcat([zero(H[i]) for l=1:j], gens(H[i]), [zero(H[i]) for l=1+j+ngens(H[i]):ngens(G)])))
  projection_dictionary[i] = pro
  return pro
end
    
##################################################
# Tensor
##################################################

@doc Markdown.doc"""
    tensor_product(G::FreeMod...; task::Symbol = :none)

Given free modules $G_i$ compute the tensor product $G_1\otimes \cdots \otimes G_n$.
If `task` is set to ":map", a map $\phi$ is returned that
maps tuples in $G_1 \times \cdots \times G_n$ to pure tensors
$g_1 \otimes \cdots \otimes g_n$. The map admits a preimage as well.
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
    tensor_product(G::ModuleFP...; task::Symbol = :none)

Given modules $G_i$ compute the tensor product $G_1\otimes \cdots \otimes G_n$.
If `task` is set to ":map", a map $\phi$ is returned that
maps tuples in $G_1 \times \cdots \times G_n$ to pure tensors
$g_1 \otimes \cdots \otimes g_n$. The map admits a preimage as well.
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

  generating_tensors = map(mF, map(tuple -> map(x -> parent(x) isa FreeMod ? x : x.repres, tuple), corresponding_tuples))
  s, emb = sub(F, generating_tensors, :with_morphism)
  #s, emb = sub(F, vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(gens(x, ambient_free_module(x)) for x = G))]), :with_morphism)
  q = vcat([vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(i == j ? rels(G[i]) : gens(ambient_free_module(G[i])) for i=1:length(G)))]) for j=1:length(G)]...) 
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
  tensor_chain = Hecke.map_type(C)[]
  tensor_modules = [tensor_product(P, domain(C.maps[1]), task=:cache_morphism)[1]]
  tensor_modules = vcat(tensor_modules, [tensor_product(P, codomain(f), task=:cache_morphism)[1] for f in C.maps])

  for i=1:length(C)
    A = tensor_modules[i]
    B = tensor_modules[i+1]

    push!(tensor_chain, hom_tensor(A,B,[identity_map(P), map(C,i)]))
  end

  return Hecke.ChainComplex(ModuleFP, tensor_chain)
end

@doc Markdown.doc"""
    tensor_product(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)

Apply $\bullet \otimes P$ to `C`.
"""
function tensor_product(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
  tensor_chain = Hecke.map_type(C)[]
  tensor_modules = [tensor_product(domain(C.maps[1]), P, task=:cache_morphism)[1]]
  tensor_modules = vcat(tensor_modules, [tensor_product(codomain(f), P, task=:cache_morphism)[1] for f in C.maps])

  for i=1:length(C)
    A = tensor_modules[i]
    B = tensor_modules[i+1]

    push!(tensor_chain, hom_tensor(A,B,[map(C,i), identity_map(P)]))
  end

  return Hecke.ChainComplex(ModuleFP, tensor_chain)
end

@doc Markdown.doc"""
    tor(M::ModuleFP, N::ModuleFP, i::Int)

Compute $\text{Tor}_i(M,N)$.
"""
function tor(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M)[1:end-2]
  lifted_resolution = tensor_product(free_res, N) #TODO only three homs are necessary
  return homology(lifted_resolution,length(lifted_resolution)-i)
end

#TODO, mF
#  (hom lift) => hom and tensor functor
#  filtrations
#  more constructors
#################################################
#
#################################################
@doc Markdown.doc"""
    lift_homomorphism_contravariant(Hom_MP::ModuleFP, Hom_NP::ModuleFP, phi::ModuleMap)

Let `Hom_MP` $= \text{Hom}(M,P)$, `Hom_NP` $= \text{Hom}(N,P)$ and `phi` $= \phi : N \to M$ a morphism.
Compute $\phi^{\ast} : \text{Hom}(M,P) \to \text{Hom}(N,P)$.
"""
function lift_homomorphism_contravariant(Hom_MP::ModuleFP, Hom_NP::ModuleFP, phi::ModuleMap)
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
  
  phi_lifted = hom(Hom_MP, Hom_NP, Vector{ModuleFPElem}([module_elem(Hom_NP, phi*homomorphism(f)) for f in gens(Hom_MP)]))
  return phi_lifted
end

@doc Markdown.doc"""
    lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, phi::ModuleMap)

Let `Hom_PM` $= \text{Hom}(P,M)$, `Hom_PN` $= \text{Hom}(P,N)$ and `phi` $= \phi : M \to N$ a morphism.
Compute $\phi_{\ast} : \text{Hom}(P,M) \to \text{Hom}(P,N)$.
"""
function lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, phi::ModuleMap)
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
    return hom(Hom_PM, Hom_PN, Vector{ModuleFPElem}([zero(Hom_PN) for _=1:ngens(Hom_PM)]))
  end
  phi_lifted = hom(Hom_PM, Hom_PN, Vector{ModuleFPElem}([module_elem(Hom_PN, homomorphism(f)*phi) for f in gens(Hom_PM)]))
  return phi_lifted
end

@doc Markdown.doc"""
    hom(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})

Apply $\text{Hom}(P,-)$ to `C`. Return the lifted chain complex.
"""
function hom(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})
  hom_chain = Hecke.map_type(C)[]
  hom_modules = [hom(P, domain(C.maps[1]))]
  hom_modules = vcat(hom_modules, [hom(P, codomain(f)) for f = C.maps])

  for i=1:length(C)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    push!(hom_chain, lift_homomorphism_covariant(A,B,map(C,i)))
  end
  return Hecke.ChainComplex(ModuleFP, hom_chain)
end

@doc Markdown.doc"""
    hom(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)

Apply $\text{Hom}(-,P)$ to `C`. Return the lifted chain complex.
"""
function hom(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
  hom_chain = Hecke.map_type(C)[]
  hom_modules = [hom(domain(C.maps[1]),P)]
  hom_modules = vcat(hom_modules, [hom(codomain(f), P) for f = C.maps])

  for i=1:length(C)
    A = hom_modules[i][1]
    B = hom_modules[i+1][1]

    push!(hom_chain, lift_homomorphism_contravariant(B,A,map(C,i)))
  end
  return Hecke.ChainComplex(ModuleFP, reverse(hom_chain))
end

#############################
@doc Markdown.doc"""
    homology(C::Hecke.ChainComplex{ModuleFP})

Compute all homology groups of `C`.
"""
function homology(C::Hecke.ChainComplex{ModuleFP})
  H = SubQuo[]
  for i=1:length(C)-1
    push!(H, quo(kernel(C.maps[i+1])[1], image(C.maps[i])[1]))
  end
  return H
end

@doc Markdown.doc"""
    homology(C::Hecke.ChainComplex{ModuleFP}, i::Int)

Compute the `i`-th homology of `C`.
"""
function homology(C::Hecke.ChainComplex{ModuleFP}, i::Int)
  @assert length(C) > 0 #TODO we need actually only the base ring
  if i == 0
    return kernel(map(C,1))[1]
  elseif i == length(C)
    return quo(obj(C,i),image(map(C,i))[1])
  elseif i < 0 || i > length(C)
    return FreeMod(base_ring(obj(C,1)),0)
  else
    return quo(kernel(map(C,i+1))[1], image(map(C,i))[1])
  end
end

#############################
# Ext
#############################
@doc Markdown.doc"""
    ext(M::ModuleFP, N::ModuleFP, i::Int)

Compute $\text{Ext}^i(M,N)$.
"""
function ext(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M)[1:end-2]
  lifted_resolution = hom(free_res, N) #TODO only three homs are necessary
  return homology(lifted_resolution,i)
end

#############################
# TODO ?
#############################
@doc Markdown.doc"""
    map_canonically(M::SubQuo, v::SubQuoElem)

Map the element `v` to an element of the module `M` using cached 
canonical homomorphisms between the parent Module of `v` and `M`.
"""
function map_canonically(M::SubQuo, v::SubQuoElem)
  N = parent(v)
  if N===M
    return v
  end

  # Breadth-First Search to find path to N:
  parent_hom = IdDict{SubQuo,ModuleMap}()
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
  result = v
  A = N
  while A !== M
    H = parent_hom[A]
    result = H(result)
    A = codomain(H)
  end
  return result
end

function all_canonical_maps(M::SubQuo, N::SubQuo)
  # from N to M

  all_paths = []

  function helper_dfs!(U::SubQuo, D::SubQuo, visited::Vector{<:ModuleMap}, path::Vector)
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

  helper_dfs!(N, M, Vector{ModuleMap}(), [])

  morphisms = Vector{ModuleMap}()
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
    register_morphism!(f::ModuleMap)

Cache the morphism `f` in the corresponding caches of the domain and codomain of `f`.
"""
function register_morphism!(f::ModuleMap)
  push!(domain(f).outgoing_morphisms, f)
  push!(codomain(f).incoming_morphisms, f)
end

#############################
#TODO move to Hecke
#  re-evaluate and use or not
#=
function differential(c::Hecke.ChainComplex, i::Int)
  return map(c,length(c)-i)
end

function module_in_complex(c::Hecke.ChainComplex, i::Int)
  return obj(c,length(c)-i)
end

getindex(c::Hecke.ChainComplex, i::Int) = module_in_complex(c,i)
=#

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

#=function default_ordering(R::MPolyRing)
  dp = Singular.ordering_dp(nvars(R))
  return Singular.sordering([dp.data[1], Singular.ordering_C().data[1]])
end=#

function default_ordering(F::FreeMod)
  if iszero(F)
    return degrevlex(gens(base_ring(F)))*Oscar.ModuleOrdering(F, Oscar.Orderings.ModOrdering(Vector{Int}(), :lex))
  end
  return degrevlex(gens(base_ring(F)))*lex(gens(F))
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
  p = _reduce(Singular.Module(base_ring(b), a), b)[1]
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
  return hom(F, G, Vector{FreeModElem}([i in indices ? G[findfirst(x->x==i,indices)] : zero(G) for i=1:ngens(F)]))
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
    if !isdefined(M, :quo)
      return false
    end
    singular_assure(M.quo.gens)
    reduced_standard_unit_vector = _reduce(M.quo.gens.SF(M.F[i]), singular_generators(std_basis(M.quo)))
    return iszero(reduced_standard_unit_vector)
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
    K = matrix_kernel(A)
    for i=1:size(K)[1], j=1:max_index
      if isunit(K[i,j])
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
    injection = FreeModuleHom(zero_module, M, [])
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
    isinjective(f::ModuleMap)

Test if `f` is injective.
"""
function isinjective(f::ModuleMap)
  return iszero(kernel(f)[1])
end

@doc Markdown.doc"""
    issurjective(f::ModuleMap)

Test if `f` is surjective.
"""
function issurjective(f::ModuleMap)
  return image(f)[1] == codomain(f)
end

@doc Markdown.doc"""
    isbijective(f::ModuleMap)

Test if `f` is bijective.
"""
function isbijective(f::ModuleMap)
  return isinjective(f) && issurjective(f)
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

  t2 = size(g2)[1]
  n = m
  m = s0*t1 + s1*t2
  rho::MatrixElem{T} = zero_matrix(R, m,n)
  for j=1:s0*t1
    b_vector = zero_matrix(R, 1,s0*t1)
    b_vector[1,j] = R(1)
    C = copy_and_reshape(b_vector, s0, t1)
    res1 = C*g1
    res2 = f1*C
    rho[j,1:length(res1)] = copy_and_reshape(res1, 1,length(res1))
    rho[j,length(res1)+1:end] = copy_and_reshape(res2, 1,length(res2))
  end
  for j=s0*t1+1:m
    b_vector = zero_matrix(R, 1,m-s0*t1)
    b_vector[1,j-s0*t1] = R(1)
    D = copy_and_reshape(b_vector, s1, t2)
    res2 = - D*g2
    rho[j,s0*t0+1:end] = copy_and_reshape(res2, 1,length(res2))
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
    to_subquotient_elem = function(H::ModuleMap)
      m = length(matrix(H))
      v = copy_and_reshape(matrix(H),1,m)
      tmp_sq = SubQuo(generator_matrix(SQ.sub)[:,1:m], generator_matrix(SQ.quo)[:,1:m])
      v = FreeModElem(sparse_row(v), FreeMod(R, length(v)))
      coeffs = coordinates(v, tmp_sq)
      return p(SubQuoElem(coeffs, SQ))
    end

    to_hom_map = MapFromFunc(to_homomorphism, to_subquotient_elem, SQ2, Hecke.MapParent(M, N, "homomorphisms"))
    set_attribute!(SQ2, :hom => (M, N), :module_to_hom_map => to_hom_map)

    return SQ2, to_hom_map
  else
    to_subquotient_elem = function(H::ModuleMap)
      m = length(matrix(H))
      v = copy_and_reshape(matrix(H),1,m)
      tmp_sq = SubQuo(generator_matrix(SQ.sub)[:,1:m], generator_matrix(SQ.quo)[:,1:m])
      v = FreeModElem(sparse_row(v), FreeMod(R, length(v)))
      coeffs = coordinates(v, tmp_sq)
      return SubQuoElem(coeffs, SQ)
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
