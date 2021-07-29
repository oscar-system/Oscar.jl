export presentation

# TODO replace asserts by error messages?

@doc Markdown.doc"""
  ModuleFP{T}

The abstract supertype of all modules. Here, all modules are finitely presented.
The type variable `T` refers to the type of the elements of the base ring.
"""
abstract type ModuleFP{T} end

@doc Markdown.doc"""
  CRing

Type union for commutative rings.
"""
const CRing = Union{MPolyRing, MPolyQuo{<:Oscar.MPolyElem}, MPolyRing_dec, MPolyQuo{<:Oscar.MPolyElem_dec}}
const CRingElem = Union{MPolyElem, MPolyQuoElem{<:Oscar.MPolyElem}, MPolyElem_dec, MPolyQuoElem{<:Oscar.MPolyElem_dec}}
#TODO: "fix" to allow QuoElem s as well...
# this requires
#  re-typeing of FreeModule
#  typing of BiModArray
# ... and all the rest.
# parametrization has to be by elem_type(coeff_ring) and not, like currently, the bottom coeff ring
# Also: qring is a Singular native. So it needs to be added to the ring creation

abstract type ModuleFPHom end
@doc Markdown.doc"""
  ModuleMap{T1, T2}

The abstract supertype of module morphisms.
`T1` and `T2` are the types of domain and codomain respectively.
"""
abstract type ModuleMap{T1, T2} <: Map{T1, T2, Hecke.HeckeMap, ModuleFPHom} end


@doc Markdown.doc"""
  FreeMod{T <: CRingElem} <: ModuleFP{T}

The type of free modules.
Free modules are determined by their base ring, the rank and the names of 
the (standard) generators.
Moreover, canonical in- and outgoing morphisms are stored if the corresponding
option is set in suitable functions.
`FreeMod{T}` is a subtype of `ModuleFP{T}`.
"""
mutable struct FreeMod{T <: CRingElem} <: ModuleFP{T}
  R::CRing
  n::Int
  S::Vector{Symbol}

  ingoing_morphisms::Vector{<:ModuleMap}
  outgoing_morphisms::Vector{<:ModuleMap}

  AbstractAlgebra.@declare_other

  function FreeMod{T}(n::Int,R::CRing,S::Vector{Symbol}) where T <: CRingElem
    r = new{elem_type(R)}()
    r.n = n
    r.R = R
    r.S = S

    r.ingoing_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
end

@doc Markdown.doc"""
  FreeMod(R::CRing, n::Int, name::String = "e"; cached::Bool = false)

Construct a free module over the ring `R` with rank `n`.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for the unit vectors.
"""
function FreeMod(R::CRing, n::Int, name::String = "e"; cached::Bool = false) # TODO cached?
  return FreeMod{elem_type(R)}(n, R, [Symbol("$name[$i]") for i=1:n])
end
@doc Markdown.doc"""
  free_module(R::CRing, n::Int, name::String = "e"; cached::Bool = false)

Construct a free module over the ring `R` with rank `n`.
Additionally one can provide names for the generators. If one does 
not provide names for the generators, the standard names e_i are used for the unit vectors.
"""
# TODO is this function neccessary?
free_module(R::CRing, n::Int, name::String = "e"; cached::Bool = false) = FreeMod(R, n, name, cached = cached)

#=XXX this cannot be as it is inherently ambigous
  - FreeModule(R, n)
  - direct sum of rings, ie. a ring
  - set of n-th powers of R
thus the "category" needs to be set explicitly

^(R::Ring_dec, n::Int) = FreeModule(R, n)
=#

function AbstractAlgebra.extra_name(F::FreeMod)
  return nothing
end

function (F::FreeMod)()
  return FreeModuleElem(sparse_row(base_ring(F)), F)
end

function show(io::IO, F::FreeMod)
  @show_name(io, F)
  @show_special(io, F)

  print(io, "Free module of rank $(F.n) over ")
  print(IOContext(io, :compact =>true), F.R)
#=
  i = 1
  while i < dim(F)
    d = F.d[i]
    j = 1
    while i+j <= dim(F) && d == F.d[i+j]
      j += 1
    end
    print(IOContext(io, :compact => true), F.R, "^$j")
    print(IOContext(io, :compact => true), "(", -d, ")")
    if i+j < dim(F)
      print(io, " + ")
    end
    i += j
  end
  =#
end

@doc Markdown.doc"""
  rank(F::FreeMod)
  ngens(F::FreeMod)
  dim(F::FreeMod)

Return the rank of `F`.
"""
dim(F::FreeMod) = F.n
rank(F::FreeMod) = F.n
ngens(F::FreeMod) = dim(F)

@doc Markdown.doc"""
  ==(F::FreeMod, G::FreeMod)

Equality check between `F` and `G`. `F` and `G` are equal iff their 
base rings, ranks and (basis) element names are equal.
"""
# two free modules are equal if the rank and the ring are
function ==(F::FreeMod, G::FreeMod)
  # TODO it this enough or e.g. stored morphisms also be considered?
  return F.R == G.R && rank(F) == rank(G) && F.S == G.S
end

@doc Markdown.doc"""
  iszero(F::FreeMod)

Check if `F` is the zero module, that is if `F` has rank 0.
"""
function iszero(F::FreeMod)
  return rank(F) == 0
end

@doc Markdown.doc"""
  FreeModuleElem{T}

The type of free module elements. A free module element
is determined by a sparse row (`SRow`) which contains the coords 
wrt to the standard basis and the parent free module.

# Example
```jldoctest
julia> R, (x,y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = FreeMod(R,2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y over Rational Field

julia> v = FreeModuleElem(sparse_row(R, [(1,x),(2,y)]),F)
(x)*e[1] + (y)*e[2]

julia> v == x*F[1]+y*F[2]
true

```
"""
struct FreeModuleElem{T}
  coords::SRow{T} # also usable via coeffs()
  parent::FreeMod{T}
end

function in(v::FreeModuleElem, M::ModuleFP)
  return parent(v) === M
end

@doc Markdown.doc"""
  coords(v::FreeModuleElem)

Return the entries (wrt to the standard basis) of `v` as a sparse row.
"""
function coords(v::FreeModuleElem)
  return v.coords
end

@doc Markdown.doc"""
  coeffs(v::FreeModuleElem)

Return the entries (wrt to the standard basis) of `v` as a sparse row.
"""
function coeffs(v::FreeModuleElem)
  return coords(v)
end

@doc Markdown.doc"""
  repres(v::FreeModuleElem)

Return just `v`. This function exists for compatiblity (with subquotient elements) reasons.
"""
function repres(v::FreeModuleElem)
  return v
end

function getindex(v::FreeModuleElem, i::Int)
  if isempty(v.coords)
    return zero(base_ring(v.parent))
  end
  return v.coords[i]
end

elem_type(::Type{FreeMod{T}}) where {T} = FreeModuleElem{T}
parent_type(::Type{FreeModuleElem{T}}) where {T} = FreeMod{T}
elem_type(::FreeMod{T}) where {T} = FreeModuleElem{T}
parent_type(::FreeModuleElem{T}) where {T} = FreeMod{T}

function show(io::IO, e::FreeModuleElem)
  if length(e.coords) == 0
    print(io, 0)
    return
  end
  i = 1
  while i <= length(e.coords)
    print(io, "(", e.coords.values[i], ")*", e.parent.S[e.coords.pos[i]])
    if i < length(e.coords)
      print(io, " + ")
    end
    i += 1
  end
end

@doc Markdown.doc"""
  basis(F::FreeMod)

Return the standard basis of `F`.
"""
function basis(F::FreeMod)
  bas = elem_type(F)[]
  for i=1:dim(F)
    s = Hecke.sparse_row(F.R, [(i, F.R(1))])
    push!(bas, FreeModuleElem(s, F))
  end
  return bas
end
@doc Markdown.doc"""
  gens(F::FreeMod)

Return the (canonical) generators of the free module `F`.
"""
gens(F::FreeMod) = basis(F)

@doc Markdown.doc"""
  gen(F::FreeMod, i::Int)

Return the `i`th generator of `F`, that is the `i`th unit vector.
"""
function gen(F::FreeMod, i::Int)
  @assert 0 < i <= ngens(F)
  s = Hecke.sparse_row(F.R, [(i, F.R(1))])
  return FreeModuleElem(s, F)
end

function Base.getindex(F::FreeMod, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

@doc Markdown.doc"""
  base_ring(F::FreeMod)

Return the base ring of the free module `F`.
"""
base_ring(F::FreeMod) = F.R

#TODO: Parent - checks everywhere!!!

# the negative of a free module element
-(a::FreeModuleElem) = FreeModuleElem(-a.coords, a.parent)

# Addition of free module elements
function +(a::FreeModuleElem, b::FreeModuleElem)
   check_parent(a, b)
   return FreeModuleElem(a.coords+b.coords, a.parent)
end

# Subtraction of free module elements
function -(a::FreeModuleElem, b::FreeModuleElem)
    check_parent(a,b)
    return FreeModuleElem(a.coords-b.coords, a.parent)
end

# Equality of free module elements
function ==(a::FreeModuleElem, b::FreeModuleElem) 
    check_parent(a,b)
    return a.coords == b.coords
end

# scalar multiplication with polynomials, integers
function *(a::MPolyElem_dec, b::FreeModuleElem)
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
  end
  return FreeModuleElem(a*b.coords, b.parent)
end
function *(a::MPolyElem, b::FreeModuleElem) 
  if parent(a) !== base_ring(parent(b))
    error("elements not compatible")
  end
  FreeModuleElem(a*b.coords, b.parent)
end
*(a::Int, b::FreeModuleElem) = FreeModuleElem(a*b.coords, b.parent)
*(a::Integer, b::FreeModuleElem) = FreeModuleElem(b.parent.R(a)*b.coords, b.parent)
*(a::fmpq, b::FreeModuleElem) = FreeModuleElem(b.parent.R(a)*b.coords, b.parent)

@doc Markdown.doc"""
  zero(F::FreeMod)

Return the zero element of the free module `F`.
"""
zero(F::FreeMod) = FreeModuleElem(sparse_row(F.R, Tuple{Int, elem_type(F.R)}[]), F)

@doc Markdown.doc"""
  parent(a::FreeModuleElem)

Return the free module where `a` lives in.
"""
parent(a::FreeModuleElem) = a.parent

@doc Markdown.doc"""
  iszero(a::FreeModuleElem)

Check whether the free module element `a` is zero.
"""
iszero(a::FreeModuleElem) = Hecke.iszero(a.coords)

# data structure for a generating systems for submodules
# contains structures for the generators, the corresponding module on the Singular side, 
# the embedding free module, the embedding free module on the Singular side
# subquotients will be built from a tuple of submodules which again are given by 
# generating sets. In this way, the Singular stuff is hidden on the higher structures
# and all the conversion is taken care of here
# a module generating system is generated from an array of free module elements
# the fields are called O,S,F,SF rename?
#
# The same could be done rather on the level of vectors, that might be preferable if 
# performance is ok.
#
mutable struct ModuleGens{T}
  O::Vector{FreeModuleElem{T}}
  S::Singular.smodule
  F::FreeMod{T}
  SF::Singular.FreeMod

  # ModuleGens from an Array of Oscar free module elements, specifying the free module 
  # and Singular free module, only useful indirectly
  function ModuleGens{T}(O::Vector{<:FreeModuleElem}, F::FreeMod{T}, SF::Singular.FreeMod) where {T}
    r = new{T}()
    r.O = O
    r.SF = SF
    r.F = F
    return r
  end

  # ModuleGens from a Singular submodule
  function ModuleGens{S}(F::FreeMod{S}, s::Singular.smodule) where {S} # FreeMod is neccessary due to type S
    r = new{S}()
    r.F = F
    if Singular.ngens(s) == 0
      r.SF = Singular.FreeModule(base_ring(s), 0)
    else
      r.SF = parent(s[1])
    end
    r.S = s
    return r
  end
end

@doc Markdown.doc"""
  ModuleGens(O::Vector{<:FreeModuleElem}, F::FreeMod{T}, SF::Singular.FreeMod) where T

Construct `ModuleGens` from an array of Oscar free module elements, specifying the Oscar free module 
and Singular free module. 
This function is only useful indirectly.
"""
ModuleGens(O::Vector{<:FreeModuleElem}, F::FreeMod{T}, SF::Singular.FreeMod) where T = ModuleGens{T}(O, F, SF)

@doc Markdown.doc"""
  ModuleGens(F::FreeMod{S}, s::Singular.smodule) where {S}

Construct `ModuleGens` from a given Singular submodule.
"""
ModuleGens(F::FreeMod{S}, s::Singular.smodule) where {S} = ModuleGens{S}(F, s)

@doc Markdown.doc"""
  ModuleGens(O::Vector{<:FreeModuleElem})

Construct `ModuleGens` from an array of Oscar free module elements.
Note: The array must not be empty.
"""
# TODO Empty generating set
function ModuleGens(O::Vector{<:FreeModuleElem})
  @assert length(O) > 0
  SF = singular_module(parent(O[1]))
  return ModuleGens(O, SF)
end

@doc Markdown.doc"""
  ModuleGens(O::Vector{<:FreeModuleElem}, F::FreeMod{T}) where {T}

Construct `ModuleGens` from an array of Oscar free module elements, specifying the Oscar free module.
Note: The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModuleElem}, F::FreeMod{T}) where {T}
  SF = singular_module(F)
  return ModuleGens{T}(O, F, SF)
end

@doc Markdown.doc"""
  ModuleGens(O::Vector{<:FreeModuleElem}, SF::Singular.FreeMod)

Construct `ModuleGens` from an array of Oscar free module elements, specifying the Singular free module.
Note: The array might be empty.
"""
function ModuleGens(O::Vector{<:FreeModuleElem}, SF::Singular.FreeMod)
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

# TODO remove output saying defined on the Singular side?
function show(io::IO, F::ModuleGens)
  println(io, "Array of length ", length(F))
  for i=1:length(F)
    if isassigned(F.O, i)
      println(io, i, " -> ", F.O[i])
    end
  end
  if isdefined(F, :S)
    println(io, "defined on the Singular side")
  end
end

@doc Markdown.doc"""
  length(F::ModuleGens)

Return the number of elements of the module generating set.
Note: This is not the length in the mathematical sense!
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
    F.O[i] = convert(F.F, singular_generators(F)[i])
  end
  return oscar_generators(F)[i]
end

# i-th entry of module generating set on Singular side
# Todo: clean up, convert or assure
function getindex(F::ModuleGens, ::Val{:S}, i::Int)
  if !isdefined(F, :S)
    F.S = Singular.Module(base_ring(F.SF), [convert(F.SF,x) for x = oscar_generators(F)]...)
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
    F.O = [convert(F.F, singular_generators(F)[i]) for i=1:Singular.ngens(singular_generators(F))]
  end
end

@doc Markdown.doc"""
  singular_assure(F::ModuleGens)

If fields of `F` from the Singular side are not defined, they
are computed, given the OSCAR side.
"""
function singular_assure(F::ModuleGens)
  if !isdefined(F, :S)
    if length(F) == 0
      singular_ring = base_ring(F.SF)
      F.S = Singular.Module(singular_ring, Singular.vector(singular_ring, singular_ring(0)))
      return 
    end
    F.S = Singular.Module(base_ring(F.SF), [convert(F.SF,x) for x = oscar_generators(F)]...)
    return
  end
  #F[Val(:S), 1]
end

# i-th entry of module generating set (taken from Oscar side)
getindex(F::ModuleGens, i::Int) = getindex(F, Val(:O), i)

@doc Markdown.doc"""
  singular_module(F::FreeMod)

Create a Singular module from an OSCAR free module.
"""
function singular_module(F::FreeMod)
  Sx = singular_ring(base_ring(F))
  return Singular.FreeModule(Sx, dim(F))
end

@doc Markdown.doc"""
  convert(SF::Singular.FreeMod, m::FreeModuleElem)

Convert an OSCAR free module element to the Singular side.
"""
function convert(SF::Singular.FreeMod, m::FreeModuleElem)
  g = Singular.gens(SF)
  e = SF()
  Sx = base_ring(SF)
  for (p,v) = m.coords
    e += Sx(v)*g[p]
  end
  return e
end

@doc Markdown.doc"""
  convert(F::FreeMod, s::Singular.svector)

Convert a Singular vector to a free module element on the OSCAR side.
"""
function convert(F::FreeMod, s::Singular.svector)
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
  return FreeModuleElem(sparse_row(base_ring(F), pv), F)
end

@doc Markdown.doc"""
  FreeModuleHom{T1, T2} <: ModuleMap{T1, T2} 

Data structure for morphisms where the domain is a free module (`FreeMod`).
`T1` and `T2` are the types of domain and codomain respectively.
`FreeModuleHom` is a subtype of `ModuleMap`.
When computed, the corresponding matrix (via `matrix()`) and inverse isomorphism
(in case there exists one) (via `inv()`) are cached.
"""
# data structure for homomorphisms of free modules
mutable struct FreeModuleHom{T1, T2} <: ModuleMap{T1, T2} 
  matrix::MatElem
  header::MapHeader
  inverse_isomorphism::ModuleMap
  Hecke.@declare_other

  # generate homomorphism of free modules from F to G where the Array a contains the images of
  # the generators of F
  function FreeModuleHom{T,S}(F::FreeMod{T}, G::S, a::Vector{<:Any}) where {T, S}
    @assert all(x->parent(x) === G, a)
    @assert length(a) == ngens(F)
    r = new{typeof(F), typeof(G)}()
    function im_func(x::FreeModuleElem)
      b = zero(G)
      for (i,v) = x.coords
        b += v*a[i]
      end
      return b
    end
    function pr_func(x)
      @assert parent(x) === G
      c = coordinates(repres(x), sub(G, a))
      return FreeModuleElem(c, F)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)

    return r
  end

  function FreeModuleHom{T,S}(F::FreeMod{T}, G::S, mat::MatElem{T}) where {T,S}
    @assert nrows(mat) == ngens(F)
    @assert ncols(mat) == ngens(G)
    if typeof(G) <: FreeMod
      hom = FreeModuleHom(F, G, [FreeModuleElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)])
    else
      hom = FreeModuleHom(F, G, [SubQuoElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)])
    end
    hom.matrix = mat
    return hom
  end
end

@doc Markdown.doc"""
  FreeModuleHom(F::FreeMod{T}, G::S, a::Vector{<:Any}) where {T, S}

Construct the morphism $F → G$ where `F[i]` is mapped to `a[i]`.
In particular, `ngens(F) == length(a)` must hold.
"""
FreeModuleHom(F::FreeMod{T}, G::S, a::Vector{<:Any}) where {T, S} = FreeModuleHom{T,S}(F, G, a)

@doc Markdown.doc"""
  FreeModuleHom(F::FreeMod{T}, G::S, mat::MatElem{T}) where {T,S}

Construct the morphism $F → G$ corresponding to the matrix `mat`.
"""
FreeModuleHom(F::FreeMod{T}, G::S, mat::MatElem{T}) where {T,S} = FreeModuleHom{T,S}(F, G, mat)

#=function Base.getproperty(f::FreeModuleHom, s::Symbol)
  if s == :matrix
    error("deprecated")
    if !isdefined(f, s)
      D = domain(f)
      C = codomain(f)
      R = base_ring(D)
      matrix = zero_matrix(R, D.n, ngens(C))
      for i=1:D.n
        image_of_gen = f(D[i])
        for j=1:ngens(C)
          matrix[i,j] = image_of_gen[j]
        end
      end
      setfield!(f, s, matrix)
    end
  end
  return getfield(f, s)
end=#

@doc Markdown.doc"""
  matrix(f::FreeModuleHom)

Compute the matrix corresponding to the morphism `f`.
The matrix is cached.
"""
function matrix(f::FreeModuleHom)
  if !isdefined(f, :matrix)
    D = domain(f)
    C = codomain(f)
    R = base_ring(D)
    matrix = zero_matrix(R, D.n, ngens(C))
    for i=1:D.n
      image_of_gen = f(D[i])
      for j=1:ngens(C)
        matrix[i,j] = image_of_gen[j]
      end
    end
    setfield!(f, :matrix, matrix)
  end
  return f.matrix
end

(h::FreeModuleHom)(a::FreeModuleElem) = image(h, a)

@doc Markdown.doc"""
  hom(F::FreeMod{T}, G, a)

Return the morphism $F → G$ where `F[i]` is mapped to `a[i]`.
"""
hom(F::FreeMod{T}, G, a) where {T} = FreeModuleHom(F, G, a)

@doc Markdown.doc"""
  identity_map(M::ModuleFP)

Return the identity map $id_M$.
"""
function identity_map(M::ModuleFP)
  return hom(M, M, gens(M))
end

@doc Markdown.doc"""
  SubModuleOfFreeModule{T} <: ModuleFP{T}

Data structure for submodules of free modules. `SubModuleOfFreeModule` shouldn't be
used by the end user.
When computed, a standard basis (computed via `std_basis()`) and generating matrix (that is the rows of the matrix
generate the submodule) (computed via `matrix()`) are cached.
"""
mutable struct SubModuleOfFreeModule{T} <: ModuleFP{T}
  F::FreeMod{T}
  gens::ModuleGens
  std_basis::ModuleGens
  matrix::MatElem

  function SubModuleOfFreeModule{R}(F::FreeMod{R}, gens::Vector{<:FreeModuleElem}) where {R}
    @assert all(x -> parent(x) === F, gens)
    r = new{R}()
    r.F = F
    r.gens = ModuleGens(gens, F)
    return r
  end

  function SubModuleOfFreeModule{R}(F::FreeMod{R}, singular_module::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    r.gens = ModuleGens(F, singular_module)
    if singular_module.isGB
      r.std_basis = r.gens
    end
    return r
  end
  
  function SubModuleOfFreeModule{R}(F::FreeMod{R}, gens::ModuleGens) where {R}
    r = new{R}()
    r.F = F
    r.gens = gens
    if singular_generators(gens).isGB
      r.std_basis = r.gens
    end
    return r
  end

  function SubModuleOfFreeModule{L}(F::FreeMod{L}, A::MatElem{L}) where {L}
    r = new{L}()
    #R = base_ring(A)
    #F = FreeMod(R, ncols(A))
    r.F = F
    O = [FreeModuleElem(sparse_row(A[i,:]), F) for i in 1:nrows(A)]
    r.gens = ModuleGens(O, F)
    r.matrix = A
    return r
  end
end

@doc Markdown.doc"""
  SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModuleElem}) where {R}

Construct the submodule of `F` generated by the elements of `gens` (the elements of 
`gens` must live in `F`).
"""
SubModuleOfFreeModule(F::FreeMod{R}, gens::Vector{<:FreeModuleElem}) where {R} = SubModuleOfFreeModule{R}(F, gens)

SubModuleOfFreeModule(F::FreeMod{R}, singular_module::Singular.smodule) where {R} = SubModuleOfFreeModule{R}(F, singular_module)

@doc Markdown.doc"""
  SubModuleOfFreeModule(F::FreeMod{R}, gens::ModuleGens{R}) where {R}

Construct the submodule of `F` generated by `gens`.
"""
SubModuleOfFreeModule(F::FreeMod{R}, gens::ModuleGens{R}) where {R} = SubModuleOfFreeModule{R}(F, gens)

@doc Markdown.doc"""
  SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}) where {L}

Construct the submodule generated by the rows of `A`. The embedding free
module is `F`. In particular, `rank(F) == ncols(A)` must hold.
"""
SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}) where {L} = SubModuleOfFreeModule{L}(F, A)

@doc Markdown.doc"""
  SubModuleOfFreeModule(A::MatElem{L}) where {L} 

Construct the submodule generated by the rows of `A`.
Note: The embedding free module of the submodule is constructed by the function and therefore
not compatible with free modules that are defined by the user or other functions. 
If you need compatiblity use `SubModuleOfFreeModule(F::FreeMod{L}, A::MatElem{L}) where {L}`.
"""
function SubModuleOfFreeModule(A::MatElem{L}) where {L} 
  R = base_ring(A)
  F = FreeMod(R, ncols(A))
  return SubModuleOfFreeModule{L}(F, A)
end

function Base.getindex(M::SubModuleOfFreeModule, i::Int)
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
  std_basis(submod::SubModuleOfFreeModule)

Compute a standard basis of `submod`.
"""
function std_basis(submod::SubModuleOfFreeModule)
  if !isdefined(submod, :std_basis)
    submod.std_basis = groebner_basis(submod.gens)
  end
  return submod.std_basis
end

@doc Markdown.doc"""
  matrix(submod::SubModuleOfFreeModule)

Return the generators of `submod` in matrix-form, that is the rows of the 
matrix generate `submod`.
"""
function matrix(submod::SubModuleOfFreeModule)
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
  isgenerated_by_unit_vectors(M::SubModuleOfFreeModule)

Check if `M` is generated by the unit vectors $e_1,\dots,e_n$ (in this order). 
This function is only used for `show`-purposes.
"""
function isgenerated_by_unit_vectors(M::SubModuleOfFreeModule)
  m,n = size(matrix(M))
  if m != n 
    return false
  end
  return matrix(M) == identity_matrix(base_ring(M),n)
end

function show(io::IO, M::SubModuleOfFreeModule)
  if ngens(M) == 1
    println(io, "Submodule with ", ngens(M), " generator")
  else
    println(io, "Submodule with ", ngens(M), " generators")
  end
  for i=1:ngens(M)
    if isassigned(M.gens.O, i)
      println(io, i, " -> ", M[i])
    end
  end
  if isdefined(M.gens, :S)
    println(io, "defined on the Singular side")
  end
end

function length(M::SubModuleOfFreeModule)
  #TODO replace all calls by ngens
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

Return the generators of `M` as an array of `FreeModuleElem`s.
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
  #return M.gens[Val(:O), i]
end

@doc Markdown.doc"""
  sum(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Compute $M+N$.
"""
function sum(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  @assert M.F === N.F
  return SubModuleOfFreeModule(M.F, vcat(collect(M.gens), collect(N.gens)))
end

function issubset(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  @assert M.F === N.F
  M_mod_N = _reduce(singular_generators(std_basis(M)), singular_generators(std_basis(N)))
  return iszero(M_mod_N)
end

@doc Markdown.doc"""
  ==(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)

Check for equality. For two submodules of free modules to be equal their embedding 
free modules must be identical (`===`) and the generators must generate equal submodules.
"""
function ==(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule)
  @assert M.F === N.F
  M_mod_N = _reduce(singular_generators(std_basis(M)), singular_generators(std_basis(N)))
  N_mod_M = _reduce(singular_generators(std_basis(N)), singular_generators(std_basis(M)))
  return iszero(M_mod_N) && iszero(N_mod_M)
end

#+(M::SubModuleOfFreeModule, N::SubModuleOfFreeModule) = sum(M, N)

@doc Markdown.doc"""
  SubQuo{T} <: ModuleFP{T}

The type of subquotient modules.
A subquotient module $M$ is a module where $M = A + B / B$ where $A$ and $B$ are 
submodules of a free module.
$A$, $B$ and $A+B$ (they have type `SubModuleOfFreeModule`) as well as the embedding
free module are stored. 
One can construct ordinary submodules of free modules by not giving $B$.
Moreover, canonical in- and outgoing morphisms are stored if the corresponding
option is set in suitable functions.
`SubQuo{T}` is a subtype of `ModuleFP{T}`.
"""
mutable struct SubQuo{T} <: ModuleFP{T}
  #meant to represent sub+ quo mod quo - as lazy as possible
  F::FreeMod{T}
  sub::SubModuleOfFreeModule
  quo::SubModuleOfFreeModule
  sum::SubModuleOfFreeModule

  ingoing_morphisms::Vector{<:ModuleMap}
  outgoing_morphisms::Vector{<:ModuleMap} # TODO is it possible to make ModuleMap to SubQuoHom?

  AbstractAlgebra.@declare_other

  function SubQuo{R}(sub::SubModuleOfFreeModule{R}) where {R}
    r = new{R}()
    r.F = sub.F
    r.sub = sub
    r.sum = r.sub

    r.ingoing_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  function SubQuo{R}(sub::SubModuleOfFreeModule{R}, quo::SubModuleOfFreeModule{R}) where {R}
    @assert sub.F === quo.F
    r = new{R}()
    r.F = sub.F
    r.sub = sub
    r.quo = quo
    r.sum = sum(r.sub, r.quo)

    r.ingoing_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  function SubQuo{R}(F::FreeMod{R}, O::Vector{<:FreeModuleElem}) where {R}
    r = new{R}()
    r.F = F
    #r.sub = ModuleGens(O, F, singular_module(F))
    r.sub = SubModuleOfFreeModule(F, O)
    r.sum = r.sub

    r.ingoing_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  function SubQuo{L}(S::SubQuo{L}, O::Vector{<:FreeModuleElem}) where {L} #TODO to be replaced by quo
    r = new{L}()
    r.F = S.F
    r.sub = S.sub
    #r.quo = ModuleGens(O, S.F, S.sub.SF)
    O_as_submodule = SubModuleOfFreeModule(S.F, O)
    r.quo = isdefined(S,:quo) ? sum(S.quo,O_as_submodule) : O_as_submodule
    #r.sum = ModuleGens(vcat(collect(r.sub), collect(r.quo)), S.F, S.sub.SF)
    r.sum = sum(r.sub, r.quo)

    r.ingoing_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  #=function SubQuo(S::SubQuo, O::Array{<:SubQuoElem, 1})
    @assert all(x->x.parent === S, O)
    r = SubQuo(S.F, [x.repres for x in O])
    r.quo = S.quo
    r.sum = sum(r.sub, r.quo)
    return r
  end=#
  function SubQuo{R}(F::FreeMod{R}, s::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    #r.sub = ModuleGens(F, s)
    r.sub = SubModuleOfFreeModule(F, s)
    r.sum = r.sub

    r.ingoing_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  function SubQuo{R}(F::FreeMod{R}, s::Singular.smodule, t::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    #r.sub = ModuleGens(F, s)
    r.sub = SubModuleOfFreeModule(F, s)
    #r.quo = ModuleGens(F, t)
    r.quo = SubModuleOfFreeModule(F, t)
    #r.sum = ModuleGens(vcat(collect(r.sub), collect(r.quo)))
    r.sum = sum(r.sub, r.quo)

    r.ingoing_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
end

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
  SubQuo(F::FreeMod{R}, O::Vector{<:FreeModuleElem}) where {R}

Construct the module generated by the elements of `O` as a subquotient.
The elements of `O` must live in `F`.

# Example
```jldoctest
julia> R, (x,y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = FreeMod(R,2)
Free module of rank 2 over Multivariate Polynomial Ring in x, y over Rational Field

julia> O = [x*F[1]+F[2],y*F[2]]
2-element Vector{Oscar.FreeModuleElem{fmpq_mpoly}}:
 (x)*e[1] + (1)*e[2]
 (y)*e[2]

julia> M = SubQuo(F, O)
Submodule with 2 generators
1 -> (x)*e[1] + (1)*e[2]
2 -> (y)*e[2]

represented as subquotient with no relations.


```
"""
SubQuo(F::FreeMod{R}, O::Vector{<:FreeModuleElem}) where {R} = SubQuo{R}(F, O)

@doc Markdown.doc"""
  SubQuo(S::SubQuo{L}, O::Vector{<:FreeModuleElem}) where {L}

Construct a subquotient where the generators are those of `S` and the relations are 
the union of `O` and the relations of `S`.
"""
SubQuo(S::SubQuo{L}, O::Vector{<:FreeModuleElem}) where {L} = SubQuo{L}(S, O)

SubQuo(F::FreeMod{R}, s::Singular.smodule) where {R} = SubQuo{R}(F, s)

SubQuo(F::FreeMod{R}, s::Singular.smodule, t::Singular.smodule) where {R} = SubQuo{R}(F, s, t)

@doc Markdown.doc"""
  SubQuo(F::FreeMod{R}, A::MatElem{R}, B::MatElem{R}) where {R}

Let $A$ be an $n_1×m$-matrix, $B$ an $n_2×m$-matrix over the ring $R$ and 
$F = R^m$. Return the subquotient $(\im(A) + \im(B))/ \im(B)$, where the 
embedding free module is $F$.
"""
function SubQuo(F::FreeMod{R}, A::MatElem{R}, B::MatElem{R}) where {R}
  @assert ncols(A) == ncols(B) == rank(F)
  return SubQuo(SubModuleOfFreeModule(F, A), SubModuleOfFreeModule(F, B))
end

@doc Markdown.doc"""
  SubQuo(A::MatElem{R}, B::MatElem{R}) where {R}

Let $A$ be an $n_1×m$-matrix, $B$ an $n_2×m$-matrix over the ring $R$. 
Then return the subquotient $(\im(A) + \im(B))/ \im(B)$. 
Note: The embedding free module of the subquotient is constructed by the function and therefore
not compatible with free modules $R^m$ that are defined by the user or other functions. 
If you need compatiblity use `SubQuo(F::FreeMod{R}, A::MatElem{R}, B::MatElem{R})`.
"""
function SubQuo(A::MatElem{R}, B::MatElem{R}) where {R}
  @assert ncols(A) == ncols(B)
  S = base_ring(A)
  F = FreeMod(S, ncols(A))
  return SubQuo(SubModuleOfFreeModule(F, A), SubModuleOfFreeModule(F, B))
end

function show(io::IO, SQ::SubQuo)
  @show_name(io, SQ)
  @show_special(io, SQ)

  if isdefined(SQ, :quo)
    println(io, "Subquotient of ", SQ.sub, "by ", SQ.quo)
  else
    #println(io, "Subquotient by ", SQ.sub)
    println(io, SQ.sub)
    println("represented as subquotient with no relations.")
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
    if isgenerated_by_unit_vectors(SQ.sub)
      println("Cokernel of ", matrix(SQ.quo))
    else
      println("Subquotient with of image of")
      display(matrix(SQ.sub))
      println("by image of")
      display(matrix(SQ.quo))
      #println("Subquotient with of image of ", SQ.sub.matrix, " by image of ", SQ.quo.matrix)
    end
  else
    println("Image of ", matrix(SQ.sub))
  end
end

@doc Markdown.doc"""
  cokernel(f::FreeModuleHom)

Let $f : F → G$ be a morphism between free modules. Return the cokernel of $f$ as a subquotient.
"""
function cokernel(f::FreeModuleHom)
  @assert typeof(codomain(f)) <: FreeMod
  return quo(codomain(f), image(f)[1])
end

@doc Markdown.doc"""
  cokernel(F::FreeMod{R}, A::MatElem{R}) where R

Let $F = R^m$ and $A$ an $n×m$-matrix. Return the subquotient $F / \im(A)$.
"""
function cokernel(F::FreeMod{R}, A::MatElem{R}) where R
  return cokernel(matrix_to_map(F,A))
end

@doc Markdown.doc"""
  cokernel(A::MatElem)

Let $A$ be an $n×m$-matrix over the ring $R$. Then return the subquotient $R^m / \im(A)$. 
Note: The free module $R^m$ is constructed by the function and therefore not compatible with 
free modules $R^m$ that are defined by the user or other functions. If you need compatiblity
use `cokernel(F::FreeMod{R}, A::MatElem{R})`.
"""
function cokernel(A::MatElem)
  return cokernel(matrix_to_map(A))
end

@doc Markdown.doc"""
  ==(M::SubQuo{T}, N::SubQuo{T}) where {T}

Check for equality. For two subquotients to be equal their embedding free modules 
must be identical (`===`) and the generators respectively (generators + relations) must 
generate equal submodules.
"""
function ==(M::SubQuo{T}, N::SubQuo{T}) where {T}
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

Compute $M+N$ along with the inclusion morphisms $M → M+N$ and $N → M+N$.
"""
function sum(M::SubQuo{T},N::SubQuo{T}) where T
  #TODO use SubModuleOfFreeModule instead of matrices
  gm1,gm2 = size(matrix(M.sub))
  gn1,gn2 = size(matrix(N.sub))
  R = base_ring(M)

  n_rel = isdefined(N, :quo) ? matrix(N.quo) : zero_matrix(R, 1, gn2)
  m_rel = isdefined(M, :quo) ? matrix(M.quo) : zero_matrix(R, 1, gm2)

  if (m_rel == n_rel) || Set([m_rel[i,:] for i=1:size(m_rel)[1]])==Set([n_rel[j,:] for j=1:size(n_rel)[1]]) || M.quo == N.quo
    SQ = SubQuo(vcat(matrix(M.sub), matrix(N.sub)), m_rel)

    # injection maps:
    M_mat = hcat(identity_matrix(R,gm1), zero_matrix(R,gm1,gn1))
    iM = SubQuoHom(M,SQ,M_mat)

    N_mat = hcat(zero_matrix(R,gn1,gm1), identity_matrix(R,gn1))
    iN = SubQuoHom(N,SQ,N_mat)

    register_morphism!(iM)
    register_morphism!(iN)

    return SQ, iM, iN
  end
  throw(ArgumentError("img(M.relations) != img(N.relations)"))
end

@doc Markdown.doc"""
  Base.:+(M::SubQuo{T},N::SubQuo{T}) where T

Compute $M+N$.
"""
function Base.:+(M::SubQuo{T},N::SubQuo{T}) where T
  return sum(M,N)[1]
end

@doc Markdown.doc"""
  Base.:intersect(M::SubQuo{T}, N::SubQuo{T}) where T

Compute the intersection $M ∩ N$ along with the inclusion morphisms $M∩N → M$ and $M∩N → N$.
"""
function Base.:intersect(M::SubQuo{T}, N::SubQuo{T}) where T
  #TODO allow task as argument?
  @assert free_module(M) === free_module(N)
  n_rel = matrix(N.quo)
  m_rel = matrix(M.quo)

  if (m_rel == n_rel) || Set([m_rel[i,:] for i=1:size(m_rel)[1]])==Set([n_rel[j,:] for j=1:size(n_rel)[1]]) || M.quo == N.quo
    n = size(matrix(N.sub))
    m = size(matrix(M.sub))
    if n[2]!=m[2]
      throw(DimensionMismatch("Matrices have different number of columns"))
    end

    global_module_matrix = vcat(matrix(M.sub), matrix(N.sub), matrix(M.quo))

    CD = matrix_kernel(global_module_matrix)

    C = CD[:,1:m[1]]
    D = CD[:,(m[1]+1):(m[1]+n[1])]
    new_gen = C*matrix(M.sub)
    SQ = SubQuo(free_module(M),new_gen, m_rel)

    M_hom = SubQuoHom(SQ,M,C)
    N_hom = SubQuoHom(SQ,N,D)
    register_morphism!(M_hom)
    register_morphism!(N_hom)

    return SQ,M_hom,N_hom
  end
  throw(ArgumentError("img(M.relations) != img(N.relations)"))
end

@doc Markdown.doc"""
  SubQuoElem{T}

An element of a subquotient module $A+B / B$ is given by a sparse row (mathematically speaking just 
a vector) $v$. The element is then just $u := \sum_i v[i]\cdot A[i]$ (where $A[i]$ are the generators of $A$).
The representative $u$ is also stored (along with the parent module where the 
element lives in).
"""
struct SubQuoElem{T} # this needs to be redone TODO
  coeffs::SRow{T}
  repres::FreeModuleElem{T}
  parent::SubQuo

  function SubQuoElem{R}(v::SRow{R}, SQ::SubQuo) where {R}
    @assert length(v) <= ngens(SQ.sub)
    if isempty(v)
      r = new{R}(v, zero(SQ.F), SQ)
      return r
    end
    r = new{R}(v, Base.sum([v[i]*SQ.sub[i] for i=1:ngens(SQ.sub)]), SQ)
    #r.coeffs = v
    #r.parent = SQ
    #r.repres = sum([v[i]*SQ.F[i] for i=1:ngens(SQ.F)]...)
    return r
  end

  function SubQuoElem{R}(a::FreeModuleElem{R}, SQ::SubQuo) where {R}
    @assert a.parent === SQ.F
    r = new{R}(coordinates(a,SQ), a, SQ)
    #r.parent = SQ
    #r.repres = a
    #r.v = coordinates(a, SQ)
    return r
  end
end

@doc Markdown.doc"""
  SubQuoElem(v::SRow{R}, SQ::SubQuo) where {R}

Return the element $\sum_i v[i]*SQ[i]$.
"""
SubQuoElem(v::SRow{R}, SQ::SubQuo) where {R} = SubQuoElem{R}(v, SQ)

@doc Markdown.doc"""
  SubQuoElem(a::FreeModuleElem{R}, SQ::SubQuo) where {R}

Construct an element $v ∈ SQ$ that is represented by $a$.
"""
SubQuoElem(a::FreeModuleElem{R}, SQ::SubQuo) where {R} = SubQuoElem{R}(a, SQ)

elem_type(::SubQuo{T}) where {T} = SubQuoElem{T}
parent_type(::SubQuoElem{T}) where {T} = SubQuo{T}
elem_type(::Type{SubQuo{T}}) where {T} = SubQuoElem{T}
parent_type(::Type{SubQuoElem{T}}) where {T} = SubQuo{T}

function in(v::SubQuoElem, M::ModuleFP)
  return parent(v) === M
end

@doc Markdown.doc"""
  getindex(v::SubQuoElem, i::Int)

Let $v ∈ M$ with $v = \sum_i a[i]*M[i]$. Return $a[i]$
"""
function getindex(v::SubQuoElem, i::Int)
  if isempty(coeffs(v))
    return zero(base_ring(v.parent))
  end
  return coeffs(v)[i]
end

@doc Markdown.doc"""
  coeffs(v::SubQuoElem)

Let $v ∈ M$. Return an `SRow` `a` such that $\sum_i a[i]*M[i] = v$.
"""
function coeffs(v::SubQuoElem)
  return v.coeffs
end

@doc Markdown.doc"""
  repres(v::SubQuoElem)

Return a free module element that is a representative of `v`.
"""
function repres(v::SubQuoElem)
  return v.repres
end

@doc Markdown.doc"""
  groebner_basis(F::ModuleGens)

Return a Gröbner basis of `F` as an object of type `ModuleGens`.
"""
function groebner_basis(F::ModuleGens)
  singular_assure(F)
  if singular_generators(F).isGB
    return F
  end
  return ModuleGens(F.F, Singular.std(singular_generators(F)))
end

function show(io::IO, b::SubQuoElem)
  print(io, b.repres)
end

@doc Markdown.doc"""
  parent(b::SubQuoElem)

Let $b ∈ M$. Return $M$.
"""
parent(b::SubQuoElem) = b.parent

@doc Markdown.doc"""
  (R::SubQuo)(a::FreeModuleElem; check::Bool = true)

Return `a` as an element of `R`. If `check == true` (default) it is checked that `a` represents
indeed an element of `R`.
"""
function (R::SubQuo)(a::FreeModuleElem; check::Bool = true)
  if check
    b = convert(R.sum.gens.SF, a)
    c = _reduce(b, singular_generators(std_basis(R.sum)))
    iszero(c) || error("not in the module")
  end
  return SubQuoElem(a, R)
end

@doc Markdown.doc"""
  (R::SubQuo)(a::SRow)

Return the subquotient element $\sum_i R[i]*a[i] ∈ R$.
"""
function (R::SubQuo)(a::SRow)
  return SubQuoElem(a, R)
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

Let $v ∈ G$ with $v$ the `i`th generator of $G$. Return `i`.
"""
function index_of_gen(v::SubQuoElem)
  @assert length(coeffs(v).pos) == 1
  @assert isone(coeffs(v).values[1])
  return coeffs(v).pos[1]
end

# function to check whether a free module element is in a particular free module
function check_parent(a::Union{FreeModuleElem,SubQuoElem}, b::Union{FreeModuleElem,SubQuoElem})
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
*(a::Int, b::SubQuoElem) = SubQuoElem(a*coeffs(b), b.parent)
*(a::Integer, b::SubQuoElem) = SubQuoElem(a*coeffs(b), b.parent)
*(a::fmpq, b::SubQuoElem) = SubQuoElem(a*coeffs(b), b.parent)
==(a::SubQuoElem, b::SubQuoElem) = iszero(a-b)

@doc Markdown.doc"""
  sub(F::FreeMod, O::Vector{<:FreeModuleElem}, task::Symbol = :none)

Return `S` as a submodule of `F`, where `S` is generated by `O`.
`S` is a represented as a subquotient module.
The elements of `O` must live in `F`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical injection morphism
$S → F$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function sub(F::FreeMod, O::Vector{<:FreeModuleElem}, task::Symbol = :none)
  s = SubQuo(F, O)
  if task == :none || task == :module
    return s
  else
    emb = hom(s, F, O)
    task == :store && register_morphism!(emb)
    task == :morphism && return emb
    return s, emb
  end
end

@doc Markdown.doc"""
  sub(F::FreeMod, O::Vector{<:SubQuoElem}, task::Symbol = :none)

Return `S` as a submodule of `F`, where `S` is generated by `O`.
The embedding module of the parent of the elements of `O` must be `F`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical injection morphism
$S → F$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function sub(F::FreeMod, O::Vector{<:SubQuoElem}, task::Symbol = :none)
  s = SubQuo(F, [x.repres for x = O])
  return sub(F, s, task)
  #=if task == :none
    return s
  else
    emb = hom(s, F, [x.repres for x = O])
  end=#
end

@doc Markdown.doc"""
  sub(F::FreeMod, s::SubQuo, task::Symbol = :none)

Return `s` as a submodule of `F`, that is the embedding free module of `s` must 
be `F` and `s` has no relations.
If `task` is set to `:with_morphism` or to `:both` return also the canonical injection morphism
$s → F$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function sub(F::FreeMod, s::SubQuo, task::Symbol = :none)
  @assert !isdefined(s, :quo)
  @assert s.F === F
  if task == :none || task == :module
    return s
  else
    emb = hom(s, F, [FreeModuleElem(x.repres.coords, F) for x in gens(s)])
    task == :store && register_morphism!(emb)
    task == :morphism && return emb 
    return s, emb
  end
end

@doc Markdown.doc"""
  sub(S::SubQuo, O::Vector{<:SubQuoElem}, task::Symbol = :none, check = true)

Compute a subquotient $T ≤ S$, where $T$ is generated by $O$. 
The elements of `O` must live in `S`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical injection morphism
$T → S$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
If `check` is set to `false` then it is not checked that the elements of `O` live in `S`.
"""
function sub(S::SubQuo, O::Vector{<:SubQuoElem}, task::Symbol = :none, check = true)
  if check
    @assert all(x -> x.parent === S, O)
  end
  t = SubQuo(S.F, [x.repres for x in O])
  if isdefined(S, :quo)
    t.quo = S.quo
    t.sum = sum(t.sub, t.quo)
  end
  if task == :none || task == :module
    return t
  else
    emb = hom(t, S, O)
    task == :store && register_morphism!(emb)
    task == :morphism && return emb 
    return t, emb
  end

  #=t = sub(S.F, O, task)
  if task != :none
    t,emb = t
  end
  if isdefined(S, :quo)
    s = quo(t, collect(S.quo.gens))
    if task == :none
      return s
    else
      emb2 = hom(s, S, [SubQuoElem(x.repres, S) for x in gens(s)])
    end
  else
    return t
  end=#
end

@doc Markdown.doc"""
  quo(F::FreeMod, O::Vector{<:FreeModuleElem}, task::Symbol = :none)

Compute $F / T$, where $T$ is generated by $O$.
The elements of `O` must live in `F`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$F → F/T$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function quo(F::FreeMod, O::Vector{<:FreeModuleElem}, task::Symbol = :none)
  S = SubQuo(F, basis(F))
  Q = SubQuo(S, O)

  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
  quo(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :none) where T

Compute $F / T$, where $T$ is generated by $O$.
The embedding free module of the parent of the elements of `O` must be `F`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$F → F/T$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function quo(F::FreeMod{T}, O::Vector{<:SubQuoElem{T}}, task::Symbol = :none) where T
  S = SubQuo(F, basis(F))
  Q = SubQuo{T}(S, [x.repres for x = O])
  
  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
  quo(F::SubQuo, O::Vector{<:FreeModuleElem}, task::Symbol = :none)

Compute $F / T$, where $T$ is generated by $O$.
The elements of `O` must be elements of the embedding free module of `S`.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$F → F/T$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function quo(F::SubQuo, O::Vector{<:FreeModuleElem}, task::Symbol = :none)
  if length(O) > 0
    @assert parent(O[1]) === F.F
  end
  if isdefined(F, :quo)
    #F.sub[Val(:S), 1]
    #[F.quo.gens[Val(:O), i] for i = 1:length(F.quo.gens.O)] 
    oscar_assure(F.quo.gens)
    s = Singular.Module(base_ring(F.quo.gens.SF), [convert(F.quo.gens.SF, x) for x = [O; oscar_generators(F.quo.gens)]]...)
    Q = SubQuo(F.F, singular_generators(F.sub.gens), s)
    return return_quo_wrt_task(F, Q, task)
  end
  Q = SubQuo(F, O)
  return return_quo_wrt_task(F, Q, task)
end

@doc Markdown.doc"""
  quo(S::SubQuo, O::Vector{<:SubQuoElem}, task::Symbol = :none)

Compute $S / T$, where $T$ is generated by $O$.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$S → S/T$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function quo(S::SubQuo, O::Vector{<:SubQuoElem}, task::Symbol = :none)
  return quo(S, [x.repres for x = O], task)
end

@doc Markdown.doc"""
  quo(S::SubQuo, T::SubQuo, task::Symbol = :none)

Compute $S / T$.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$S → S/T$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function quo(S::SubQuo, T::SubQuo, task::Symbol = :none)
#  @assert !isdefined(T, :quo)
  # TODO @assert S.quo == T.quo or too expensive?
  Q = SubQuo(S, oscar_generators(T.sum.gens))
  return return_quo_wrt_task(S, Q, task)
end

@doc Markdown.doc"""
  quo(F::FreeMod{R}, T::SubQuo{R}, task::Symbol = :none) where R

Compute $F / T$.
If `task` is set to `:with_morphism` or to `:both` return also the canonical projection morphism
$F → F/T$.
If `task` is set to `:store` the morphism is also cached.
If `task` is set to `:morphism` return only the morphism.
"""
function quo(F::FreeMod{R}, T::SubQuo{R}, task::Symbol = :none) where R
  @assert !isdefined(T, :quo)
  return quo(F, gens(T), task)
end

@doc Markdown.doc"""
  return_quo_wrt_task(M::ModuleFP, Q::ModuleFP, task)

This helper function returns the module `Q = M / N` for some `N` 
along with the canonical projection morphism $M → Q$ according to the 
given `task`.
"""
function return_quo_wrt_task(M::ModuleFP, Q::ModuleFP, task)
  if task == :none || task == :module
    return Q
  else
    pro = hom(M, Q, gens(Q))
    task == :store && register_morphism!(pro)
    task == :morphism && return pro
    return Q, pro
  end
end

@doc Markdown.doc"""
  syzygy_module(F::ModuleGens; sub = FreeMod(base_ring(F.F), length(oscar_generators(F))))
"""
function syzygy_module(F::ModuleGens; sub = FreeMod(base_ring(F.F), length(oscar_generators(F))))
  F[Val(:S), 1] #to force the existence of F.S
  s = Singular.syz(singular_generators(F)) # TODO syz is sometimes too slow, example [8*x^2*y^2*z^2 + 13*x*y*z^2 12*x^2 + 7*y^2*z; 13*x*y^2 + 12*y*z^2 4*x^2*y^2*z + 8*x*y*z; 9*x*y^2 + 4*z 12*x^2*y*z^2 + 9*x*y^2*z]
  #=if sub !== 0
    G = sub
  else
    G = FreeMod(base_ring(F.F), length(oscar_generators(F)))
  end
  return SubQuo(G, s)=#
  return SubQuo(sub, s)
end

@doc Markdown.doc"""
  gens(F::SubQuo{T}) where T

Return the generators of `F`.
"""
function gens(F::SubQuo{T}) where T
  return [gen(F,i) for i=1:ngens(F)]
end

@doc Markdown.doc"""
  gen(F::SubQuo{T}, i::Int) where T

Return the `i`th generator of `F`.
"""
function gen(F::SubQuo{T}, i::Int) where T
  R = base_ring(F)
  v = sparse_row(R)
  v.pos = [i]
  v.values = [R(1)]
  return SubQuoElem{T}(v, F)
end

@doc Markdown.doc"""
  ngens(F::SubQuo)

Return the number of generators of `F`.
"""
ngens(F::SubQuo) = ngens(F.sub)

@doc Markdown.doc"""
  base_ring(SQ::SubQuo)

Let `SQ` be an `R`-module. Return `R`.
"""
base_ring(SQ::SubQuo) = base_ring(SQ.F)

@doc Markdown.doc"""
  zero(SQ::SubQuo)

Return the zero element in `SQ`.
"""
zero(SQ::SubQuo) = SubQuoElem(zero(SQ.F), SQ)

@doc Markdown.doc"""
  Base.iszero(F::SubQuo)

Check if `F` is the zero module.
"""
function Base.iszero(F::SubQuo)
  return all(iszero, gens(F))
end

@doc Markdown.doc"""
  Base.getindex(F::SubQuo, i::Int)

Return the `i`th generator of `F`.
"""
function Base.getindex(F::SubQuo, i::Int)
  i == 0 && return zero(F)
  return gen(F, i)
end

function Base.iterate(F::ModuleGens, i::Int = 1)
  if i>length(F)
    return nothing
  else
    return F[i], i+1
  end
end
Base.eltype(::ModuleGens{T}) where {T} = FreeModuleElem{T} 

#??? A scalar product....
function *(a::FreeModuleElem, b::Vector{FreeModuleElem})
  @assert dim(parent(a)) == length(b)
  s = zero(parent(a))
  for (p,v) = a.coords
    s += v*b[p]
  end
  return s
end

@doc Markdown.doc"""
  presentation(SQ::SubQuo)

Compute a free resolution of `SQ`. 
"""
function presentation(SQ::SubQuo)
  #A+B/B is generated by A and B
  #the relations are A meet B? written wrt to A
  s = syzygy_module(SQ.sum.gens)
  #TODO: wait for Hans to release Modulo(A, B) that does exactly this
  c = collect(s.sub.gens)
  R = base_ring(SQ)
  F = FreeMod(R, ngens(SQ.sub))
  q = elem_type(F)[]

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
    #=if length(b) == 0 #TODO why was this here
      continue
    end=#
    push!(q, FreeModuleElem(b, F))
  end
  #want R^a -> R^b -> SQ -> 0
  #TODO sort decoration and fix maps, same decoration should be bundled (to match pretty printing)
  G = FreeMod(R, ngens(s.sub))
  h_G_F = hom(G, F, q)
  h_F_SQ = hom(F, SQ, gens(SQ)) # DO NOT CHANGE THIS LINE, see present and preimage

  Z = FreeMod(F.R, 0)
  Hecke.set_special(Z, :name => "Zero")
  h_SQ_Z = hom(SQ, Z, [zero(Z) for i=1:ngens(SQ)])
  return Hecke.ChainComplex(Oscar.ModuleFP, Oscar.ModuleMap[h_G_F, h_F_SQ, h_SQ_Z], check = false)
end

@doc Markdown.doc"""
  presentation(F::FreeMod)

Return a free resolution of $F$.
"""
function presentation(F::FreeMod)
  Z = FreeMod(F.R, 0)
  Hecke.set_special(Z, :name => "Zero")
  return Hecke.ChainComplex(ModuleFP, ModuleMap[hom(Z, F, FreeModuleElem[]), hom(F, F, gens(F)), hom(F, Z, [zero(Z) for i=1:ngens(F)])], check = false)
end

@doc Markdown.doc"""
  present_as_cokernel(SQ::SubQuo, task::Symbol = :none)

Return a subquotient $M = R^n / im(φ) $, i.e. $M = \text{coker}(φ)$, such that
$M ≅ SQ$.
If `task` is set to `:with_morphism` or to `:both` then return also an isomorphism $M → SQ$. Calling `inv()`
on this isomorphism is cheap.
If `task` is set to `:store` then the isomorphism is cached.
If `task` is set to `:morphism` then return only the isomorphism.
"""
function present_as_cokernel(SQ::SubQuo, task::Symbol = :none)
  chainComplex = presentation(SQ)
  R_b = obj(chainComplex, 1)
  f = map(chainComplex, 1)
  g = map(chainComplex, 2)
  presentation_module = quo(R_b, image(f)[1])

  if task == :none
    return presentation_module
  end
  
  # The isomorphism is just the identity matrix
  isomorphism = hom(presentation_module, SQ, [g(x) for x in gens(R_b)])
  inverse_isomorphism = hom(SQ, presentation_module, [presentation_module[i] for i=1:ngens(SQ)])
  isomorphism.inverse_isomorphism = inverse_isomorphism

  if task == :store
    register_morphism!(isomorphism)
    register_morphism!(inverse_isomorphism)
  end
  task == :morphism && return isomorphism
  
  return presentation_module, isomorphism
end

@doc Markdown.doc"""
  is_equal_with_morphism(M::SubQuo{T}, N::SubQuo{T}, task::Symbol = :none) where {T}
If $M = N$ (mathematically, but with (possibly) different generating systems), return $\phi : M → N$ 
which is mathematically the identity. 
If `task == :inverse` also the inverse map is computed and cached (in the morphism).
If `task == :store` the inverse map is also cached in `M` and `N`.
"""
function is_equal_with_morphism(M::SubQuo{T}, N::SubQuo{T}, task::Symbol = :none) where {T}
  @assert M == N

  M_to_N = hom(M, N, [SubQuoElem(coordinates(m.repres, N), N) for m in gens(M)])

  if task == :store || task == :inverse
    N_to_M = hom(N, M, [SubQuoElem(coordinates(n.repres, M), M) for n in gens(N)])
    M_to_N.inverse_isomorphism = N_to_M
    N_to_M.inverse_isomorphism = M_to_N

    if task == :store
      register_morphism!(M_to_N) 
      register_morphism!(N_to_M)
    end
  end
  
  return M_to_N
end

mutable struct SubQuoHom{T1, T2} <: ModuleMap{T1, T2}
  matrix::MatElem
  header::Hecke.MapHeader
  im::Vector
  inverse_isomorphism::ModuleMap

  function SubQuoHom{T1,T2}(D::SubQuo, C::ModuleFP, im::Vector) where {T1,T2}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    #r = new{SubQuo, typeof(C)}()
    r = new{T1, T2}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    if C isa FreeMod
      r.im = Vector{FreeModuleElem}(im)
    elseif C isa SubQuo
      r.im = Vector{SubQuoElem}(im)
    else
      r.im = im
    end

    return r
  end
end

@doc Markdown.doc"""
  SubQuoHom(D::SubQuo, C::ModuleFP, im::Vector)

Return the morphism $D → C$ for a subquotient $D$ where `D[i]` is mapped to `im[i]`.
In particular, `length(im) == ngens(D)` must hold.
"""
SubQuoHom(D::SubQuo, C::ModuleFP, im::Vector) = SubQuoHom{typeof(D), typeof(C)}(D, C, im)

@doc Markdown.doc"""
  SubQuoHom(D::SubQuo, C::ModuleFP, mat::MatElem)

Return the morphism $D → C$ corresponding to the given matrix, where $D$ is a subquotient.
`mat` must have `ngens(D)` many rows and `ngens(C)` many columns.
"""
function SubQuoHom(D::SubQuo, C::ModuleFP, mat::MatElem)
  @assert nrows(mat) == ngens(D)
  @assert ncols(mat) == ngens(C)
  if typeof(C) <: FreeMod
    hom = SubQuoHom(D, C, [FreeModuleElem(sparse_row(mat[i,:]), C) for i=1:ngens(D)])
    return hom
  else
    hom = SubQuoHom(D, C, [SubQuoElem(sparse_row(mat[i,:]), C) for i=1:ngens(D)])
    return hom
  end
end

#=function Base.getproperty(f::SubQuoHom, s::Symbol)
  if s == :matrix
    error("deprecated")
    if !isdefined(f, s)
      D = domain(f)
      C = codomain(f)
      R = base_ring(D)
      matrix = zero_matrix(R, ngens(D), ngens(C))
      for i=1:ngens(D), j=1:ngens(C)
        matrix[i,j] = f.im[i][j]
      end
      setfield!(f, s, matrix)
    end
  end
  return getfield(f,s)
end=#

@doc Markdown.doc"""
  matrix(f::SubQuoHom)

Return the matrix corresponding to `f`. 
The matrix is cached.
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
  hom_tensor(G::ModuleFP, H::ModuleFP, A::Vector{ <: ModuleMap})

Let $G = G_1 ⊗ \cdot ⊗ G_n$, $H = H_1 ⊗ \cdot H_n$ and `A` an array of morphisms
$φ_1, \cdot, φ_n$ with $φ_i : G_i → H_i$. Return then $φ_1 ⊗ \cdot ⊗ φ_n$.
"""
function hom_tensor(G::ModuleFP, H::ModuleFP, A::Vector{ <: ModuleMap})
  tG = get_special(G, :tensor_product)
  tG === nothing && error("both modules must be tensor products")
  tH = get_special(H, :tensor_product)
  tH === nothing && error("both modules must be tensor products")
  @assert length(tG) == length(tH) == length(A)
  @assert all(i-> domain(A[i]) === tG[i] && codomain(A[i]) === tH[i], 1:length(A))
  #gens of G are G[i][j] tensor G[h][l] for i != h and all j, l
  #such a pure tensor is mapped to A[i](G[i][j]) tensor A[h](G[j][l])
  #thus need the pure map - and re-create the careful ordering of the generators as in the 
  # constructor
  #store the maps? and possibly more data, like the ordeing
  decompose_G = get_special(G, :tensor_generator_decompose_function)
  pure_H = get_special(H, :tensor_pure_function)
  function map_gen(g) # Is there something that generalizes FreeModuleElem and SubQuoElem?
    g_decomposed = decompose_G(g)
    image_as_tuple = Tuple(f(x) for (f,x) in zip(A,g_decomposed))
    res = pure_H(image_as_tuple)
    return res
  end
  #error("not done yet")
  return hom(G, H, map(map_gen, gens(G)))
end

@doc Markdown.doc"""
  hom_prod_prod(G::ModuleFP, H::ModuleFP, A::Matrix{<:ModuleMap})

Let $G = G_1 ⊕ \cdot ⊕ G_n$, $H = H_1 ⊕ \cdot H_m$ and $A = (φ_{ij})_{1≤i≤n,1≤j≤m}$ with 
$φ_{ij} : G_i → H_j$ morphisms. Return the morphism $φ : G → H$ corresponding to $A$.
"""
function hom_prod_prod(G::ModuleFP, H::ModuleFP, A::Matrix{<:ModuleMap})
  tG = get_special(G, :direct_product)
  tG === nothing && error("both modules must be direct products")
  tH = get_special(H, :direct_product)
  tH === nothing && error("both modules must be direct products")
  @assert length(tG) == size(A, 1) && length(tH) == size(A, 2)
  @assert all(ij -> domain(A[ij[1],ij[2]]) === tG[ij[1]] && codomain(A[ij[1],ij[2]]) === tH[ij[2]], Base.Iterators.ProductIterator((1:size(A, 1), 1:size(A, 2))))
  #need the canonical maps..., maybe store them as well?
  return hom(G,H,[Base.sum([Hecke.canonical_injection(H,j)(Base.sum([A[i,j](Hecke.canonical_projection(G,i)(g)) for i=1:length(tG)])) for j=1:length(tH)]) for g in gens(G)])
  #error("not done yet")
end
# hom(prod -> X), hom(x -> prod)
# if too much time: improve the hom(A, B) in case of A and/or B are products - or maybe not...
# tensor and hom functors for chain complex
# dual: ambig: hom(M, R) or hom(M, Q(R))?

@doc Markdown.doc"""
  coordinates(a::FreeModuleElem, SQ::SubQuo)::Union{Nothing,Oscar.SRow}
Compute a sparse row `r` such that `a` is a representative of `SubQuoElem(r, SQ)`.
If no such `r` exists, `nothing` is returned.
"""
function coordinates(a::FreeModuleElem, SQ::SubQuo)::Union{Nothing,Oscar.SRow}
  if iszero(a)
    return sparse_row(base_ring(parent(a)))
  end
  #[SQ.sub[Val(:O), i] for i = 1:length(SQ.sub.O)]
  oscar_assure(SQ.sub.gens)
  if isdefined(SQ, :quo)
    #[SQ.quo[Val(:O), i] for i = 1:length(SQ.quo.O)]
    oscar_assure(SQ.quo.gens)
    generators = sum(SQ.sub, SQ.quo)
  else
    generators = SQ.sub
  end
  S = singular_generators(generators.gens)
  #S = Singular.smodule{elem_type(base_ring(SQ.sub.gens.SF))}(base_ring(SQ.sub.gens.SF), [convert(SQ.sub.gens.SF, x) for x = generators]...)
  b = ModuleGens([a], SQ.sum.gens.SF)
  singular_assure(b)
  s, r = Singular.lift(S, singular_generators(b))
  if Singular.ngens(s) == 0 || iszero(s[1])
    #error("elem not in module")
    return nothing
  end
  Rx = base_ring(SQ)
  return sparse_row(Rx, s[1], 1:ngens(SQ))
end

@doc Markdown.doc"""
  represents_element(a::FreeModuleElem, SQ::SubQuo)

Check if `a` represents an element `SQ`.
"""
function represents_element(a::FreeModuleElem, SQ::SubQuo)
  return !isnothing(coordinates(a,SQ))
end

hom(D::SubQuo, C::ModuleFP, A::Vector{<:Any}) = SubQuoHom(D, C, A)

@doc Markdown.doc"""
  image(f::SubQuoHom, a::SubQuoElem)
Return $f(a)$.
"""
function image(f::SubQuoHom, a::SubQuoElem)
  # TODO matrix vector multiplication
  @assert a.parent === domain(f)
  i = zero(codomain(f))
  D = domain(f)
  b = coeffs(a)
  #b = coordinates(a.repres, D)
  for (p,v) = b
    i += v*f.im[p]
  end
  return i
end

@doc Markdown.doc"""
  image(f::SubQuoHom, a::FreeModuleElem)
Return $f(a)$. `a` must represent an element in the domain of `f`.
"""
function image(f::SubQuoHom, a::FreeModuleElem)
  return image(f, SubQuoElem(a, domain(f)))
  #=i = zero(codomain(f))
  D = domain(f)
  b = coordinates(a, D)
  for (p,v) = b
    i += v*f.im[p]
  end
  return i=#
end

@doc Markdown.doc"""
  preimage(f::SubQuoHom, a::Union{SubQuoElem,FreeModuleElem})
Compute a preimage of `a` under `f`.
"""
function preimage(f::SubQuoHom, a::Union{SubQuoElem,FreeModuleElem})
  @assert parent(a) === codomain(f)
  D = domain(f)
  i = zero(D)
  b = coordinates(typeof(a) <: FreeModuleElem ? a : a.repres, image(f)[1])
  for (p,v) = b
    i += v*gen(D, p)
  end
  return i
end

(f::SubQuoHom)(a::FreeModuleElem) = image(f, a)
(f::SubQuoHom)(a::SubQuoElem) = image(f, a)

@doc Markdown.doc"""
  iszero(a::SubQuoElem)

Return if the subquotient element `a` is zero.
"""
function iszero(a::SubQuoElem)
  C = parent(a)
  if !isdefined(C, :quo)
    return iszero(a.repres)
  end
  x = _reduce(convert(C.quo.gens.SF, a.repres), singular_generators(std_basis(C.quo)))
  return iszero(x)
end

@doc Markdown.doc"""
  hom(F::FreeMod, G::FreeMod)
Return a subquotient $S$ such that $Hom(F,G) \cong S$ along with a function 
that converts elements from $S$ into morphisms $F → G$.
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
  function im(x::FreeModuleElem)
    return hom(F, G, [FreeModuleElem(x.coords[R, (i-1)*m+1:i*m], G) for i=1:n])
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
    return FreeModuleElem(s, GH)
  end
  to_hom_map = Hecke.MapFromFunc(im, pre, GH, X)
  Hecke.set_special(GH, :show => Hecke.show_hom, :hom => (F, G), :module_to_hom_map => to_hom_map)
  return GH, to_hom_map
end

@doc Markdown.doc"""
  kernel(h::FreeModuleHom)

Compute the kernel $K$ of `h` along with the inclusion morphism 
of $K$ into the domain of `h`.
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
    k = sub(G, [FreeModuleElem(x.coords[R,1:dim(G)], G) for x = s])
  else
    #the syzygie_module creates a new free module to work in
    k = sub(G, [FreeModuleElem(x.coords, G) for x = collect(k.sub.gens)])
  end
  @assert k.F === G
  c = collect(k.sub.gens)
  return k, hom(k, parent(c[1]), c)
end

@doc Markdown.doc"""
  image(h::FreeModuleHom)

Compute the image of `h`. Return also the inclusion map into the codomain of `h`.
"""
function image(h::FreeModuleHom)
  si = [x for x = map(h, basis(domain(h))) if !iszero(x)]
  s = sub(codomain(h), si)
  return s, hom(s, codomain(h), si)
end

@doc Markdown.doc"""
  image(h::SubQuoHom)

Compute the image of `h`. Return also the inclusion map into the codomain of `h`.
"""
function image(h::SubQuoHom)
  #=if typeof(codomain(h)) <: FreeMod
    image_vector::Array{FreeModuleElem} = h.im
    s = sub(codomain(h), image_vector)
    return s, hom(s, codomain(h), image_vector)
  else
    image_vector2::Array{SubQuoElem} = h.im
    s = sub(codomain(h), image_vector2)
    return s, hom(s, codomain(h), image_vector2)
  end=#
  s = sub(codomain(h), h.im)
  return s, hom(s, codomain(h), h.im)
end

@doc Markdown.doc"""
  kernel(h::SubQuoHom)

Compute the kernel $K$ of `h` along with the inclusion morphism 
of $K$ into the domain of `h`.
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

Compute a free resolution of the free module `F`.
"""
function free_resolution(F::FreeMod)
  return presentation(F)
end

@doc Markdown.doc"""
  free_resolution(S::SubQuo, limit::Int = -1)
Compute a free resolution of `S`. If `limit != -1` the free resolution
is only computed up to the `limit`-th free module.
"""
function free_resolution(S::SubQuo, limit::Int = -1)
  p = presentation(S)
  mp = [map(p, j) for j=1:length(p)]
  while true
    k, mk = kernel(mp[1])
    nz = findall(x->!iszero(x), gens(k))
    if length(nz) == 0 
      Z = FreeMod(base_ring(S), 0)
      Hecke.set_special(Z, :name => "Zero")
      h = hom(Z, domain(mp[1]), FreeModuleElem[])
      insert!(mp, 1, h)
      break
    elseif limit != -1 && length(mp) > limit
      break
    end
    F = FreeMod(base_ring(S), length(nz))
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
  hom(M::ModuleFP, N::ModuleFP)
Return a subquotient $S$ such that $Hom(M,N) \cong S$ along with a function 
that converts elements from $S$ into morphisms $M → N$.
"""
function hom(M::ModuleFP, N::ModuleFP)
  p1 = presentation(M)
  p2 = presentation(N)
  k, mk = kernel(map(p2, 1))
  #Janko: have R^t1 -- g1 = map(p2, 0) -> R^t0 -> G
  #kernel g1: k -> R^t1
  #source: Janko's CA script: https://www.mathematik.uni-kl.de/~boehm/lehre/17_CA/ca.pdf
  F = FreeMod(base_ring(M), ngens(k))
  g2 = hom(F, codomain(mk), collect(k.sub.gens)) #not clean - but maps not (yet) working
  #step 2
  H_s0_t0, mH_s0_t0 = hom(domain(map(p1, 2)), domain(map(p2, 2)))
  H_s1_t1, mH_s1_t1 = hom(domain(map(p1, 1)), domain(map(p2, 1)))
  D, pro, emb = direct_product(H_s0_t0, H_s1_t1, task = :both)

  H_s1_t0, mH_s1_t0 = hom(domain(map(p1, 1)), domain(map(p2, 2)))

  delta = hom(D, H_s1_t0, [preimage(mH_s1_t0, map(p1, 1)*mH_s0_t0(pro[1](g))-mH_s1_t1(pro[2](g))*map(p2, 1)) for g = gens(D)])

  H_s0_t1, mH_s0_t1 = hom(domain(map(p1, 2)), domain(map(p2, 1)))
  H_s1_t2, mH_s1_t2 = hom(domain(map(p1, 1)), F)

  E, pr = direct_product(H_s0_t1, H_s1_t2, task = :prod)

  rho = hom(E, D, [emb[1](preimage(mH_s0_t0, mH_s0_t1(pr[1](g))*map(p2, 1))) + 
                   emb[2](preimage(mH_s1_t1, map(p1, 1)*mH_s0_t1(pr[1](g)) - mH_s1_t2(pr[2](g))*g2)) for g = gens(E)])
  #need quo(kern(delta), image(rho))                 
 
  kDelta = kernel(delta)

  #psi = kDelta[2]*pro[1]
  #psi = hom(kDelta[1], H_s0_t0, [psi(g) for g = gens(kDelta[1])])

  H = quo(sub(D, kDelta[1]), image(rho)[1])

  #x in ker delta: mH_s0_t0(pro[1](x)) should be a hom from M to N
  function im(x::SubQuoElem)
    @assert parent(x) === H
    #return SubQuoHom(M, N, mH_s0_t0(pro[1](x.repres)).matrix)
    return hom(M, N, [map(p2, 2)(mH_s0_t0(pro[1](x.repres))(preimage(map(p1, 2), g))) for g = gens(M)])
  end

  function pre(f::ModuleMap)
    @assert domain(f) === M
    @assert codomain(f) === N
    Rs0 = domain(map(p1, 2))
    Rt0 = domain(map(p2, 2))
    g = hom(Rs0, Rt0, [preimage(map(p2, 2), f(map(p1, 2)(g))) for g = gens(Rs0)])

    #return H(preimage(psi, (preimage(mH_s0_t0, g))).repres)
    return SubQuoElem(preimage(kDelta[2], emb[1](preimage(mH_s0_t0, g))).repres, H)
    #return SubQuoElem(emb[1](preimage(mH_s0_t0, g)), H) #???
  end
  to_hom_map = MapFromFunc(im, pre, H, Hecke.MapParent(M, N, "homomorphisms"))
  Hecke.set_special(H, :show => Hecke.show_hom, :hom => (M, N), :module_to_hom_map => to_hom_map)
  return H, to_hom_map
end

@doc Markdown.doc"""
  homomorphism(f::Union{SubQuoElem,FreeModuleElem})
If `f` is an element in a module created via `hom(M,N)` for some `M` and `N`, 
return the morphism $\phi : M → N$ that corresponds to `f`.
"""
function homomorphism(f::Union{SubQuoElem,FreeModuleElem})
  H = f.parent
  to_hom_map = get_special(H, :module_to_hom_map)
  to_hom_map === nothing && error("element doesn't live in a hom module")  
  return to_hom_map(f)
end

@doc Markdown.doc"""
  homomorphism_to_module_elem(H::ModuleFP, phi::ModuleMap)
Let `H` is created via `hom(M,N)` for some `M` and `N`. Return 
the element in `H` corresponding to `phi`.
"""
function homomorphism_to_module_elem(H::ModuleFP, phi::ModuleMap)
  to_hom_map = get_special(H, :module_to_hom_map)
  to_hom_map === nothing && error("module must be a hom module")
  map_to_hom = to_hom_map.g
  return map_to_hom(phi)
end

#TODO
#  replace the +/- for the homs by proper constructors for homs and direct sums
#  relshp to store the maps elsewhere

@doc Markdown.doc"""
  *(h::ModuleMap, g::ModuleMap)

Return the composition $g ∘ h$.
"""
function *(h::ModuleMap, g::ModuleMap)
  @assert codomain(h) === domain(g)
  return hom(domain(h), codomain(g), [g(h(x)) for x = gens(domain(h))])
end
-(h::FreeModuleHom, g::FreeModuleHom) = hom(domain(h), codomain(h), [h(x) - g(x) for x = gens(domain(h))])
+(h::FreeModuleHom, g::FreeModuleHom) = hom(domain(h), codomain(h), [h(x) + g(x) for x = gens(domain(h))])


@doc Markdown.doc"""
  restrict_codomain(H::ModuleMap, M::SubQuo)
Return, if possible, a homomorphism, which is mathematically identical to `H`,
but has codomain `M`. `M` has to be a submodule of the codomain of `H`.
"""
function restrict_codomain(H::ModuleMap, M::SubQuo)
  @assert typeof(codomain(H)) <: SubQuo
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
  i = sub(domain(H), map(m -> SubQuoElem(repres(m), domain(H)), gens(M)), :store, false)[2]
  return i*H
end

@doc Markdown.doc"""
  Base.inv(H::ModuleMap)
Compute $H^{-1}$. `H` must be bijective.
"""
function Base.inv(H::ModuleMap)
  if isdefined(H, :inverse_isomorphism)
    return H.inverse_isomorphism
  end
  @assert isbijective(H)
  N = domain(H)
  M = codomain(H)

  Hinv = hom(M,N, [preimage(H,m) for m in gens(M)])
  Hinv.inverse_isomorphism = H
  H.inverse_isomorphism = Hinv

  return Hinv
end

#=@doc Markdown.doc"""
  pseudo_inv(f::ModuleMap)
Compute $g : im(f) -> domain(f)$ such that $f ∘ g = id$.
"""
function pseudo_inv(f::ModuleMap)
  if isdefined(f, :inverse_isomorphism)
    return f.inverse_isomorphism
  end
  image_f,emb = image(f)
  f_pseudo_inv = hom(image_f, domain(f), [preimage(f,emb(g)) for g in gens(image_f)])

  return f_pseudo_inv
end=#

##################################################
# direct product
##################################################
@doc Markdown.doc"""
  direct_product(F::FreeMod{T}...; task::Symbol = :sum)

Given free modules $F_i$ compute the direct product $P := F_1\oplus \cdots \oplus F_n$.
If `task` is set to ":prod", an array of maps $\phi_1, \cdot, \phi_n$ is returned such
that $\phi_i$ is the canonical projection $P → F_i$.
If `task` is set to ":sum", an array of maps $\psi_1, \cdot, \psi_n$ is returned such 
that $\psi_i$ is the canonical injection $F_i → P$.
If `task` is set to ":both", both, the array of projections and the array of injections,
are returned (with projections first).
"""
function direct_product(F::FreeMod{T}...; task::Symbol = :sum) where {T}
  R = base_ring(F[1])
  G = FreeMod(R, Base.sum([rank(f) for f = F]))
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
  Hecke.set_special(G, :show => Hecke.show_direct_product, :direct_product => F)
  emb = []
  pro = []
  projection_dictionary = IdDict{Int,ModuleMap}()
  injection_dictionary = IdDict{Int,ModuleMap}()
  Hecke.set_special(G, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  i = 0
  for f = F
    if task in [:sum, :both]
      push!(emb, hom(f, G, [gen(G, j+i) for j=1:ngens(f)]))
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
  direct_product(G::ModuleFP...; task::Symbol = :none)

Given modules $G_i$ compute the direct product $P := G_1\oplus \cdots \oplus G_n$.
If `task` is set to ":prod", an array of maps $\phi_1, \cdot, \phi_n$ is returned such
that $\phi_i$ is the canonical projection $P → G_i$.
If `task` is set to ":sum", an array of maps $\psi_1, \cdot, \psi_n$ is returned such 
that $\psi_i$ is the canonical injection $G_i → P$.
If `task` is set to ":both", both, the array of projections and the array of injections,
are returned (with projections first).
"""
function direct_product(G::ModuleFP...; task::Symbol = :none)
  F, pro, mF = direct_product([free_module(x) for x = G]..., task = :both)
  s, emb_sF = sub(F, vcat([[mF[i](y) for y = gens(G[i], free_module(G[i]))] for i=1:length(G)]...), :both)
  q = vcat([[mF[i](y) for y = rels(G[i])] for i=1:length(G)]...)
  pro_quo = nothing
  if length(q) != 0
    s, pro_quo = quo(s, q, :both)
  end
  Hecke.set_special(s, :show => Hecke.show_direct_product, :direct_product => G)
  projection_dictionary = IdDict{Int,ModuleMap}()
  injection_dictionary = IdDict{Int,ModuleMap}()
  Hecke.set_special(s, :projection_morphisms => projection_dictionary, :injection_morphisms => injection_dictionary)
  if task == :none
    return s
  end
  if task == :prod || task != :sum
    if pro_quo === nothing
      for i=1:length(pro)
        pro[i] = hom(s, G[i], [G[i](pro[i](emb_sF(gen))) for gen in gens(s)]) # TODO distinction between pro on the left and pro on the right side!
        projection_dictionary[i] = pro[i]
      end
    else
      for i=1:length(pro)
        pro[i] = hom(s, G[i], [G[i](pro[i](emb_sF(preimage(pro_quo,gen)))) for gen in gens(s)])
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
        mF[i] = hom(G[i], s, [preimage(emb_sF, mF[i](typeof(G) <: FreeMod ? g : g.repres)) for g in gens(G[i])])
        injection_dictionary[i] = mF[i]
      end
    else
      for i=1:length(mF)
        mF[i] = hom(G[i], s, [pro_quo(preimage(emb_sF, mF[i](typeof(G) <: FreeMod ? g : g.repres))) for g in gens(G[i])])
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
⊕(M::ModuleFP...) = direct_product(M..., task = :none)


@doc Markdown.doc"""
  Hecke.canonical_injection(G::ModuleFP, i::Int)

Return the canonical injection $G_i → G$ where $G = G_1 ⊕ \cdot ⊕ G_n$.
"""
function Hecke.canonical_injection(G::ModuleFP, i::Int)
  H = Hecke.get_special(G, :direct_product)
  if H === nothing
    error("module not a direct product")
  end
  injection_dictionary = Hecke.get_special(G, :injection_morphisms)
  if haskey(injection_dictionary, i)
    return injection_dictionary[i]
  end
  0<i<= length(H) || error("index out of bound")
  j = i == 1 ? 0 : sum(ngens(H[l]) for l=1:i-1) -1
  emb = hom(H[i], G, [G[l+j] for l = 1:ngens(H[i])])
  injection_dictionary[i] = emb
  return emb
end

@doc Markdown.doc"""
  Hecke.canonical_projection(G::ModuleFP, i::Int)

Return the canonical projection $G → G_i$ where $G = G_1 ⊕ \cdot ⊕ G_n$.
"""
function Hecke.canonical_projection(G::ModuleFP, i::Int)
  H = Hecke.get_special(G, :direct_product)
  if H === nothing
    error("module not a direct product")
  end
  projection_dictionary = Hecke.get_special(G, :projection_morphisms)
  if haskey(projection_dictionary, i)
    return projection_dictionary[i]
  end
  0<i<= length(H) || error("index out of bound")
  j = i == 1 ? 0 : sum(ngens(H[l]) for l=1:i-1) 
  pro = hom(G, H[i], vcat([zero(H[i]) for l=1:j], gens(H[i]), [zero(H[i]) for l=1+j+ngens(H[i]):ngens(G)]))
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
  Hecke.set_special(F, :show => Hecke.show_tensor_product, :tensor_product => G)

  function pure(g::FreeModuleElem...)
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
    return FreeModuleElem(sparse_row(F.R, indices, zz), F)
  end
  function pure(T::Tuple)
    return pure(T...)
  end
  function inv_pure(e::FreeModuleElem)
    if length(e.coords.pos) == 0
      return Tuple(zero(g) for g = G)
    end
    @assert length(e.coords.pos) == 1
    @assert isone(e.coords.values[1])
    return Tuple(gen(G[i], t[e.coords.pos[1]][i]) for i = 1:length(G))
  end

  Hecke.set_special(F, :tensor_pure_function => pure, :tensor_generator_decompose_function => inv_pure)

  if task == :none
    return F
  end

  return F, MapFromFunc(pure, inv_pure, Hecke.TupleParent(Tuple([g[0] for g = G])), F)
end

⊗(G::ModuleFP...) = tensor_product(G..., task = :none)

@doc Markdown.doc"""
  free_module(F::FreeMod)

Just return `F`. This function exists only for compatiblity reasons.
"""
function free_module(F::FreeMod)
  return F
end

@doc Markdown.doc"""
  free_module(F::SubQuo)

Return the embedding free module of `F`.
"""
function free_module(F::SubQuo)
  return F.F
end

@doc Markdown.doc"""
  gens(F::FreeMod, G::FreeMod)

Return the generators of `F` as elements of `G` where `G === F`. This function 
exists only for compatiblity reasons.
"""
function gens(F::FreeMod, G::FreeMod)
  @assert F === G
  return gens(F)
end

@doc Markdown.doc"""
  gens(F::SubQuo, G::FreeMod)

Return the generators of `F` as elements of `G` where `G` is the embedding free module.
"""
function gens(F::SubQuo, G::FreeMod)
  @assert F.F === G
  return [FreeModuleElem(x.repres.coords, G) for x = gens(F)]
end
rels(F::FreeMod) = elem_type(F)[]
rels(F::SubQuo) = isdefined(F, :quo) ? collect(F.quo.gens) : elem_type(F.F)[]

@doc Markdown.doc"""
    tensor_product(G::ModuleFP...; task::Symbol = :none)

Given modules $G_i$ compute the tensor product $G_1\otimes \cdots \otimes G_n$.
If `task` is set to ":map", a map $\phi$ is returned that
maps tuples in $G_1 \times \cdots \times G_n$ to pure tensors
$g_1 \otimes \cdots \otimes g_n$. The map admits a preimage as well.
"""
function tensor_product(G::ModuleFP...; task::Symbol = :none)
  F, mF = tensor_product([free_module(x) for x = G]..., task = :map)
  # We want to store a dict where the keys are tuples of indices and the values
  # are the corresponding pure vectors (i.e. a tuple (2,1,5) represents the 
  # 2nd, 1st and 5th generator of the 1st, 2nd and 3rd module, which we are 
  # tensoring. The corresponding value is then G[1][2] ⊗ G[2][1] ⊗ G[2][5]). 
  corresponding_tuples_as_indices = vec([x for x = Base.Iterators.ProductIterator(Tuple(1:ngens(x) for x = G))])
  # In corresponding_tuples we store tuples of the actual generators, so in 
  # the example above we would store (G[1][2], G[2][1], G[2][5]).
  corresponding_tuples = map(index_tuple -> Tuple(map(index -> G[index][index_tuple[index]],1:length(index_tuple))), corresponding_tuples_as_indices)

  generating_tensors = map(mF, map(tuple -> map(x -> typeof(parent(x)) <: FreeMod ? x : x.repres, tuple), corresponding_tuples))
  s, emb = sub(F, generating_tensors, :map)
  #s, emb = sub(F, vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(gens(x, free_module(x)) for x = G))]), :map)
  q = vcat([vec([mF(x) for x = Base.Iterators.ProductIterator(Tuple(i == j ? rels(G[i]) : gens(free_module(G[i])) for i=1:length(G)))]) for j=1:length(G)]...) 
  local projection_map
  if length(q) != 0
    s, projection_map = quo(s, q, :map)
  end

  tuples_pure_tensors_dict = IdDict(zip(corresponding_tuples_as_indices, gens(s)))
  Hecke.set_special(s, :show => Hecke.show_tensor_product, :tensor_product => G)

  
  function pure(tuple_elems::Union{SubQuoElem,FreeModuleElem}...)
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

  Hecke.set_special(s, :tensor_pure_function => pure, :tensor_generator_decompose_function => decompose_generator)

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
Apply $P⊗-$ to `C`.
"""
function tensor_product(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})
  tensor_chain = Hecke.map_type(C)[]
  tensor_modules = [tensor_product(P, domain(C.maps[1]), task=:store)[1]]
  tensor_modules = vcat(tensor_modules, [tensor_product(P, codomain(f), task=:store)[1] for f in C.maps])

  for i=1:length(C)
    A = tensor_modules[i]
    B = tensor_modules[i+1]

    push!(tensor_chain, hom_tensor(A,B,[identity_map(P), map(C,i)]))
  end

  return Hecke.ChainComplex(ModuleFP, tensor_chain)
end

@doc Markdown.doc"""
  tensor_product(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
Apply $-⊗P$ to `C`.
"""
function tensor_product(C::Hecke.ChainComplex{ModuleFP}, P::ModuleFP)
  tensor_chain = Hecke.map_type(C)[]
  tensor_modules = [tensor_product(domain(C.maps[1]), P, task=:store)[1]]
  tensor_modules = vcat(tensor_modules, [tensor_product(codomain(f), P, task=:store)[1] for f in C.maps])

  for i=1:length(C)
    A = tensor_modules[i]
    B = tensor_modules[i+1]

    push!(tensor_chain, hom_tensor(A,B,[map(C,i), identity_map(P)]))
  end

  return Hecke.ChainComplex(ModuleFP, tensor_chain)
end

@doc Markdown.doc"""
  tor(M::ModuleFP, N::ModuleFP, i::Int)
Compute $Tor_i(M,N)$.
"""
function tor(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M)[1:end-2]
  lifted_resolution = tensor_product(free_res, N) #TODO only three homs are neccessary
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
Let $Hom_MP = Hom(M,P)$, $Hom_NP = Hom(N,P)$ and $phi = φ : N → M$ a morphism.
Compute $φ^{\ast} : Hom(M,P) → Hom(N,P)$.
"""
function lift_homomorphism_contravariant(Hom_MP::ModuleFP, Hom_NP::ModuleFP, phi::ModuleMap)
  # phi : N -> M
  M_P = get_special(Hom_MP, :hom)
  M_P === nothing && error("Both modules must be hom modules")
  N_P = get_special(Hom_NP, :hom)
  N_P === nothing && error("Both modules must be hom modules")
  
  @assert M_P[2] === N_P[2]
  M,P = M_P
  N,_ = N_P
  @assert domain(phi) === N
  @assert codomain(phi) === M
  
  phi_lifted = hom(Hom_MP, Hom_NP, [homomorphism_to_module_elem(Hom_NP, phi*homomorphism(f)) for f in gens(Hom_MP)])
  return phi_lifted
end

@doc Markdown.doc"""
  lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, phi::ModuleMap)
Let $Hom_PM = Hom(P,M)$, $Hom_PN = Hom(P,N)$ and $phi = φ : M → N$ a morphism.
Compute $φ_{\ast} : Hom(P,M) → Hom(P,N)$.
"""
function lift_homomorphism_covariant(Hom_PM::ModuleFP, Hom_PN::ModuleFP, phi::ModuleMap)
  # phi : M -> N
  P_M = get_special(Hom_PM, :hom)
  P_M === nothing && error("Both modules must be hom modules")
  P_N = get_special(Hom_PN, :hom)
  P_N === nothing && error("Both modules must be hom modules")

  @assert P_M[1] === P_N[1]
  P,M = P_M
  _,N = P_N
  @assert domain(phi) === M
  @assert codomain(phi) === N

  if iszero(Hom_PN)
    return hom(Hom_PM, Hom_PN, [zero(Hom_PN) for _=1:ngens(Hom_PM)])
  end
  phi_lifted = hom(Hom_PM, Hom_PN, [homomorphism_to_module_elem(Hom_PN, homomorphism(f)*phi) for f in gens(Hom_PM)])
  return phi_lifted
end

@doc Markdown.doc"""
  hom(P::ModuleFP, C::Hecke.ChainComplex{ModuleFP})
Apply $Hom(P,-)$ to `C`. Return the lifted chain complex.
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
Apply $Hom(-,P)$ to `C`. Return the lifted chain complex.
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
Compute $Ext^i(M,N)$.
"""
function ext(M::ModuleFP, N::ModuleFP, i::Int)
  free_res = free_resolution(M)[1:end-2]
  lifted_resolution = hom(free_res, N) #TODO only three homs are neccessary
  return homology(lifted_resolution,i)
end

#############################
# TODO ?
#############################
@doc Markdown.doc"""
  map_canonically(M::SubQuo, v::SubQuoElem)
Map the element `v` to an elemet of the module `M` using cached 
canonical homomorphisms between the parent Module of `v` and `M`.
"""
function map_canonically(M::SubQuo, v::SubQuoElem)
  N = parent(v)
  if N===M
    return v
  end

  # Breadth-First Search to find path to N:
  parent_hom = IdDict{SubQuo,ModuleMap}()
  modules = Set([M])
  found_N = false
  for A in modules
    for H in A.ingoing_morphisms
      B = domain(H)
      if B!==A # on trees "B!==A" is enough!
        if !(B in modules)
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

  function helper_dfs(U::SubQuo, D::SubQuo, visited::Set, path::Vector)
    if U === D
      push!(all_paths, path)
      return
    end
    for neighbor_morphism in U.outgoing_morphisms
      if !(neighbor_morphism in visited)

        #push!(visited, neighbor_morphism)
        #push!(path, neighbor_morphism)
        #helper_dfs(codomain(neighbor_morphism), D, visited, path)
        helper_dfs(codomain(neighbor_morphism), D, union(visited, Set([neighbor_morphism])), union(path, [neighbor_morphism]))
        #helper_dfs(codomain(neighbor_morphism), D, visited, union(path, [neighbor_morphism]))
      end
    end
  end

  helper_dfs(N, M, Set(), [])

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
  push!(codomain(f).ingoing_morphisms, f)
end

#############################
#TODO move to Hecke
#  re-evaluate and use or not
function differential(c::Hecke.ChainComplex, i::Int)
  return map(c,length(c)-i)
end

function module_in_complex(c::Hecke.ChainComplex, i::Int)
  return obj(c,length(c)-i)
end

getindex(c::Hecke.ChainComplex, i::Int) = module_in_complex(c,i)

function Base.getindex(r::Hecke.SRow, u::UnitRange)
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

function Base.getindex(r::Hecke.SRow, R::AbstractAlgebra.Ring, u::UnitRange)
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

function getindex(a::Hecke.SRow, b::AbstractArray{Int, 1})
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
Return the canonical projection from $F = R^I$ to $R^(indices)$ where $indices ⊂ I$.
"""
function projection(F::FreeMod, indices::AbstractArray)
  @assert all(x -> x <= ngens(F), indices)
  @assert length(Set(indices)) == length(indices) # unique indices
  R = base_ring(F)
  G = FreeMod(R, length(indices))
  return hom(F, G, [i in indices ? G[findfirst(x->x==i,indices)] : zero(G) for i=1:ngens(F)])
end

@doc Markdown.doc"""
    preimage(H::SubQuoHom,elems::Vector{SubQuoElem{T}}, task::Symbol = :none) where {T}
Return the preimage of the submodule generated by the Elements $elems$ under $H$
as a subquotient, as well as the injection homomorphism into the domain of $H$.
"""
function preimage(H::SubQuoHom,elems::Vector{SubQuoElem{T}}, task::Symbol = :none) where {T}
  if length(elems)==0
      throw(ArgumentError("too few arguments"))
  end
  R = base_ring(domain(H))
  row_length = ngens(codomain(H))
  submod = vcat((dense_row(coeffs(e), row_length) for e in elems)...)
  C = matrix(present_as_cokernel(codomain(H)).quo)
  A = vcat(matrix(H), C, submod)
  G = FreeMod(R, nrows(A))
  A = FreeModuleHom(G, FreeMod(R, ncols(A)), A)
  #K = kernel(A)
  K,kernel_injection = kernel(A)
  N = domain(H)
  n = ngens(N)
  generators = Vector{SubQuoElem{T}}()
  projection_map = projection(G, 1:n)
  for i=1:ngens(K)
      coeffs_for_new = projection_map(kernel_injection(K[i])).coords
      if isempty(coeffs_for_new)
        continue
      end
      new = SubQuoElem(coeffs_for_new, N)
      if !iszero(new)
          push!(generators, new)
      end
  end
  if length(generators)==0
      push!(generators,zero(N))
  end

  
  preimage, emb = sub(domain(H), generators, :map)
  preimage_pruned, prune_isomorphism, _ = simplify(preimage)
  if task != :none
    return preimage_pruned, prune_isomorphism*emb
  else
    return preimage_pruned
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
with the injection map $N --> M$ and the projection map $M --> N$.
"""
function simplify(M::SubQuo)
  function unit_vector_in_relations(i::Int, M::SubQuo)
    if !isdefined(M, :quo)
      return false
    end
    reduced_unit_vector = _reduce(convert(M.quo.gens.SF, M.F[i]), singular_generators(std_basis(M.quo)))
    return iszero(reduced_unit_vector)
  end

  function delete_rows(A::MatElem, to_delete::Vector{Int})
    Mat = A[setdiff(1:nrows(A),to_delete),:]
    return Mat
  end
  function delete_columns(A::MatElem, to_delete::Vector{Int})
    return delete_rows(A', to_delete)'
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

  M_generators = matrix(M.sub)
  M_relations = isdefined(M, :quo) ? matrix(M.quo) : zero_matrix(R, 1,ncols(M_generators))

  to_delete::Vector{Int} = []
  for i=1:size(M_relations)[2]
    if unit_vector_in_relations(i, M)
      push!(to_delete, i)
    end
  end

  new_generators = delete_columns(M_generators, to_delete)
  new_relations = delete_columns(M_relations, to_delete)

  to_delete,_,_ = rows_to_delete(vcat(new_generators, new_relations)',size(new_relations)[2])

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
      unit_vector_index = i-length(filter(x -> x < i, to_delete))
      unit_vector = [j == unit_vector_index ? R(1) : R(0) for j=1:size(projection_matrix)[2]]
      assign_row!(projection_matrix, unit_vector, i)
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

  #=if isgraded(M)
    ishomogeneous_function_lookup_table[SQ] = x -> custom_ishomogeneous(injection(x))
    degree_function_lookup_table[SQ] = x -> custom_degree(injection(x))
  end=#

  return SQ,injection,projection
end

######################################
# Not only for testing
######################################
@doc Markdown.doc"""
  matrix_to_map(F::FreeMod{T}, A::MatElem{T}) where T
Converts a given n×m-matrix into the corresponding morphism $A : R^n → F$, 
with `rank(F) == m`.
"""
function matrix_to_map(F::FreeMod{T}, A::MatElem{T}) where T
  R = base_ring(F)
  F_domain = FreeMod(R, nrows(A))

  phi = FreeModuleHom(F_domain, F, A)
  return phi
end

@doc Markdown.doc"""
  matrix_to_map(A::MatElem)
Converts a given n×m-matrix into the corresponding morphism $A : R^n → R^m$.
"""
function matrix_to_map(A::MatElem)
  R = base_ring(A)
  F_codomain = FreeMod(R, ncols(A))
  return matrix_to_map(F_codomain,A)
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

@doc Markdown.doc"""
    iswelldefined(H::ModuleMap)
Test if `H` is well-defined.
"""
function iswelldefined(H::ModuleMap)
  if typeof(H) <: FreeModuleHom
    return true
  end
  M = domain(H)
  C = matrix(present_as_cokernel(M).quo)
  n,m = size(C)
  g = M[1]
  ImH = map(x -> H(x), gens(M))
  for i=1:n
    if !iszero(Base.sum([C[i,j]*ImH[j] for j=1:m]))
      return false
    end
  end
  return true
end

#############################################
# Test hom
#############################################
function to_julia_matrix(A::Union{MatElem})
  return eltype(A)[A[i, j] for i = 1:nrows(A), j = 1:ncols(A)]
end

function Base.reshape(M::MatElem, n, m)
  julia_matrix = to_julia_matrix(M)
  julia_matrix = reshape(julia_matrix, n, m)
  R = base_ring(M)
  mat_space = MatrixSpace(R, n, m)
  return mat_space(julia_matrix)
end

function hom2(f1::MatElem{T}, g1::MatElem{T}) where T
  R = base_ring(f1)
  s1, s0 = size(f1)
  t1, t0 = size(g1)

  g2 = matrix_kernel(g1)
  n = s1*t0
  m = s0*t0 + s1*t1
  delta::MatrixElem{T} = zero_matrix(R, m,n)
  #delta::Matrix{T} = zeros(R, n, m)
  for j=1:m
    #b_vector::Vector{T} = [i == j ? R(1) : R(0) for i=1:m]
    b_vector::MatrixElem{T} = zero_matrix(R, 1,m)
    b_vector[1,j] = R(1)
    #println(typeof(Core.Array(b_vector[1,1:s0*t0])))
    A = reshape(b_vector[1,1:s0*t0], s0, t0)
    #Mat = AbstractAlgebra.MatrixSpace(R, t0, s0)
    #A = Mat(A)
    B = reshape(b_vector[1,s0*t0+1:length(b_vector)], s1, t1)
    #Mat = AbstractAlgebra.MatrixSpace(R, t1, s1)
    #B = Mat(B)
    res = f1*A - B*g1
    res = reshape(res, 1, length(res))
    for i=1:length(res)
      delta[j,i] = res[1,i]
    end
  end

  #Mat = AbstractAlgebra.MatrixSpace(R, size(delta)...)
  #display(delta)
  gamma = matrix_kernel(delta)

  t2 = size(g2)[1]
  n = m
  m = s0*t1 + s1*t2
  #rho::Matrix{T} = zeros(R, n, m)
  rho::MatrixElem{T} = zero_matrix(R, m,n)
  for j=1:m
    #b_vector::Vector{T} = [R(0) for i=1:m]
    b_vector = zero_matrix(R, 1,m)
    b_vector[1,j] = R(1)
    C = reshape(b_vector[1,1:s0*t1], s0, t1)
    #Mat = AbstractAlgebra.MatrixSpace(R, t1, s0)
    #C = Mat(C)
    D = reshape(b_vector[1,s0*t1+1:length(b_vector)], s1, t2)
    #Mat = AbstractAlgebra.MatrixSpace(R, t2, s1)
    #D = Mat(D)
    res1 = C*g1
    res1 = reshape(res1, 1,length(res1))
    res2 = f1*C - D*g2
    res2 = reshape(res2, 1,length(res2))
    for i=1:length(res1)
        rho[j,i] = res1[1,i]
    end
    for i=1:length(res2)
        #rho[:,j] = vcat(vec(res1), vec(res2))
        rho[j,i+length(res1)] = res2[1,i]
    end
  end

  M = SubQuo(gamma, rho)

  function convert_to_matrix(v::SubQuoElem{T})
    if parent(v) !== M
      throw(DomainError("v does not represent a homomorphism"))
    end
    R = base_ring(M)
    #A = AbstractAlgebra.matrix(R,map(x -> convert_poly(R,x),Array(v.singular_v)))
    #A = map(x -> convert_poly(R,x),Array(v.singular_v))
    #A = AbstractAlgebra.matrix(R, t0, s0, A[1:s0*t0])'
    A = reshape(dense_row(v.repres.coords[R, 1:s0*t0], s0*t0), s0, t0)
    #A = matrix(R, t0, s0, dense_row(v.repres.coords[R, 1:s0*t0], s0*t0))'
    #A = v.v*v.parent_sq.generators  # ineffizienter
    #A = reshape(A[1,1:s0*t0], s0, t0)
    #N = cokernel(g1)
    #for i=1:s0ems) <: Tuple ? [i for i in elems] : [elems])

    #    row = lift(N,A[i,:])
    #    for j=1:t0
    #        A[i,j]=row[1,j]
    #    end
    #end
    return A
  end

  return M, convert_to_matrix
end

# lookup-tables only used in the old code
# lookup-tables to find the correct function, that converts subquotient elements to homomorphisms or homomorphisms to subquotient elements:
#global subquotient_elem_to_homomorphism_lookup_table = IdDict{Subquotient,Function}()
#global homomorphism_to_subquotient_elem_lookup_table = IdDict{Tuple{Subquotient,Subquotient},Function}()


# since hom yields possibly wrong results, this is an alternative implementation for comparison, 
# which is known to work correctly in various examples (by comparing with Macaulay2 and by testing 
# whether the homomorphism corresponding to elements of the subquotient are well-defined
# the function also supports simplification of the subquotient, while retaining the correspondence to 
# morphisms, which is essential for checking correctness and practical use of the output
@doc Markdown.doc"""
  hom2(M::Subquotient,N::Subquotient)
> Return a subquotient S such that $Hom(M,N) \cong S$
"""
function hom2(M::SubQuo{T},N::SubQuo{T},simplify=true) where T
  f1 = matrix(present_as_cokernel(M).quo)
  g1 = matrix(present_as_cokernel(N).quo)
  #f1 = presentation_matrix(M)
  #g1 = presentation_matrix(N)
  SQ, convert_to_matrix = hom2(f1,g1)
  if simplify
    SQ2, i, p = simplify(SQ)
    to_homomorphism = function(elem::SubQuoElem{T})
      elem2 = i(elem)
      A = convert_to_matrix(elem2)
      return SubQuoHom(M,N,A)
      #return AbstractAlgebra.ModuleHomomorphism(M,N,A)
    end
    to_subquotient_elem = function(H::ModuleMap)
      m = length(matrix(H))
      v = reshape(matrix(H),1,m)
      #v = i(v)
      tmp_sq = SubQuo(matrix(SQ.sub)[:,1:m], matrix(SQ.quo)[:,1:m])
      #tmp_sq = Subquotient(SQ.generators[:,1:m], SQ.relations[:,1:m])
      #coeffs = Nemo.lift(tmp_sq, v)
      v = FreeModuleElem(sparse_row(v), FreeMod(R, length(v)))
      coeffs = coordinates(v, tmp_sq)
      return p(SubQuoElem(coeffs, SQ))
      #return p(SubquotientElem(SQ, coeffs))
    end
    #subquotient_elem_to_homomorphism_lookup_table[SQ2] = elem -> to_homomorphism(elem)
    #homomorphism_to_subquotient_elem_lookup_table[(M,N)] = H -> to_subquotient_elem(H)

    to_hom_map = MapFromFunc(to_homomorphism, to_subquotient_elem, SQ2, Hecke.MapParent(M, N, "homomorphisms"))
    Hecke.set_special(SQ2, :hom => (M, N), :module_to_hom_map => to_hom_map)

    return SQ2
  else
    to_subquotient_elem = function(H::ModuleMap)
      m = length(matrix(H))
      v = reshape(matrix(H),1,m)
      #tmp_sq = Subquotient(SQ.generators[:,1:m], SQ.relations[:,1:m])
      tmp_sq = SubQuo(matrix(SQ.sub)[:,1:m], matrix(SQ.quo)[:,1:m])
      v = FreeModuleElem(sparse_row(v), FreeMod(R, length(v)))
      coeffs = coordinates(v, tmp_sq)
      #coeffs = Nemo.lift(tmp_sq, v)
      return SubQuoElem(coeffs, SQ)
      #return SubquotientElem(SQ, coeffs)
    end
    to_homomorphism = function(elem::SubQuoElem{T})
      A = convert_to_matrix(elem)
      return SubQuoHom(M,N,A)
      #return AbstractAlgebra.ModuleHomomorphism(M,N,A)
    end
    #subquotient_elem_to_homomorphism_lookup_table[SQ] = to_homomorphism
    #homomorphism_to_subquotient_elem_lookup_table[(M,N)] = H -> to_subquotient_elem(H)

    to_hom_map = MapFromFunc(to_homomorphism, to_subquotient_elem, SQ, Hecke.MapParent(M, N, "homomorphisms"))
    Hecke.set_special(SQ, :hom => (M, N), :module_to_hom_map => to_hom_map)

    return SQ
  end
end
