export ModuleFP, ModuleFPElem, ModuleFPHom, FreeMod,
       FreeModElem, SubQuo, SubQuoElem, FreeModuleHom, SubQuoHom,
       FreeMod_dec, FreeModElem_dec, FreeModuleHom_dec, FreeResolution

@doc Markdown.doc"""
    ModuleFP{T}

The abstract supertype of all finitely presented modules.
The type variable `T` refers to the type of the elements of the base ring.
"""
abstract type ModuleFP{T} end

@doc Markdown.doc"""
    AbstractFreeMod{T} <: ModuleFP{T}

The abstract supertype of all finitely generated free modules.
"""
abstract type AbstractFreeMod{T} <: ModuleFP{T} end

@doc Markdown.doc"""
    AbstractSubQuo{T} <: ModuleFP{T}

The abstract supertype of all finitely presented subquotient modules.
"""
abstract type AbstractSubQuo{T} <: ModuleFP{T} end


@doc Markdown.doc"""
    ModuleFPElem{T} <: ModuleElem{T}

The abstract supertype of all elements of finitely presented modules.
"""
abstract type ModuleFPElem{T} <: ModuleElem{T} end

@doc Markdown.doc"""
    AbstractFreeModElem{T} <: ModuleFPElem{T}

The abstract supertype of all elements of finitely generated free modules.
"""
abstract type AbstractFreeModElem{T} <: ModuleFPElem{T} end

@doc Markdown.doc"""
    AbstractSubQuoElem{T} <: ModuleFPElem{T}

The abstract supertype of all elements of subquotient modules.
"""
abstract type AbstractSubQuoElem{T} <: ModuleFPElem{T} end

abstract type ModuleFPHomDummy end

@doc Markdown.doc"""
    ModuleFPHom{T1, T2}

The abstract supertype for morphisms of finitely presented modules over multivariate polynomial rings .
`T1` and `T2` are the types of domain and codomain respectively.
"""
abstract type ModuleFPHom{T1, T2} <: Map{T1, T2, Hecke.HeckeMap, ModuleFPHomDummy} end

parent(f::ModuleFPHom) = Hecke.MapParent(domain(f), codomain(f), "homomorphisms")

@doc Markdown.doc"""
    FreeMod{T <: RingElem} <: ModuleFP{T}

The type of free modules.
Free modules are determined by their base ring, the rank and the names of 
the (standard) generators.
Moreover, canonical incoming and outgoing morphisms are stored if the corresponding
option is set in suitable functions.
`FreeMod{T}` is a subtype of `ModuleFP{T}`.
"""
@attributes mutable struct FreeMod{T <: RingElem} <: AbstractFreeMod{T}
  R::Ring
  n::Int
  S::Vector{Symbol}

  incoming_morphisms::Vector{<:ModuleFPHom}
  outgoing_morphisms::Vector{<:ModuleFPHom}

  function FreeMod{T}(n::Int,R::Ring,S::Vector{Symbol}) where T <: RingElem
    r = new{elem_type(R)}()
    r.n = n
    r.R = R
    r.S = S

    r.incoming_morphisms = Vector{ModuleFPHom}()
    r.outgoing_morphisms = Vector{ModuleFPHom}()

    return r
  end
end

@doc Markdown.doc"""
    FreeModElem{T}

The type of free module elements. An element of a free module $F$ is given by a sparse row (`SRow`)
which specifies its coordinates with respect to the basis of standard unit vectors of $F$.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> F = free_module(R, 3)
Free module of rank 3 over Multivariate Polynomial Ring in x, y over Rational Field

julia> f = F(sparse_row(R, [(1,x),(3,y)]))
x*e[1] + y*e[3]

julia> typeof(f)
FreeModElem{fmpq_mpoly}

julia> g = x*F[1] + y*F[3]
x*e[1] + y*e[3]

julia> f == g
true
```
"""
struct FreeModElem{T} <: AbstractFreeModElem{T}
  coords::SRow{T} # also usable via coeffs()
  parent::FreeMod{T}

  function FreeModElem{T}(coords::SRow{T}, parent::FreeMod{T}) where T
    r = new{T}(coords,parent)
    return r
  end
end

@doc Markdown.doc"""
    ModuleGens{T}

Data structure for a generating systems for submodules.
Contains structures for the generators, the corresponding module on the Singular side, 
the embedding free module, the embedding free module on the Singular side.
Subquotients will be built from a tuple of submodules which again are given by 
generating sets. In this way, the Singular stuff is hidden on the higher structures
and all the conversion is taken care of here.

This data structure is also used for representing Gröbner / standard bases.
Relative Gröbner / standard bases are also supported.
"""
@attributes mutable struct ModuleGens{T} # T is the type of the elements of the ground ring.
  O::Vector{FreeModElem{T}}
  S::Singular.smodule
  F::FreeMod{T}
  SF::Singular.FreeMod

  isGB::Bool
  is_reduced::Bool
  ordering::ModuleOrdering
  quo_GB::ModuleGens{T} # Pointer to the quotient GB when having a relative GB

  function ModuleGens{T}(O::Vector{<:FreeModElem}, F::FreeMod{T}) where {T}
    r = new{T}()
    r.O = O
    r.F = F
    return r
  end

  # ModuleGens from an Array of Oscar free module elements, specifying the free module 
  # and Singular free module, only useful indirectly
  function ModuleGens{T}(O::Vector{<:FreeModElem}, F::FreeMod{T}, SF::Singular.FreeMod) where {T}
    r = new{T}()
    r.O = O
    r.SF = SF
    r.F = F
    return r
  end

  # ModuleGens from a Singular submodule
  function ModuleGens{S}(F::FreeMod{S}, s::Singular.smodule) where {S} # FreeMod is necessary due to type S
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
    SubModuleOfFreeModule{T} <: ModuleFP{T}

Data structure for submodules of free modules. `SubModuleOfFreeModule` shouldn't be
used by the end user.
When computed, a standard basis (computed via `standard_basis()`) and generating matrix (that is the rows of the matrix
generate the submodule) (computed via `generator_matrix()`) are cached.
"""
@attributes mutable struct SubModuleOfFreeModule{T} <: ModuleFP{T}
  F::FreeMod{T}
  gens::ModuleGens{T}
  groebner_basis::Dict{ModuleOrdering, ModuleGens{T}}
  default_ordering::ModuleOrdering
  matrix::MatElem

  function SubModuleOfFreeModule{R}(F::FreeMod{R}) where {R}
    # this does not construct a valid SubModuleOfFreeModule
    r = new{R}()
    r.F = F
    r.groebner_basis = Dict()
    return r
  end
end


@doc Markdown.doc"""
    SubQuo{T} <: ModuleFP{T}

The type of subquotient modules.
A subquotient module $M$ is a module where $M = A + B / B$ where $A$ and $B$ are 
submodules of a free module.
$A$, $B$ and $A+B$ (they have type `SubModuleOfFreeModule`) as well as the embedding
free module are stored. 
One can construct ordinary submodules of free modules by not giving $B$.
Moreover, canonical incoming and outgoing morphisms are stored if the corresponding
option is set in suitable functions.
`SubQuo{T}` is a subtype of `ModuleFP{T}`.
"""
@attributes mutable struct SubQuo{T} <: AbstractSubQuo{T}
  #meant to represent sub+ quo mod quo - as lazy as possible
  F::FreeMod{T}
  sub::SubModuleOfFreeModule
  quo::SubModuleOfFreeModule
  sum::SubModuleOfFreeModule

  groebner_basis::Dict{ModuleOrdering, ModuleGens{T}}

  incoming_morphisms::Vector{<:ModuleFPHom}
  outgoing_morphisms::Vector{<:ModuleFPHom} # TODO is it possible to make ModuleFPHom to SubQuoHom?

  function SubQuo{R}(F::FreeMod{R}) where {R}
    # this does not construct a valid subquotient
    r = new{R}()
    r.F = F

    r.groebner_basis = Dict()
    r.incoming_morphisms = Vector{ModuleFPHom}()
    r.outgoing_morphisms = Vector{ModuleFPHom}()

    return r
  end
end


@doc Markdown.doc"""
    SubQuoElem{T}

The type of subquotient elements. An element $f$ of a subquotient $M$ over the ring $R$
is given by a sparse row (`SRow`) which specifies the coefficients of an $R$-linear 
combination of the generators of $M$ which defines $f$.

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

julia> f = SubQuoElem(sparse_row(R, [(1,z),(2,one(R))]),M)
(x*z + y)*e[1]

julia> g = z*M[1] + one(R)*M[2]
(x*z + y)*e[1]

julia> typeof(g)
SubQuoElem{fmpq_mpoly}

julia> f == g
true
```
"""
struct SubQuoElem{T} <: AbstractSubQuoElem{T} 
  coeffs::SRow{T}
  repres::FreeModElem{T}
  parent::SubQuo

  function SubQuoElem{R}(v::SRow{R}, SQ::SubQuo) where {R}
    @assert length(v) <= ngens(SQ.sub)
    if isempty(v)
      r = new{R}(v, zero(SQ.F), SQ)
      return r
    end
    r = new{R}(v, sum([v[i]*SQ.sub[i] for i=1:ngens(SQ.sub)]), SQ)
    return r
  end

  function SubQuoElem{R}(a::FreeModElem{R}, SQ::SubQuo) where {R}
    @assert a.parent === SQ.F
    r = new{R}(coordinates(a,SQ), a, SQ)
    return r
  end
end


mutable struct SubQuoHom{
    T1<:AbstractSubQuo, 
    T2<:ModuleFP, 
    RingMapType<:Any
  } <: ModuleFPHom{T1, T2}
  matrix::MatElem
  header::Hecke.MapHeader
  im::Vector
  inverse_isomorphism::ModuleFPHom
  ring_map::RingMapType

  # Constructors for maps without change of base ring
  function SubQuoHom{T1,T2,RingMapType}(D::SubQuo, C::FreeMod, im::Vector) where {T1,T2,RingMapType}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, Nothing}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{FreeModElem}(im)
    return r
  end

  function SubQuoHom{T1,T2,RingMapType}(D::SubQuo, C::SubQuo, im::Vector) where {T1,T2,RingMapType}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, Nothing}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{SubQuoElem}(im)
    return r
  end

  function SubQuoHom{T1,T2,RingMapType}(D::SubQuo, C::ModuleFP, im::Vector) where {T1,T2,RingMapType}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, Nothing}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = im
    return r
  end

  # Constructors for maps with change of base ring
  function SubQuoHom{T1,T2,RingMapType}(
      D::SubQuo, 
      C::FreeMod, 
      im::Vector, 
      h::RingMapType
    ) where {T1,T2,RingMapType}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, RingMapType}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{FreeModElem}(im)
    r.ring_map = h
    return r
  end

  function SubQuoHom{T1,T2,RingMapType}(
      D::SubQuo, 
      C::SubQuo, 
      im::Vector, 
      h::RingMapType
    ) where {T1,T2,RingMapType}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, RingMapType}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{SubQuoElem}(im)
    r.ring_map = h
    return r
  end

  function SubQuoHom{T1,T2,RingMapType}(
      D::SubQuo, 
      C::ModuleFP, 
      im::Vector, 
      h::RingMapType
    ) where {T1,T2,RingMapType}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, RingMapType}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = im
    r.ring_map = h
    return r
  end

end


###############################################################################
# Graded modules
###############################################################################
const CRing_dec = Union{MPolyRing_dec, MPolyQuo{<:Oscar.MPolyElem_dec}}
const CRingElem_dec = Union{MPolyElem_dec, MPolyQuoElem{<:Oscar.MPolyElem_dec}}
#TODO: other name for CRing_dec -> which?

@doc Markdown.doc"""
    FreeMod_dec{T <: CRingElem_dec} <: ModuleFP_dec{T}

The type of decorated (graded or filtered) free modules.
Decorated free modules are determined by their base ring, the rank,
the grading or filtration and the names of the (standard) generators.
Moreover, canonical incoming and outgoing morphisms are stored if the corresponding
option is set in suitable functions.
`FreeMod_dec{T}` is a subtype of `ModuleFP{T}`.
"""
@attributes mutable struct FreeMod_dec{T <: CRingElem_dec} <: AbstractFreeMod{T}
  F::FreeMod{T}
  d::Vector{GrpAbFinGenElem}

  function FreeMod_dec{T}(F::FreeMod, d::Vector{GrpAbFinGenElem}) where T <: CRingElem_dec
    @assert length(d) == rank(F)
    r = new{elem_type(base_ring(F))}(F, d)
    return r
  end

  function FreeMod_dec{T}(R::CRing_dec,S::Vector{Symbol},d::Vector{GrpAbFinGenElem}) where T <: CRingElem_dec
    r = new{elem_type(R)}()
    r.F = FreeMod{T}(length(d),R,S)
    r.d = d
    return r
  end
end

@doc Markdown.doc"""
    FreeModElem_dec{T}

The type of decorated free module elements. An element of a decorated free module $F$ is 
given by a sparse row (`SRow`) which specifies its coordinates with respect to the basis
of standard unit vectors of $F$.
"""
struct FreeModElem_dec{T} <: AbstractFreeModElem{T}
  coords::SRow{T} # also usable via coeffs()
  parent::FreeMod_dec{T}

  function FreeModElem_dec{T}(coords::SRow{T}, parent::FreeMod_dec{T}) where T
    r = new{T}(coords,parent)
    return r
  end
end





const ModuleFP_dec{T} = Union{FreeMod_dec{T}} # SubQuo_dec{T} will be included
const ModuleFPElem_dec{T} = Union{FreeModElem_dec{T}} # SubQuoElem_dec{T} will be included


@doc Markdown.doc"""
    FreeModuleHom{T1, T2, RingMapType} <: ModuleFPHom{T1, T2} 

Data structure for morphisms where the domain is a free module (`FreeMod`).
`T1` and `T2` are the types of domain and codomain respectively.
`FreeModuleHom` is a subtype of `ModuleFPHom`.
When computed, the corresponding matrix (via `matrix()`) and inverse isomorphism
(in case there exists one) (via `inv()`) are cached.
"""
@attributes mutable struct FreeModuleHom{
    T1 <: AbstractFreeMod,
    T2 <: ModuleFP,
    RingMapType <: Any} <: ModuleFPHom{T1, T2} 
  header::MapHeader
  ring_map::RingMapType
  
  matrix::MatElem
  inverse_isomorphism::ModuleFPHom

  # generate homomorphism of free modules from F to G where the vector a contains the images of
  # the generators of F
  function FreeModuleHom(
      F::AbstractFreeMod, G::S, a::Vector{ModuleElemType}
    ) where {S<:ModuleFP, ModuleElemType<:ModuleFPElem}
    @assert all(x->parent(x) === G, a)
    @assert length(a) == ngens(F)
    r = new{typeof(F), typeof(G), Nothing}()
    function im_func(x::AbstractFreeModElem)
      b = zero(G)
      for (i,v) = x.coords
        b += v*a[i]
      end
      return b
    end
    function pr_func(x)
      @assert parent(x) === G
      c = coordinates(repres(x), sub(G, a, :module))
      return FreeModElem(c, F)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)
    return r
  end

  function FreeModuleHom(
      F::AbstractFreeMod, G::T2, a::Vector{ModuleElemType}, h::RingMapType
    ) where {T2, ModuleElemType<:ModuleFPElem, RingMapType}
    @assert all(x->parent(x) === G, a)
    @assert length(a) == ngens(F)
    @assert h(one(base_ring(F))) == one(base_ring(G))
    r = new{typeof(F), T2, RingMapType}()
    function im_func(x::AbstractFreeModElem)
      b = zero(G)
      for (i,v) = x.coords
        b += h(v)*a[i]
      end
      return b
    end
    r.header = MapHeader{typeof(F), T2}(F, G, im_func)
    r.ring_map = h
    return r
  end

end

# Further constructors taking matrices as input
function FreeModuleHom(
    F::AbstractFreeMod{T}, G::S, mat::MatElem{T}
  ) where {T<:RingElem,S<:AbstractFreeMod}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  hom = FreeModuleHom(F, G, [FreeModElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)])
  hom.matrix = mat
  return hom
end

function FreeModuleHom(
    F::AbstractFreeMod{T}, G::S, mat::MatElem{T}
  ) where {T<:RingElem, S<:ModuleFP}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  hom = FreeModuleHom(F, G, [SubQuoElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)])
  hom.matrix = mat
  return hom
end

function FreeModuleHom(
    F::AbstractFreeMod, G::S, mat::MatElem, h::RingMapType
  ) where {S<:AbstractFreeMod, RingMapType}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  @assert base_ring(mat) === base_ring(G)
  hom = FreeModuleHom(F, G, [FreeModElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)], h)
  hom.matrix = mat
  return hom
end

function FreeModuleHom(
    F::AbstractFreeMod, G::S, mat::MatElem, h::RingMapType
  ) where {S<:ModuleFP, RingMapType}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  @assert base_ring(mat) === base_ring(G)
  hom = FreeModuleHom(F, G, [SubQuoElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)], h)
  hom.matrix = mat
  return hom
end

struct FreeModuleHom_dec{
    T1 <: AbstractFreeMod,
    T2 <: ModuleFP,
    RingMapType <: Any} <: ModuleFPHom{T1, T2}
  f::FreeModuleHom{T1,T2, RingMapType}
  header::MapHeader
  # TODO degree and homogeneity

  function FreeModuleHom_dec(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, a::Vector) where {T}
    f = FreeModuleHom(F,G,a)
    r = new{typeof(F), typeof(G), Nothing}(f, f.header)
    return r
  end

  function FreeModuleHom_dec(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, mat::MatElem{T}) where {T}
    f = FreeModuleHom(F,G,mat)
    r = new{typeof(F), typeof(G), Nothing}(f, f.header)
    return r
  end

  function FreeModuleHom_dec(F::FreeMod_dec, G::ModuleFP_dec, a::Vector, h::RingMapType) where {RingMapType}
    f = FreeModuleHom(F,G,a,h)
    r = new{typeof(F), typeof(G), RingMapType}(f, f.header)
    return r
  end

  function FreeModuleHom_dec(F::FreeMod_dec, G::ModuleFP_dec{T}, mat::MatElem{T}, h::RingMapType) where {T, RingMapType}
    f = FreeModuleHom(F,G,mat,h)
    r = new{typeof(F), typeof(G), RingMapType}(f, f.header)
    return r
  end
end

@doc Markdown.doc"""
    FreeResolution{T}

Data structure for free resolutions.
"""
mutable struct FreeResolution{T}
    C::Hecke.ChainComplex
    complete::Bool

    function FreeResolution(C::Hecke.ChainComplex{T}) where {T}
        FR = new{T}()
        FR.C = C

        return FR
    end
end

Base.getindex(FR::FreeResolution, i::Int) = FR.C[i]

function Base.show(io::IO, FR::FreeResolution)
    C   = FR.C
    rk  = Dict{Int, String}()
    rng = Hecke.map_range(C)
    len = 0

    println(io, "")
    print(io, "rank   | ")
    # get names
    for i = rng
        if i != last(rng)
            rk[i] =  "$(rank(C[i]))"
            len   += length(rk[i]) + 2
            continue
        end
    end
    for i = rng
        if i != last(rng)
            print(io, rk[i], "  ")
        end
    end
    println(io, "")
    print(io, "-------|")
    print(io, repeat("-", len))
    println(io, "")
    print(io, "degree | ")
    for i = rng
        if i != last(rng)
            print(io, i, repeat(" ", length(rk[i]) + 1))
            continue
        end
    end
end
