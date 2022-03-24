export ModuleFP, ModuleFPElem, ModuleFPHom, ModuleMap, FreeMod,
       FreeModElem, SubQuo, SubQuoElem, FreeModuleHom, SubQuoHom,
       FreeMod_dec, FreeModElem_dec, FreeModuleHom_dec

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

abstract type ModuleFPHom end

@doc Markdown.doc"""
    ModuleMap{T1, T2}

The abstract supertype of module morphisms.
`T1` and `T2` are the types of domain and codomain respectively.
"""
abstract type ModuleMap{T1, T2} <: Map{T1, T2, Hecke.HeckeMap, ModuleFPHom} end

parent(f::ModuleMap) = Hecke.MapParent(domain(f), codomain(f), "homomorphisms")

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

  incoming_morphisms::Vector{<:ModuleMap}
  outgoing_morphisms::Vector{<:ModuleMap}

  function FreeMod{T}(n::Int,R::Ring,S::Vector{Symbol}) where T <: RingElem
    r = new{elem_type(R)}()
    r.n = n
    r.R = R
    r.S = S

    r.incoming_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

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
mutable struct ModuleGens{T} # T is the type of the elements of the ground ring.
  O::Vector{FreeModElem{T}}
  S::Singular.smodule
  F::FreeMod{T}
  SF::Singular.FreeMod

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

  function ModuleGens{T}(O::Vector{FreeModElemType}, F::FreeMod{T}) where {T, FreeModElemType<:FreeModElem}
    r = new{T}()
    r.O = O
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


###############################################################################
# Gröbner basis
###############################################################################
mutable struct ModuleGB{T}
  groebner_basis::ModuleGens{T}
  quo_groebner_basis::ModuleGB{T}
  reduced_groebner_basis::ModuleGens{T}

  leading_monomials::ModuleGens{T}
  ordering::ModuleOrdering
  #ordering::Singular.sordering # Remove

  #function ModuleGB{R}(gb::ModuleGens{R}, lm::ModuleGens{R}, ord::Singular.sordering) where R
  function ModuleGB{R}(gb::ModuleGens{R}, lm::ModuleGens{R}, ord::ModuleOrdering) where R
    r = new{R}()
    r.groebner_basis = gb
    r.leading_monomials = lm
    r.ordering = ord

    return r
  end

  #function ModuleGB{R}(gb::ModuleGens{R}, qgb::ModuleGB{R}, lm::ModuleGens{R}, ord::Singular.sordering) where R # Remove
  function ModuleGB{R}(gb::ModuleGens{R}, qgb::ModuleGB{R}, lm::ModuleGens{R}, ord::ModuleOrdering) where R
    r = new{R}()
    r.groebner_basis = gb
    r.quo_groebner_basis = qgb
    r.leading_monomials = lm
    r.ordering = ord

    return r
  end
end

@doc Markdown.doc"""
    SubModuleOfFreeModule{T} <: ModuleFP{T}

Data structure for submodules of free modules. `SubModuleOfFreeModule` shouldn't be
used by the end user.
When computed, a standard basis (computed via `std_basis()`) and generating matrix (that is the rows of the matrix
generate the submodule) (computed via `generator_matrix()`) are cached.
"""
mutable struct SubModuleOfFreeModule{T} <: ModuleFP{T}
  F::FreeMod{T}
  gens::ModuleGens{T}
  #groebner_basis::Dict{Singular.sordering, ModuleGB{T}} # Remove
  groebner_basis::Dict{ModuleOrdering, ModuleGB{T}}
  #default_ordering::Singular.sordering # Remove
  default_ordering::ModuleOrdering
  matrix::MatElem

  function SubModuleOfFreeModule{R}(F::FreeMod{R}, gens::Vector{<:FreeModElem}, 
                                       default_ordering::ModuleOrdering = default_ordering(base_ring(F))) where {R}
  #                                    default_ordering::Singular.sordering = default_ordering(base_ring(F))) where {R} # Remove
    @assert all(x -> parent(x) === F, gens)
    # @assert exactly_one_module_ordering(default_ordering) # Remove
    # @assert ordering_compatible_with_ring(base_ring(F), default_ordering) # Remove

    r = new{R}()
    r.F = F
    r.gens = ModuleGens(gens, F, default_ordering)
    r.default_ordering = default_ordering
    r.groebner_basis = Dict()
    return r
  end

  function SubModuleOfFreeModule{R}(F::FreeMod{R}, singular_module::Singular.smodule,
                                       default_ordering::ModuleOrdering = default_ordering(base_ring(F))) where {R}
  #                                    default_ordering::Singular.sordering = default_ordering(base_ring(F))) where {R} # Remove
    # @assert exactly_one_module_ordering(default_ordering) # Remove
    # @assert ordering_compatible_with_ring(base_ring(F), default_ordering) # Remove
    
    r = new{R}()
    r.F = F
    r.gens = ModuleGens(F, singular_module)
    r.default_ordering = default_ordering
    r.groebner_basis = Dict()
    return r
  end
  
  function SubModuleOfFreeModule{R}(F::FreeMod{R}, gens::ModuleGens,
                                       default_ordering::ModuleOrdering = default_ordering(base_ring(F))) where {R}
#                                      default_ordering::Singular.sordering = default_ordering(base_ring(F))) where {R} # Remove
    # @assert exactly_one_module_ordering(default_ordering) # Remove
    # @assert ordering_compatible_with_ring(base_ring(F), default_ordering) # Remove

    r = new{R}()
    r.F = F
    r.gens = gens
    r.default_ordering = default_ordering
    r.groebner_basis = Dict()
    return r
  end

  function SubModuleOfFreeModule{L}(F::FreeMod{L}, A::MatElem{L},
                                       default_ordering::ModuleOrdering = default_ordering(base_ring(F))) where {L}
#                                      default_ordering::Singular.sordering = default_ordering(base_ring(F))) where {L} # Remove
    # @assert exactly_one_module_ordering(default_ordering) # Remove
    # @assert ordering_compatible_with_ring(base_ring(F), default_ordering) # Remove

    r = new{L}()
    r.F = F
    O = [FreeModElem(sparse_row(A[i,:]), F) for i in 1:nrows(A)]
    r.gens = ModuleGens(O, F, default_ordering)
    r.matrix = A
    r.default_ordering = default_ordering
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

  #groebner_basis::Dict{Singular.sordering, ModuleGB{T}} # Remove
  groebner_basis::Dict{ModuleOrdering, ModuleGB{T}}

  incoming_morphisms::Vector{<:ModuleMap}
  outgoing_morphisms::Vector{<:ModuleMap} # TODO is it possible to make ModuleMap to SubQuoHom?

  function SubQuo{R}(sub::SubModuleOfFreeModule{R}) where {R}
    r = new{R}()
    r.F = sub.F
    r.sub = sub
    r.sum = r.sub

    r.groebner_basis = Dict()

    r.incoming_morphisms = Vector{ModuleMap}()
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

    r.groebner_basis = Dict()

    r.incoming_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  function SubQuo{R}(F::FreeMod{R}, O::Vector{<:FreeModElem}) where {R}
    r = new{R}()
    r.F = F
    r.sub = SubModuleOfFreeModule(F, O)
    r.sum = r.sub

    r.groebner_basis = Dict()

    r.incoming_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  function SubQuo{L}(S::SubQuo{L}, O::Vector{<:FreeModElem}) where {L} #TODO to be replaced by quo
    r = new{L}()
    r.F = S.F
    r.sub = S.sub
    O_as_submodule = SubModuleOfFreeModule(S.F, O)
    r.quo = isdefined(S,:quo) ? sum(S.quo,O_as_submodule) : O_as_submodule
    r.sum = sum(r.sub, r.quo)

    r.groebner_basis = Dict()

    r.incoming_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  #=function SubQuo(S::SubQuo, O::Vector{<:SubQuoElem})
    @assert all(x->x.parent === S, O)
    r = SubQuo(S.F, [x.repres for x in O])
    r.quo = S.quo
    r.sum = sum(r.sub, r.quo)
    return r
  end=#
  function SubQuo{R}(F::FreeMod{R}, s::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    r.sub = SubModuleOfFreeModule(F, s)
    r.sum = r.sub

    r.groebner_basis = Dict()

    r.incoming_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

    return r
  end
  function SubQuo{R}(F::FreeMod{R}, s::Singular.smodule, t::Singular.smodule) where {R}
    r = new{R}()
    r.F = F
    r.sub = SubModuleOfFreeModule(F, s)
    r.quo = SubModuleOfFreeModule(F, t)
    r.sum = sum(r.sub, r.quo)

    r.groebner_basis = Dict()

    r.incoming_morphisms = Vector{ModuleMap}()
    r.outgoing_morphisms = Vector{ModuleMap}()

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


mutable struct SubQuoHom{T1, T2} <: ModuleMap{T1, T2}
  matrix::MatElem
  header::Hecke.MapHeader
  im::Vector
  inverse_isomorphism::ModuleMap

  function SubQuoHom{T1,T2}(D::SubQuo, C::ModuleFP, im::Vector) where {T1,T2}
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2}()
    r.header = Hecke.MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    if C isa FreeMod
      r.im = Vector{FreeModElem}(im)
    elseif C isa SubQuo
      r.im = Vector{SubQuoElem}(im)
    else
      r.im = im
    end

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
    FreeModuleHom{T1, T2} <: ModuleMap{T1, T2} 

Data structure for morphisms where the domain is a free module (`FreeMod`).
`T1` and `T2` are the types of domain and codomain respectively.
`FreeModuleHom` is a subtype of `ModuleMap`.
When computed, the corresponding matrix (via `matrix()`) and inverse isomorphism
(in case there exists one) (via `inv()`) are cached.
"""
@attributes mutable struct FreeModuleHom{T1, T2} <: ModuleMap{T1, T2} 
  matrix::MatElem
  header::MapHeader
  inverse_isomorphism::ModuleMap

  # generate homomorphism of free modules from F to G where the vector a contains the images of
  # the generators of F
  function FreeModuleHom{T,S}(F::AbstractFreeMod{T}, G::S, a::Vector) where {T, S}
    @assert all(x->parent(x) === G, a)
    @assert length(a) == ngens(F)
    r = new{typeof(F), typeof(G)}()
    function im_func(x::AbstractFreeModElem)
      b = zero(G)
      for (i,v) = x.coords
        b += v*a[i]
      end
      return b
    end
    function pr_func(x)
      @assert parent(x) === G
      c = coordinates(repres(x), sub(G, a))
      return FreeModElem(c, F)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)

    return r
  end

  function FreeModuleHom{T,S}(F::AbstractFreeMod{T}, G::S, mat::MatElem{T}) where {T,S}
    @assert nrows(mat) == ngens(F)
    @assert ncols(mat) == ngens(G)
    if typeof(G) <: AbstractFreeMod
      hom = FreeModuleHom(F, G, [FreeModElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)])
    else
      hom = FreeModuleHom(F, G, [SubQuoElem(sparse_row(mat[i,:]), G) for i=1:ngens(F)])
    end
    hom.matrix = mat
    return hom
  end
end

struct FreeModuleHom_dec{T1, T2} <: ModuleMap{T1, T2}
  f::FreeModuleHom{T1,T2}
  header::MapHeader
  # TODO degree and homogeneity

  function FreeModuleHom_dec{T}(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, a::Vector) where {T}
    f = FreeModuleHom(F,G,a)
    r = new{typeof(F), typeof(G)}(f, f.header)
    return r
  end

  function FreeModuleHom_dec{T}(F::FreeMod_dec{T}, G::ModuleFP_dec{T}, mat::MatElem{T}) where {T}
    f = FreeModuleHom(F,G,mat)
    r = new{typeof(F), typeof(G)}(f, f.header)
    return r
  end
end
