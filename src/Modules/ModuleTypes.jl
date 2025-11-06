import AbstractAlgebra.WeakKeyIdDict

# The following type union gathers all types for elements for which 
# we expect the `ModuleFP` framework to be functional and/or working.
# It can gradually be extended, but should not be considered to be 
# visible or accessible to the outside world. 
const AdmissibleModuleFPRingElem = Union{RingElem, PBWAlgQuoElem, PBWAlgElem}
const AdmissibleModuleFPRing = Union{Ring, PBWAlgQuo, PBWAlgRing}

@doc raw"""
    ModuleFP{T}

The abstract supertype of all finitely presented modules.
The type variable `T` refers to the type of the elements of the base ring.
"""
abstract type ModuleFP{T <: AdmissibleModuleFPRingElem} <: AbstractAlgebra.Module{T} end

@doc raw"""
    AbstractFreeMod{T} <: ModuleFP{T}

The abstract supertype of all finitely generated free modules.
"""
abstract type AbstractFreeMod{T <: AdmissibleModuleFPRingElem} <: ModuleFP{T} end

@doc raw"""
    AbstractSubQuo{T} <: ModuleFP{T}

The abstract supertype of all finitely presented subquotient modules.
"""
abstract type AbstractSubQuo{T <: AdmissibleModuleFPRingElem} <: ModuleFP{T} end


@doc raw"""
    ModuleFPElem{T} <: ModuleElem{T}

The abstract supertype of all elements of finitely presented modules.
"""
abstract type ModuleFPElem{T <: AdmissibleModuleFPRingElem} <: ModuleElem{T} end

@doc raw"""
    AbstractFreeModElem{T} <: ModuleFPElem{T}

The abstract supertype of all elements of finitely generated free modules.
"""
abstract type AbstractFreeModElem{T <: AdmissibleModuleFPRingElem} <: ModuleFPElem{T} end

@doc raw"""
    AbstractSubQuoElem{T} <: ModuleFPElem{T}

The abstract supertype of all elements of subquotient modules.
"""
abstract type AbstractSubQuoElem{T <: AdmissibleModuleFPRingElem} <: ModuleFPElem{T} end

abstract type ModuleFPHomDummy end

@doc raw"""
    ModuleFPHom{T1, T2, RingMapType}

The abstract supertype for morphisms of finitely presented modules over multivariate polynomial rings .
`T1` and `T2` are the types of domain and codomain respectively.
`RingMapType` is a type for a homomorphism of rings ``f : R → S`` whenever the 
`base_ring` ``R`` of the domain is different from the `base_ring` ``S`` of the codomain 
and the codomain is considered as an ``R``-module via ``f``. 
In case there is no base change, this parameter is set to `Nothing`.
"""
abstract type ModuleFPHom{T1, T2, RingMapType} <: Map{T1, T2, Hecke.HeckeMap, ModuleFPHomDummy} end

parent(f::ModuleFPHom) = Hecke.MapParent(domain(f), codomain(f), "homomorphisms")

@doc raw"""
    FreeMod{T <: RingElem} <: AbstractFreeMod{T}

The type of free modules.
Free modules are determined by their base ring, the rank and the names of 
the (standard) generators.
Moreover, canonical incoming and outgoing morphisms are stored if the corresponding
option is set in suitable functions.
`FreeMod{T}` is a subtype of `AbstractFreeMod{T}`.
"""
@attributes mutable struct FreeMod{T <: AdmissibleModuleFPRingElem} <: AbstractFreeMod{T}
  R::NCRing
  n::Int
  S::Union{Function, Vector{Symbol}} # The symbols for printing. This is either a 
                                     # ready-made list, or a function to provide them 
                                     # in a lazy way. 
  d::Union{Vector{FinGenAbGroupElem}, Nothing}
  default_ordering::ModuleOrdering

  # We register the incoming and outgoing natural morphisms.
  # This must be done in a way that objects can be collected by the 
  # garbage collector. In particular, we can not store the actual 
  # map as the value for a specific key (domain or codomain depending 
  # on whether the map is incoming or outgoing), because then the 
  # value has a reference to the key and thus the pair will never be 
  # deleted. 
  #
  # Instead, we store a sparse matrix which allows us to reconstruct 
  # the map and potentially a change of rings. This allows us to 
  # reconstruct the map on request (which should be of relatively 
  # low cost).
  incoming::WeakKeyIdDict{<:ModuleFP, <:Tuple{<:SMat, <:Any}}
  outgoing::WeakKeyIdDict{<:ModuleFP, <:Tuple{<:SMat, <:Any}}

  function FreeMod{T}(n::Int, R::AdmissibleModuleFPRing, S::Vector{Symbol}) where T <: AdmissibleModuleFPRingElem
    r = new{elem_type(R)}()
    r.n = n
    r.R = R
    r.S = S
    r.d = nothing

    r.incoming = WeakKeyIdDict{ModuleFP, Tuple{SMat, Any}}()
    r.outgoing = WeakKeyIdDict{ModuleFP, Tuple{SMat, Any}}()
    return r
  end
  
  function FreeMod{T}(n::Int, R::AdmissibleModuleFPRing, symbol_fun::Function) where T <: AdmissibleModuleFPRingElem
    r = new{elem_type(R)}(R, n, symbol_fun, nothing)

    r.incoming = WeakKeyIdDict{ModuleFP, Tuple{SMat, Any}}()
    r.outgoing = WeakKeyIdDict{ModuleFP, Tuple{SMat, Any}}()
    return r
  end
end

@doc raw"""
    FreeModElem{T <: AdmissibleModuleFPRingElem} <: AbstractFreeModElem{T}

The type of free module elements. An element of a free module $F$ is given by a sparse row (`SRow`)
which specifies its coordinates with respect to the basis of standard unit vectors of $F$.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> F = free_module(R, 3)
Free module of rank 3 over R

julia> f = F(sparse_row(R, [(1,x),(3,y)]))
x*e[1] + y*e[3]

julia> typeof(f)
FreeModElem{QQMPolyRingElem}

julia> g = x*F[1] + y*F[3]
x*e[1] + y*e[3]

julia> f == g
true
```
"""
mutable struct FreeModElem{T <: AdmissibleModuleFPRingElem} <: AbstractFreeModElem{T}
  coords::SRow{T} # also usable via coeffs()
  parent::FreeMod{T}
  d::Union{FinGenAbGroupElem, Nothing}

  function FreeModElem{T}(coords::SRow{T}, parent::FreeMod{T}) where T
      r = new{T}(coords, parent, nothing)
      return r
  end
end

@doc raw"""
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
@attributes mutable struct ModuleGens{T <: AdmissibleModuleFPRingElem} # T is the type of the elements of the ground ring.
  O::Vector{FreeModElem{T}}
  S::Singular.smodule
  F::FreeMod{T}
  SF::Singular.FreeMod

  isGB::Bool
  is_reduced::Bool
  ordering::ModuleOrdering
  quo_GB::ModuleGens{T} # Pointer to the quotient GB when having a relative GB

  function ModuleGens{T}(O::Vector{<:FreeModElem}, F::FreeMod{T}) where {T <: AdmissibleModuleFPRingElem}
    r = new{T}()
    r.O = O
    r.F = F
    return r
  end

  # ModuleGens from an Array of Oscar free module elements, specifying the free module 
  # and Singular free module, only useful indirectly
  function ModuleGens{T}(O::Vector{<:FreeModElem}, F::FreeMod{T}, SF::Singular.FreeMod) where {T <: AdmissibleModuleFPRingElem}
    r = new{T}()
    r.O = O
    r.SF = SF
    r.F = F
    return r
  end

  # ModuleGens from a Singular submodule
  function ModuleGens{S}(F::FreeMod{S}, s::Singular.smodule) where {S <: AdmissibleModuleFPRingElem} # FreeMod is necessary due to type S
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

@doc raw"""
    SubModuleOfFreeModule{T <: AdmissibleModuleFPRingElem} <: ModuleFP{T}

Data structure for submodules of free modules. `SubModuleOfFreeModule` shouldn't be
used by the end user.
When computed, a standard basis (computed via `standard_basis()`) and generating matrix (that is the rows of the matrix
generate the submodule) (computed via `generator_matrix()`) are cached.
"""
@attributes mutable struct SubModuleOfFreeModule{T <: AdmissibleModuleFPRingElem} <: ModuleFP{T}
  F::FreeMod{T}
  groebner_basis::Dict{ModuleOrdering, ModuleGens{T}}
  gens::ModuleGens{T}
  default_ordering::ModuleOrdering{FreeMod{T}}
  any_gb::ModuleGens{T} # A field to store the first groebner basis ever computed.
                        # Lookups in the above dictionary is tentatively expensive. 
                        # So this field stores any gb for cases where the actual 
                        # ordering does not matter. Then this field here can be used. 
  any_gb_with_transition::ModuleGens{T} # The same but for one with transition matrix
  matrix::MatElem

  function SubModuleOfFreeModule{R}(F::FreeMod{R}) where {R}
    # this does not construct a valid SubModuleOfFreeModule
    return new{R}(F, Dict{ModuleOrdering, ModuleGens{R}}())
  end
end


@doc raw"""
    SubquoModule{T <: AdmissibleModuleFPRingElem} <: AbstractSubQuo{T}

The type of subquotient modules.
A subquotient module $M$ is a module where $M = A + B / B$ where $A$ and $B$ are 
submodules of a free module.
$A$, $B$ and $A+B$ (they have type `SubModuleOfFreeModule`) as well as the embedding
free module are stored. 
One can construct ordinary submodules of free modules by not giving $B$.
Moreover, canonical incoming and outgoing morphisms are stored if the corresponding
option is set in suitable functions.
`SubquoModule{T}` is a subtype of `ModuleFP{T}`.
"""
@attributes mutable struct SubquoModule{T <: AdmissibleModuleFPRingElem} <: AbstractSubQuo{T}
  #meant to represent sub+ quo mod quo - as lazy as possible
  F::FreeMod{T}

  groebner_basis::Dict{ModuleOrdering, ModuleGens{T}}

  incoming::WeakKeyIdDict{<:ModuleFP, <:Tuple{<:SMat, <:Any}}
  outgoing::WeakKeyIdDict{<:ModuleFP, <:Tuple{<:SMat, <:Any}}

  sub::SubModuleOfFreeModule{T}
  quo::SubModuleOfFreeModule{T}
  sum::SubModuleOfFreeModule{T}

  function SubquoModule{R}(F::FreeMod{R}) where {R}
    # this does not construct a valid subquotient
    return new{R}(F,
                  Dict(),  # groebner_basis
                  WeakKeyIdDict{ModuleFP, Tuple{SMat, Any}}(), # incoming
                  WeakKeyIdDict{ModuleFP, Tuple{SMat, Any}}(), # outgoing
                  )
  end
end


@doc raw"""
    SubquoModuleElem{T <: AdmissibleModuleFPRingElem} <: AbstractSubQuoElem{T} 

The type of subquotient elements. An element $f$ of a subquotient $M$ over the ring $R$
is given by a sparse row (`SRow`) which specifies the coefficients of an $R$-linear 
combination of the generators of $M$ which defines $f$.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[x, y, z])

julia> A = R[x; y]
[x]
[y]

julia> B = R[x^2; x*y; y^2; z^4]
[x^2]
[x*y]
[y^2]
[z^4]

julia> M = SubquoModule(A, B)
Subquotient of submodule with 2 generators
  1: x*e[1]
  2: y*e[1]
by submodule with 4 generators
  1: x^2*e[1]
  2: x*y*e[1]
  3: y^2*e[1]
  4: z^4*e[1]

julia> f = SubquoModuleElem(sparse_row(R, [(1,z),(2,one(R))]),M)
(x*z + y)*e[1]

julia> g = z*M[1] + one(R)*M[2]
(x*z + y)*e[1]

julia> typeof(g)
SubquoModuleElem{QQMPolyRingElem}

julia> f == g
true
```
"""
mutable struct SubquoModuleElem{T <: AdmissibleModuleFPRingElem} <: AbstractSubQuoElem{T} 
  parent::SubquoModule{T}
  coeffs::SRow{T} # Need not be defined! Use `coordinates` as getter
  repres::FreeModElem{T} # Need not be defined! Use `repres` as getter
  is_reduced::Bool # `false` by default. Will be set by `simplify` and `simplify!`.

  function SubquoModuleElem{R}(v::SRow{R}, SQ::SubquoModule) where {R <: AdmissibleModuleFPRingElem}
    @assert length(v) <= ngens(SQ.sub)
    if isempty(v)
      r = new{R}(SQ, v, zero(SQ.F))
      return r
    end
    r = new{R}(SQ, v)
    #, sum(a*SQ.sub[i] for (i, a) in v; init = zero(SQ.sub)), SQ)
    r.is_reduced = false
    return r
  end

  function SubquoModuleElem{R}(a::FreeModElem{R}, SQ::SubquoModule; is_reduced::Bool=false) where {R <: AdmissibleModuleFPRingElem}
    @assert a.parent === SQ.F
    r = new{R}(SQ)
    r.repres = a
    r.is_reduced = is_reduced
    return r
  end
end


mutable struct SubQuoHom{
    T1<:AbstractSubQuo, 
    T2<:ModuleFP, 
    RingMapType<:Any
  } <: ModuleFPHom{T1, T2, RingMapType}
  matrix::MatElem
  header::MapHeader{T1,T2}
  im::Vector # The images of the generators; use `images_of_generators` as a getter.
  inverse_isomorphism::ModuleFPHom
  ring_map::RingMapType
  d::FinGenAbGroupElem
  generators_map_to_generators::Union{Bool, Nothing} # A flag to allow for shortcut in evaluation;
                                                     # value `nothing` by default and to be set manually.

  # Constructors for maps without change of base ring
  function SubQuoHom{T1,T2,RingMapType}(D::SubquoModule, C::FreeMod, im::Vector;
                                        check::Bool=true
                                       ) where {T1,T2,RingMapType}
    ###@assert is_graded(D) == is_graded(C)
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)
    @check true # Will throw if checks are not supposed to be enabled

    r = new{T1, T2, Nothing}()
    r.header = MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{elem_type(C)}(im)
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

  function SubQuoHom{T1,T2,RingMapType}(D::SubquoModule, C::SubquoModule, im::Vector;
                                        check::Bool=true
                                       ) where {T1,T2,RingMapType}
    ###@assert is_graded(D) == is_graded(C)
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, Nothing}()
    r.header = MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{elem_type(C)}(im)
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

  function SubQuoHom{T1,T2,RingMapType}(D::SubquoModule, C::ModuleFP, im::Vector;
                                        check::Bool=true
                                       ) where {T1,T2,RingMapType}
    ###@assert is_graded(D) == is_graded(C)
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, Nothing}()
    r.header = MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{elem_type(C)}(im)
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

  # Constructors for maps with change of base ring
  function SubQuoHom{T1,T2,RingMapType}(
      D::SubquoModule, 
      C::FreeMod, 
      im::Vector, 
      h::RingMapType;
      check::Bool=true
    ) where {T1,T2,RingMapType}
    ###@assert is_graded(D) == is_graded(C)
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, RingMapType}()
    r.header = MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{elem_type(C)}(im)
    r.ring_map = h
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

  function SubQuoHom{T1,T2,RingMapType}(
      D::SubquoModule, 
      C::SubquoModule, 
      im::Vector, 
      h::RingMapType;
      check::Bool=true
    ) where {T1,T2,RingMapType}
    ###@assert is_graded(D) == is_graded(C)
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, RingMapType}()
    r.header = MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{elem_type(C)}(im)
    r.ring_map = h
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

  function SubQuoHom{T1,T2,RingMapType}(
      D::SubquoModule, 
      C::ModuleFP, 
      im::Vector, 
      h::RingMapType;
      check::Bool=true
    ) where {T1,T2,RingMapType}
    ###@assert is_graded(D) == is_graded(C)
    @assert length(im) == ngens(D)
    @assert all(x-> parent(x) === C, im)

    r = new{T1, T2, RingMapType}()
    r.header = MapHeader(D, C)
    r.header.image = x->image(r, x)
    r.header.preimage = x->preimage(r, x)
    r.im = Vector{elem_type(C)}(im)
    r.ring_map = h
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

end


###############################################################################
# Graded modules
###############################################################################
const CRing_dec = Union{MPolyDecRing, MPolyQuoRing{<:Oscar.MPolyDecRingElem}}
const CRingElem_dec = Union{MPolyDecRingElem, MPolyQuoRingElem{<:Oscar.MPolyDecRingElem}}
#TODO: other name for CRing_dec -> which?

@doc raw"""
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
  d::Vector{FinGenAbGroupElem}

  function FreeMod_dec{T}(F::FreeMod, d::Vector{FinGenAbGroupElem}) where T <: CRingElem_dec
    @assert length(d) == rank(F)
    r = new{elem_type(base_ring(F))}(F, d)
    return r
  end

  function FreeMod_dec{T}(R::CRing_dec,S::Vector{Symbol},d::Vector{FinGenAbGroupElem}) where T <: CRingElem_dec
    r = new{elem_type(R)}()
    r.F = FreeMod{T}(length(d),R,S)
    r.d = d
    return r
  end
end

@doc raw"""
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





const ModuleFP_dec{T} = Union{FreeMod_dec{T}} # SubquoDecModule{T} will be included
const ModuleFPElem_dec{T} = Union{FreeModElem_dec{T}} # SubquoDecModuleElem{T} will be included


@doc raw"""
    FreeModuleHom{T1, T2, RingMapType} <: ModuleFPHom{T1, T2, RingMapType} 

Data structure for morphisms where the domain is a free module (`FreeMod`).
`T1` and `T2` are the types of domain and codomain respectively.
`FreeModuleHom` is a subtype of `ModuleFPHom`.
When computed, the corresponding matrix (via `matrix()`) and inverse isomorphism
(in case there exists one) (via `inv()`) are cached.
"""
@attributes mutable struct FreeModuleHom{
    T1 <: AbstractFreeMod,
    T2 <: ModuleFP,
    RingMapType <: Any} <: ModuleFPHom{T1, T2, RingMapType} 
  header::MapHeader{T1, T2}
  ring_map::RingMapType
  d::FinGenAbGroupElem
  imgs_of_gens::Vector # stored here for easy evaluation; use `images_of_generators` as getter
  
  matrix::MatElem
  inverse_isomorphism::ModuleFPHom
  generators_map_to_generators::Union{Bool, Nothing} # A flag to allow for shortcut in evaluation
                                                     # of the map; `nothing` by default and to be 
                                                     # set manually.

  # generate homomorphism of free modules from F to G where the vector a contains the images of
  # the generators of F
  function FreeModuleHom(
      F::AbstractFreeMod, G::S, a::Vector{ModuleElemType};
      check::Bool=true
    ) where {S<:ModuleFP, ModuleElemType<:ModuleFPElem}
    ###@assert is_graded(F) == is_graded(G)
    @assert all(x->parent(x) === G, a)
    @assert length(a) == ngens(F)
    r = new{typeof(F), typeof(G), Nothing}()
    a=Vector{elem_type(G)}(a)
    image_module = sub_object(G, a)
    function im_func(x::AbstractFreeModElem)
     # The lines below were an attempt to speed up mapping. 
     # However, it turns out that checking the equality is more 
     # expensive in average than the gain for actual mappings. 
     # Apparently, maps are likely to be used just once, or only 
     # few times. 
     # But the flag can (and probably should) be set by the constructors
     # of maps whenever applicable.
     #if r.generators_map_to_generators === nothing
     #  r.generators_map_to_generators = images_of_generators(r) == gens(codomain(r))
     #end
      r.generators_map_to_generators === true && return codomain(r)(coordinates(x))
      # The following code is a slight manual speedup of the following original line:
      # return sum(b*a[i] for (i, b) in coordinates(x); init=zero(codomain(r)))
      cod_ring = base_ring(G)
      res_coord = sparse_row(cod_ring)
      for (i, b) in coordinates(x)
        Hecke.add_scaled_row!(coordinates(a[i]), res_coord, b)
      end
      return G(res_coord)
    end
    function pr_func(x)
      @assert parent(x) === G
      r.generators_map_to_generators === true && return FreeModElem(coordinates(simplify!(x)), F)
      #c = coordinates(repres(simplify!(x)), image_module)
      c = coordinates(repres(x), image_module)
      return FreeModElem(c, F)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)
    r.imgs_of_gens = Vector{elem_type(G)}(a)
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

  # We need to introduce a separate constructor for decorated modules here
  # because some of the generic functionality (sub_object) is not there yet. 
  function FreeModuleHom(
      F::FreeMod_dec, G::S, a::Vector{ModuleElemType};
      check::Bool=true
    ) where {S<:ModuleFP, ModuleElemType<:ModuleFPElem}
    ###@assert is_graded(F) == is_graded(G)
    @assert all(x->parent(x) === G, a)
    @assert length(a) == ngens(F)
    r = new{typeof(F), typeof(G), Nothing}()
    a=Vector{elem_type(G)}(a)
    function im_func(x::AbstractFreeModElem)
     # The lines below were an attempt to speed up mapping. 
     # However, it turns out that checking the equality is more 
     # expensive in average than the gain for actual mappings. 
     # Apparently, maps are likely to be used just once, or only 
     # few times. 
     # But the flag can (and probably should) be set by the constructors
     # of maps whenever applicable.
     #if r.generators_map_to_generators === nothing
     #  r.generators_map_to_generators = images_of_generators(r) == gens(codomain(r))
     #end
      r.generators_map_to_generators === true && return codomain(r)(coordinates(x))
      return sum(b*a[i] for (i, b) in coordinates(x); init=zero(codomain(r)))
    end
    function pr_func(x)
      @assert parent(x) === G
      c = coordinates(repres(simplify!(x)), sub_object(G, a))
      return FreeModElem(c, F)
    end
    r.header = MapHeader{typeof(F), typeof(G)}(F, G, im_func, pr_func)
    r.imgs_of_gens = Vector{elem_type(G)}(a)
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

  function FreeModuleHom(
      F::AbstractFreeMod, G::T2, a::Vector{ModuleElemType}, h::RingMapType;
      check::Bool=true
    ) where {T2, ModuleElemType<:ModuleFPElem, RingMapType}
    ###@assert is_graded(F) == is_graded(G)
    @assert all(x->parent(x) === G, a)
    @assert length(a) == ngens(F)
    @assert h(one(base_ring(F))) == one(base_ring(G))
    r = new{typeof(F), T2, RingMapType}()
    a=Vector{elem_type(G)}(a)
    image_module = sub_object(G, a)
    function im_func(x::AbstractFreeModElem)
      iszero(x) && return zero(codomain(r))
      # See the above comment
     #if r.generators_map_to_generators === nothing
     #  r.generators_map_to_generators = images_of_generators(r) == gens(codomain(r))
     #end
      r.generators_map_to_generators === true && return codomain(r)(map_entries(h, coordinates(x)))
      # in-place version of 
      # return sum(h(b)*a[i] for (i, b) in coordinates(x); init=zero(codomain(r)))
      S = base_ring(G)
      pre_res = sparse_row(S)
      for (i, b) in coordinates(x)
        pre_res = Hecke.add_scaled_row!(coordinates(a[i]), pre_res, h(b))
      end
      return FreeModElem(pre_res, G)
    end
    function pr_func(x)
      @assert parent(x) === G
      #r.generators_map_to_generators === true && return FreeModElem(map_entries(x->preimage(h, x), coordinates(simplify!(x))), F)
      c = coordinates(repres(x), image_module)
      cc = map_entries(x->preimage(h, x), c)
      return FreeModElem(cc, F)
    end
    r.header = MapHeader{typeof(F), T2}(F, G, im_func, pr_func)
    r.ring_map = h
    r.imgs_of_gens = Vector{elem_type(G)}(a)
    r.generators_map_to_generators = nothing
    return set_grading(r; check)
  end

end

# Further constructors taking matrices as input
function FreeModuleHom(
    F::AbstractFreeMod{T}, G::S, mat::MatElem{T};
    check::Bool=true
  ) where {T<:RingElem,S<:AbstractFreeMod}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  hom = FreeModuleHom(F, G, [FreeModElem(sparse_row(mat[i:i,:]), G) for i=1:ngens(F)])
  hom.matrix = mat
  return hom
end

function FreeModuleHom(
    F::AbstractFreeMod{T}, G::S, mat::MatElem{T};
    check::Bool=true
  ) where {T<:RingElem, S<:ModuleFP}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  hom = FreeModuleHom(F, G, [SubquoModuleElem(sparse_row(mat[i:i,:]), G) for i=1:ngens(F)]; check)
  hom.matrix = mat
  return hom
end

function FreeModuleHom(
    F::AbstractFreeMod, G::S, mat::MatElem, h::RingMapType;
    check::Bool=true
  ) where {S<:AbstractFreeMod, RingMapType}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  @assert base_ring(mat) === base_ring(G)
  hom = FreeModuleHom(F, G, [FreeModElem(sparse_row(mat[i:i,:]), G) for i=1:ngens(F)], h; check)
  hom.matrix = mat
  return hom
end

function FreeModuleHom(
    F::AbstractFreeMod, G::S, mat::MatElem, h::RingMapType;
    check::Bool=true
  ) where {S<:ModuleFP, RingMapType}
  @assert nrows(mat) == ngens(F)
  @assert ncols(mat) == ngens(G)
  @assert base_ring(mat) === base_ring(G)
  hom = FreeModuleHom(F, G, [SubquoModuleElem(sparse_row(mat[i:i,:]), G) for i=1:ngens(F)], h; check)
  hom.matrix = mat
  return hom
end

struct FreeModuleHom_dec{
    T1 <: AbstractFreeMod,
    T2 <: ModuleFP,
    RingMapType <: Any} <: ModuleFPHom{T1, T2, RingMapType}
  f::FreeModuleHom{T1,T2, RingMapType}
  header::MapHeader{T1,T2}
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

@doc raw"""
    FreeResolution{T}

Data structure for free resolutions.
"""
mutable struct FreeResolution{T}
    C::Hecke.ComplexOfMorphisms

    function FreeResolution(C::Hecke.ComplexOfMorphisms{T}) where {T}
        FR = new{T}()
        FR.C = C

        return FR
    end
end

Base.getindex(FR::FreeResolution, i::Int) = FR.C[i]

mutable struct BettiTable
  B::Dict{Tuple{Int, Any}, Int}
  project::Union{FinGenAbGroupElem, Nothing}
  reverse_direction::Bool
  function BettiTable(B::Dict{Tuple{Int, Any}, Int}; project::Union{FinGenAbGroupElem, Nothing}=nothing, reverse_direction::Bool=false)
      return new(B, project, reverse_direction)
  end
end
