###############################################################################
#
#  Elevations context
#
###############################################################################

@doc raw"""
A structure encoding a computational context to enumerate elevations of a given
pair `(L, f)`.

Given a sorted list of integers `L` and an integer `d`, we call a `d`-elevation
of `L` any set of indices for which the sum of the corresponding elements in `L`
is `d`. The set of all `d`-elevations of `L` is called the `d`-elevator of `L`.
For any `d`-elevations `e` of `L`, `d` is called the degree of `e`.

If `(L, f)` consists of a list `L` of objects of the same type and a $\mathbb Z$-
valued function `f` such that `f(L)` is a sorted list of integers, the `d`-elevator
of `(L, f)` is defined to be the `d`-elevator of `f(L)`.
"""
@attributes mutable struct ElevCtx{T, U}
  L::Vector{T}
  d::Int
  it::SubObjectIterator{PointVector{ZZRingElem}}
  f::U
  bounds::Tuple{Vector{Int}, Vector{Int}}

  function ElevCtx(L::Vector{T}, d::Int, it::SubObjectIterator{PointVector{ZZRingElem}}, f::U, bounds::Tuple{Vector{Int}, Vector{Int}}) where {T, U}
    z = new{T, U}()
    z.L = L
    z.d = d
    z.it = it
    z.f = f
    z.bounds = bounds
    return z
  end
end

###############################################################################
#
#  Representations
#
###############################################################################

### Representation ring

@doc raw"""
A structure encoding a representation ring of a finite group over a splitting field
of characteristic zero.

Here a representation ring is a "parent" for linear representations, used to store
the underlying group, a splitting field of characteristic zero, a set of generators
for the group, the character table, the irreducible characters and, iteratively,
representations affording some of the irreducible characters (only the one needed so far).
"""
@attributes mutable struct RepRing{S, T}
  field::S
  group::T
  gens::Vector{GAPGroupElem{T}}
  ct::Oscar.GAPGroupCharacterTable
  irr::Vector{Oscar.GAPGroupClassFunction}
  
  function RepRing{S, T}(F::S, E::T) where {S, T}
    ct = character_table(E)
    Irr = collect(ct)
    H = gens(E)
    RR = new{S, T}(F, E, H, ct, Irr)
  end

end

elem_type(::Type{RepRing{S, T}}) where {S, T} = LinRep{S, T, Oscar.elem_type(S)}

###############################################################################
### Linear representations of finite groups

@doc raw"""
A structure encoding a linear representation of a finite group.

Here a linear representation is an object of a representation ring, which is used as a parent
via caching. We also store the character. Representations are represented by a homomorphism
from the underlying group of their representation ring, to a matrix group whose matrices
have entries in the base field of the representation ring.
"""
@attributes mutable struct LinRep{S, T, U}
  rep_ring::RepRing{S, T}
  f::GAPGroupHomomorphism{T, MatrixGroup{U, AbstractAlgebra.Generic.MatSpaceElem{U}}}
  char::Oscar.GAPGroupClassFunction

  function LinRep{S, T, U}(RR::RepRing{S, T}, f::GAPGroupHomomorphism{T, MatrixGroup{U, AbstractAlgebra.Generic.MatSpaceElem{U}}}, char::Oscar.GAPGroupClassFunction) where {S, T, U}
    z = new{S, T, U}()
    z.rep_ring = RR
    z.f = f
    z.char = char
    return z
  end
end

###############################################################################
### Projective representation of finite groups

@doc raw"""
A structure encoding a projective representation of a finite group.

Here a projective representation of a group is represented by a linear lift along a
Schur cover of the group. We store the underlying Schur cover and the linear lift.
"""
@attributes mutable struct ProjRep{S, T, U, V}
  LR::LinRep{S, T, U}
  p::V

  function ProjRep{S, T, U, V}(LR::LinRep{S, T, U}, p::V) where {S, T, U, V}
    z = new{S, T, U, V}()
    z.LR = LR
    z.p = p
    return z
  end
end

###############################################################################
#
# Symmetric grassmannians
#
###############################################################################

@doc raw"""
Abstract supertype for different kinds of symmetric Grassmannians:

- `IsotGrass` for parametrising submodules of an isotypical group algebra module
  affording a given character;
- `CharGrass` for parametrising submodules of a group algebra module affording a
  given character;
- `DetGrass` for parametrising submodules of a group algebra module with a given
  determinant character;
- `InvGrass` for parametrising submodules of a group algebra module of a given
  dimension.
"""
abstract type SymGrass{S, T, U} end

@doc raw"""
A structure encoding a symmetric Grassmannian parametrising submodules of
a fixed isotypical group algebra module affording a given character.
"""
@attributes mutable struct IsotGrass{S, T, U} <: SymGrass{S, T, U}
  chi::Oscar.GAPGroupClassFunction
  rep_mod::LinRep{S, T, U}
  t::Int
  vs_struct::MapFromFunc{AbstractAlgebra.Generic.FreeModule{U}, AbstractAlgebra.Generic.MatSpace{U}}
  
  function IsotGrass(chi::Oscar.GAPGroupClassFunction, rep::LinRep{S, T, U}, dim::Int, vs::MapFromFunc{AbstractAlgebra.Generic.FreeModule{U}, AbstractAlgebra.Generic.MatSpace{U}}) where {S, T, U}
    z = new{S, T, U}()
    z.chi = chi
    z.rep_mod = rep
    z.t = dim
    z.vs_struct = vs
    return z
  end
end

@doc raw"""
A structure encoding a symmetric Grassmannian parametrising submodules of
a fixed group algebra module affording a given character.
"""
@attributes mutable struct CharGrass{S, T, U} <: SymGrass{S, T, U}
  dec::Vector{Oscar.GAPGroupClassFunction}
  constituents::Vector{IsotGrass{S, T, U}}

  function CharGrass(prod::Vector{IsotGrass{S, T, U}}) where {S, T, U}
    z = new{S, T, U}()
    dec = [M.chi for M in prod]
    z.dec = dec
    z.constituents = prod
    return z
  end
end

@doc raw"""
A structure encoding a symmetric Grassmannian parametrising submodules of a
fixed group algebra module having a given determinant character.
"""
@attributes mutable struct DetGrass{S, T, U} <: SymGrass{S, T, U}
  rep::LinRep{S, T, U}
  det_char::Oscar.GAPGroupClassFunction
  irr_comp::Vector{CharGrass{S, T, U}}
  d::Int

  function DetGrass(rep::LinRep{S, T, U}, dc::Oscar.GAPGroupClassFunction, ic::Vector{CharGrass{S, T, U}}, d::Int) where {S, T, U}
    return new{S, T, U}(rep, dc, ic, d)
  end
end

@doc raw"""
A structure encoding a symmetric Grassmannian parametrising submodules of a
fixed group algebra module of a given dimension.
"""
@attributes mutable struct InvGrass{S, T, U} <: SymGrass{S, T, U}
  rep::LinRep{S, T, U}
  irr_comp::Vector{CharGrass{S, T, U}}
  d::Int
  
  function InvGrass(r::LinRep{S, T, U}, ic::Vector{CharGrass{S, T, U}}, d::Int) where {S, T, U}
    return new{S, T, U}(r, ic, d)
  end
end

###############################################################################
#
# Symmetric Complete Intersections
#
###############################################################################

@doc raw"""
A structure encoding a parameter space for ideals generated by homogeneous
polynomials of the same degree, and fixed under a linear action of a group
on the projective space of coordinates.

We represent such a space by the underlying parameter space of group algebra modules
(generating the invariant ideals), together with an injection of the corresponding
homogeneous component of the polynomial algebra (seen as an abstract vector space).
"""
@attributes mutable struct SymInter{S, T, U, V}
  prep::ProjRep{S, T, U, V}
  para::CharGrass{S, T, U}
  j::MapFromFunc{AbstractAlgebra.Generic.FreeModule{U}, MPolyDecRing{U, AbstractAlgebra.Generic.MPolyRing{U}}}
  
  function SymInter(prep::ProjRep{S, T, U, V}, para::CharGrass{S, T, U}, j::MapFromFunc) where {S, T, U, V}
    z = new{S, T, U, V}(prep, para, j)
    return z
  end
end

