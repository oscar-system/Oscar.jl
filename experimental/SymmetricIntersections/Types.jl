import Oscar: elem_type

"""
We work always with finite groups and splitting fields of characteristic zero.
So by default, we take `QQAb` which is algebraically close of characteristic
zero.
"""

export CharacterGrassmannian,
       DeterminantGrassmannian,
       ElevCtx,
       InvariantGrassmannian,
       IsotypicalGrassmannian,
       LinRep,
       ParameterMap,
       ProjRep,
       RepRing,
       SymmetricGrassmannian,
       SymmetricIntersections

###############################################################################
#
#  Elevations context
#
###############################################################################

@attributes mutable struct ElevCtx{T, U}
  L::Vector{T}
  d::Int
  it::SubObjectIterator{PointVector{fmpz}}
  f::U
  bounds::Tuple{Vector{Int}, Vector{Int}}

  function ElevCtx(L::Vector{T}, d::Int, it::SubObjectIterator{PointVector{fmpz}}, f::U, bounds::Tuple{Vector{Int}, Vector{Int}}) where {T, U}
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

"""
Here a representation ring is a "parent" for linear representation, used to store
the underlying group, a splitting field (`QQAb` by default), a set of generators
for the group, the character table, the irreducible characters and, iteratively,
representations affording some of the irreducible characters (only the one needed
so far). There could be less accessors since irreducible characters store the character
tale which store the group and the field, by it is easily accessed like this.
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

elem_type(RR::RepRing{S, T}) where {S, T} = LinRep{S, T, Oscar.elem_type(S)}

elem_type(::Type{RepRing{S, T}}) where {S, T} = LinRep{S, T, Oscar.elem_type(S)}

###############################################################################
### Linear representations of finite groups

"""
Linear representations: an object of a representation ring, which is used as a parent
via caching. We also store the character. Representations are represented by an homomorphism
from the cached underlying group of their representation ring, to a matrix group whose matrices
have entries in the cached splitting field of the representation ring.
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

"""
Projective representations of a group are represented by a linear lift to a
Schur cover of the group. We store the underlying group and the linear lift.
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
export ProjRep

###############################################################################
#
# Symmetric grassmannians
#
###############################################################################

abstract type SymmetricGrassmannian{S, T, U} end

@attributes mutable struct IsotypicalGrassmannian{S, T, U} <: SymmetricGrassmannian{S, T, U}
  chi::Oscar.GAPGroupClassFunction
  rep_mod::LinRep{S, T, U}
  t::Int
  vs_struct::MapFromFunc{AbstractAlgebra.Generic.FreeModule{U}, AbstractAlgebra.Generic.MatSpace{U}}
  
  function IsotypicalGrassmannian(chi::Oscar.GAPGroupClassFunction, rep::LinRep{S, T, U}, dim::Int, vs::MapFromFunc{AbstractAlgebra.Generic.FreeModule{U}, AbstractAlgebra.Generic.MatSpace{U}}) where {S, T, U}
    z = new{S, T, U}()
    z.chi = chi
    z.rep_mod = rep
    z.t = dim
    z.vs_struct = vs
    return z
  end
end

@attributes mutable struct CharacterGrassmannian{S, T, U} <: SymmetricGrassmannian{S, T, U}
  dec::Vector{Oscar.GAPGroupClassFunction}
  constituents::Vector{IsotypicalGrassmannian{S, T, U}}

  function CharacterGrassmannian(prod::Vector{IsotypicalGrassmannian{S, T, U}}) where {S, T, U}
    z = new{S, T, U}()
    dec = [M.chi for M in prod]
    z.dec = dec
    z.constituents = prod
    return z
  end
end

@attributes mutable struct DeterminantGrassmannian{S, T, U} <: SymmetricGrassmannian{S, T, U}
  rep::LinRep{S, T, U}
  det_char::Oscar.GAPGroupClassFunction
  irr_comp::Vector{CharacterGrassmannian{S, T, U}}
  d::Int

  function DeterminantGrassmannian(rep::LinRep{S, T, U}, dc::Oscar.GAPGroupClassFunction, ic::Vector{CharacterGrassmannian{S, T, U}}, d::Int) where {S, T, U}
    return new{S, T, U}(rep, dc, ic, d)
  end
end

@attributes mutable struct InvariantGrassmannian{S, T, U} <: SymmetricGrassmannian{S, T, U}
  rep::LinRep{S, T, U}
  irr_comp::Vector{CharacterGrassmannian{S, T, U}}
  d::Int
  
  function InvariantGrassmannian(r::LinRep{S, T, U}, ic::Vector{CharacterGrassmannian{S, T, U}}, d::Int) where {S, T, U}
    return new{S, T, U}(r, ic, d)
  end
end

###############################################################################
#
# Symmetric Complete Intersections
#
###############################################################################

@attributes mutable struct SymmetricIntersections{S, T, U, V}
  prep::ProjRep{S, T, U, V}
  para::CharacterGrassmannian{S, T, U}
  j::MapFromFunc{AbstractAlgebra.Generic.FreeModule{U}, MPolyRing_dec{U, AbstractAlgebra.Generic.MPolyRing{U}}}
  
  function SymmetricIntersections(prep::ProjRep{S, T, U, V}, para::CharacterGrassmannian{S, T, U}, j::MapFromFunc) where {S, T, U, V}
    z = new{S, T, U, V}(prep, para, j)
    return z
  end
end

