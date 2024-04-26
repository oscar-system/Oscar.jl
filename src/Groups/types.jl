@doc raw"""
    GAPGroup <: AbstractAlgebra.Group

Each object of the abstract type `GAPGroup` stores a group object from
the GAP system,
and thus can delegate questions about this object to GAP.

For expert usage, you can extract the underlying GAP object via `GapObj`,
i.e., if `G` is a `GAPGroup`, then `GapObj(G)` is the `GapObj` underlying `G`.

Concrete subtypes of `GAPGroup` are `PermGroup`, `FPGroup`, `PcGroup`,
`SubPcGroup`, and `MatrixGroup`.
"""
abstract type GAPGroup <: AbstractAlgebra.Group end

## `GapGroup` to GAP group
GAP.julia_to_gap(obj::GAPGroup) = obj.X

@doc raw"""
    GAPGroupElem <: AbstractAlgebra.GroupElem

Each object of the abstract type `GAPGroupElem` stores a group element
object from the GAP system,
and thus can delegate questions about this object to GAP.

For expert usage, you can extract the underlying GAP object via `GapObj`,
i.e., if `g` is a `GAPGroupElem`, then `GapObj(g)` is the `GapObj` underlying `g`.
"""
abstract type GAPGroupElem{T<:GAPGroup} <: AbstractAlgebra.GroupElem end

## `GapGroupElem` to GAP group element
GAP.julia_to_gap(obj::GAPGroupElem) = obj.X

@doc raw"""
    BasicGAPGroupElem{T<:GAPGroup} <: GAPGroupElem{T}

The type `BasicGAPGroupElem` gathers all types of group elements
described *only* by an underlying GAP object.

If $x$ is an element of the group `G` of type `T`,
then the type of $x$ is `BasicGAPGroupElem{T}`.
"""
struct BasicGAPGroupElem{T<:GAPGroup} <: GAPGroupElem{T}
   parent::T
   X::GapObj
end

function Base.deepcopy_internal(x::BasicGAPGroupElem, dict::IdDict)
  X = Base.deepcopy_internal(GapObj(x), dict)
  return BasicGAPGroupElem(x.parent, X)
end

Base.hash(x::GAPGroup, h::UInt) = h # FIXME
Base.hash(x::GAPGroupElem, h::UInt) = h # FIXME


"""
    PermGroup

Groups of permutations.
Every group of this type is a subgroup of Sym(n) for some n.

# Examples
- `symmetric_group(n::Int)`: the symmetric group Sym(n)
- `alternating_group(n::Int)`: the alternating group Alt(n)
- subgroups of Sym(n)
- `dihedral_group(PermGroup, n::Int)`:
  the dihedral group of order `n` as a group of permutations.
  Same holds replacing `dihedral_group` by `quaternion_group`

If `G` is a permutation group and `x` is a permutation,
`G(x)` returns a permutation `x` with parent `G`;
an exception is thrown if `x` does not embed into `G`.
```jldoctest
julia> G=symmetric_group(5)
Sym(5)

julia> x=cperm([1,2,3])
(1,2,3)

julia> parent(x)
Sym(3)

julia> y=G(x)
(1,2,3)

julia> parent(y)
Sym(5)
```

If `G` is a permutation group and `x` is a vector of integers,
`G(x)` returns a [`PermGroupElem`](@ref) with parent `G`;
an exception is thrown if the element does not embed into `G`.

# Examples
```jldoctest
julia> G = symmetric_group(6)
Sym(6)

julia> x = G([2,4,6,1,3,5])
(1,2,4)(3,6,5)

julia> parent(x)
Sym(6)
```
"""
@attributes mutable struct PermGroup <: GAPGroup
   X::GapObj
   deg::Int64       # G < Sym(deg)
   
   function PermGroup(G::GapObj)
     @assert GAPWrap.IsPermGroup(G)
     n = GAPWrap.LargestMovedPoint(G)::Int
     if n == 0
       # We support only positive degrees.
       # (`symmetric_group(0)` yields an error,
       # and `symmetric_group(1)` yields a GAP group with `n == 0`.)
       n = 1
     end
     z = new(G, n)
     return z
   end
   
   function PermGroup(G::GapObj, deg::Int)
     @assert GAPWrap.IsPermGroup(G) && deg > 0 && deg >= GAPWrap.LargestMovedPoint(G)::Int
     z = new(G, deg)
     return z
   end
end

permutation_group(G::GapObj) = PermGroup(G)
permutation_group(G::GapObj, deg::Int) = PermGroup(G, deg)

"""
    PermGroupElem

Element of a group of permutations.
It is displayed as product of disjoint cycles.
# Assumptions:
- for `x`,`y` in Sym(n), the product `xy` is read from left to right;
- for `x` in Sym(n) and `i` in {1,...,n}, `i^x` and `x(i)` return the image of `i` under the action of `x`.
"""
const PermGroupElem = BasicGAPGroupElem{PermGroup}


"""
    PcGroup

Polycyclic group, a group that is defined by a finite presentation
of a special kind, a so-called polycyclic presentation.
Contrary to arbitrary finitely presented groups
(see [Finitely presented groups](@ref)),
this presentation allows for efficient computations with the group elements.

For a group `G` of type `PcGroup`, the elements in `gens(G)` satisfy the
relators of the underlying presentation.

Functions that compute subgroups of `G` return groups of type `SubPcGroup`.

# Examples
- `cyclic_group(n::Int)`: cyclic group of order `n`
- `abelian_group(PcGroup, v::Vector{Int})`:
  direct product of cyclic groups of the orders
  `v[1]`, `v[2]`, ..., `v[length(v)]`
"""
@attributes mutable struct PcGroup <: GAPGroup
  X::GapObj

  function PcGroup(G::GapObj)
    _is_full_pc_group(G) && return new(G)
    # Switch to a full pcp or pc group.
    if GAP.Globals.IsPcpGroup(G)::Bool
      return new(GAP.Globals.PcpGroupByPcp(GAP.Globals.Pcp(G)::GapObj)::GapObj)
    elseif GAPWrap.IsPcGroup(G)::Bool
      return new(GAP.Globals.PcGroupWithPcgs(GAP.Globals.Pcgs(G)::GapObj)::GapObj)
    end
    throw(ArgumentError("G must be in IsPcGroup or IsPcpGroup"))
  end
end

pc_group(G::GapObj) = PcGroup(G)

# Return `true` if the generators of `G` fit to those of its pc presentation.
function _is_full_pc_group(G::GapObj)
  return GAP.Globals.GroupGeneratorsDefinePresentation(G)::Bool
end
#T the code for PcpGroups is too expensive!


"""
    PcGroupElem

Element of a polycyclic group.

The generators of a polycyclic group are displayed as `f1`, `f2`, `f3`, etc.,
and every element of a polycyclic group is displayed as product of the
generators.

# Examples

```jldoctest
julia> G = abelian_group(PcGroup, [2, 4]);

julia> G[1], G[2]
(f1, f2)

julia> G[2]*G[1]
f1*f2
```

Note that this does not define Julia variables named `f1`, `f2`, etc.!
To get the generators of the group `G`, use `gens(G)`;
for convenience they can also be accessed as `G[1]`, `G[2]`,
as shown in Section [Elements of groups](@ref elements_of_groups).
"""
const PcGroupElem = BasicGAPGroupElem{PcGroup}


"""
    SubPcGroup

Subgroup of a polycyclic group,
a group that is defined by generators that are elements of a group `G`
of type [`PcGroup`](@ref).
The arithmetic operations with elements are thus performed using the
polycyclic presentation of `G`.

Operations for computing subgroups of a group of type `PcGroup` or
`SubPcGroup`, such as `derived_subgroup` and `sylow_subgroup`,
return groups of type `SubPcGroup`.
"""
@attributes mutable struct SubPcGroup <: GAPGroup
  X::GapObj
  full_group::PcGroup
#T better create an embedding!

  function SubPcGroup(G::GapObj)
    @assert GAPWrap.IsPcGroup(G) || GAP.Globals.IsPcpGroup(G)::Bool
    if GAPWrap.IsPcGroup(G)
      full = GAP.Globals.GroupOfPcgs(GAP.Globals.FamilyPcgs(G)::GapObj)::GapObj
#T use GAPWrap!
    else
      full = GAP.Globals.PcpGroupByCollectorNC(GAP.Globals.Collector(G)::GapObj)::GapObj
    end
    z = new(G, PcGroup(full))
    return z
  end
end

sub_pc_group(G::GapObj) = SubPcGroup(G)


"""
    SubPcGroupElem

Element of a subgroup of a polycyclic group.

# Examples

```jldoctest

... hier! ...

```
"""
const SubPcGroupElem = BasicGAPGroupElem{SubPcGroup}


"""
    FPGroup

Finitely presented group.
Such groups can be constructed a factors of free groups,
see [`free_group`](@ref).
"""
@attributes mutable struct FPGroup <: GAPGroup
  X::GapObj
  
  function FPGroup(G::GapObj)
    @assert GAPWrap.IsSubgroupFpGroup(G)
    z = new(G)
    return z
  end
end

fp_group(G::GapObj) = FPGroup(G)

"""
TODO: document this
"""
const FPGroupElem = BasicGAPGroupElem{FPGroup}

abstract type AbstractMatrixGroupElem <: GAPGroupElem{GAPGroup} end

# NOTE: always defined are deg, ring and at least one between { X, gens, descr }
"""
    MatrixGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup

Type of groups `G` of `n x n` matrices over the ring `R`, where `n = degree(G)` and `R = base_ring(G)`.
"""
@attributes mutable struct MatrixGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup
   deg::Int
   ring::Ring
   X::GapObj
   gens::Vector{<:AbstractMatrixGroupElem}
   descr::Symbol                       # e.g. GL, SL, symbols for isometry groups
   ring_iso::MapFromFunc # Isomorphism from the Oscar base ring to the GAP base ring

   function MatrixGroup{RE,T}(F::Ring, m::Int) where {RE,T}
     G = new{RE, T}()
     G.deg = m
     G.ring = F
     return G
   end
end

# NOTE: at least one of the fields :elm and :X must always defined, but not necessarily both of them.
"""
    MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}} <: AbstractMatrixGroupElem

Elements of a group of type `MatrixGroup{RE<:RingElem, T<:MatElem{RE}}`
"""
mutable struct MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}} <: AbstractMatrixGroupElem
   parent::MatrixGroup{RE, T}
   elm::T                         # Oscar matrix
   X::GapObj                     # GAP matrix. If x isa MatrixGroupElem, then x.X = map_entries(x.parent.ring_iso, x.elm)

   # full constructor
   MatrixGroupElem{RE,T}(G::MatrixGroup{RE,T}, x::T, x_gap::GapObj) where {RE, T} = new{RE,T}(G, x, x_gap)

   # constructor which leaves `X` undefined
   MatrixGroupElem{RE,T}(G::MatrixGroup{RE,T}, x::T) where {RE, T} = new{RE,T}(G, x)

   # constructor which leaves `elm` undefined
   function MatrixGroupElem{RE,T}(G::MatrixGroup{RE,T}, x_gap::GapObj) where {RE, T}
      z = new{RE,T}(G)
      z.X = x_gap
      return z
   end
end

################################################################################
#
# Construct an Oscar group wrapping the GAP group `obj`
# *and* compatible with a given Oscar group `G`.

const _sub_types = Dict{Type, Type}()

_sub_types[PcGroup] = SubPcGroup

function sub_type(T::Type)
  haskey(_sub_types, T) && return _sub_types[T]
  return T
end

sub_type(G::GAPGroup) = sub_type(typeof(G))

# _oscar_group is used to create the subgroup of `G`
# that is described by the GAP group `obj`;
# default: ignore `G`
function _oscar_group(obj::GapObj, G::GAPGroup)
  S = sub_type(G)(obj)
  @assert GAP.Globals.FamilyObj(S.X) === GAP.Globals.FamilyObj(G.X)
  return S
end
#T better rename to _oscar_subgroup?

# `PermGroup`: set the degree of `G`
function _oscar_group(obj::GapObj, G::PermGroup)
  n = GAPWrap.LargestMovedPoint(obj)
  N = degree(G)
  n <= N || error("requested degree ($N) is smaller than the largest moved point ($n)")
  return permutation_group(obj, N)
end

# `MatrixGroup`: set dimension and ring of `G`
function _oscar_group(obj::GapObj, G::MatrixGroup)
  d = GAP.Globals.DimensionOfMatrixGroup(obj)
  d == G.deg || error("requested dimension of matrices ($(G.deg)) does not match the given matrix dimension ($d)")

  R = G.ring
  iso = _ring_iso(G)
  GAPWrap.IsSubset(codomain(iso), GAP.Globals.FieldOfMatrixGroup(obj)) || error("matrix entries are not in the requested ring ($(codomain(iso)))")

  M = matrix_group(R, d)
  M.X = obj
  M.ring = R
  M.ring_iso = iso
  return M
end


################################################################################
#
# "Coerce" an Oscar group `G` to one that is compatible with
# the given Oscar group `S`.
#T what does compatible mean?
compatible_group(G::T, S::T) where T <: GAPGroup = _oscar_group(GapObj(G), S)


################################################################################

abstract type GSet{T} end


################################################################################
#
#   Conjugacy Classes
#
################################################################################

"""
    GroupConjClass{T, S}

It can be either the conjugacy class of an element or of a subgroup of type `S`
in a group `G` of type `T`.
"""
abstract type GroupConjClass{T, S} <: GSet{T} end


################################################################################
#
#  Group Homomorphism
#
################################################################################

abstract type GAPMap <: SetMap end

struct GAPGroupHomomorphism{S<: GAPGroup, T<: GAPGroup} <: Map{S,T,GAPMap,GAPGroupHomomorphism{S,T}}
   domain::S
   codomain::T
   map::GapObj

   function GAPGroupHomomorphism(G::S, H::T, mp::GapObj) where {S<: GAPGroup, T<: GAPGroup}
     return new{S, T}(G, H, mp)
   end
end

GapObj(f::GAPGroupHomomorphism) = f.map


"""
    AutomorphismGroup{T} <: GAPGroup

Group of automorphisms over a group of type `T`. It can be defined via the function `automorphism_group`
"""
@attributes mutable struct AutomorphismGroup{T} <: GAPGroup
  X::GapObj
  G::T

  function AutomorphismGroup{T}(G::GapObj, H::T) where T
    @assert GAPWrap.IsGroupOfAutomorphisms(G)
    z = new{T}(G, H)
    return z
  end
end

function AutomorphismGroup(G::GapObj, H::T) where T
  return AutomorphismGroup{T}(G, H)
end

(aut::AutomorphismGroup{T} where T)(x::GapObj) = group_element(aut,x)

const AutomorphismGroupElem{T} = BasicGAPGroupElem{AutomorphismGroup{T}} where T

function Base.show(io::IO, AGE::AutomorphismGroupElem{FinGenAbGroup}) 
    print(io, "Automorphism of ", FinGenAbGroup, " with matrix representation ", matrix(AGE))
end

################################################################################
#
#  Composite Groups
#
################################################################################


"""
    DirectProductGroup

Either direct product of two or more groups of any type, or subgroup of a direct product of groups.
"""
@attributes mutable struct DirectProductGroup <: GAPGroup
  X::GapObj
  L::Vector{<:GAPGroup}   # list of groups
  Xfull::GapObj      # direct product of the GAP groups of L
  isfull::Bool     # true if G is direct product of the groups of L, false if it is a proper subgroup

  function DirectProductGroup(X::GapObj, L::Vector{<:GAPGroup}, Xfull::GapObj, isfull::Bool)
    return new(X, L, Xfull, isfull)
  end
end


"""
    SemidirectProductGroup{S,T}

Semidirect product of two groups of type `S` and `T` respectively, or
subgroup of a semidirect product of groups.
"""
@attributes mutable struct SemidirectProductGroup{S<:GAPGroup, T<:GAPGroup} <: GAPGroup
  X::GapObj
  N::S              # normal subgroup
  H::T              # group acting on N
  f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}        # action of H on N
  Xfull::GapObj         # full semidirect product: X is a subgroup of Xfull.
  isfull::Bool     # true if X==Xfull

  function SemidirectProductGroup{S, T}(X::GapObj, N::S, H::T, f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}, Xfull::GapObj, isfull::Bool) where {S<:GAPGroup, T<:GAPGroup}
    return new{S, T}(X, N, H, f, Xfull, isfull)
  end
end

"""
    WreathProductGroup

Wreath product of a group `G` and a group of permutations `H`, or a generic
group `H` together with the homomorphism `a` from `H` to a permutation
group.
"""
@attributes mutable struct WreathProductGroup <: GAPGroup
  X::GapObj
  G::GAPGroup
  H::GAPGroup
  a::GAPGroupHomomorphism   # morphism from H to the permutation group
  Xfull::GapObj            # if H does not move all the points, this is the wreath product of (G, Sym(degree(H))
  isfull::Bool             # true if Xfull == X

  function WreathProductGroup(X::GapObj, G::GAPGroup, H::GAPGroup, a::GAPGroupHomomorphism, Xfull::GapObj, isfull::Bool)
    return new(X, G, H, a, Xfull, isfull)
  end
end


"""
    elem_type(::Type{T}) where T <: GAPGroup
    elem_type(::T) where T <: GAPGroup

`elem_type` maps (the type of) a group to the type of its elements.
For now, a group of type `T` has elements of type `BasicGAPGroupElem{T}`.
So we provide it mostly for consistency with other parts of OSCAR.
In the future, a more elaborate setup for group element types
might also be needed.
"""
elem_type(::Type{T}) where T <: GAPGroup = BasicGAPGroupElem{T}

Base.eltype(::Type{T}) where T <: GAPGroup = BasicGAPGroupElem{T}

# `parent_type` is defined and documented in AbstractAlgebra.
parent_type(::Type{BasicGAPGroupElem{T}}) where T <: GAPGroup = T

#
# The array _gap_group_types contains pairs (X,Y) where
# X is a GAP filter such as IsPermGroup, and Y is a corresponding
# Julia type such as `PermGroup`.
#
 
const _gap_group_types = Tuple{GAP.GapObj, Type}[]

# important:
# The function returned by _get_type is not allowed to
# create a new GAP group (pc group or fp group)
function _get_type(G::GapObj)
  for pair in _gap_group_types
    if pair[1](G)
      if pair[2] == MatrixGroup
#T HACK: We need more information in the case of matrix groups.
#T (Usually we should not need to guess the Oscar side of a GAP group.)
        return function(dom::GAP.GapObj)
                 deg = GAP.Globals.DimensionOfMatrixGroup(dom)
                 iso = iso_gap_oscar(GAP.Globals.FieldOfMatrixGroup(dom))
                 ring = codomain(iso)
                 matgrp = matrix_group(ring, deg)
                 matgrp.ring_iso = inv(iso)
                 matgrp.X = dom
                 return matgrp
               end
      elseif pair[2] == AutomorphismGroup
        return function(A::GAP.GapObj)
                 actdom_gap = GAP.Globals.AutomorphismDomain(A)
                 actdom_oscar = _get_type(actdom_gap)(actdom_gap)
                 return AutomorphismGroup(A, actdom_oscar)
               end
      else
        # The result will be used as a constructor for an Oscar group,
        # in order to wrap `G`.
        # Thus we have to switch to the appropriate subgroup type
        # if and only if `G` is not the full group.
        if pair[2] === PcGroup
          return _is_full_pc_group(G) ? pair[2] : SubPcGroup
        end
        return pair[2]
      end
    end
  end
  error("Not a known type of group")
end

# Check the compatibility of two groups in the sense that an element in the
# first group can be multiplied with an element in the second.
#T To which group does the product belong?
#T In which functions is this used?
# The group *types* can be different,

# and check their *contexts*
# (in case of f.p. and p.c. groups and matrix groups)

# The underlying GAP groups must belong to the same family.
function _check_compatible(G1::GAPGroup, G2::GAPGroup; error::Bool = true)
  GAPWrap.FamilyObj(G1.X) === GAPWrap.FamilyObj(G2.X) && return true
  error && throw(ArgumentError("G1 and G2 are not compatible"))
  return false
end

# Any two permutation groups are compatible.
_check_compatible(G1::PermGroup, G2::PermGroup; error::Bool = true) = true

# The groups must have the same dimension and the same base ring.
function _check_compatible(G1::MatrixGroup, G2::MatrixGroup; error::Bool = true)
  base_ring(G1) == base_ring(G2) && degree(G1) == degree(G2) && return true
  error && throw(ArgumentError("G1 and G2 must have same base_ring and degree"))
  return false
end

#T FinGenAbGroup: how do they behave?
