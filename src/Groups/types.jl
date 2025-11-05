@doc raw"""
    GAPGroup <: Group

Each object of the abstract type `GAPGroup` stores a group object from
the GAP system,
and thus can delegate questions about this object to GAP.

For expert usage, you can extract the underlying GAP object via `GapObj`,
i.e., if `G` is a `GAPGroup`, then `GapObj(G)` is the `GapObj` underlying `G`.

Concrete subtypes of `GAPGroup` are `PermGroup`, `FPGroup`, `SubFPGroup`,
`PcGroup`, `SubPcGroup`, and `MatrixGroup`.
"""
abstract type GAPGroup <: Group end

## `GapGroup` to underlying GAP group
GAP.@install GapObj(obj::GAPGroup) = obj.X

@doc raw"""
    GAPGroupElem <: GroupElem

Each object of the abstract type `GAPGroupElem` stores a group element
object from the GAP system,
and thus can delegate questions about this object to GAP.

For expert usage, you can extract the underlying GAP object via `GapObj`,
i.e., if `g` is a `GAPGroupElem`, then `GapObj(g)` is the `GapObj` underlying `g`.
"""
abstract type GAPGroupElem{T<:GAPGroup} <: GroupElem end

## `GapGroupElem` to GAP group element
GAP.@install GapObj(obj::GAPGroupElem) = obj.X

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
Symmetric group of degree 5

julia> x=cperm([1,2,3])
(1,2,3)

julia> parent(x)
Symmetric group of degree 3

julia> y=G(x)
(1,2,3)

julia> parent(y)
Symmetric group of degree 5
```

If `G` is a permutation group and `x` is a vector of integers,
`G(x)` returns a [`PermGroupElem`](@ref) with parent `G`;
an exception is thrown if the element does not embed into `G`.

# Examples
```jldoctest
julia> G = symmetric_group(6)
Symmetric group of degree 6

julia> x = G([2,4,6,1,3,5])
(1,2,4)(3,6,5)

julia> parent(x)
Symmetric group of degree 6
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

function Base.hash(x::PermGroupElem, h::UInt)
  return UInt(GAPWrap.HashPermutation(GapObj(x), GapInt(h)))
end


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
  as_sub_pc_group::GAPGroup  # we cannot prescribe `SubPcGroup`

  function PcGroup(G::GapObj)
    # The constructor is not allowed to replace the given GAP group.
    # (The function `pc_group` may do this.)
    @assert _is_full_pc_group(G)
    return new(G)
  end
end

function pc_group(G::GapObj)
  _is_full_pc_group(G) && return PcGroup(G)

  # Switch to a full pcp or pc group.
  if GAPWrap.IsPcpGroup(G)
    return PcGroup(GAP.Globals.PcpGroupByPcp(GAP.Globals.Pcp(G)::GapObj)::GapObj)
  elseif GAPWrap.IsPcGroup(G)
    return PcGroup(GAP.Globals.PcGroupWithPcgs(GAP.Globals.Pcgs(G)::GapObj)::GapObj)
  end
  throw(ArgumentError("G must be in IsPcGroup or IsPcpGroup"))
end

# Return `true` if the generators of `G` fit to those of its pc presentation.
function _is_full_pc_group(G::GapObj)
  GAPWrap.IsPcpGroup(G) || GAPWrap.IsPcGroup(G) || return false
  return GAP.Globals.GroupGeneratorsDefinePresentation(G)::Bool
end
#T the code for PcpGroups is too expensive!

function as_sub_pc_group(G::PcGroup)
  if !isdefined(G, :as_sub_pc_group)
    G.as_sub_pc_group = sub_pc_group(G)
  end
  return G.as_sub_pc_group::SubPcGroup
end

"""
    PcGroupElem

Element of a polycyclic group.

The generators of a polycyclic group are displayed as `f1`, `f2`, `f3`, etc.,
and every element of a polycyclic group is displayed as product of the
generators.

# Examples

```jldoctest
julia> G = abelian_group(PcGroup, [2, 3]);

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
    @assert GAPWrap.IsPcGroup(G) || GAPWrap.IsPcpGroup(G)
    if GAPWrap.IsPcGroup(G)
      full = GAPWrap.GroupOfPcgs(GAPWrap.FamilyPcgs(G))
    else
      full = GAPWrap.PcpGroupByCollectorNC(GAPWrap.Collector(G))
    end
    z = new(G, PcGroup(full))
    return z
  end
end

sub_pc_group(G::GapObj) = SubPcGroup(G)


"""
    SubPcGroupElem

Element of a subgroup of a polycyclic group.

The elements are displayed in the same way as the elements of full
polycyclic groups, see [`PcGroupElem`](@ref).

# Examples

```jldoctest
julia> G = abelian_group(SubPcGroup, [4, 2]);

julia> G[1], G[2]
(f1, f3)

julia> G[2]*G[1]
f1*f3
```
"""
const SubPcGroupElem = BasicGAPGroupElem{SubPcGroup}

function Base.hash(x::Union{PcGroupElem,SubPcGroupElem}, h::UInt)
  return hash(letters(x), hash(parent(x), h))
end


"""
    FPGroup

Finitely presented group.
Such groups can be constructed a factors of free groups,
see [`free_group`](@ref).

For a group `G` of type `FPGroup`, the elements in `gens(G)` satisfy the
relators of the underlying presentation.

Functions that compute subgroups of `G` return groups of type `SubFPGroup`.
"""
@attributes mutable struct FPGroup <: GAPGroup
  X::GapObj
  as_sub_fp_group::GAPGroup  # we cannot prescribe `SubFPGroup`

  function FPGroup(G::GapObj)
    # Accept only full f.p. groups.
    @assert GAPWrap.IsFpGroup(G)
    return new(G)
  end
end

function fp_group(G::GapObj)
  _is_full_fp_group(G) && return FPGroup(G)

  # Switch to a full fp group.
  @req GAPWrap.IsSubgroupFpGroup(G) "G must be in IsSubgroupFpGroup"

  f = GAP.Globals.IsomorphismFpGroup(G)::GapObj
  return FPGroup(GAPWrap.Range(f))
end

# Return `true` if the generators of `G` fit
# to those of its underlying presentation.
_is_full_fp_group(G::GapObj) = GAPWrap.IsFpGroup(G)

function as_sub_fp_group(G::FPGroup)
  if !isdefined(G, :as_sub_fp_group)
    G.as_sub_fp_group = sub_fp_group(G)
  end
  return G.as_sub_fp_group::SubFPGroup
end

"""
    FPGroupElem

Element of a finitely presented group.

The generators of a finitely presented group are displayed as
`f1`, `f2`, `f3`, etc.,
and every element of a finitely presented group is displayed as product of the
generators.
"""
const FPGroupElem = BasicGAPGroupElem{FPGroup}

"""
    SubFPGroup

Subgroup of a finitely presented group,
a group that is defined by generators that are elements of a group `G`
of type [`FPGroup`](@ref).

Operations for computing subgroups of a group of type `FPGroup` or
`SubFPGroup`, such as `derived_subgroup` and `sylow_subgroup`,
return groups of type `SubFPGroup`.

Note that functions such as [`relators`](@ref) do not make sense for proper
subgroups of a finitely presented group.
"""
@attributes mutable struct SubFPGroup <: GAPGroup
  X::GapObj
  full_group::FPGroup
#T better create an embedding!

  function SubFPGroup(G::GapObj)
    @assert GAPWrap.IsSubgroupFpGroup(G)
    full = GAP.getbangproperty(GAPWrap.FamilyObj(G), :wholeGroup)::GapObj
    z = new(G, FPGroup(full))
    return z
  end
end

sub_fp_group(G::GapObj) = SubFPGroup(G)


"""
    SubFPGroupElem

Element of a subgroup of a finitely presented group.

The elements are displayed in the same way as the elements of full
finitely presented groups, see [`FPGroupElem`](@ref).
"""
const SubFPGroupElem = BasicGAPGroupElem{SubFPGroup}

function Base.hash(x::Union{FPGroupElem,SubFPGroupElem}, h::UInt)
  return hash(letters(x), hash(parent(x), h))
end


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

Type of elements of a group of type `MatrixGroup{RE, T}`.
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

function Base.hash(x::MatrixGroupElem, h::UInt)
  return hash(matrix(x), hash(parent(x), h))
end

################################################################################
#
# Construct an Oscar group wrapping the GAP group `obj`
# *and* compatible with a given Oscar group `G`.

sub_type(T::Type) = T
sub_type(::Type{PcGroup}) = SubPcGroup
sub_type(::Type{FPGroup}) = SubFPGroup
sub_type(G::GAPGroup) = sub_type(typeof(G))

# `_oscar_subgroup(obj, G)` is used to create the subgroup of `G`
# that is described by the GAP group `obj`;
# default: ignore `G`
function _oscar_subgroup(obj::GapObj, G::GAPGroup)
  S = sub_type(G)(obj)
  @assert GAP.Globals.FamilyObj(GapObj(S)) === GAP.Globals.FamilyObj(GapObj(G))
  return S
end

# `PermGroup`: set the degree of `G`
function _oscar_subgroup(obj::GapObj, G::PermGroup)
  n = GAPWrap.LargestMovedPoint(obj)
  N = degree(G)
  n <= N || error("requested degree ($N) is smaller than the largest moved point ($n)")
  return permutation_group(obj, N)
end

# `MatrixGroup`: set dimension and ring of `G`
function _oscar_subgroup(obj::GapObj, G::MatrixGroup)
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

abstract type GSet{T,S} end


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
abstract type GroupConjClass{T,S} <: GSet{T,S} end


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
  is_known_to_be_full::Bool

  function AutomorphismGroup{T}(G::GapObj, H::T, full::Bool = false) where T
    @assert GAPWrap.IsGroupOfAutomorphisms(G)
    z = new{T}(G, H, full)
    return z
  end
end

function AutomorphismGroup(G::GapObj, H::T, full::Bool = false) where T
  return AutomorphismGroup{T}(G, H, full)
end

(aut::AutomorphismGroup{T} where T)(x::GapObj) = group_element(aut,x)

const AutomorphismGroupElem{T} = BasicGAPGroupElem{AutomorphismGroup{T}} where T

function Base.show(io::IO, ::MIME"text/plain", f::AutomorphismGroupElem{FinGenAbGroup})
  D = domain(parent(f))
  io = pretty(io)
  println(io, "Automorphism of")
  println(io, Indent(), Lowercase(), D, Dedent())
  println(io, "with matrix representation")
  print(io, Indent())
  show(io, MIME"text/plain"(), matrix(f))
  print(io, Dedent())
end

function Base.show(io::IO, f::AutomorphismGroupElem{FinGenAbGroup})
  print(io, matrix(f))
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
const _gap_group_types = Tuple{GapObj, Type}[]

# `_oscar_group(G)` wraps the GAP group `G` into a suitable Oscar group `OG`,
# such that `GapObj(OG)` is equal to `G`.
# The function is not allowed to create an independent new GAP group object;
# this would be fatal at least for pc groups and fp groups.
function _oscar_group(G::GapObj)
  for pair in _gap_group_types
    if pair[1](G)
      if pair[2] == MatrixGroup
#T HACK: We need more information in the case of matrix groups.
#T (Usually we should not need to guess the Oscar side of a GAP group.)
        deg = GAP.Globals.DimensionOfMatrixGroup(G)
        iso = iso_gap_oscar(GAP.Globals.FieldOfMatrixGroup(G))
        ring = codomain(iso)
        matgrp = matrix_group(ring, deg)
        matgrp.ring_iso = inv(iso)
        set_attribute!(ring, :iso_oscar_gap, matgrp.ring_iso)
        matgrp.X = G
        return matgrp
      elseif pair[2] == AutomorphismGroup
        actdom_gap = GAP.Globals.AutomorphismDomain(G)
        actdom_oscar = _oscar_group(actdom_gap)
        return AutomorphismGroup(G, actdom_oscar)
      elseif pair[2] === PcGroup && !_is_full_pc_group(G)
        # We have to switch to the appropriate subgroup type
        # if and only if `G` is not the full group.
        return SubPcGroup(G)
      elseif pair[2] === FPGroup && !_is_full_fp_group(G)
        # We have to switch to the appropriate subgroup type
        # if and only if `G` is not the full group.
        return SubFPGroup(G)
      else
        return pair[2](G)
      end
    end
  end
  error("Not a known type of group")
end


# Check the compatibility of two groups in the sense that an element in the
# first group can be multiplied with an element in the second.
#T To which group does the product belong?
#T In which functions is this used?
# The group *types* can be different.

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

#TODO FinGenAbGroup: How do they behave?
#     And does this question arise, since embeddings are used throughout?
