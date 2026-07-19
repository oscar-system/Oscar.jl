@doc raw"""
    GAPGroup <: Group

Each object of the abstract type `GAPGroup` stores a group object from
the GAP system,
and thus can delegate questions about this object to GAP.

For expert usage, you can extract the underlying GAP object via `GapObj`,
i.e., if `G` is a `GAPGroup`, then `GapObj(G)` is the `GapObj` underlying `G`.

Concrete subtypes of `GAPGroup` are `PermGroup`, `FPGroup`, `SubFPGroup`,
`PcGroup`, `SubPcGroup`, and `MatGroup`.
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
     n = GAPWrap.LargestMovedPoint(G)
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
     @assert GAPWrap.IsPermGroup(G) && deg > 0 && deg >= GAPWrap.LargestMovedPoint(G)
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
    G.as_sub_pc_group = sub_pc_group(GapObj(G), G, check = false)
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
    # Create the `PcGroup` on the Oscar side anew.
    @assert GAPWrap.IsPcGroup(G) || GAPWrap.IsPcpGroup(G)
    if GAPWrap.IsPcGroup(G)
      full = GAPWrap.GroupOfPcgs(GAPWrap.FamilyPcgs(G))
    else
      full = GAPWrap.PcpGroupByCollectorNC(GAPWrap.Collector(G))
    end
    z = new(G, PcGroup(full))
    return z
  end

  function SubPcGroup(G::GapObj, F::PcGroup; check::Bool=true)
    # Keep object identity of the available `PcGroup` on the Oscar side.
    if check
      @assert (GAPWrap.IsPcGroup(G) || GAPWrap.IsPcpGroup(G)) &&
              GAPWrap.FamilyObj(G) === GAPWrap.FamilyObj(F.X)
    end
    z = new(G, F)
    return z
  end
end

sub_pc_group(G::GapObj) = SubPcGroup(G)
sub_pc_group(G::GapObj, F::PcGroup; check::Bool=true) = SubPcGroup(G, F; check=check)


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


"""
    FPGroup

Finitely presented group.
Such groups can be constructed as factors of free groups,
see [`free_group`](@ref).

For a group `G` of type `FPGroup`, the elements in `gens(G)` satisfy the
relators of the underlying presentation.

Functions that compute subgroups of `G` return groups of type `SubFPGroup`.
"""
@attributes mutable struct FPGroup <: GAPGroup
  X::GapObj
  as_sub_fp_group::GAPGroup  # we cannot prescribe `SubFPGroup`
  free_group::GAPGroup  # we cannot prescribe `FPGroup`

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

function fp_group(G::GapObj, F::FPGroup)
  # We want to store `F` as the `free_group` value of the result.
  # For that, `G` must be a full f.p. group on the GAP side
  # (otherwise prescribing a known free group on the Oscar side
  # does not make sense)
  # and `GapObj(F)` must be identical with the free group on the GAP side
  # that is stored in `G`
  gapF = GapObj(F)
  @assert GAPWrap.IsFpGroup(G)
  @assert gapF === GAPWrap.FreeGroupOfFpGroup(G)
  FG = FPGroup(G)
  FG.free_group = F
  return FG
end

# Return `true` if the generators of `G` fit
# to those of its underlying presentation.
_is_full_fp_group(G::GapObj) = GAPWrap.IsFpGroup(G)

function as_sub_fp_group(G::FPGroup)
  if !isdefined(G, :as_sub_fp_group)
    G.as_sub_fp_group = sub_fp_group(GapObj(G), G, check = false)
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
    # Create the `FPGroup` on the Oscar side anew.
    @assert GAPWrap.IsSubgroupFpGroup(G)
    full = GAP.getbangproperty(GAPWrap.FamilyObj(G), :wholeGroup)::GapObj
    z = new(G, FPGroup(full))
    return z
  end

  function SubFPGroup(G::GapObj, F::FPGroup; check::Bool=true)
    # Keep object identity of the available `FPGroup` on the Oscar side.
    if check
      @assert GAPWrap.IsSubgroupFpGroup(G) &&
              GAPWrap.FamilyObj(G) === GAPWrap.FamilyObj(F.X)
    end
    z = new(G, F)
    return z
  end
end

sub_fp_group(G::GapObj) = SubFPGroup(G)
sub_fp_group(G::GapObj, F::FPGroup; check::Bool=true) = SubFPGroup(G, F; check=check)


"""
    SubFPGroupElem

Element of a subgroup of a finitely presented group.

The elements are displayed in the same way as the elements of full
finitely presented groups, see [`FPGroupElem`](@ref).
"""
const SubFPGroupElem = BasicGAPGroupElem{SubFPGroup}


abstract type AbstractMatGroupElem <: GAPGroupElem{GAPGroup} end

# NOTE: always defined are deg, ring and at least one between { X, gens, descr }
"""
    MatGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup

Type of groups `G` of `n x n` matrices over the ring `R`, where `n = degree(G)` and `R = base_ring(G)`.
"""
@attributes mutable struct MatGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup
   deg::Int
   ring::Ring
   X::GapObj
   gens::Vector{<:AbstractMatGroupElem}
   descr::Symbol                       # e.g. GL, SL, symbols for isometry groups
   ring_iso::MapFromFunc # Isomorphism from the Oscar base ring to the GAP base ring

   function MatGroup{RE,T}(F::Ring, m::Int) where {RE,T}
     G = new{RE, T}()
     G.deg = m
     G.ring = F
     return G
   end
end

# Construct an Oscar matrix group from a GAP matrix group,
# guess the ring over which the Oscar group will live.
# (This is usually not what we want to do.)
function matrix_group(G::GapObj)
  @req GAPWrap.IsMatrixGroup(G) "the GAP group must be a matrix group"
  iso = inv(iso_gap_oscar(GAPWrap.FieldOfMatrixGroup(G)))
  set_attribute!(domain(iso), :iso_oscar_gap, iso)
  return matrix_group(iso, G; check = false)
end

# Construct an Oscar matrix group from a GAP matrix group.
# In order to construct the group over a given ring,
# prescribe the `ring_iso` of the desired matrix group.
function matrix_group(ring_iso::Map, G::GapObj; check::Bool=true)
  if check
    @req GAPWrap.IsMatrixGroup(G) "the GAP group must be a matrix group"
    C = codomain(ring_iso)
    @req isa(C, GapObj) "the codomain of the given map must be a GapObj"
    @req GAPWrap.IsSubset(C, GAPWrap.FieldOfMatrixGroup(G)) "the codomain of the given map must contain the matrix entries"
  end
  deg = GAPWrap.DimensionOfMatrixGroup(G)
  ring = domain(ring_iso)
  matgrp = matrix_group(ring, deg)
  matgrp.ring_iso = ring_iso
  matgrp.X = G
  return matgrp
end


# NOTE: at least one of the fields :elm and :X must always defined, but not necessarily both of them.
"""
    MatGroupElem{RE<:RingElem, T<:MatElem{RE}} <: AbstractMatGroupElem

Type of elements of a group of type `MatGroup{RE, T}`.
"""
mutable struct MatGroupElem{RE<:RingElem, T<:MatElem{RE}} <: AbstractMatGroupElem
   parent::MatGroup{RE, T}
   elm::T                         # Oscar matrix
   X::GapObj                     # GAP matrix. If x isa MatGroupElem, then x.X = map_entries(x.parent.ring_iso, x.elm)

   # full constructor
   MatGroupElem{RE,T}(G::MatGroup{RE,T}, x::T, x_gap::GapObj) where {RE, T} = new{RE,T}(G, x, x_gap)

   # constructor which leaves `X` undefined
   MatGroupElem{RE,T}(G::MatGroup{RE,T}, x::T) where {RE, T} = new{RE,T}(G, x)

   # constructor which leaves `elm` undefined
   function MatGroupElem{RE,T}(G::MatGroup{RE,T}, x_gap::GapObj) where {RE, T}
      z = new{RE,T}(G)
      z.X = x_gap
      return z
   end
end


################################################################################
#
# Construct an Oscar group wrapping the GAP group `obj`
# *and* compatible with a given Oscar group `G`.

sub_type(T::Type) = T
sub_type(::Type{PcGroup}) = SubPcGroup
sub_type(::Type{FPGroup}) = SubFPGroup
sub_type(G::GAPGroup) = sub_type(typeof(G))

# `_oscar_subgroup(obj, G; check::Bool = true)` is used to create
# the subgroup of `G` that is described by the GAP group `obj`;
# default: ignore `G` and `check`
function _oscar_subgroup(obj::GapObj, G::GAPGroup; check::Bool = true)
  S = sub_type(G)(obj)
  @assert GAPWrap.FamilyObj(GapObj(S)) === GAPWrap.FamilyObj(GapObj(G))
  return S
end

# `PermGroup`: set the degree of `G`, ignore `check`
function _oscar_subgroup(obj::GapObj, G::PermGroup; check::Bool = true)
  n = GAPWrap.LargestMovedPoint(obj)
  N = degree(G)
  n <= N || error("requested degree ($N) is smaller than the largest moved point ($n)")
  return permutation_group(obj, N)
end

# `MatGroup`: set dimension and ring of `G`
function _oscar_subgroup(obj::GapObj, G::MatGroup; check::Bool = true)
  d = GAPWrap.DimensionOfMatrixGroup(obj)
  @req G.deg == d "requested dimension of matrices ($(G.deg)) does not match the given matrix dimension ($d)"
  return matrix_group(_ring_iso(G), obj; check = check)
end

# `FPGroup`: keep the object identity of `G`
function _oscar_subgroup(obj::GapObj, G::FPGroup; check::Bool = true)
  return SubFPGroup(obj, G; check = check)
end

# `PcGroup`: keep the object identity of `G`
function _oscar_subgroup(obj::GapObj, G::PcGroup; check::Bool = true)
  return SubPcGroup(obj, G; check = check)
end

################################################################################
#
#  Cosets
#
################################################################################

# T=type of the group, S=type of the element
@doc raw"""
    GroupCoset{TG <: GAPGroup, TH <: GAPGroup, S <: GAPGroupElem}

Type of right and left cosets of subgroups in groups.

For an element $g$ in a group $G$, and a subgroup $H$ of $G$,
the set $Hg = \{ hg; h \in H \}$ is a right coset of $H$ in $G$,
and the set $gH = \{ gh; h \in H \}$ is a left coset of $H$ in $G$.

- [`group(C::GroupCoset)`](@ref) returns $G$.

- [`acting_group(C::GroupCoset)`](@ref) returns $H$.

- [`representative(C::GroupCoset)`](@ref) returns an element
  (the same element for each call) of `C`.

- [`is_right(C::GroupCoset)`](@ref) and [`is_left(C::GroupCoset)`](@ref)
  return whether `C` is a right or left coset, respectively.

Two cosets are equal if and only if they are both left or right, respectively,
and they contain the same elements.
"""
struct GroupCoset{TG <: GAPGroup, TH <: GAPGroup, S <: GAPGroupElem}
  G::TG                   # big group containing the subgroup and the element
  H::TH                   # subgroup (may have a different type)
  repr::S                 # element
  side::Symbol            # says if the coset is left or right
  X::Ref{GapObj}          # GapObj(H*repr)

  function GroupCoset(G::TG, H::TH, representative::S, side::Symbol) where {TG <: GAPGroup, TH <: GAPGroup, S <:GAPGroupElem}
    return new{TG, TH, S}(G, H, representative, side, Ref{GapObj}())
  end
end

@doc raw"""
    SubgroupTransversal{T<: GAPGroup, S<: GAPGroup, E<: GAPGroupElem}

Type of left/right transversals of subgroups in groups.

For a group $G$ and a subgroup $H$ of $G$, $T$ is a right
(resp. left) transversal for $H$ in $G$ if $T$ contains
precisely one element of each right (resp. left) cosets of $H$ in $G$.

Objects of this type are created by [`right_transversal`](@ref) and
[`left_transversal`](@ref).

- [`group(T::SubgroupTransversal)`](@ref) returns $G$.

- [`subgroup(T::SubgroupTransversal)`](@ref) returns $H$.

# Note for developers

The elements are encoded via a right transversal object in GAP.
(Note that GAP does not support left transversals.)
"""
struct SubgroupTransversal{T<: GAPGroup, S<: GAPGroup, E<: GAPGroupElem} <: AbstractVector{E}
  G::T                    # big group containing the subgroup
  H::S                    # subgroup
  side::Symbol            # says if the transversal is left or right
  X::GapObj               # underlying *right* transversal in GAP
end

@doc raw"""
    GroupDoubleCoset{T<: Group, S <: GAPGroupElem}

Type of double cosets of subgroups in groups.

For an element $g$ in a group $G$, and two subgroups $H$, $K$ of $G$,
the set $HgK = \{ hgk; h \in H, k \in K \}$ is a $H-K$-double coset in $G$.

- [`group(C::GroupDoubleCoset)`](@ref) returns $G$.

- [`left_acting_group(C::GroupDoubleCoset)`](@ref) returns $H$.

- [`right_acting_group(C::GroupDoubleCoset)`](@ref) returns $K$.

- [`representative(C::GroupDoubleCoset)`](@ref) returns an element
  (the same element for each call) of `C`.

Two double cosets are equal if and only if they contain the same elements.
"""
struct GroupDoubleCoset{T <: GAPGroup, S <: GAPGroupElem}
  # T=type of the group, S=type of the element
  G::T
  H::GAPGroup
  K::GAPGroup
  repr::S
  X::Ref{GapObj}
  size::Ref{ZZRingElem}
  right_coset_reps::Ref{Dict{GAPGroupElem, Tuple{GAPGroupElem, GAPGroupElem}}}

  function GroupDoubleCoset(G::T, H::GAPGroup, K::GAPGroup, representative::S) where {T<: GAPGroup, S<:GAPGroupElem}
    return new{T, S}(G, H, K, representative, Ref{GapObj}(), Ref{ZZRingElem}(),
                     Ref{Dict{GAPGroupElem, Tuple{GAPGroupElem, GAPGroupElem}}}())
  end
end

################################################################################
#
#  G-Sets
#
################################################################################

abstract type GSet{T,S} end

"""
    GSetByElements{T,S} <: GSet{T,S}

Objects of this type represent G-sets that are willing to write down
orbits and elements lists as vectors.
These G-sets are created by default by [`gset`](@ref).

The fields are
- the group that acts, of type `T`,
- the Julia function (for example `on_tuples`) that describes the action,
- the seeds (something iterable of eltype `S`) whose closure under the action is the G-set
- the dictionary used to store attributes (orbits, elements, ...).
"""
@attributes mutable struct GSetByElements{T,S} <: GSet{T,S}
  group::T
  action_function::Function
  seeds

  function GSetByElements(G::T, fun::Function, seeds; closed::Bool = false, check::Bool = true) where {T<:Union{Group, FinGenAbGroup}}
    @req !isempty(seeds) "seeds for G-set must be nonempty"
    check && @req hasmethod(fun, (typeof(first(seeds)), elem_type(T))) "action function does not fit to seeds"
    Omega = new{T,eltype(seeds)}(G, fun, seeds, Dict{Symbol,Any}())
    closed && set_attribute!(Omega, :elements => unique!(collect(seeds)))
    return Omega
  end
end
#TODO: How can I specify that `seeds` should be an iterable object?

##  wrapper objects for elements of G-sets,
##  with fields `gset` (the G-set) and `objects` (the unwrapped object)
##
##  These objects are optional ("syntactic sugar"), they can be used to
##  - apply group elements via `^`,
##    not via the action function stored in the G-set,
##  - write something like `orbit(omega)`, `stabilizer(omega)`.
struct ElementOfGSet{T, S, G <: GSet{T, S}}
  gset::G
  obj::S
end

@doc raw"""
    GSetBySubgroupTransversal{T, S, E} <: GSet{T}

Objects of this type represent G-sets that describe the left or right cosets
of a subgroup $H$ in a group $G$.
The group $G$ acts on the G-set by multiplication from the right or (after
taking inverses) from the left.
These G-sets store just transversals,
see [`right_transversal`](@ref) and [`left_transversal`](@ref).
The construction of explicit right or left cosets is not necessary in order
to compute the permutation action of elements of $G$ on the cosets.

The fields are
- the group that acts, of type `T`, with elements of type `E`,
- the subgroup whose cosets are the elements, of type `S`,
- the side from which the group acts (`:right` or `:left`),
- the (left or right) transversal, of type `SubgroupTransversal{T, S, E}`,
- the dictionary used to store attributes (orbits, elements, ...).
"""
@attributes mutable struct GSetBySubgroupTransversal{T, S, E} <: GSet{T,GroupCoset{T, S, E}}
  group::T
  subgroup::S
  side::Symbol
  transversal::SubgroupTransversal{T, S, E}

  function GSetBySubgroupTransversal(G::T, H::S, side::Symbol; check::Bool = true) where {T<:GAPGroup, S<:GAPGroup}
    check && @req is_subgroup(H, G)[1] "H must be a subgroup of G"
    E = eltype(G)
    if side == :right
      tr = right_transversal(G, H)
    elseif side == :left
      tr = left_transversal(G, H)
    else
      throw(ArgumentError("side must be :right or :left"))
    end
    return new{T, S, E}(G, H, side, tr, Dict{Symbol,Any}())
  end
end

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

@attributes mutable struct GAPGroupConjClass{T<:GAPGroup, S<:Union{GAPGroupElem,GAPGroup}} <: GroupConjClass{T, S}
  X::T
  repr::S
  CC::GapObj

  function GAPGroupConjClass(G::T, obj::S, C::GapObj) where T<:GAPGroup where S<:Union{GAPGroupElem, GAPGroup}
    return new{T, S}(G, obj, C, Dict{Symbol,Any}())
  end
end

struct FinGenAbGroupConjClass{T<:FinGenAbGroup, S<:Union{FinGenAbGroupElem,FinGenAbGroup}} <: GroupConjClass{T, S}
  X::T
  repr::S
end

################################################################################
#
#  Group Homomorphism
#
#  We support the following types of homomorphism objects.
#
################################################################################

abstract type GAPMap <: SetMap end

"""
    GAPGroupHomomorphism{S<:GAPGroup, T<:GAPGroup}

An object of the type `GAPGroupHomomorphism` stores a homomorphism object
as its `GapObj` value, and delegates all computations to it.
"""
struct GAPGroupHomomorphism{S<: GAPGroup, T<: GAPGroup} <: Map{S,T,GAPMap,GAPGroupHomomorphism{S,T}}
   domain::S
   codomain::T
   map::GapObj

   function GAPGroupHomomorphism(G::S, H::T, mp::GapObj) where {S<: GAPGroup, T<: GAPGroup}
     return new{S, T}(G, H, mp)
   end
end

"""
    GAPGroupEmbedding{S<:GAPGroup, T<:GAPGroup}

An object `emb` of the type `GAPGroupEmbedding` knows that
the corresponding map on the GAP side is the identity map
on the GAP object corresponding to `domain(emb)`.

A corresponding homomorphism `GapObj(emb)` on the GAP side
is created only on demand, in many situations this is not needed.

Many natural embeddings between `GapGroup` objects are realized as
`GAPGroupEmbedding` objects.
"""
struct GAPGroupEmbedding{S<: GAPGroup, T<: GAPGroup} <: Map{S,T,GAPMap,GAPGroupEmbedding{S,T}}
   domain::S
   codomain::T
   map::Ref{GapObj} # may be missing

   function GAPGroupEmbedding(G::S, H::T) where {S<: GAPGroup, T<: GAPGroup}
     return new{S, T}(G, H, Ref{GapObj}())
   end
end

GapObj(f::GAPGroupHomomorphism) = f.map

function GapObj(f::GAPGroupEmbedding)
  isassigned(f.map) && return f.map[]::GapObj
  mp = GAPWrap.GroupHomomorphismByFunction(GapObj(f.domain), GapObj(f.codomain),
           GAP.Globals.IdFunc, false, GAP.Globals.IdFunc)
  f.map[] = mp
  return mp
end


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

mutable struct GroupIsomorphismFromFunc{R, T} <: Map{R, T, Hecke.HeckeMap, MapFromFunc}
  map::MapFromFunc{R, T}
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
      if pair[2] == MatGroup
        return matrix_group(G)
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
function _check_compatible(G1::MatGroup, G2::MatGroup; error::Bool = true)
  base_ring(G1) == base_ring(G2) && degree(G1) == degree(G2) && return true
  error && throw(ArgumentError("G1 and G2 must have same base_ring and degree"))
  return false
end

#TODO FinGenAbGroup: How do they behave?
#     And does this question arise, since embeddings are used throughout?

################################################################################
#
#  Character Tables
#
################################################################################

abstract type GroupCharacterTable end

"""
    GAPGroupCharacterTable <: GroupCharacterTable

This is the type of (ordinary or Brauer) character tables that can delegate
tasks to an underlying character table object in the GAP system
(field `GAPTable`).

The value of the field `characteristic` determines whether the table
is an ordinary one (value `0`) or a `p`-modular one (value `p`).

A group can (but need not) be stored in the field `group`.
If it is available then also the field `isomorphism` is available,
its value is a bijective map from the `group` value to a group in GAP.

Objects of type `GAPGroupCharacterTable` support [`get_attribute`](@ref),
for example in order to store the already computed `p`-modular tables
in an ordinary table, and to store the corresponding ordinary table
in a `p`-modular table.
"""
@attributes mutable struct GAPGroupCharacterTable <: GroupCharacterTable
  GAPTable::GapObj  # the GAP character table object
  characteristic::T where T <: IntegerUnion
  group::Union{GAPGroup, FinGenAbGroup}    # the underlying group, if any
  isomorphism::Map  # isomorphism from `group` to a group in GAP

  function GAPGroupCharacterTable(G::Union{GAPGroup, FinGenAbGroup}, tab::GapObj, iso::Map, char::T) where T <: IntegerUnion
    return new(tab, char, G, iso)
  end

  function GAPGroupCharacterTable(tab::GapObj, char::T) where T <: IntegerUnion
    # group and isomorphism are left undefined
    return new(tab, char)
  end
end

abstract type GroupClassFunction end

struct GAPGroupClassFunction <: GroupClassFunction
  table::GAPGroupCharacterTable
  values::GapObj
end

"""
    GAPGroupCharacterTableRational <: GroupCharacterTable

This is the type of ordinary *rational* character tables
that can delegate tasks to the underlying character table
(field `character_table`).

As a collection, a rational character table stores the *rational irreducible*
characters of the underlying character table `t`, that is, the Galois sums of
the irreducible characters of `t` (field `irr`).

The norms of these characters (the lengths of the Galois orbits)
are stored in the field `norms`.
"""
@attributes mutable struct GAPGroupCharacterTableRational <: GroupCharacterTable
  character_table::GAPGroupCharacterTable
  irr::Vector{GAPGroupClassFunction}
  norms::Vector{Int}

  function GAPGroupCharacterTableRational(tbl::GAPGroupCharacterTable)
    # irr and norms are left undefined
    return new(tbl)
  end
end

################################################################################
#
#  Tables of marks
#
################################################################################

abstract type GroupTableOfMarks end

"""
    GAPGroupTableOfMarks <: GroupTableOfMarks

This is the type of tables of marks that can delegate
tasks to an underlying table of marks object in the GAP system
(field `GAPTable`).

A group can (but need not) be stored in the field `group`.
If it is available then also the field `isomorphism` is available,
its value is a bijective map from the `group` value to a group in GAP.

Objects of type `GAPGroupTableOfMarks` support [`get_attribute`](@ref).
"""
@attributes mutable struct GAPGroupTableOfMarks <: GroupTableOfMarks
  GAPTable::GapObj  # the GAP table of marks object
  group::Union{GAPGroup, FinGenAbGroup}    # the underlying group, if any
  isomorphism::Map  # isomorphism from `group` to a group in GAP

  function GAPGroupTableOfMarks(G::Union{GAPGroup, FinGenAbGroup}, tom::GapObj, iso::Map)
    @req GAPWrap.IsTableOfMarks(tom) "tom must be a GAP table of marks"
    return new(tom, G, iso)
  end

  function GAPGroupTableOfMarks(tom::GapObj)
    # group and isomorphism are left undefined
    @req GAPWrap.IsTableOfMarks(tom) "tom must be a GAP table of marks"
    return new(tom)
  end
end

################################################################################
#
#  Marks vectors
#
################################################################################

#  In order to describe Burnside rings over the integers,
#  we implement marks vectors as wrapped GAP vectors,
#  with parent the table of marks in question.
#  Note that marks vectors can be shorter than the number of columns
#  of the table of marks,
#  meaning that the values at larger positions are zero.

abstract type GroupMarksVector end

struct GAPGroupMarksVector <: GroupMarksVector
  table::GAPGroupTableOfMarks
  values::GapObj
end

################################################################################
#
#  Cycles
#
################################################################################

struct CycleType <: AbstractVector{Pair{Int64, Int64}}
  # pairs 'cycle length => number of times it occurs'
  # so 'n => 1' is a single n-cycle and  '1 => n' is the identity on n points
  s::Vector{Pair{Int, Int}}

  # take a vector of cycle lengths
  function CycleType(c::Vector{Int})
    s = Vector{Pair{Int, Int}}()
    for i = c
      _push_cycle!(s, i)
    end
    sort!(s; by=first)
    return new(s)
  end
  function CycleType(v::Vector{Pair{Int, Int}}; sorted::Bool = false)
    sorted && return new(v)
    return new(sort(v; by=first))
    #TODO: check that each cycle length is specified at most once?
  end
end

################################################################################
#
#  Subspaces iterator
#
################################################################################

# For an iterator of `k`-dimensional subspaces in an `n`-dimensional space
# over the field with `q` elements,
# the state is `(comb_iter, choice, prod_iter, values)`
# where
# - `comb_iter` is the iterator `combinations(1:n, k)`,
# - `choice` is the current state of `comb_iter`,
#   a vector that denotes the `k` pivot column positions in the matrix,
# - `prod_iter` is an iterator of vectors, where the current vector
#   corresponds to those entries in the not necessarily zero positions
#   of the matrix that are not pivots, and
# - `values` is the current state ot `prod_iter`.

struct SubspacesIterator{T <: FinFieldElem}
  V::AbstractAlgebra.Generic.AbstractAlgebra.Generic.FreeModule{T}
  k::Int

  function SubspacesIterator{T}(V, k) where T
    @req 0 <= k <= vector_space_dim(V) "wrong dimension"
    return new(V, k)
  end
end

################################################################################
#
#  Words iterator (for printing and showing character tables)
#
################################################################################

# Utility:
# Create strings in length-lexicographical ordering w.r.t. the
# alphabet 'alphabet'.
# (If `alphabet` is `"ABCDEFGHIJKLMNOPQRSTUVWXYZ"` then the strings
# have the form `"A", "B", ..., "Z", "AA", ...`.)
mutable struct WordsIterator
  alphabet::String
end

################################################################################
#
#  Collector
#
################################################################################

# Create an Oscar collector object, its type parameter `T` describes
# the type of integers that occur as exponents.

abstract type Collector{T} end

mutable struct GAP_Collector{T} <: Collector{T}
  ngens::Int
  relorders::Vector{T}
  powers::Vector{Vector{Pair{Int, T}}}
  conjugates::Matrix{Vector{Pair{Int, T}}}
  X::GapObj   # a collector in GAP, if defined
  F::FPGroup  # the free group in Oscar that belongs to `X`, if defined

  function GAP_Collector{T}(n::Int) where T <: IntegerUnion
    return new(n, zeros(T, n), # relative orders undefined
               repeat([Pair{Int, T}[]], n), # default powers are identity
               Matrix{Vector{Pair{Int, T}}}(undef, n, n)) # conjugates undefined
  end
end

################################################################################
#
#  Recognition
#
################################################################################

@attributes mutable struct GroupRecognitionTree{T <: GAPGroup}
  input_group::T
  gap_tree::GapObj

  function GroupRecognitionTree{T}(G::GAPGroup, gap_tree::GapObj) where T
    @req G isa PermGroup || G isa MatGroup "only matrix and permutation groups are supported"
    res = new{T}(G, gap_tree)
    return res
  end
end
