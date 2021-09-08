# further possible functions: similar, literal_pow

import Base: ^, Base.Vector

using Random

export GroupConjClass

export
    comm,
    comm!,
    complement_system, hascomplement_system, setcomplement_system,
    conjugacy_class,
    conjugacy_classes_maximal_subgroups,
    conjugacy_classes_subgroups,
    conjugacy_classes,
    core,
    coset_decomposition,
    cperm,
    degree,
    describe,
    div_left,
    div_left!,
    div_right,
    div_right!,
    elements,
    hasexponent, setexponent,
    fitting_subgroup, hasfitting_subgroup, setfitting_subgroup,
    frattini_subgroup, hasfrattini_subgroup, setfrattini_subgroup,
    gap_perm, # HACK
    gen,
    gens,
    hall_subgroup,
    hall_subgroups_representatives,
    hall_system, hashall_system, sethall_system,
    inv!,
    isalmostsimple, hasisalmostsimple, setisalmostsimple,
    isconjugate,
    isfinite, hasisfinite, setisfinite,
    isfinite_order,
    isperfect, hasisperfect, setisperfect,
    ispgroup,
    issimple, hasissimple, setissimple,
    moved_points, hasmoved_points, setmoved_points,
    mul,
    mul!,
    nilpotency_class, hasnilpotency_class, setnilpotency_class,
    ngens,
    normal_closure,
    normalizer,
    number_conjugacy_classes, hasnumber_conjugacy_classes, setnumber_conjugacy_classes,
    number_moved_points, hasnumber_moved_points, setnumber_moved_points,
    one!,
    order, hasorder, setorder,
    pcore,
    radical_subgroup, hasradical_subgroup, setradical_subgroup,
    rand_pseudo,
    relators,
    representative_action,
    right_coset,
    right_cosets ,
    right_transversal,
    socle, hassocle, setsocle,
    sylow_subgroup,
    sylow_system, hassylow_system, setsylow_system


# TODO: as soon as GAP packages like `polycyclic` or `rcwa` are loaded,
# the custom group types and isos they define should be added to the arrays
# _gap_group_types resp. _iso_function


group_element(G::T, x::GapObj) where T <: GAPGroup = BasicGAPGroupElem{T}(G, x)

function elements(G::T) where T <: GAPGroup
  els = GAP.gap_to_julia(Vector{GapObj},GAP.Globals.Elements(G.X))
  elems = Vector{elem_type(G)}(undef, length(els))
  i = 1
  for x in els
    elems[i] = group_element(G, x)
    i += 1
  end
  return elems
end

function parent(x::GAPGroupElem)
  return x.parent
end

import Base.isfinite

@gapattribute isfinite(G::GAPGroup) = GAP.Globals.IsFinite(G.X)::Bool
"""
    isfinite(G::GAPGroup) -> Bool

Return `true` if `G` is finite, and `false` otherwise.

# Examples
```jldoctest
julia> isfinite(symmetric_group(5))
true

julia> isfinite(free_group(2))
false

```
""" isfinite(G::GAPGroup)

Base.isfinite(G::PermGroup) = true

Base.isfinite(G::PcGroup) = true

"""
    isfinite_order(g::GAPGroupElem) -> Bool

Return `true` if `g` has finite order, and `false` otherwise.

# Examples
```jldoctest
julia> isfinite_order(gen(symmetric_group(5), 1))
true

julia> isfinite_order(gen(free_group(2), 1))
false

```
"""
isfinite_order(x::GAPGroupElem) = GAP.Globals.IsInt(GAP.Globals.Order(x.X))::Bool

"""
    degree(G::PermGroup) -> Int

Return the degree of `G` as a permutation group, that is,
an integer `n` that is stored in `G`, with the following meaning.

- `G` embeds into `symmetric_group(n)`.
- Two permutation groups of different degrees are regarded as not equal,
  even if they contain the same permutations.
- Subgroups constructed with `derived_subgroup`, `sylow_subgroup`, etc.,
  get the same degree as the given group.
- The range `1:degree(G)` is used as the default set of points on which
  `G` and its element acts.

!!! note
    The degree of a group of permutations is not necessarily equal to the largest moved point of the group `G`. For example, the trivial subgroup of `symmetric_group(n)` has degree `n` even though it fixes `n`.

# Examples
```jldoctest
julia> degree(symmetric_group(4))
4

julia> t4 = trivial_subgroup(symmetric_group(4))[1];

julia> degree(t4)
4

julia> t4 == trivial_subgroup(symmetric_group(5))[1]
false

julia> show(Vector(gen(symmetric_group(4), 2)))
[2, 1, 3, 4]
julia> show(Vector(gen(symmetric_group(5), 2)))
[2, 1, 3, 4, 5]
```
"""
degree(x::PermGroup) = x.deg

@gapattribute moved_points(x::Union{PermGroupElem,PermGroup}) = [y for y in GAP.gap_to_julia(GAP.Globals.MovedPoints(x.X))]
# This is more efficient than `Vector{Int}(GAP.Globals.MovedPoints(x.X))`.
# (And note that `GAP.gap_to_julia(GAP.Globals.MovedPoints(G.X))` may return
# a range.)
"""
    moved_points(x::PermGroupElem)
    moved_points(G::PermGroup)

Return the vector of those points in `1:degree(x)` or `1:degree(G)`,
respectively, that are not mapped to themselves under the action `^`.

# Examples
```jldoctest
julia> g = symmetric_group(4);  s = sylow_subgroup(g, 3)[1];

julia> length(moved_points(s))
3

julia> length(moved_points(gen(s, 1)))
3

```
""" moved_points

@gapattribute number_moved_points(x::Union{PermGroupElem,PermGroup}) = fmpz(GAP.Globals.NrMovedPoints(x.X))::fmpz
"""
    number_moved_points(::Type{T} = fmpz, x::PermGroupElem) where T <: Union{Integer, fmpz}
    number_moved_points(::Type{T} = fmpz, G::PermGroup) where T <: Union{Integer, fmpz}

Return the number of those points in `1:degree(x)` or `1:degree(G)`,
respectively, that are not mapped to themselves under the action `^`,
as an instance of `T`.

# Examples
```jldoctest
julia> g = symmetric_group(4);  s = sylow_subgroup(g, 3)[1];

julia> number_moved_points(s)
3

julia> number_moved_points(Int, s)
3

julia> number_moved_points(gen(s, 1))
3

```
""" number_moved_points

number_moved_points(::Type{T}, x::Union{PermGroupElem,PermGroup}) where T <: Union{Base.Integer, fmpz} = T(GAP.Globals.NrMovedPoints(x.X))::T


"""
    order(::Type{T} = fmpz, x::Union{GAPGroupElem, GAPGroup}) where T <: Union{Integer, fmpz}

Return the order of `x`, as an instance of `T`.

For a group element `x` in the group `G`, the order of `x` is the smallest
positive integer `n` such that `x^n` is the identity of `G`.
For a group `x`, the order of `x` is the number of elements in `x`.

An exception is thrown if the order of `x` is infinite,
use [`isfinite`](@ref) in order to check for finiteness.
"""
function order(::Type{T}, x::Union{GAPGroupElem, GAPGroup}) where T <: Union{Base.Integer, fmpz}
   ord = GAP.Globals.Order(x.X)
   if ord === GAP.Globals.infinity
      error("order() not supported for infinite groups, use isfinite()")
   end
   return T(ord)
end

order(x::Union{GAPGroupElem, GAPGroup}) = order(fmpz, x)

@gapwrap hasorder(G::GAPGroup) = GAP.Globals.HasSize(G.X)
@gapwrap setorder(G::GAPGroup, val::T) where T<:Union{Base.Integer,fmpz} = GAP.Globals.SetSize(G.X, GapObj(val))

import Base.exponent

@gapattribute exponent(x::GAPGroup) = fmpz(GAP.Globals.Exponent(x.X))
@doc Markdown.doc"""
    exponent(::Type{T} = fmpz, G::GAPGroup) where T <: Union{Integer, fmpz}

Return the exponent of `G`, as an instance of `T`,
i. e., the smallest positive integer $e$ such that
$g^e$ is the identity of `G` for every $g$ in `G`.
""" exponent(x::GAPGroup)

Base.exponent(::Type{T}, G::GAPGroup) where T <: Union{Base.Integer, fmpz} = T(GAP.Globals.Exponent(G.X))

"""
    rand(rng::Random.AbstractRNG = Random.GLOBAL_RNG, G::Group)

Return a random element of `G`, using the random number generator `rng`.
"""
Base.rand(G::GAPGroup) = Base.rand(Random.GLOBAL_RNG, G)

function Base.rand(rng::Random.AbstractRNG, G::GAPGroup)
   s = GAP.Globals.Random(GAP.wrap_rng(rng), G.X)
   return group_element(G, s)
end

"""
    rand_pseudo(G::Group)

Return a pseudo random element of `G`.  This works faster than `rand`,
but the returned elements are not necessarily uniformly distributed.
"""
function rand_pseudo(G::GAPGroup)
   s = GAP.Globals.PseudoRandom(G.X)
   return group_element(G,s)
end


function _maxgroup(x::T, y::T) where T <: GAPGroup
   # A typical situation should be that the two groups are identical,
   # but GAP's `IsSubset` check is not as cheap as one wants;
   # there is an `IsSubset` method that checks for identity,
   # but it is not always the first choice.
   if x.X === y.X
     return x
   elseif GAP.Globals.IsSubset(x.X, y.X)
     return x
   elseif GAP.Globals.IsSubset(y.X, x.X)
     return y
   else
     error("Not yet implemented")
   end
end

#We need a lattice of groups to implement this properly
function _prod(x::T, y::T) where T <: GAPGroupElem
  G = _maxgroup(parent(x), parent(y))
  return group_element(G, x.X*y.X)
end

Base.:*(x::GAPGroupElem, y::GAPGroupElem) = _prod(x, y)

==(x::PermGroup, y::PermGroup) = x.deg == y.deg && x.X == y.X

==(x::GAPGroup, y::GAPGroup) = x.X == y.X

==(x::PermGroupElem, y::PermGroupElem) = degree(parent(x)) == degree(parent(y)) && x.X == y.X

==(x::T, y::T) where T <: BasicGAPGroupElem = x.X == y.X

"""
    one(G::GAPGroup) -> elem_type(G)

Return the identity of the group `G`.
"""
Base.one(x::GAPGroup) = group_element(x, GAP.Globals.Identity(x.X))

"""
    one(x::GAPGroupElem{T}) -> GAPGroupElem{T}

Return the identity of the parent group of `x`.
"""
Base.one(x::GAPGroupElem) = one(parent(x))
one!(x::GAPGroupElem) = one(parent(x))

Base.show(io::IO, x::GAPGroupElem) = print(io, String(GAP.Globals.StringViewObj(x.X)))
Base.show(io::IO, x::GAPGroup) = print(io, String(GAP.Globals.StringViewObj(x.X)))

Base.isone(x::GAPGroupElem) = GAP.Globals.IsOne(x.X)

Base.inv(x::GAPGroupElem) = group_element(parent(x), GAP.Globals.Inverse(x.X))

inv!(out::GAPGroupElem, x::GAPGroupElem) = inv(x)  #if needed later

Base.:^(x::GAPGroupElem, y::Int) = group_element(parent(x), x.X ^ y)

Base.:^(x::GAPGroupElem, y::fmpz) = Hecke._generic_power(x, y) # TODO: perhaps  let GAP handle this; also handle arbitrary Integer subtypes?

Base.:<(x::PermGroupElem, y::PermGroupElem) = x.X < y.X

Base.:/(x::GAPGroupElem, y::GAPGroupElem) = x*y^-1

mul(x::GAPGroupElem, y::GAPGroupElem) = x*y
mul!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x*y

div_right(x::GAPGroupElem, y::GAPGroupElem) = x*inv(y)
div_left(x::GAPGroupElem, y::GAPGroupElem) = inv(y)*x
div_right!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x*inv(y)
div_left!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = inv(y)*x

Base.conj(x::GAPGroupElem, y::GAPGroupElem) = x^y
Base.conj!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x^y

"""
    comm(x::GAPGroupElem, y::GAPGroupElem)

Return the commutator of `x` and `y`,
which is defined as `x^-1*y^-1*x*y`,
and usually denoted as `[x,y]` in the literature.
"""
comm(x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y
comm!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y

Base.IteratorSize(::Type{<:GAPGroup}) = Base.SizeUnknown()
Base.IteratorSize(::Type{PermGroup}) = Base.HasLength()

function Base.iterate(G::GAPGroup)
  L=GAP.Globals.Iterator(G.X)
  i = GAP.Globals.NextIterator(L)
  return group_element(G, i), L
end

function Base.iterate(G::GAPGroup, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return group_element(G, i), state
end

# need this function just for the iterator
Base.length(x::GAPGroup)::Int = order(x)

"""
    Base.in(g::GAPGroupElem, G::GAPGroup)

Return whether `g` is an element of `G`.
The parent of `g` need not be equal to `G`.
"""
Base.in(g::GAPGroupElem, G::GAPGroup) = GAP.Globals.in(g.X, G.X)

# FIXME: clashes with AbstractAlgebra.perm method
#function perm(L::AbstractVector{<:Base.Integer})
#   return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.julia_to_gap(L)))
#end
# FIXME: use name gap_perm for now
@doc Markdown.doc"""
    gap_perm(L::AbstractVector{<:Base.Integer})

Return the permutation $x$ which maps every $i$ from `1` to $n$` = length(L)`
to `L`$[i]$.
The parent of $x$ is set to [`symmetric_group`](@ref)$(n)$.
An exception is thrown if `L` does not contain every integer from 1 to $n$
exactly once.

# Examples
```jldoctest
julia> gap_perm([2,4,6,1,3,5])
(1,2,4)(3,6,5)
```
"""
function gap_perm(L::AbstractVector{<:Base.Integer})
  return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.GapObj(L)))
end

gap_perm(L::AbstractVector{<:fmpz}) = gap_perm([Int(y) for y in L])

@doc Markdown.doc"""
    perm(G::PermGroup, L::AbstractVector{<:Integer})
    (G::PermGroup)(L::AbstractVector{<:Integer})

Return the permutation $x$ which maps every `i` from 1 to $n$` = length(L)`
to `L`$[i]$.
The parent of $x$ is `G`.
An exception is thrown if $x$ is not contained in `G`
or `L` does not contain every integer from 1 to $n$ exactly once.

For [`gap_perm`](@ref),
the parent group of $x$ is set to [`symmetric_group`](@ref)$(n)$.

# Examples
```jldoctest
julia> perm(symmetric_group(6),[2,4,6,1,3,5])
(1,2,4)(3,6,5)
```
"""
function perm(g::PermGroup, L::AbstractVector{<:Base.Integer})
   x = GAP.Globals.PermList(GAP.julia_to_gap(L))
   if length(L) <= degree(g) && GAP.Globals.IN(x,g.X) 
     return PermGroupElem(g, x)
   end
   throw(ArgumentError("the element does not embed in the group"))
end

perm(g::PermGroup, L::AbstractVector{<:fmpz}) = perm(g, [Int(y) for y in L])

function (g::PermGroup)(L::AbstractVector{<:Base.Integer})
   x = GAP.Globals.PermList(GAP.julia_to_gap(L))
   if length(L) <= degree(g) && GAP.Globals.IN(x,g.X) 
     return PermGroupElem(g, x)
   end
   throw(ArgumentError("the element does not embed in the group"))
end

(g::PermGroup)(L::AbstractVector{<:fmpz}) = g([Int(y) for y in L])

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarly disjoint)
@doc Markdown.doc"""
    cperm(L::AbstractVector{<:T}...) where T <: Union{Integer, fmpz}
    cperm(G::PermGroup, L::AbstractVector{<:T}...)

For given lists $[a_1, a_2, \ldots, a_n], [b_1, b_2, \ldots , b_m], \ldots$
of positive integers, return the
permutation $x = (a_1, a_2, \ldots, a_n) * (b_1, b_2, \ldots, b_m) * \ldots$.
Arrays of the form `[n, n+1, ..., n+k]` can be replaced by `n:n+k`.

The parent of $x$ is `G`.
If `G` is not specified then the parent of $x$ is set to
[`symmetric_group`](@ref)$(n)$,
where $n$ is the largest integer that occurs in an entry of `L`.

An exception is thrown if $x$ is not contained in `G`
or one of the given vectors is empty or contains duplicates.

# Examples
```jldoctest
julia> cperm([1,2,3],4:7)
(1,2,3)(4,5,6,7)

julia> cperm([1,2],[2,3])
(1,3,2)

julia> p = cperm([1,2,3],[7])
(1,2,3)

julia> degree(parent(p))
7

```

At the moment, the input vectors of the function `cperm` need not be disjoint.

!!! warning
    If the function `perm` is evaluated in a vector of integers
    without specifying the group `G`,
    then the returned value is an element of the AbstractAlgebra.jl type
    `Perm{Int}`.
    For this reason, if one wants a permutation of type
    `GAPGroupElem{PermGroup}` without specifying a parent,
    one has to use the function `gap_perm`.
"""
function cperm(L::AbstractVector{T}...) where T <: Union{Base.Integer, fmpz}
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([PermGroupElem(symmetric_group(maximum(y)), GAP.Globals.CycleFromList(GAP.julia_to_gap([Int(k) for k in y]))) for y in L])
#TODO: better create the product of GAP permutations?
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of vectors (not necessarly disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::PermGroup,L::AbstractVector{T}...) where T <: Union{Base.Integer, fmpz}
   if length(L)==0
      return one(g)
   else
      x=prod(y -> GAP.Globals.CycleFromList(GAP.julia_to_gap([Int(k) for k in y])), L)
      if length(L) <= degree(g) && GAP.Globals.IN(x,g.X)
         return PermGroupElem(g, x)
      else
         throw(ArgumentError("the element does not embed in the group"))
      end
   end
end

@deprecate listperm(x::PermGroupElem) Vector(x)

"""
    Vector{T}(x::PermGroupElem, n::Int = x.parent.deg) where T <: Union{Integer, fmpz}
    Vector(x::PermGroupElem, n::Int = x.parent.deg)

Return the list of length `n` that contains `x(i)` at position `i`. If not specified, `T` is set as `Int`.

# Examples
```jldoctest
julia> pi = cperm(1:3)
(1,2,3)
julia> Vector(pi)
3-element Vector{Int64}:
 2
 3
 1
julia> Vector(pi, 2)
2-element Vector{Int64}:
 2
 3
julia> Vector(pi, 4)
4-element Vector{Int64}:
 2
 3
 1
 4
julia> Vector{fmpz}(pi, 2)
2-element Vector{fmpz}:
 2
 3

```
"""
Base.Vector{T}(x::PermGroupElem, n::Int = x.parent.deg) where T <: Union{Base.Integer, fmpz} = T[x(i) for i in 1:n]
Base.Vector(x::PermGroupElem, n::Int = x.parent.deg) = Vector{Int}(x,n)

"""
    gens(G::Group)

Return a vector of generators of the group `G`.
To get the `i`-th generator,
use `G[i]` or `gen(G,i)` (see [`gen`](@ref)) instead of `gens(G)[i]`,
as that is more efficient.

# Examples
```jldoctest
julia> g = symmetric_group(5);  gens(g)
2-element Vector{PermGroupElem}:
 (1,2,3,4,5)
 (1,2)

julia> g[2]
(1,2)

```

!!! note
    The output of `gens(G)` is not, in general, the minimal list of generators for `G`.
"""
function gens(G::GAPGroup)
   L = GAP.Globals.GeneratorsOfGroup(G.X)
   res = Vector{elem_type(G)}(undef, length(L))
   for i = 1:length(res)
     res[i] = group_element(G, L[i])
   end
   return res
end

"""
    gen(G::GAPGroup, i::Integer)

Return the `i`-th element of the vector `gens(G)`.
This is equivalent to `G[i]`, and returns `gens(G)[i]`
but may be more efficient than the latter.

An exception is thrown if `i` is larger than the length of `gens(G)`.
"""
function gen(G::GAPGroup, i::Int)
   L = GAP.Globals.GeneratorsOfGroup(G.X)
   @assert length(L) >= i "The number of generators is lower than the given index"
   return group_element(G, L[i])
end
Base.getindex(G::GAPGroup, i::Int) = gen(G, i)

"""
    ngens(G::GAPGroup) -> Int

Return the length of the vector [`gens`](@ref)`(G)`.

!!! warning "WARNING:" 
    this is *NOT*, in general, the minimum number of generators for G.
"""
ngens(G::GAPGroup) = length(GAP.Globals.GeneratorsOfGroup(G.X))


Base.sign(x::PermGroupElem) = GAP.Globals.SignPerm(x.X)

Base.isless(x::PermGroupElem, y::PermGroupElem) = x<y

#embedding of a permutation in permutation group
function (G::PermGroup)(x::PermGroupElem)
   if !GAP.Globals.IN(x.X,G.X)
      throw(ArgumentError("the element does not embed in the group"))
   end
   return group_element(G, x.X)
end

#evaluation function
function (x::PermGroupElem)(n::T) where T <: Union{Base.Integer,fmpz}
   return T(GAP.Globals.OnPoints(GAP.GapObj(n), x.X))
end

(x::PermGroupElem)(n::Int) = GAP.Globals.OnPoints(n,x.X)

^(n::T, x::PermGroupElem) where T <: Union{Base.Integer,fmpz} = T(GAP.Globals.OnPoints(GAP.GapObj(n), x.X))

^(n::Int, x::PermGroupElem) = GAP.Globals.OnPoints(n,x.X)

################################################################################
#
#   Conjugacy Classes
#
################################################################################

"""
    GroupConjClass

It could be either the conjugacy class of an element or of a subgroup in a group G. It is displayed as
```
     cc = x ^ G
```
where G is a group and x = `representative`(`cc`) is either an element or a subgroup of G.
"""
struct GroupConjClass{T<:GAPGroup, S<:Union{GAPGroupElem,GAPGroup}}
   X::T
   repr::S
   CC::GapObj
end

Base.eltype(::Type{GroupConjClass{T,S}}) where {T,S} = S
Base.hash(x::GroupConjClass, h::UInt) = h # FIXME

function Base.show(io::IO, x::GroupConjClass)
  print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.repr.X)),
            " ^ ",
            GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X.X)))
end

function _conjugacy_class(G, g, cc::GapObj)         # function for assignment
  return GroupConjClass{typeof(G), typeof(g)}(G, g, cc)
end

==(a::GroupConjClass{T, S}, b::GroupConjClass{T, S}) where S where T = a.CC == b.CC 

Base.length(C::GroupConjClass) = GAP.Globals.Size(C.CC)

representative(C::GroupConjClass) = C.repr

@gapattribute number_conjugacy_classes(G::GAPGroup) = GAP.Globals.NrConjugacyClasses(G.X)

# START elements conjugation

"""
    conjugacy_class(G::Group, g::GAPGroupElem) -> GroupConjClass

Return the conjugacy class `cc` of `g` in `G`, where `g` = `representative`(`cc`).
"""
function conjugacy_class(G::GAPGroup, g::GAPGroupElem)
   return _conjugacy_class(G, g, GAP.Globals.ConjugacyClass(G.X,g.X))
end

function Base.rand(C::GroupConjClass{S,T}) where S where T<:GAPGroupElem
   return Base.rand(Random.GLOBAL_RNG, C)
end

function Base.rand(rng::Random.AbstractRNG, C::GroupConjClass{S,T}) where S where T<:GAPGroupElem
   return group_element(C.X, GAP.Globals.Random(GAP.wrap_rng(rng), C.CC))
end

@deprecate elements(C::GroupConjClass) collect(C)

"""
    conjugacy_classes(G::Group)

Return the vector of all conjugacy classes of elements in G.
It is guaranteed that the class of the identity is in the first position.
"""
function conjugacy_classes(G::GAPGroup)
   L=GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClasses(G.X))
   return GroupConjClass{typeof(G), elem_type(G)}[ _conjugacy_class(G,group_element(G,GAP.Globals.Representative(cc)),cc) for cc in L]
end

Base.:^(x::T, y::T) where T <: GAPGroupElem = group_element(_maxgroup(parent(x), parent(y)), x.X ^ y.X)

@doc Markdown.doc"""
    isconjugate(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem)

Return whether `x` and `y` are conjugate elements in `G`,
i. e., there is an element $z$ in `G` such that `x^`$z$ equals `y`.
"""
isconjugate(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem) = GAP.Globals.IsConjugate(G.X,x.X,y.X)

"""
    representative_action(G::Group, x::GAPGroupElem, y::GAPGroupElem)

If `x` and `y` are conjugate in `G`,
return `(true, z)`, where `x^z == y` holds;
otherwise, return `(false, nothing)`.
"""
function representative_action(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem)
   conj = GAP.Globals.RepresentativeAction(G.X, x.X, y.X)
   if conj != GAP.Globals.fail
      return true, group_element(G, conj)
   else
      return false, nothing
   end
end
# END elements conjugation 

# START subgroups conjugation
"""
    conjugacy_class(G::T, H::T) where T<:Group -> GroupConjClass

Return the subgroup conjugacy class `cc` of `H` in `G`, where `H` = `representative`(`cc`).
"""
function conjugacy_class(G::T, g::T) where T<:GAPGroup
   return _conjugacy_class(G, g, GAP.Globals.ConjugacyClassSubgroups(G.X,g.X))
end

function Base.rand(C::GroupConjClass{S,T}) where S where T<:GAPGroup
   return Base.rand(Random.GLOBAL_RNG, C)
end

function Base.rand(rng::Random.AbstractRNG, C::GroupConjClass{S,T}) where S where T<:GAPGroup
   return T(GAP.Globals.Random(GAP.wrap_rng(rng), C.CC))
end

"""
    conjugacy_classes_subgroups(G::Group)

Return the vector of all conjugacy classes of subgroups of G.
"""
function conjugacy_classes_subgroups(G::GAPGroup)
   L=GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClassesSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,typeof(G)(GAP.Globals.Representative(cc)),cc) for cc in L]
end

"""
    conjugacy_classes_maximal_subgroups(G::Group)

Return the vector of all conjugacy classes of maximal subgroups of G.
"""
function conjugacy_classes_maximal_subgroups(G::GAPGroup)
  L = GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClassesMaximalSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,typeof(G)(GAP.Globals.Representative(cc)),cc) for cc in L]
end

Base.:^(H::GAPGroup, y::GAPGroupElem) = typeof(H)(H.X ^ y.X)

function conjugate_subgroup(G::T, x::GAPGroupElem) where T<:GAPGroup
  return T(GAP.Globals.ConjugateSubgroup(G.X,x.X))
end

"""
    isconjugate(G::GAPGroup, H::GAPGroup, K::GAPGroup)

Return whether `H` and `K` are conjugate subgroups in `G`.
"""
isconjugate(G::GAPGroup, H::GAPGroup, K::GAPGroup) = GAP.Globals.IsConjugate(G.X,H.X,K.X)

"""
    representative_action(G::Group, H::Group, K::Group)

If `H` and `K` are conjugate subgroups in `G`, return `true, z`
where `H^z = K`; otherwise, return `false, nothing`.
```
"""
function representative_action(G::GAPGroup, H::GAPGroup, K::GAPGroup)
   conj = GAP.Globals.RepresentativeAction(G.X, H.X, K.X)
   if conj != GAP.Globals.fail
      return true, group_element(G, conj)
   else
      return false, nothing
   end
end

# END subgroups conjugation


# START iterator
Base.IteratorSize(::Type{<:GroupConjClass}) = Base.SizeUnknown()

Base.iterate(cc::GroupConjClass) = iterate(cc, GAP.Globals.Iterator(cc.CC))

function Base.iterate(cc::GroupConjClass{S,T}, state::GapObj) where {S,T}
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  if T <: GAPGroupElem
     return group_element(cc.X, i), state
  else
     return _as_subgroup(cc.X, i)[1], state
  end
end

################################################################################
#
# Normal Structure
#
################################################################################

"""
    normalizer(G::Group, H::Group)

Return `N, f`, where `N` is the normalizer of `H` in `G`,
i. e., the largest subgroup of `G` in which `H` is normal,
and `f` is the embedding morphism of `N` into `G`.
"""
normalizer(G::T, H::T) where T<:GAPGroup = _as_subgroup(G, GAP.Globals.Normalizer(G.X,H.X))

"""
    normalizer(G::Group, x::GAPGroupElem)

Return `N, f`, where `N` is the normalizer of the cyclic subgroup generated
by `x` in `G` and `f` is the embedding morphism of `N` into `G`.
"""
normalizer(G::GAPGroup, x::GAPGroupElem) = _as_subgroup(G, GAP.Globals.Normalizer(G.X,x.X))

"""
    core(G::Group, H::Group)

Return `C, f`, where `C` is the normal core of `H` in `G`,
that is, the largest normal subgroup of `G` that is contained in `H`,
and `f` is the embedding morphism of `C` into `G`.
"""
core(G::T, H::T) where T<:GAPGroup = _as_subgroup(G, GAP.Globals.Core(G.X,H.X))

"""
    normal_closure(G::Group, H::Group)

Return `N, f`, where `N` is the normal closure of `H` in `G`,
that is, the smallest normal subgroup of `G` that contains `H`,
and `f` is the embedding morphism of `N` into `G`.

Note that `H` must be a subgroup of `G`.
"""
normal_closure(G::T, H::T) where T<:GAPGroup = _as_subgroup(G, GAP.Globals.NormalClosure(G.X,H.X))

# Note:
# GAP admits `NormalClosure` also when `H` is not a subgroup of `G`,
# and in this case the result is not contained in `G`.
# (We should test whether `H` is a subgroup of `G`,
# but then the user should have the possibility to omit this check.)

"""
    pcore(G::Group, p::Int64)

Return `C, f`, where `C` is the `p`-core
(i.e. the largest normal `p`-subgroup) of `G`
and `f` is the embedding morphism of `C` into `G`.
"""
function pcore(G::GAPGroup, p::Int64)
   if !isprime(p)
      throw(ArgumentError("p is not a prime"))
   end
   return _as_subgroup(G, GAP.Globals.PCore(G.X,p))
end



################################################################################
#
# Specific Subgroups
#
################################################################################

# commutator_subgroup(G::T, H::T) where T<:GAPGroup = T(GAP.Globals.CommutatorSubgroup(G.X,H.X))
# In the literature, the name commutator subgroup is often used as a synonym
# of derived subgroup.
# GAP defines `CommutatorSubgroup( G, H )` for arbitrary groups `G`, `H` in
# the same family; thus the result is in general not contained in `G` or `H`,
# and we do not know into which group the result should be embedded.
# Is this function useful at all?
# (The name `Commutator*Subgroup*` is irritating, isn't it?)

@gapattribute fitting_subgroup(G::GAPGroup) = _as_subgroup(G, GAP.Globals.FittingSubgroup(G.X))
"""
    fitting_subgroup(G::GAPGroup)

Return the Fitting subgroup of `G`, i.e.,
the largest nilpotent normal subgroup of `G`.
""" fitting_subgroup

@gapattribute frattini_subgroup(G::GAPGroup) = _as_subgroup(G, GAP.Globals.FrattiniSubgroup(G.X))
"""
    frattini_subgroup(G::GAPGroup)

Return the Frattini subgroup of `G`, i.e.,
the intersection of all maximal subgroups of `G`.
""" frattini_subgroup

@gapattribute radical_subgroup(G::GAPGroup) = _as_subgroup(G, GAP.Globals.RadicalGroup(G.X))
"""
    radical_subgroup(G::GAPGroup)

Return the solvable radical of `G`, i.e.,
the largest solvable normal subgroup of `G`.
""" radical_subgroup
#T wrong name, already in GAP!

@gapattribute socle(G::GAPGroup) = _as_subgroup(G, GAP.Globals.Socle(G.X))
"""
    socle(G::GAPGroup)

Return the socle of `G`, i.e.,
the subgroup generated by all minimal normal subgroups of `G`,
see [`minimal_normal_subgroups`](@ref).
""" socle


################################################################################
#
# Sylow & Hall Subgroups
#
################################################################################

"""
    sylow_subgroup(G::Group, p::Int64)

Return a Sylow `p`-subgroup of the finite group `G`, for a prime `p`.
This is a subgroup of `p`-power order in `G`
whose index in `G` is coprime to `p`.

# Examples
```jldoctest
julia> g = symmetric_group(4); order(g)
24

julia> s = sylow_subgroup(g, 2); order(s[1])
8

julia> s = sylow_subgroup(g, 3); order(s[1])
3

```
"""
function sylow_subgroup(G::GAPGroup, p::Int64)
   if !isprime(p)
      throw(ArgumentError("p is not a prime"))
   end
   return _as_subgroup(G,GAP.Globals.SylowSubgroup(G.X,p))
end

# no longer documented, better use `hall_subgroups_representatives`
function hall_subgroup(G::GAPGroup, P::AbstractVector{<:Base.Integer})
   P = unique(P)
   for p in P
      if !isprime(p)
         throw(ArgumentError("The integers must be prime"))
      end
   end
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroup(G,GAP.Globals.HallSubgroup(G.X,GAP.julia_to_gap(P)))
end

"""
    hall_subgroups_representatives(G::Group, P::Vector{Int})

Return a vector that contains representatives of conjugacy classes of
Hall `P`-subgroups of the finite group `G`, for a vector `P` of primes.
A Hall `P`-subgroup of `G` is a subgroup the order of which is only divisible
by primes in `P` and whose index in `G` is coprime to all primes in `P`.

For solvable `G`, Hall `P`-subgroups exist and are unique up to conjugacy.
For nonsolvable `G`, Hall `P`-subgroups may not exist or may not be unique
up to conjugacy.

# Examples
```jldoctest
julia> g = dihedral_group(30);

julia> h = hall_subgroups_representatives(g, [2, 3]);

julia> (length(h), order(h[1]))
(1, 6)

julia> g = GL(3, 2)
GL(3,2)

julia> h = hall_subgroups_representatives(g, [2, 3]);

julia> (length(h), order(h[1]))
(2, 24)

julia> h = hall_subgroups_representatives(g, [2, 7]); length(h)
0

```
"""
function hall_subgroups_representatives(G::GAPGroup, P::AbstractVector{<:Base.Integer})
   P = unique(P)
   for p in P
      if !isprime(p)
         throw(ArgumentError("The integers must be prime"))
      end
   end
   res_gap = GAP.Globals.HallSubgroup(G.X, GAP.julia_to_gap(P))
   if res_gap == GAP.Globals.fail
     return typeof(G)[]
   elseif GAP.Globals.IsList(res_gap)
     return _as_subgroups(G, res_gap)
   else
     return [_as_subgroup_bare(G, res_gap)]
   end
end

@gapattribute function sylow_system(G::GAPGroup)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroups(G, GAP.Globals.SylowSystem(G.X))
end
@doc Markdown.doc"""
    sylow_system(G::Group)

Return a vector of Sylow $p$-subgroups of the finite group `G`,
where $p$ runs over the prime factors of the order of `G`,
such that every two such subgroups commute with each other (as subgroups).

Sylow systems exist only for solvable groups,
an exception is thrown if `G` is not solvable.
""" sylow_system

@gapattribute function complement_system(G::GAPGroup)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroups(G, GAP.Globals.ComplementSystem(G.X))
end
@doc Markdown.doc"""
    complement_system(G::Group)

Return a vector of $p'$-Hall subgroups of the finite group `G`,
where $p$ runs over the prime factors of the order of `G`.

Complement systems exist only for solvable groups,
an exception is thrown if `G` is not solvable.
""" complement_system

@gapattribute function hall_system(G::GAPGroup)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroups(G, GAP.Globals.HallSystem(G.X))
end
@doc Markdown.doc"""
    hall_system(G::Group)

Return a vector of $P$-Hall subgroups of the finite group `G`,
where $P$ runs over the subsets of prime factors of the order of `G`.

Hall systems exist only for solvable groups,
an exception is thrown if `G` is not solvable.
""" hall_system


################################################################################
#
# Some Properties
#
################################################################################

@gapattribute isperfect(G::GAPGroup) = GAP.Globals.IsPerfectGroup(G.X)::Bool
"""
    isperfect(G)

Return whether `G` is a perfect group, i.e., equal to its derived subgroup.
""" isperfect

@gapattribute issimple(G::GAPGroup) = GAP.Globals.IsSimpleGroup(G.X)::Bool
"""
    issimple(G)

Return whether `G` is a simple group, i.e.,
`G` is not trivial and has no non-trivial normal subgroups.
""" issimple

@gapattribute isalmostsimple(G::GAPGroup) = GAP.Globals.IsAlmostSimpleGroup(G.X)::Bool
@doc Markdown.doc"""
    isalmostsimple(G)

Return whether `G` is an almost simple group,
i. e., `G` is isomorphic to a group $H$ with the property
$S \leq H \leq Aut(S)$, for some non-abelian simple group $S$.
""" isalmostsimple

"""
    ispgroup(G)

Return `(true, nothing)` if `G` is the trivial group,
`(true, p)` if the order of every element in `G` is a power of a prime `p`,
and `(false, nothing)` otherwise.

For finite groups `G`, the first return value is `true` if and only if
the order of `G` is a prime power.
"""
function ispgroup(G::GAPGroup)
   if GAP.Globals.IsPGroup(G.X)
      p = GAP.Globals.PrimePGroup(G.X)
      if p != GAP.Globals.fail
         return true, p
      else
         return true, nothing
      end
   end
   return false, nothing
end

@doc Markdown.doc"""
    relators(G::FPGroup)

Return a vector of relators for the finitely presented group, i. e.,
elements $[x_1, x_2, \ldots, x_n]$ in $F =$ `free_group(ngens(G))` such that
`G` is isomorphic with $F/[x_1, x_2, \ldots, x_n]$.
"""
function relators(G::FPGroup)
   L=GAP.Globals.RelatorsOfFpGroup(G.X)
   F=free_group(G)
   return [group_element(F,L[i]) for i in 1:length(L)]
end

@gapattribute function nilpotency_class(G::GAPGroup)
   @assert isnilpotent(G) "The group is not nilpotent."
   return GAP.Globals.NilpotencyClassOfGroup(G.X)
end

@doc Markdown.doc"""
    nilpotency_class(G::GAPGroup) -> Int

Return the nilpotency class of `G`, i.e.,
the smallest integer $n$ such that `G` has a central series of length $n$.

An exception is thrown if `G` is not nilpotent.
""" nilpotency_class(G::GAPGroup)

#
describe(G::GAPGroup) = String(GAP.Globals.StructureDescription(G.X))
