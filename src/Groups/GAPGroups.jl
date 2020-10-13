# further possible functions: similar, literal_pow, parent_type

import Base: ^, Base.Vector

export GroupConjClass

export
    comm,
    comm!,
    complement_system,
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
    fitting_subgroup,
    frattini_subgroup,
    gap_perm, # HACK
    gen,
    gens,
    hall_subgroup,
    hall_system,
    inv!,
    isalmostsimple,
    isconjugate,
    isfinite,
    isperfect,
    ispgroup,
    issimple,
    listperm,
    mul,
    mul!,
    nilpotency_class,
    ngens,
    normal_closure,
    normaliser,
    normalizer,
    number_conjugacy_classes,
    one!,
    order,
    pcore,
    radical_subgroup,
    rand_pseudo,
    relators,
    right_coset,
    right_cosets ,
    right_transversal,
    socle,
    sylow_subgroup,
    sylow_system


# TODO: as soon as GAP packages like `polycyclic` or `rcwa` are loaded,
# the custom group types and isos they define should be added to the arrays
# _gap_group_types resp. _iso_function


function group_element(G::T, x::GapObj) where T <: GAPGroup
  if T<:MatrixGroup return MatrixGroupElem(G,x)
  else return BasicGAPGroupElem{T}(G, x)
  end
end

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

function Base.isfinite(G::PermGroup)
  return true
end

function Base.isfinite(G::PcGroup)
  return true
end

function Base.isfinite(G::GAPGroup)
  return GAP.Globals.IsFinite(G.X)
end

"""
    degree(G::PermGroup)

Return the degree as permutation group, that is the integer `n` such that `G < Sym(n)`.

!!! warning "Note"
    This is *not* the smallest `k` such that `G` embeds in `Sym(k)`.
"""
function degree(x::PermGroup)
   return x.deg
end

function order(x::Union{GAPGroupElem, GAPGroup})
   return GAP.gap_to_julia(GAP.Globals.Order(x.X))
end

function order(::Type{T}, x::Union{GAPGroupElem, GAPGroup}) where T<:Number
   return GAP.gap_to_julia(T, GAP.Globals.Order(x.X))
end

"""
    exponent(G::Group)

Return the exponent of `G`, i.e. the smallest positive integer `e` such that `g`^`e`=1 for every `g` in `G`.
"""
Base.exponent(x::GAPGroup) = GAP.Globals.Exponent(x.X)

"""
    rand(G::Group)

Return a random element of the group `G`.
"""
function Base.rand(x::GAPGroup)
   s = GAP.Globals.Random(x.X)
   return group_element(x, s)
end

"""
    rand_pseudo(G::Group)

Return a random element of the group `G`. It works faster than `rand`, but the elements are not necessarily equally distributed.
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

function ==(x::PermGroup, y::PermGroup)
   return x.X == y.X && x.deg == y.deg
end

function ==(x::GAPGroup, y::GAPGroup)
   return x.X == y.X
end

function ==(x::PermGroupElem, y::PermGroupElem)
   return x.X == y.X && degree(parent(x))==degree(parent(y))
end

function ==(x::T, y::T) where T <: BasicGAPGroupElem
   return x.X == y.X
end

"""
    one(G::Group) -> x::GAPGroupElem{typeof(G)}

Return the identity of the group `G`.
"""
Base.one(x::GAPGroup) = group_element(x, GAP.Globals.Identity(x.X))

"""
    one(x::GAPGroupElement{T}) -> x::GAPGroupElem{T}

Return the identity of the parent group of x.
"""
Base.one(x::GAPGroupElem) = one(parent(x))
one!(x::GAPGroupElem) = one(parent(x))

Base.show(io::IO, x::GAPGroupElem) =  print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X)))
Base.show(io::IO, x::GAPGroup) = print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X)))

Base.isone(x::GAPGroupElem) = x == one(parent(x))

Base.inv(x::GAPGroupElem) = group_element(parent(x), GAP.Globals.Inverse(x.X))

inv!(out::GAPGroupElem, x::GAPGroupElem) = inv(x)  #if needed later

Base.:^(x::GAPGroupElem, y::Int) = group_element(parent(x), x.X ^ y)

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

Return `[x,y]` = `x^-1*y^-1*x*y`.
"""
comm(x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y
comm!(out::GAPGroupElem, x::GAPGroupElem, y::GAPGroupElem) = x^-1*x^y

Base.IteratorSize(::Type{<:GAPGroup}) = Base.SizeUnknown()
Base.IteratorSize(::Type{PermGroup}) = Base.HasLength()

function Base.iterate(G::GAPGroup)
  L=GAP.Globals.Iterator(G.X)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
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

Base.in(g::GAPGroupElem, G::GAPGroup) = GAP.Globals.in(g.X, G.X)

# FIXME: clashes with AbstractAlgebra.perm method
#function perm(L::AbstractVector{<:Base.Integer})
#   return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.julia_to_gap(L)))
#end
# FIXME: use name gap_perm for now
"""
    gap_perm(L::AbstractVector{<:Integer})

Return the permutation `x` which maps every `i` from `1` to `length(L)` to `L[i]`. `L` must contain every integer from 1 to `length(L)` exactly, otherwise an exception is thrown.
The parent of `x` is set as Sym(`n`).
```jldoctest
julia> gap_perm([2,4,6,1,3,5])
(1,2,4)(3,6,5)
```
"""
function gap_perm(L::AbstractVector{<:Base.Integer})
  return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.julia_to_gap(L)))
end


"""
    perm(G::PermGroup, L::AbstractVector{<:Integer})
    (G::PermGroup)(L::AbstractVector{<:Integer})

Return the permutation `x` which maps every `i` from `1` to `length(L)` to `L[i]`. `L` must contain every integer from 1 to `length(L)` exactly, otherwise an exception is thrown.
The parent of `x` is `G`. If `x` is not contained in `G`, an ERROR is returned. For `gap_perm`, the parent group of `x` is set as Sym(`n`), where `n` is the largest moved point of `x`. Example:
```jldoctest
julia> perm(symmetric_group(6),[2,4,6,1,3,5])
(1,2,4)(3,6,5)
```
"""
function perm(g::PermGroup, L::AbstractVector{<:Base.Integer})
   x = GAP.Globals.PermList(GAP.julia_to_gap(L))
   if GAP.Globals.IN(x,g.X) 
     return PermGroupElem(g, x)
   end
   throw(ArgumentError("the element does not embed in the group"))
end

function (g::PermGroup)(L::AbstractVector{<:Base.Integer})
   x = GAP.Globals.PermList(GAP.julia_to_gap(L))
   if GAP.Globals.IN(x,g.X) 
     return PermGroupElem(g, x)
   end
   throw(ArgumentError("the element does not embed in the group"))
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
"""
    cperm(L::AbstractVector{<:Integer}...)
    cperm(G::PermGroup, L::AbstractVector{<:Integer}...)

For given lists of positive integers `[a_1, a_2, ..., a_n],[b_1, b_2, ... , b_m], ...` return the
permutation `x = (a_1,a_2,...,a_n)(b_1,b_2,...,b_m)...`. The array `[n,n+1,...,n+k]` can be replaced by `n:n+k`.
  
If a list is empty or contains duplicates, it fails.
The parent of `x` is `G`. If `x` is not contained in `G`, an ERROR is returned. If `G` is not specified, then the parent of `x` is set as Sym(`n`), where `n` is the largest moved point of `x`. Example:
```jldoctest;
julia> cperm([1,2,3],4:7)
(1,2,3)(4,5,6,7)

julia> cperm([1,2],[2,3])
(1,3,2)
```
"""
function cperm(L::AbstractVector{<:Base.Integer}...)
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([PermGroupElem(symmetric_group(maximum(y)), GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y)))) for y in L])
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::PermGroup,L::AbstractVector{<:Base.Integer}...)
   if length(L)==0
      return one(g)
   else
      x=GAP.Globals.Product(GAP.julia_to_gap([GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y))) for y in L]))
      if GAP.Globals.IN(x,g.X) return PermGroupElem(g, x)
      else throw(ArgumentError("the element does not embed in the group"))
      end
   end
end

"""
    listperm(x::PermGroupElem)

Return the list L defined by L = [ `x`(i) for i in 1:n ], where `n` is the degree of `parent(x)`.
"""
function listperm(x::PermGroupElem)
   return [x(i) for i in 1:x.parent.deg]
end
#TODO: Perhaps omit `listperm` and use just `Vector`?

"""
    Vector{T}(x::PermGroupElem, n::Int = x.parent.deg) where {T}

Return the list of length `n` that contains `x(i)` at position `i`.
"""
Base.Vector{T}(x::PermGroupElem, n::Int = x.parent.deg) where {T} = T[x(i) for i in 1:n]

"""
    gens(G::Group)

Return an array of generators of the group `G`. To get a specific generator,
use `G[i]` or `gen(G,i)` instead of `gens(G)[i]`, as that is more efficient.
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
    gen(G::Group, i::Integer)

Return the `i`-th element of the array gens(`G`). This is equivalent to `G[i]`, and returns `gens(G)[i]` but may be more efficient than the latter. If `i` is greater than the length of gens(`G`), an ERROR is thrown.
"""
function gen(G::GAPGroup, i::Int)
   L = GAP.Globals.GeneratorsOfGroup(G.X)
   @assert length(L) >= i "The number of generators is lower than the given index"
   return group_element(G, L[i])
end

"""
    ngens(G::Group) -> Int

Return the length of the array gens(G).

!!! warning "WARNING:" 
    this is *NOT*, in general, the minimum number of generators for G.
"""
ngens(G::GAPGroup) = length(GAP.Globals.GeneratorsOfGroup(G.X))


Base.getindex(G::GAPGroup, i::Int) = gen(G, i)
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
function (x::PermGroupElem)(n::Int)
   return GAP.Globals.OnPoints(n,x.X)
end

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

number_conjugacy_classes(G::GAPGroup) = GAP.Globals.NrConjugacyClasses(G.X)

# START elements conjugation

"""
    conjugacy_class(G::Group, g::GAPGroupElem) -> GroupConjClass

Return the conjugacy class `cc` of `g` in `G`, where `g` = `representative`(`cc`).
"""
function conjugacy_class(G::GAPGroup, g::GAPGroupElem)
   return _conjugacy_class(G, g, GAP.Globals.ConjugacyClass(G.X,g.X))
end

function Base.rand(C::GroupConjClass{S,T}) where S where T<:GAPGroupElem
   return group_element(C.X, GAP.Globals.Random(C.CC))
end

"""
    elements(C::GroupConjClass)

Return the array of the elements in C.
"""
function elements(C::GroupConjClass{S, T}) where S where T<:GAPGroupElem
   L=GAP.Globals.AsList(C.CC)
   l = Vector{T}(undef, length(L))
   for i in 1:length(l)
      l[i] = group_element(C.X,L[i])
   end
   return l
end

"""
    conjugacy_classes(G::Group)

Return the array of all conjugacy classes of elements in G. It is guaranteed that the class of the identity is in the first position.
"""
function conjugacy_classes(G::GAPGroup)
   L=GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClasses(G.X))
   return GroupConjClass{typeof(G), elem_type(G)}[ _conjugacy_class(G,group_element(G,GAP.Globals.Representative(cc)),cc) for cc in L]
end

Base.:^(x::GAPGroupElem, y::GAPGroupElem) = group_element(_maxgroup(parent(x), parent(y)), x.X ^ y.X)

"""
    isconjugate(G::Group, x::GAPGroupElem, y::GAPGroupElem)

If `x`,`y` are conjugate in `G`, return 
```
true, z
```
where `x^z=y`; otherwise, return
```
false, nothing
```
"""
function isconjugate(G::GAPGroup, x::GAPGroupElem, y::GAPGroupElem)
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
   return T(GAP.Globals.Random(C.CC))
end

function elements(C::GroupConjClass{S, T}) where S where T<:GAPGroup
   L=GAP.Globals.AsList(C.CC)
   l = Vector{T}(undef, length(L))
   for i in 1:length(l)
      l[i] = _as_subgroup(C.X, L[i])[1]
   end
   return l
end

"""
    conjugacy_classes_subgroups(G::Group)

Return the array of all conjugacy classes of subgroups of G.
"""
function conjugacy_classes_subgroups(G::GAPGroup)
   L=GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClassesSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,typeof(G)(GAP.Globals.Representative(cc)),cc) for cc in L]
end

"""
    conjugacy_classes_maximal_subgroups(G::Group)

Return the array of all conjugacy classes of maximal subgroups of G.
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
    isconjugate(G::Group, H::Group, K::Group)

If `H`,`K` are conjugate subgroups in `G`, return 
```
true, z
```
where `H^z=K`; otherwise, return
```
false, nothing
```
"""
function isconjugate(G::GAPGroup, H::GAPGroup, K::GAPGroup)
   conj = GAP.Globals.RepresentativeAction(G.X, H.X, K.X)
   if conj != GAP.Globals.fail
      return true, group_element(G, conj)
   else
      return false, nothing
   end
end

# END subgroups conjugation


################################################################################
#
# Normal Structure
#
################################################################################

"""
    normalizer(G::Group, H::Group)
    normaliser(G::Group, H::Group)

Return `N,f`, where `N` is the normalizer of `H` in `G` and `f` is the embedding morphism of `N` into `G`.
"""
normalizer(G::T, H::T) where T<:GAPGroup = _as_subgroup(G, GAP.Globals.Normalizer(G.X,H.X))

"""
    normalizer(G::Group, x::GAPGroupElem)
    normaliser(G::Group, x::GAPGroupElem)

Return `N,f`, where `N` is the normalizer of <`x`> in `G` and `f` is the embedding morphism of `N` into `G`.
"""
normalizer(G::GAPGroup, x::GAPGroupElem) = _as_subgroup(G, GAP.Globals.Normalizer(G.X,x.X))

normaliser = normalizer

"""
    core(G::Group, H::Group)

Return `C,f`, where `C` is the normal core of `H` in `G`,
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

Return `C,f`, where `C` is the `p`-core (i.e. the largest normal `p`-subgroup) of `G` and `f` is the embedding morphism of `C` into `G`.
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
# we don't know how G,H embed into [G,H]

"""
    fitting_subgroup(G::GAPGroup)

Return the Fitting subgroup of `G`, the largest nilpotent normal subgroup of `G`.
"""
fitting_subgroup(G::GAPGroup) = _as_subgroup(G, GAP.Globals.FittingSubgroup(G.X))

"""
    frattini_subgroup(G::GAPGroup)

Return the Frattini subgroup of `G`, the intersection of all maximal subgroups of `G`.
"""
frattini_subgroup(G::GAPGroup) = _as_subgroup(G, GAP.Globals.FrattiniSubgroup(G.X))

"""
    radical_subgroup(G::GAPGroup)

Return the radical subgroup of `G`, the largest solvable normal subgroup of `G`.
"""
radical_subgroup(G::GAPGroup) = _as_subgroup(G, GAP.Globals.RadicalGroup(G.X))

"""
    socle(G::GAPGroup)

Return the socle of `G`, the subgroup generated by all minimal normal subgroups of `G`.
"""
socle(G::GAPGroup) = _as_subgroup(G, GAP.Globals.Socle(G.X))


################################################################################
#
# Sylow & Hall Subgroups
#
################################################################################

"""
    sylow_subgroup(G::Group, p::Int64)

Return a Sylow `p`-subgroup of `G`.
"""
function sylow_subgroup(G::GAPGroup, p::Int64)
   if !isprime(p)
      throw(ArgumentError("p is not a prime"))
   end
   return _as_subgroup(G,GAP.Globals.SylowSubgroup(G.X,p))
end

"""
    hall_subgroup(G::Group, P::Array{Int64})

Return a Hall `P`-subgroup of `G`. It works only if `G` is solvable.
"""
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
    sylow_system(G::Group)

Return an array of Sylow ``p``-subgroups of `G`, where ``p`` runs over the prime factors of |`G`|, such that every two such subgroups commute each other (as subgroups). It works only if `G` is solvable.
"""
function sylow_system(G::GAPGroup)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroups(G, GAP.Globals.SylowSystem(G.X))
end

"""
    complement_system(G::Group)

Return an array of ``p'``-Hall subgroups of `G`, where ``p`` runs over the prime factors of |`G`|. It works only if `G` is solvable.
"""
function complement_system(G::GAPGroup)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroups(G, GAP.Globals.ComplementSystem(G.X))
end

"""
    hall_system(G::Group)

Return an array of ``P``-Hall subgroups of `G`, where ``P`` runs over the subsets of prime factors of |`G`|. It works only if `G` is solvable.
"""
function hall_system(G::GAPGroup)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroups(G, GAP.Globals.HallSystem(G.X))
end


################################################################################
#
# Some Properties
#
################################################################################

"""
    isperfect(G)

Return whether `G` is a perfect group, i.e., equal to its derived subgroup.
"""
isperfect(G::GAPGroup) = GAP.Globals.IsPerfectGroup(G.X)

"""
    issimple(G)

Return whether `G` is a simple group, i.e., it has not non-trivial normal subgroups.
"""
issimple(G::GAPGroup) = GAP.Globals.IsSimpleGroup(G.X)

"""
    isalmostsimple(G)

Return whether `G` is an almost simple group, i.e., if `S` < `G` < `Aut(S)` for some non-abelian simple group `S`.
"""
isalmostsimple(G::GAPGroup) = GAP.Globals.IsAlmostSimpleGroup(G.X)

"""
    ispgroup(G)

Return (``true``,``p``) if |`G`| is a non-trivial ``p``-power, (``false``,nothing) otherwise.
"""
function ispgroup(G::GAPGroup)
   if GAP.Globals.IsPGroup(G.X)
      p = GAP.Globals.PrimePGroup(G.X)
      if p != GAP.Globals.fail
         return true, p
      end
   end
   return false, nothing
end

function relators(G::FPGroup)
   L=GAP.Globals.RelatorsOfFpGroup(G.X)
   F=free_group(G)
   return [group_element(F,L[i]) for i in 1:length(L)]
end

function nilpotency_class(G::GAPGroup)
   @assert isnilpotent(G) "The group is not nilpotent."
   return GAP.Globals.NilpotencyClassOfGroup(G.X)
end

#
describe(G::GAPGroup) = GAP.gap_to_julia(GAP.Globals.StructureDescription(G.X))
