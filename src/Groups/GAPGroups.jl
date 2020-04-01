# further possible functions: similar, iterate, literal_pow, parent_type

export order, perm, cperm, isfinite, gens, ngens, comm, comm!, inv!, rand_pseudo, one!, div_right,
div_left, div_right!, div_left!, deg, mul, mul!, listperm, degree, elements, right_coset, coset_decomposition,
right_cosets , right_transversal, conjugacy_class, conjugacy_classes, number_conjugacy_classes, isconjugate
     #conj!, conj


# TODO: as soon as GAP packages like `polycyclic` or `rcwa` are loaded,
# the custom group types and isos they define should be added to the arrays
# _gap_group_types resp. _iso_function


function group_element(G::T, x::GapObj) where T <: Group
  return GroupElem{T}(G, x)
end

function elements(G::T) where T <: Group
  els = GAP.gap_to_julia(GAP.Globals.Elements(G.X))
  elems = Vector{elem_type(G)}(undef, length(els))
  i = 1
  for x in els
    elems[i] = group_element(G, x)
    i += 1
  end
  return elems
end

# to be fixed later
Base.:hash(x::GroupElem) = 0

function parent(x::GroupElem)
  return x.parent
end

function Base.isfinite(G::PermGroup)
  return true
end

function Base.isfinite(G::PcGroup)
  return true
end

function Base.isfinite(G::Group)
  return GAP.Globals.IsFinite(G.X)
end

function degree(x::PermGroup)
   return x.deg
end

function order(x::Group)
   return GAP.gap_to_julia(GAP.Globals.Size(x.X))
end

function order(x::GroupElem)
   return GAP.gap_to_julia(GAP.Globals.Order(x.X))
end

function order(::Type{T}, x::Union{GroupElem, Group}) where T<:Number
   return T(order(x))
end

Base.:exponent(x::Group) = GAP.Globals.Exponent(x.X)

Base.:length(x::Group) = order(x)

function rand(x::Group)
   s = GAP.Globals.Random(x.X)
   return group_element(x, s)
end

# one of the following should be non-parametric
rand_pseudo(G::Group) = rand(G)


function _maxgroup(x::T, y::T) where T <: Group
   if x == y
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
function _prod(x::T, y::T) where T <: GroupElem
  G = _maxgroup(parent(x), parent(y))
  return group_element(G, x.X*y.X)
end

Base.:*(x::GroupElem, y::GroupElem) = _prod(x, y)

function ==(x::Group, y::Group)
   return x.X == y.X
end

function ==(x::T, y::T) where T <: GroupElem
   return x.X == y.X
end

Base.:one(x::Group) = group_element(x, GAP.Globals.Identity(x.X))
Base.:one(x::GroupElem) = one(parent(x))
one!(x::GroupElem) = one(parent(x))

Base.:show(io::IO, x::GroupElem) =  print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.X)))
Base.:show(io::IO, x::Group) = print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.X)))

Base.:isone(x::GroupElem) = x == one(parent(x))

Base.:inv(x::GroupElem) = group_element(parent(x), GAP.Globals.Inverse(x.X))

inv!(out::GroupElem, x::GroupElem) = inv(x)  #if needed later

Base.:^(x::GroupElem, y::Int) = group_element(parent(x), x.X ^ y)

Base.:<(x::PermGroupElem, y::PermGroupElem) = x.X < y.X

Base.:/(x::GroupElem, y::GroupElem) = x*y^-1

mul(x::GroupElem, y::GroupElem) = x*y
mul!(out::GroupElem, x::GroupElem, y::GroupElem) = x*y

div_right(x::GroupElem, y::GroupElem) = x*inv(y)
div_left(x::GroupElem, y::GroupElem) = inv(y)*x
div_right!(out::GroupElem, x::GroupElem, y::GroupElem) = x*inv(y)
div_left!(out::GroupElem, x::GroupElem, y::GroupElem) = inv(y)*x

conj(x::GroupElem, y::GroupElem) = x^y
conj!(out::GroupElem, x::GroupElem, y::GroupElem) = x^y

comm(x::GroupElem, y::GroupElem) = x^-1*x^y
comm!(out::GroupElem, x::GroupElem, y::GroupElem) = x^-1*x^y

function iterate(G::Group)
  L=GAP.Globals.Iterator(G.X)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  return group_element(G, i), L
end

function iterate(G::Group, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return group_element(G, i), state
end

function perm(L::Array{Int64,1})
   return PermGroupElem(symmetric_group(length(L)), GAP.Globals.PermList(GAP.julia_to_gap(L)))
end

function perm(g::PermGroup, L::Array{Int64,1})
   x = GAP.Globals.PermList(GAP.julia_to_gap(L))
   if GAP.Globals.IN(x,g.X) 
     return PermGroupElem(g, GAP.Globals.PermList(GAP.julia_to_gap(L)))
   end
   throw(ArgumentError("the element does not embed in the group"))
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
function cperm(L::Union{Array{Int64,1},UnitRange{Int64}}...)
   if length(L)==0
      return one(symmetric_group(1))
   else
      return prod([PermGroupElem(symmetric_group(maximum(y)), GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y)))) for y in L])
   end
end

# cperm stays for "cycle permutation", but we can change name if we want
# takes as input a list of arrays (not necessarly disjoint)
# WARNING: we allow e.g. PermList([2,3,1,4,5,6]) in Sym(3)
function cperm(g::PermGroup,L::Union{Array{Int64,1},UnitRange{Int64}}...)
   if length(L)==0
      return one(g)
   else
      x=GAP.Globals.Product(GAP.julia_to_gap([GAP.Globals.CycleFromList(GAP.julia_to_gap(collect(y))) for y in L]))
      if GAP.Globals.IN(x,g.X) return PermGroupElem(g, x)
      else throw(ArgumentError("the element does not embed in the group"))
      end
   end
end

function listperm(x::PermGroupElem)
   return [x(i) for i in 1:x.parent.deg]
end

function gens(G::Group)
   L = GAP.Globals.GeneratorsOfGroup(G.X)
   res = Vector{elem_type(G)}(undef, length(L))
   for i = 1:length(res)
     res[i] = group_element(G, L[i])
   end
   return res
end

function gens(G::Group, i::Int)
   L = GAP.Globals.GeneratorsOfGroup(G.X)
   @assert length(L) >= i "The number of generators is lower than the given index"
   return group_element(G, L[i])
end

ngens(G::Group) = length(GAP.Globals.GeneratorsOfGroup(G.X))


Base.getindex(G::Group, i::Int) = gens(G, i)
Base.:sign(x::PermGroupElem) = GAP.Globals.SignPerm(x.X)

Base.:isless(x::PermGroupElem, y::PermGroupElem) = x<y

#embedding of a permutation in permutation group
function (G::PermGroup)(x::PermGroupElem)
   if !GAP.Globals.IN(x.X,G.X)
      throw(ArgumentError("the element does not embed in the group"))
   end
   return group_element(G, x.X)
end

#evaluation function
function (x::PermGroupElem)(n)
   return GAP.Globals.OnPoints(n,x.X)
end

function conjugate_subgroup(G::T, x::GroupElem) where T<:Group
  return T(GAP.Globals.ConjugateSubgroup(G.X,x.X))
end

################################################################################
#
#   Conjugacy Classes
#
################################################################################

struct GroupConjClass{T<:Group, S<:Union{GroupElem,Group}}
   X::T
   repr::S
   CC::GapObj
end

Base.:show(io::IO, x::GroupConjClass) = print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.repr))*" ^ "*GAP.gap_to_julia(GAP.Globals.StringView(x.X)))

function _conjugacy_class(G, g, cc::GapObj)         # function for assignment
  return GroupConjClass{typeof(G), typeof(g)}(G, g, cc)
end

==(a::GroupConjClass{T, S}, b::GroupConjClass{T, S}) where S where T = a.CC == b.CC 

order(C::GroupConjClass) = GAP.Globals.Size(C.CC)

representative(C::GroupConjClass) = C.repr

number_conjugacy_classes(G::Group) = GAP.Globals.NrConjugacyClasses(G.X)

# START elements conjugation
function conjugacy_class(G::Group, g::GroupElem)
   return _conjugacy_class(G, g, GAP.Globals.ConjugacyClass(G.X,g.X))
end

function rand(C::GroupConjClass{S,T}) where S where T<:GroupElem
   return group_element(C.X, GAP.Globals.Random(C.CC))
end

function elements(C::GroupConjClass{S, T}) where S where T<:GroupElem
   L=GAP.gap_to_julia(GAP.Globals.AsList(C.CC))
   return T[group_element(C.X,x) for x in L]
end

function conjugacy_classes(G::Group)
   L=GAP.gap_to_julia(GAP.Globals.ConjugacyClasses(G.X))
   return GroupConjClass{typeof(G), elem_type(G)}[ _conjugacy_class(G,group_element(G,GAP.Globals.Representative(cc)),cc) for cc in L]
end

Base.:^(x::GroupElem, y::GroupElem) = group_element(_maxgroup(parent(x), parent(y)), x.X ^ y.X)

function isconjugate(G::Group, x::GroupElem, y::GroupElem)
   if GAP.Globals.IsConjugate(G.X, x.X, y.X)
      return true, group_element(G,GAP.Globals.RepresentativeAction(G.X, x.X, y.X))
   else
      return false, nothing
   end
end
# END elements conjugation 

# START subgroups conjugation
function conjugacy_class(G::T, g::T) where T<:Group
   return _conjugacy_class{typeof(G), typeof(g)}(G, g, GAP.Globals.ConjugacyClassSubgroups(G.X,g.X))
end

function rand(C::GroupConjClass{S,T}) where S where T<:Group
   return T(GAP.Globals.Random(C.CC))
end

function elements(C::GroupConjClass{S, T}) where S where T<:Group
   L=GAP.gap_to_julia(GAP.Globals.AsList(C.CC))
   return T[T(x) for x in L]
end

function conjugacy_classes_subgroups(G::Group)
   L=GAP.gap_to_julia(GAP.Globals.ConjugacyClassesSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,typeof(G)(GAP.Globals.Representative(cc)),cc) for cc in L]
end

function conjugacy_classes_maximal_subgroups(G::Group)
  lS = GAP.Globals.ConjugacyClassesMaximalSubgroups(G.X)
  res = Vector{GroupConjClass{typeof(G), elem_type(G)}}(undef, length(lS))
  for i = 1:length(res)
    res[i] = _conjugacy_class(G, typeof(G)(GAP.Globals.Representative(cc)), ls[i])
  end
  return res
end

Base.:^(H::Group, y::GroupElem) = typeof(H)(H.X ^ y.X)

function is_conjugate(G::Group, H::Group, K::Group)
   if GAP.Globals.IsConjugate(G.X, H.X, K.X)
      return true, group_element(G,GAP.Globals.RepresentativeAction(G.X, H.X, K.X))
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

normalizer(G::T, H::T) where T<:Group = _as_subgroup(GAP.Globals.Normalizer(G.X,H.X),G)

normalizer(G::Group, x::GroupElem) = _as_subgroup(GAP.Globals.Normalizer(G.X,x.X),G)

normaliser = normalizer

core(G::T, H::T) where T<:Group = _as_subgroup(GAP.Globals.Core(G.X,H.X),G)

# normal_closure(G::T, H::T) where T<:Group = T(GAP.Globals.NormalClosure(G.X,H.X))
# we don't know how G,H embeds into N_C(G,H) 

function pcore(G::Group, p::Int64)
   if !GAP.Globals.IsPrime(p)
      throw(ArgumentError("p is not a prime"))
   end
   return _as_subgroup(GAP.Globals.PCore(G.X,p),G)
end



################################################################################
#
# Specific Subgroups
#
################################################################################

# commutator_subgroup(G::T, H::T) where T<:Group = T(GAP.Globals.CommutatorSubgroup(G.X,H.X))
# we don't know how G,H embed into [G,H]

fitting_subgroup(G::Group) = _as_subgroup(GAP.Globals.FittingSubgroup(G.X),G)

frattini_subgroup(G::Group) = _as_subgroup(GAP.Globals.FrattiniSubgroup(G.X),G)

radical_subgroup(G::Group) = _as_subgroup(GAP.Globals.RadicalGroup(G.X),G)

socle(G::Group) = _as_subgroup(GAP.Globals.Socle(G.X),G)


################################################################################
#
# Sylow & Hall Subgroups
#
################################################################################


function sylow_subgroup(G::Group, p::Int64)
   if !GAP.Globals.IsPrime(p)
      throw(ArgumentError("p is not a prime"))
   end
   return _as_subgroup(GAP.Globals.SylowSubgroup(G.X),G)
end

function hall_subgroup(G::Group, P::Array{Int64})
   noprime=false
   for p in P
      if !GAP.Globals.IsPrime(p)
         noprime=true
         break
      end
   end
   if length(P) != length(Set(P)) noprime=true end        # avoid repetitions in P
   if noprime throw(ArgumentError("The integers must be prime")) end
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   return _as_subgroup(GAP.Globals.HallSubgroup(G.X,GAP.julia_to_gap(P)),G)
end

function sylow_system(G::Group)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   L=GAP.gap_to_julia(GAP.Globals.SylowSystem(G.X))
   return typeof(G)[_as_subgroup(G,x) for x in L]
end

function complement_system(G::Group)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   L=GAP.gap_to_julia(GAP.Globals.ComplementSystem(G.X))
   return typeof(G)[_as_subgroup(G,x) for x in L]
end

function hall_system(G::Group)
   if !issolvable(G) throw(ArgumentError("The group is not solvable")) end
   L=GAP.gap_to_julia(GAP.Globals.HallSystem(G.X))
   return typeof(G)[_as_subgroup(G,x) for x in L]
end



################################################################################
#
# Some Properties
#
################################################################################

"""
    isperfect(G)

Test whether `G` is a perfect group, i.e., equal to its derived subgroup.
"""
isperfect(G::Group) = GAP.Globals.IsPerfectGroup(G.X)

issimple(G::Group) = GAP.Globals.IsSimpleGroup(G.X)

isalmostsimple(G::Group) = GAP.Globals.IsAlmostSimpleGroup(G.X)

function ispgroup(G::Group)
   if GAP.Globals.IsPGroup(G.X)
      return true, PrimePGroup(G.X)
   else
      return false, nothing
   end
end

function relators(G::FPGroup)
   L=GAP.gap_to_julia(GAP.Globals.RelatorsOfFpGroup(G.X))
   F=free_group(G)
   return [group_element(F,x) for x in L]
end

