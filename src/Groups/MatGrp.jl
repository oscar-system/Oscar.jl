import AbstractAlgebra: MatElem, MatSpace, parent_type, Ring, RingElem
import Hecke: base_ring, det, fq_nmod, FqNmodFiniteField, tr, trace
import GAP: FFE

export
    general_linear_group,
    MatrixGroup,
    MatrixGroupElem,
    omega_group,
    orthogonal_group,
    special_linear_group,
    special_orthogonal_group,
    special_unitary_group,
    symplectic_group,
    unitary_group,
    GL, GO, GU, SL, SO, Sp, SU



mutable struct MatrixGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup
   X::GapObj
   gens::Vector{<:T}
   ring::Ring
   deg::Int
   descr::Symbol                       # e.g. GL
   ring_iso::GenRingIso
   mat_iso::GenMatIso

   function MatrixGroup(m::Int64, F::S) where S<:Ring
      K = elem_type(S)
      r = new{K, MatElem{K}}()
      r.deg = m
      r.ring = F
      return r
   end

   function MatrixGroup(m::Int64, F::S, G_gap::GapObj) where S<:Ring
      K = elem_type(S)
      r = new{K, MatElem{K}}()
      r.deg = m
      r.ring = F
      r.X = G_gap
      return r
   end
end

mutable struct MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}} <: GAPGroupElem{MatrixGroup}
   parent::MatrixGroup{RE, T}
   elm::T
   X::GapObj

   function MatrixGroupElem(G::MatrixGroup{RE,MatElem{RE}}, x::MatElem{RE}) where RE <: RingElem
      z = new{RE, MatElem{RE}}()
      z.parent = G
      z.elm = x
      return z
   end

   function MatrixGroupElem(G::MatrixGroup{RE,MatElem{RE}}, x::MatElem{RE}, x_gap::GapObj) where RE <: RingElem
      z = new{RE, MatElem{RE}}()
      z.parent = G
      z.elm = x
      z.X = x_gap
      return z
   end

   function MatrixGroupElem(G::MatrixGroup{RE,MatElem{RE}}, x_gap::GapObj) where RE <: RingElem
      z = new{RE, MatElem{RE}}()
      z.parent = G
      z.X = x_gap
      return z
   end
end

ring_elem_type(::Type{MatrixGroup{S,T}}) where {S,T} = S
mat_elem_type(::Type{MatrixGroup{S,T}}) where {S,T} = T

########################################################################
#
# Basic
#
########################################################################

function Base.show(io::IO, x::MatrixGroup)
   if isdefined(x, :descr)
      if x.descr==:GU || x.descr==:SU
         print(io, string(x.descr), "(",x.deg,",",characteristic(x.ring)^(div(degree(x.ring),2)),")")
      else
         print(io, string(x.descr), "(",x.deg,",",order(x.ring),")")
      end
   else
      print(io, "Matrix group of degree ", x.deg, " over ", x.ring)
   end
end

function Base.show(io::IO, x::MatrixGroupElem)          #TODO: is this correct?
   if isdefined(x, :elm)
      display(x.elm)
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X)))
   end
end

group_element(G::MatrixGroup, x::GapObj) = MatrixGroupElem(G,x)

function Base.getproperty(G::MatrixGroup, sym::Symbol)

   if sym === :mat_iso || sym === :ring_iso
      if isdefined(G,sym) return getfield(G,sym) end
      fm = gen_mat_iso(G.deg, G.ring)
      setfield!(G,:ring_iso,fm.fr)
      setfield!(G,:mat_iso,fm)

   elseif sym === :X
      if isdefined(G,:X) return getfield(G,:X) end
      if isdefined(G,:descr)
         if !isdefined(G,:ring_iso)
            fm = gen_mat_iso(G.deg, G.ring)
            setfield!(G,:ring_iso,fm.fr)
            setfield!(G,:mat_iso,fm)
         end
         if G.descr==:GL setfield!(G,:X,GAP.Globals.GL(G.deg,G.ring_iso.codomain))
         elseif G.descr==:SL setfield!(G,:X,GAP.Globals.SL(G.deg,G.ring_iso.codomain))
         elseif G.descr==:Sp setfield!(G,:X,GAP.Globals.Sp(G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("GO+") setfield!(G,:X,GAP.Globals.GO(1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("SO+") setfield!(G,:X,GAP.Globals.SO(1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("Omega+") setfield!(G,:X,GAP.Globals.Omega(1,G.deg,Int(order(G.ring))))
         elseif G.descr==Symbol("GO-") setfield!(G,:X,GAP.Globals.GO(-1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("SO-") setfield!(G,:X,GAP.Globals.SO(-1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("Omega-") setfield!(G,:X,GAP.Globals.Omega(-1,G.deg,Int(order(G.ring))))
         elseif G.descr==:GO setfield!(G,:X,GAP.Globals.GO(0,G.deg,G.ring_iso.codomain))
         elseif G.descr==:SO setfield!(G,:X,GAP.Globals.SO(0,G.deg,G.ring_iso.codomain))
         elseif G.descr==:Omega setfield!(G,:X,GAP.Globals.Omega(0,G.deg,Int(order(G.ring))))
         elseif G.descr==:GU setfield!(G,:X,GAP.Globals.GU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) )))
         elseif G.descr==:SU setfield!(G,:X,GAP.Globals.SU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) )))
         end
      elseif isdefined(G,:gens)
         if !isdefined(G,:ring_iso)
            fm = gen_mat_iso(G.deg, G.ring)
            setfield!(G,:ring_iso,fm.fr)
            setfield!(G,:mat_iso,fm)
         end
         V = GAP.julia_to_gap([G.mat_iso(g) for g in G.gens])
         setfield!(G,:X,GAP.Globals.Group(V))
      end

   elseif sym === :gens
      if isdefined(G, :gens) return getfield(G, :gens) end
      if !isdefined(G, :X)
         if !isdefined(G,:ring_iso)
            fm = gen_mat_iso(G.deg, G.ring)
            setfield!(G,:ring_iso,fm.fr)
            setfield!(G,:mat_iso,fm)
         end
         if G.descr==:GL setfield!(G,:X,GAP.Globals.GL(G.deg,G.ring_iso.codomain))
         elseif G.descr==:SL setfield!(G,:X,GAP.Globals.SL(G.deg,G.ring_iso.codomain))
         elseif G.descr==:Sp setfield!(G,:X,GAP.Globals.Sp(G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("GO+") setfield!(G,:X,GAP.Globals.GO(1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("SO+") setfield!(G,:X,GAP.Globals.SO(1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("Omega+") setfield!(G,:X,GAP.Globals.Omega(1,G.deg,Int(order(G.ring))))
         elseif G.descr==Symbol("GO-") setfield!(G,:X,GAP.Globals.GO(-1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("SO-") setfield!(G,:X,GAP.Globals.SO(-1,G.deg,G.ring_iso.codomain))
         elseif G.descr==Symbol("Omega-") setfield!(G,:X,GAP.Globals.Omega(-1,G.deg,Int(order(G.ring))))
         elseif G.descr==:GO setfield!(G,:X,GAP.Globals.GO(0,G.deg,G.ring_iso.codomain))
         elseif G.descr==:SO setfield!(G,:X,GAP.Globals.SO(0,G.deg,G.ring_iso.codomain))
         elseif G.descr==:Omega setfield!(G,:X,GAP.Globals.Omega(0,G.deg,Int(order(G.ring))))
         elseif G.descr==:GU setfield!(G,:X,GAP.Globals.GU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) )))
         elseif G.descr==:SU setfield!(G,:X,GAP.Globals.SU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) )))
         end
      end
      L = GAP.Globals.GeneratorsOfGroup(getfield(G,:X))
      setfield!(G,:gens,[G.mat_iso(L[i]) for i in 1:length(L)])

   else
      return getfield(G, sym)

   end

end

# add field :X
function Base.getproperty(x::MatrixGroupElem, sym::Symbol)
   if sym === :X && !isdefined(x,:X)
      if !isdefined(x.parent,:ring_iso)
         fm = gen_mat_iso(x.parent.deg, x.parent.ring)
         setfield!(x.parent,:ring_iso,fm.fr)
         setfield!(x.parent,:mat_iso,fm)
      end
      setfield!(x, :X, x.parent.mat_iso(x.elm))
   elseif sym == :elm && !isdefined(x, :elm)                        # assuming that x.X, hence x.parent.ring_iso, are defined
      setfield!(x, :elm, x.mat_iso(x.X))
   end
   return getfield(x,sym)
end

Base.IteratorSize(::Type{<:MatrixGroup}) = Base.SizeUnknown()

function Base.iterate(G::MatrixGroup)
  L=GAP.Globals.Iterator(G.X)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  return MatrixGroupElem(G, G.mat_iso(i),i), L
end

function Base.iterate(G::MatrixGroup, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return MatrixGroupElem(G, G.mat_iso(i),i), state
end

########################################################################
#
# Membership
#
########################################################################

# this saves the value of x.X
# x_gap = x.X if this is already known, x_gap = 0 otherwise
function lies_in(x::MatElem, G::MatrixGroup, x_gap)
   if base_ring(x)!=G.ring || nrows(x)!=G.deg return false, Nothing end
   if isone(x) return true, Nothing end
   if isdefined(G,:gens) && x in G.gens
      return true, Nothing
   else
      if x_gap==0 x_gap = G.mat_iso(x) end
     # x_gap !=0 || x_gap = G.mat_iso(x)
      return GAP.Globals.in(x_gap,G.X), x_gap
   end
end

Base.in(x::MatElem, G::MatrixGroup) = lies_in(x,G,0)[1]

function Base.in(x::MatrixGroupElem, G::MatrixGroup)
   if isdefined(x,:X) return lies_in(x.elm,G,x.X)[1]
   else return lies_in(x.elm,G,0)[1]
   end
end

# embedding an element of type MatElem into a group G
function (G::MatrixGroup)(x::MatElem)
   vero, x_gap = lies_in(x,G,0)
   vero || throw(ArgumentError("Element not in the group"))
   if x_gap==Nothing return MatrixGroupElem(G,x)
   else return MatrixGroupElem(G,x,x_gap)
   end
end


########################################################################
#
# Methods on elements
#
########################################################################

# TODO: we are not currently keeping track of the parent
function ==(x::MatrixGroupElem{S,T},y::MatrixGroupElem{S,T}) where {S,T}
   if isdefined(x,:X) && isdefined(y,:X) return x.X==y.X
   else return x.elm==y.elm
   end
end

# TODO: we wish to multiply / conjugate also matrices with different parents ?
function _prod(x::MatrixGroupElem,y::MatrixGroupElem)
   if x.parent==y.parent
      if isdefined(x,:X) && isdefined(y,:X) && !(isdefined(x,:elm) && isdefined(y,:elm))
         return MatrixGroupElem(x.parent, x.X*y.X)
      else
         return MatrixGroupElem(x.parent, x.elm*y.elm)
      end
   else
      throw(ArgumentError("Matrices not in the same group"))
   end
end
# Base.:* is defined in src/Groups/GAPGroups.jl

Base.:^(x::MatrixGroupElem, n::Int) = MatrixGroupElem(x.parent, x.elm^n)

Base.isone(x::MatrixGroupElem) = isone(x.elm)

Base.inv(x::MatrixGroupElem) = MatrixGroupElem(x.parent, inv(x.elm))

# TODO: we wish to multiply / conjugate also matrices with different parents ?
function Base.:^(x::MatrixGroupElem, y::MatrixGroupElem)
   if x.parent==y.parent
      if isdefined(x,:X) && isdefined(y,:X) && !(isdefined(x,:elm) && isdefined(y,:elm))
         return MatrixGroupElem(x.parent, inv(y.X)*x.X*y.X)
      else
         return MatrixGroupElem(x.parent,inv(y.elm)*x.elm*y.elm)
      end
   else
      throw(ArgumentError("Matrices not in the same group"))
   end
end

comm(x::MatrixGroupElem, y::MatrixGroupElem) = inv(x)*conj(x,y)

det(x::MatrixGroupElem) = det(x.elm)
base_ring(x::MatrixGroupElem) = x.parent.ring

parent(x::MatrixGroupElem) = x.parent

Base.getindex(x::MatrixGroupElem, i::Int, j::Int) = x.elm[i,j]

trace(x::MatrixGroupElem) = trace(x.elm)
tr(x::MatrixGroupElem) = tr(x.elm)

order(x::MatrixGroupElem) = GAP.Globals.Order(x.X)

########################################################################
#
# Methods on groups
#
########################################################################

base_ring(G::MatrixGroup) = G.ring

Base.one(G::MatrixGroup) = MatrixGroupElem(G, identity_matrix(G.ring, G.deg))

function Base.rand(G::MatrixGroup)
   x_gap = GAP.Globals.Random(G.X)
   x_oscar = G.mat_iso(x_gap)
   x = MatrixGroupElem(G,x_oscar,x_gap)
   return x
end

function gens(G::MatrixGroup)
   V = [MatrixGroupElem(G,g) for g in G.gens]
   return V
end

gens(G::MatrixGroup, i::Int) = gens(G)[i]

Base.getindex(G::MatrixGroup, i::Int) = gens(G, i)

function ngens(G::MatrixGroup)
   return length(G.gens)
end

function order(G::MatrixGroup)
   return GAP.Globals.Order(G.X)
end

########################################################################
#
# Constructors
#
########################################################################

function general_linear_group(n::Int, F::Ring)
   G = MatrixGroup(n,F)
   G.descr = :GL
   return G
end

function general_linear_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return general_linear_group(n, GF(b,a)[1])
end

function special_linear_group(n::Int, F::Ring)
   G = MatrixGroup(n,F)
   G.descr = :SL
   return G
end

function special_linear_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return special_linear_group(n, GF(b,a)[1])
end

function symplectic_group(n::Int, F::Ring)
   iseven(n) || throw(ArgumentError("The dimension must be even"))
   G = MatrixGroup(n,F)
   G.descr = :Sp
   return G
end

function symplectic_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return symplectic_group(n, GF(b,a)[1])
end

function orthogonal_group(e::Int, n::Int, F::Ring)
   if e==1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("GO+")
   elseif e==-1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("GO-")
   elseif e==0
      isodd(n) || throw(ArgumentError("The dimension must be odd"))
      G = MatrixGroup(n,F)
      G.descr = :GO
   else
      throw(ArgumentError("Invalid description of orthogonal group"))
   end
   return G
end

function orthogonal_group(e::Int, n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return orthogonal_group(e, n, GF(b,a)[1])
end

orthogonal_group(n::Int, F::Ring) = orthogonal_group(0,n,F)
orthogonal_group(n::Int, q::Int) = orthogonal_group(0,n,q)

function special_orthogonal_group(e::Int, n::Int, F::Ring)
   if e==1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("SO+")
   elseif e==-1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("SO-")
   elseif e==0
      isodd(n) || throw(ArgumentError("The dimension must be odd"))
      G = MatrixGroup(n,F)
      G.descr = :SO
   else
      throw(ArgumentError("Invalid description of orthogonal group"))
   end
   return G
end

function special_orthogonal_group(e::Int, n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return special_orthogonal_group(e, n, GF(b,a)[1])
end

special_orthogonal_group(n::Int, F::Ring) = special_orthogonal_group(0,n,F)
special_orthogonal_group(n::Int, q::Int) = special_orthogonal_group(0,n,q)

function omega_group(e::Int, n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   if e==1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,GF(b,a)[1])
      G.descr = Symbol("Omega+")
   elseif e==-1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,GF(b,a)[1])
      G.descr = Symbol("Omega-")
   elseif e==0
      isodd(n) || throw(ArgumentError("The dimension must be odd"))
      G = MatrixGroup(n,GF(b,a)[1])
      G.descr = :Omega
   else
      throw(ArgumentError("Invalid description of orthogonal group"))
   end
   return G
end

omega_group(n::Int, q::Int) = omega_group(0,n,q)

function unitary_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   G = MatrixGroup(n,GF(b,2*a)[1])
   G.descr = :GU
   return G
end

function special_unitary_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   G = MatrixGroup(n,GF(b,2*a)[1])
   G.descr = :SU
   return G
end

const GL = general_linear_group
const SL = special_linear_group
const Sp = symplectic_group
const GO = orthogonal_group
const SO = special_orthogonal_group
const GU = unitary_group
const SU = special_unitary_group


########################################################################
#
# Subgroup
#
########################################################################

function _as_subgroup(G::T, H::GapObj) where T <: MatrixGroup
  return MatrixGroup(G.deg,G.ring,H)
end

# convert a GAP list of subgroups into a vector of Julia groups objects
function _as_subgroups(G::T, subs::GapObj) where T <: MatrixGroup
  res = Vector{T}(undef, length(subs))
  for i = 1:length(res)
    res[i] = _as_subgroup(G, subs[i])
  end
  return res
end

function intersect(V::T...) where T<:MatrixGroup
   L = GAP.julia_to_gap([G.X for G in V])
   K = GAP.Globals.Intersection(L)
   return _as_subgroup(V[1], K)
end

function intersect(V::AbstractVector{T}) where T<:MatrixGroup
   L = GAP.julia_to_gap([G.X for G in V])
   K = GAP.Globals.Intersection(L)
   return _as_subgroup(V[1], K)
end


########################################################################
#
# Conjugation
#
########################################################################

function Base.show(io::IO, x::GroupConjClass{T,S}) where T <: MatrixGroup where S <: MatrixGroupElem
  show(x.repr)
  print(" ^ ")
  show(x.X)
end

function Base.show(io::IO, x::GroupConjClass{T,S}) where T <: MatrixGroup where S <: MatrixGroup
  show(x.repr)
  print(" ^ ")
  show(x.X)
end

function Base.:^(H::MatrixGroup, y::MatrixGroupElem)
   if isdefined(H,:gens)
      K = MatrixGroup(H.deg, H.ring)
      K.gens = [inv(y.elm)*x*y.elm for x in H.gens]
      return K
   else
      return MatrixGroup(H.deg,H.ring,H.X^y.X)
   end
end

function conjugacy_classes_subgroups(G::MatrixGroup)
   L=GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClassesSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,MatrixGroup(G.deg,G.ring,GAP.Globals.Representative(cc)),cc) for cc in L]
end

function conjugacy_classes_maximal_subgroups(G::MatrixGroup)
  L = GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClassesMaximalSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,MatrixGroup(G.deg,G.ring,GAP.Globals.Representative(cc)),cc) for cc in L]
end

function Base.rand(C::GroupConjClass{S,T}) where S<:MatrixGroup where T<:MatrixGroup
   return MatrixGroup(C.X.deg,C.X.ring,GAP.Globals.Random(C.CC))
end

function elements(C::GroupConjClass{S, T}) where S where T<:MatrixGroup
   L=GAP.Globals.AsList(C.CC)
   l = Vector{T}(undef, length(L))
   for i in 1:length(l)
      l[i] = _as_subgroup(C.X, L[i])
   end
   return l
end





#=

########################################################################
#
# conversions     (this should be no longer necessary, and it's also wrong)
#
########################################################################
#TODO: support other field types

function FieldGapToOscar(F::GapObj)
   p = GAP.Globals.Characteristic(F)
   q = GAP.Globals.Size(F)
   d = Base.Integer(log(p,q))

   return GF(p,d)
end

function FieldOscarToGap(F::FqNmodFiniteField)
   p = Int64(characteristic(F))
   d = Int64(degree(F))

   return GAP.Globals.GF(p,d)
end

function FieldElemGapToOscar(x::FFE, Fgap::GapObj, F::FqNmodFiniteField)
   q = GAP.Globals.Size(Fgap)
   z = gen(F)
   if GAP.Globals.IsZero(x)
      return 0*z
   else
      d = Integer(GAP.Globals.LogFFE(x, GAP.Globals.Z(q)))
      return z^d
   end
end

function FieldElemGapToOscar(x::FFE, F::GapObj)
   q = GAP.Globals.Size(F)
   K,z = FieldGapToOscar(F)
   if GAP.Globals.IsZero(x)
      return 0*z
   else
      d = Integer(GAP.Globals.LogFFE(x, GAP.Globals.Z(q)))
      return z^d
   end
end

function FieldElemOscarToGap(x::fq_nmod, F::FqNmodFiniteField)
   q = Int64(order(F))
   d = degree(F)
   z = gen(F)
   v = [Int64(coeff(x,i)) for i in 0:d]
   y = sum([GAP.Globals.Z(q)^i*v[i+1] for i in 0:d])
   
   return y
end

function MatGapToOscar(x::GapObj, Fgap::GapObj, F::FqNmodFiniteField)
   n = GAP.Globals.Size(x)
   Arr = [GAP.gap_to_julia(x[i]) for i in 1:n]
   L = [FieldElemGapToOscar(Arr[i][j],Fgap,F) for i in 1:n for j in 1:n]

   return matrix(F,n,n,L)
end

function MatGapToOscar(x::GapObj, F::GapObj)
   n = GAP.Globals.Size(x)
   Arr = [GAP.gap_to_julia(x[i]) for i in 1:n]
   L = [FieldElemGapToOscar(Arr[i][j],F) for i in 1:n for j in 1:n]

   return matrix(FieldGapToOscar(F)[1],n,n,L)
end

function MatOscarToGap(x::MatElem, F::FqNmodFiniteField)
   r = nrows(x)
   c = ncols(x)
   S = Vector{GapObj}(undef,r)
   for i in 1:r
      S[i] = GAP.julia_to_gap([FieldElemOscarToGap(x[i,j],F) for j in 1:c])
   end

   return GAP.julia_to_gap(S)
end

=#
