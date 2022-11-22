export
    general_linear_group,
    mat_elem_type,
    matrix_group,
    MatrixGroup,
    MatrixGroupElem,
    omega_group,
    orthogonal_group,
    ring_elem_type,
    special_linear_group,
    special_orthogonal_group,
    special_unitary_group,
    symplectic_group,
    unitary_group,
    GL, GO, GU, SL, SO, Sp, SU



abstract type AbstractMatrixGroupElem <: GAPGroupElem{GAPGroup} end

# NOTE: always defined are deg, ring and at least one between { X, gens, descr }
"""
    MatrixGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup

Type of groups `G` of `n x n` matrices over the ring `R`, where `n = degree(G)` and `R = base_ring(G)`.

At the moment, only rings of type `FqNmodFiniteField` are supported.
"""
@attributes mutable struct MatrixGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup
   deg::Int
   ring::Ring
   X::GapObj
   gens::Vector{<:AbstractMatrixGroupElem}
   descr::Symbol                       # e.g. GL, SL, symbols for isometry groups
   ring_iso::MapFromFunc # Isomorphism from the Oscar base ring to the GAP base ring

   MatrixGroup{RE,T}(m::Int, F::Ring) where {RE,T} = new{RE,T}(m,F)

end

MatrixGroup(m::Int, F::Ring) = MatrixGroup{elem_type(F), dense_matrix_type(elem_type(F))}(m,F)

# build a MatrixGroup given a list of generators, given as array of either MatrixGroupElem or AbstractAlgebra matrices
# WARNING: if the elements of V have type MatElem, it does not check whether the determinant is nonzero
function MatrixGroup{RE,S}(m::Int, F::Ring, V::AbstractVector{T}) where {RE,S} where T<:Union{MatElem,AbstractMatrixGroupElem}
   G = MatrixGroup(m,F)

   L = Vector{elem_type(G)}(undef, length(V))
   for i in 1:length(V)
      @assert base_ring(V[i])==F "The elements must have the same base ring"
      @assert nrows(V[i])==m "The elements must have the same dimension"
      if T<:MatElem
         L[i] = MatrixGroupElem(G,V[i])
      else
# TODO: this part of code from here
         if isdefined(V[i],:elm)
            if isdefined(V[i],:X)
               L[i] = MatrixGroupElem(G,V[i].elm,V[i].X)
            else
               L[i] = MatrixGroupElem(G,V[i].elm)
            end
         else
            L[i] = MatrixGroupElem(G,V[i].X)
         end
# to here

# can be replaced by the following line once deepcopy works for GAP objects
# > L[i] = deepcopy(V[i]); L[i].parent = G;
      end
   end
   G.gens = L

   return G
end

MatrixGroup(m::Int, F::Ring, V::AbstractVector{T}) where T<:Union{MatElem,AbstractMatrixGroupElem} = MatrixGroup{elem_type(F), dense_matrix_type(elem_type(F))}(m,F,V)


# `MatrixGroup`: compare types, dimensions, and coefficient rings
function check_parent(G::T, g::GAPGroupElem) where T <: MatrixGroup
  P = g.parent
  return T === typeof(P) && G.deg == P.deg && G.ring == P.ring
end

# `MatrixGroup`: set dimension and ring of `G`
function _oscar_group(obj::GapObj, G::MatrixGroup)
  d = GAP.Globals.DimensionOfMatrixGroup(obj)
  d == G.deg || error("requested dimension of matrices ($(G.deg)) does not match the given matrix dimension ($d)")

  R = G.ring
  iso = G.ring_iso
  GAPWrap.IsSubset(codomain(iso), GAP.Globals.FieldOfMatrixGroup(obj)) || error("matrix entries are not in the requested ring ($(codomain(iso)))")

  M = MatrixGroup(d, R)
  M.X = obj
  M.ring = R
  M.ring_iso = iso
  return M
end


function _as_subgroup_bare(G::T, H::GapObj) where {T <: MatrixGroup}
  H1 = T(G.deg,G.ring)
  H1.ring_iso = G.ring_iso
  H1.X = H
  return H1
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

MatrixGroupElem(G::MatrixGroup{RE,T}, x::T, x_gap::GapObj) where {RE,T} = MatrixGroupElem{RE,T}(G, x, x_gap)

MatrixGroupElem(G::MatrixGroup{RE,T}, x::T) where {RE, T} = MatrixGroupElem{RE,T}(G,x)
MatrixGroupElem(G::MatrixGroup{RE,T}, x_gap::GapObj) where {RE, T} = MatrixGroupElem{RE,T}(G,x_gap)

ring_elem_type(::Type{MatrixGroup{S,T}}) where {S,T} = S
mat_elem_type(::Type{MatrixGroup{S,T}}) where {S,T} = T
_gap_filter(::Type{<:MatrixGroup}) = GAP.Globals.IsMatrixGroup

elem_type(::Type{MatrixGroup{S,T}}) where {S,T} = MatrixGroupElem{S,T}
elem_type(::MatrixGroup{S,T}) where {S,T} = MatrixGroupElem{S,T}
Base.eltype(::Type{MatrixGroup{S,T}}) where {S,T} = MatrixGroupElem{S,T}

# `parent_type` is defined and documented in AbstractAlgebra.
parent_type(::Type{T}) where T<:MatrixGroupElem{RE,S} where {RE,S} = MatrixGroup{RE,S}
parent_type(::T) where T<:MatrixGroupElem{RE,S} where {RE,S} = MatrixGroup{RE,S}


function Base.deepcopy_internal(x::MatrixGroupElem, dict::IdDict)
  if isdefined(x, :X)
    X = Base.deepcopy_internal(x.X, dict)
    if isdefined(x, :elm)
      elm = Base.deepcopy_internal(x.elm, dict)
      return MatrixGroupElem(x.parent, elm, X)
    else
      return MatrixGroupElem(x.parent, X)
    end
  elseif isdefined(x, :elm)
    elm = Base.deepcopy_internal(x.elm, dict)
    return MatrixGroupElem(x.parent, elm)
  end
  error("$x has neither :X nor :elm")
end


########################################################################
#
# Basic
#
########################################################################

function Base.show(io::IO, x::MatrixGroup)
   AbstractAlgebra.@show_name(io, x)
   AbstractAlgebra.@show_special(io,x)

   if isdefined(x, :descr)
      if x.descr==:GU || x.descr==:SU
         print(io, string(x.descr), "(",x.deg,",",characteristic(x.ring)^(div(degree(x.ring),2)),")")
      else
         print(io, string(x.descr), "(",x.deg,",",order(x.ring),")")
      end
   else
      print(io, "Matrix group of degree ", x.deg, " over ")
      show(IOContext(io, :compact => true), x.ring)
   end
end

Base.show(io::IO, x::MatrixGroupElem) = show(io, x.elm)
Base.show(io::IO, mi::MIME"text/plain", x::MatrixGroupElem) = show(io, mi, x.elm)

group_element(G::MatrixGroup, x::GapObj) = MatrixGroupElem(G,x)

function assign_from_description(G::MatrixGroup)
   F = codomain(G.ring_iso)
   if G.descr==:GL G.X=GAP.Globals.GL(G.deg, F)
   elseif G.descr==:SL G.X=GAP.Globals.SL(G.deg, F)
   elseif G.descr==:Sp G.X=GAP.Globals.Sp(G.deg, F)
   elseif G.descr==Symbol("GO+") G.X=GAP.Globals.GO(1, G.deg, F)
   elseif G.descr==Symbol("SO+") G.X=GAP.Globals.SO(1, G.deg, F)
   elseif G.descr==Symbol("Omega+")
      # FIXME/TODO: Work around GAP issue <https://github.com/gap-system/gap/issues/500>
      # using the following inefficient code. In the future, we should use appropriate
      # generators for Omega (e.g. by applying a form change matrix to the Omega
      # generators returned by GAP).
      # This also compensates for the fact that Omega(-1,2,q) is not supported in GAP.
      L = GAP.Globals.SubgroupsOfIndexTwo(GAP.Globals.SO(1, G.deg, F))
      if G.deg==4 && order(G.ring)==2  # this is the only case SO(n,q) has more than one subgroup of index 2
         for y in L
            _ranks = [GAP.Globals.Rank(u) for u in GAP.Globals.GeneratorsOfGroup(y)]
            if all(r->iseven(r),_ranks)
               G.X=y
               break
            end
         end
      else
         @assert length(L) == 1
         G.X=L[1]
      end
   elseif G.descr==Symbol("GO-") G.X=GAP.Globals.GO(-1, G.deg, F)
   elseif G.descr==Symbol("SO-") G.X=GAP.Globals.SO(-1, G.deg, F)
   elseif G.descr==Symbol("Omega-") G.X=GAP.Globals.SubgroupsOfIndexTwo(GAP.Globals.SO(-1, G.deg, F))[1]
   elseif G.descr==:GO G.X=GAP.Globals.GO(0, G.deg, F)
   elseif G.descr==:SO G.X=GAP.Globals.SO(0, G.deg, F)
   elseif G.descr==:Omega
     # For even q or d = 1, \Omega(d,q) is equal to SO(d,q).
     # Otherwise, \Omega(d,q) has index 2 in SO(d,q).
     # Here d is odd, and we do not get here if d == 1 holds
     # because `omega_group` delegates to `SO` in this case.
     @assert G.deg > 1
     if iseven(GAPWrap.Size(F))
       G.X = GAP.Globals.SO(0, G.deg, F)
     else
       L = GAP.Globals.SubgroupsOfIndexTwo(GAP.Globals.SO(0, G.deg, F))
       @assert length(L) == 1
       G.X = L[1]
     end
   elseif G.descr==:GU G.X=GAP.Globals.GU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) ))
   elseif G.descr==:SU G.X=GAP.Globals.SU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) ))
   else error("unsupported description")
   end
end

# return the G.sym if isdefined(G, :sym); otherwise, the field :sym is computed and set using information from other defined fields
function Base.getproperty(G::MatrixGroup, sym::Symbol)

   isdefined(G,sym) && return getfield(G,sym)

   if sym === :ring_iso
      G.ring_iso = iso_oscar_gap(G.ring)

   elseif sym === :X
      if isdefined(G,:descr)
         assign_from_description(G)
      elseif isdefined(G,:gens)
         V = GAP.GapObj([g.X for g in gens(G)])
         G.X = isempty(V) ? GAP.Globals.Group(V, one(G).X) : GAP.Globals.Group(V)
      else
         error("Cannot determine underlying GAP object")
      end
   end

   return getfield(G, sym)

end


function Base.getproperty(x::MatrixGroupElem, sym::Symbol)

   isdefined(x,sym) && return getfield(x,sym)

   if sym === :X
      x.X = map_entries(x.parent.ring_iso, x.elm)
   elseif sym == :elm
      x.elm = preimage_matrix(x.parent.ring_iso, x.X)
   end
   return getfield(x,sym)
end

Base.IteratorSize(::Type{<:MatrixGroup}) = Base.SizeUnknown()

function Base.iterate(G::MatrixGroup)
  L=GAP.Globals.Iterator(G.X)::GapObj
  @assert ! GAPWrap.IsDoneIterator(L)
  i = GAPWrap.NextIterator(L)
  return MatrixGroupElem(G, i), L
end

function Base.iterate(G::MatrixGroup, state::GapObj)
  if GAPWrap.IsDoneIterator(state)
    return nothing
  end
  i = GAPWrap.NextIterator(state)
  return MatrixGroupElem(G, i), state
end

########################################################################
#
# Membership
#
########################################################################


function ==(G::MatrixGroup,H::MatrixGroup)
   G.deg==H.deg || return false
   G.ring==H.ring || return false
   if isdefined(G, :descr) && isdefined(H, :descr)
      return G.descr == H.descr
   end
   if isdefined(G, :gens) && isdefined(H, :gens)
      gens(G)==gens(H) && return true
   end
   return G.X==H.X
end


# this saves the value of x.X
# x_gap = x.X if this is already known, x_gap = nothing otherwise
function lies_in(x::MatElem, G::MatrixGroup, x_gap)
   if base_ring(x)!=G.ring || nrows(x)!=G.deg return false, x_gap end
   if isone(x) return true, x_gap end
   if isdefined(G,:gens)
      for g in gens(G)
         if x==g.elm
            return true, x_gap
         end
      end
   end
   if isdefined(G,:descr) && G.descr==:GL
      return det(x)!=0, x_gap
   elseif isdefined(G,:descr) && G.descr==:SL
      return det(x)==1, x_gap
   else
      if x_gap==nothing x_gap = map_entries(G.ring_iso, x) end
     # x_gap !=nothing || x_gap = map_entries(G.ring_iso, x)
      return (x_gap in G.X), x_gap
   end
end

Base.in(x::MatElem, G::MatrixGroup) = lies_in(x,G,nothing)[1]

function Base.in(x::MatrixGroupElem, G::MatrixGroup)
   isdefined(x,:X) && return lies_in(x.elm,G,x.X)[1]
   _is_true, x_gap = lies_in(x.elm,G,nothing)
   if x_gap !=nothing x.X = x_gap end
   return _is_true
end

# embedding an element of type MatElem into a group G
# if check=false, there are no checks on the condition `x in G`
function (G::MatrixGroup)(x::MatElem; check::Bool=true)
   if check
      _is_true, x_gap = lies_in(x,G,nothing)
      _is_true || throw(ArgumentError("Element not in the group"))
      x_gap != nothing && return MatrixGroupElem(G,x,x_gap)
   end
   return MatrixGroupElem(G,x)
end

# embedding an element of type MatrixGroupElem into a group G
# if check=false, there are no checks on the condition `x in G`
function (G::MatrixGroup)(x::MatrixGroupElem; check::Bool=true)
   if !check
      z = x
      z.parent = G
      return z
   end
   if isdefined(x,:X)
      if isdefined(x,:elm)
         _is_true = lies_in(x.elm,G,x.X)[1]
         _is_true || throw(ArgumentError("Element not in the group"))
         return MatrixGroupElem(G,x.elm,x.X)
      else
         x.X in G.X || throw(ArgumentError("Element not in the group"))
         return MatrixGroupElem(G,x.X)
      end
   else
      _is_true, x_gap = lies_in(x.elm,G,nothing)
      _is_true || throw(ArgumentError("Element not in the group"))
      if x_gap==nothing return MatrixGroupElem(G,x.elm)
      else return MatrixGroupElem(G,x.elm,x_gap)
      end
   end
end

# embedding a n x n array into a group G
function (G::MatrixGroup)(L::AbstractVecOrMat; check::Bool=true)
   x = matrix(G.ring, G.deg, G.deg, L)
   return G(x; check=check)
end


########################################################################
#
# Methods on elements
#
########################################################################

# we are not currently keeping track of the parent; two elements coincide iff their matrices coincide
function ==(x::MatrixGroupElem{S,T},y::MatrixGroupElem{S,T}) where {S,T}
   if isdefined(x,:X) && isdefined(y,:X) return x.X==y.X
   else return x.elm==y.elm
   end
end

function _common_parent_group(x::T, y::T) where T <: MatrixGroup
   @assert x.deg == y.deg
   @assert x.ring === y.ring
   x === y && return x
   return GL(x.deg, x.ring)::T
end

# Base.:* is defined in src/Groups/GAPGroups.jl, and it calls the function _prod below
# if the parents are different, the parent of the output product is set as GL(n,q)
function _prod(x::T,y::T) where {T <: MatrixGroupElem}
   G = _common_parent_group(parent(x), parent(y))

   # if the underlying GAP matrices are both defined, but not both Oscar matrices,
   # then use the GAP matrices.
   # Otherwise, use the Oscar matrices, which if necessary are implicitly computed
   # by the Base.getproperty(::MatrixGroupElem, ::Symbol) method .
   if isdefined(x,:X) && isdefined(y,:X) && !(isdefined(x,:elm) && isdefined(y,:elm))
      return T(G, x.X*y.X)
   else
      return T(G, x.elm*y.elm)
   end
end

Base.:*(x::MatrixGroupElem{RE, T}, y::T) where RE where T = x.elm*y
Base.:*(x::T, y::MatrixGroupElem{RE, T}) where RE where T = x*y.elm

Base.:^(x::MatrixGroupElem, n::Int) = MatrixGroupElem(x.parent, x.elm^n)

Base.isone(x::MatrixGroupElem) = isone(x.elm)

Base.inv(x::MatrixGroupElem) = MatrixGroupElem(x.parent, inv(x.elm))

# if the parents are different, the parent of the output is set as GL(n,q)
function Base.:^(x::MatrixGroupElem, y::MatrixGroupElem)
   G = x.parent==y.parent ? x.parent : GL(x.parent.deg, x.parent.ring)
   if isdefined(x,:X) && isdefined(y,:X) && !(isdefined(x,:elm) && isdefined(y,:elm))
      return MatrixGroupElem(G, inv(y.X)*x.X*y.X)
   else
      return MatrixGroupElem(G,inv(y.elm)*x.elm*y.elm)
   end
end

comm(x::MatrixGroupElem, y::MatrixGroupElem) = inv(x)*conj(x,y)

"""
    det(x::MatrixGroupElem)
    
Return the determinant of the underlying matrix of `x`.
"""
det(x::MatrixGroupElem) = det(matrix(x))

"""
    base_ring(x::MatrixGroupElem)

Return the base ring of the underlying matrix of `x`.
"""
base_ring(x::MatrixGroupElem) = x.parent.ring

parent(x::MatrixGroupElem) = x.parent

"""
    matrix(x::MatrixGroupElem)

Return the underlying matrix of `x`.
"""
matrix(x::MatrixGroupElem) = x.elm

Base.getindex(x::MatrixGroupElem, i::Int, j::Int) = x.elm[i,j]

"""
    nrows(x::MatrixGroupElem)

Return the number of rows of the underlying matrix of `x`.
"""
nrows(x::MatrixGroupElem) = nrows(matrix(x))

"""
    ncols(x::MatrixGroupElem)

Return the number of columns of the underlying matrix of `x`.
"""
ncols(x::MatrixGroupElem) = ncols(matrix(x))


"""
    tr(x::MatrixGroupElem)

Return the trace of the underlying matrix of `x`.
"""
tr(x::MatrixGroupElem) = tr(matrix(x))

#FIXME for the following functions, the output may not belong to the parent group of x
#=
frobenius(x::MatrixGroupElem, n::Int) = MatrixGroupElem(x.parent, matrix(x.parent.ring, x.parent.deg, x.parent.deg, [frobenius(y,n) for y in x.elm]))
frobenius(x::MatrixGroupElem) = frobenius(x,1)

transpose(x::MatrixGroupElem) = MatrixGroupElem(x.parent, transpose(x.elm))
=#

########################################################################
#
# Methods on groups
#
########################################################################

"""
    base_ring(G::MatrixGroup)

Return the base ring of the matrix group `G`.
"""
base_ring(G::MatrixGroup) = G.ring

"""
    degree(G::MatrixGroup)

Return the degree of the matrix group `G`, i.e. the number of rows of its matrices.
"""
degree(G::MatrixGroup) = G.deg

Base.one(G::MatrixGroup) = MatrixGroupElem(G, identity_matrix(G.ring, G.deg))

function Base.rand(rng::Random.AbstractRNG, G::MatrixGroup)
   x_gap = GAP.Globals.Random(GAP.wrap_rng(rng), G.X)::GapObj
   return MatrixGroupElem(G, x_gap)
end

function gens(G::MatrixGroup)
   if !isdefined(G,:gens)
      L = GAP.Globals.GeneratorsOfGroup(G.X)::GapObj
      G.gens = [MatrixGroupElem(G, a) for a in L]
   end
   return G.gens
end

gen(G::MatrixGroup, i::Int) = gens(G)[i]

ngens(G::MatrixGroup) = length(gens(G))


compute_order(G::GAPGroup) = fmpz(GAPWrap.Order(G.X))

function compute_order(G::MatrixGroup{T}) where {T <: Union{nf_elem, fmpq}}
  #=
    - For a matrix group G over the Rationals or over a number field,
    the GAP group G.X does usually not store the flag `IsHandledByNiceMonomorphism`.
    - If we know a reasonable ("nice") faithful permutation action of `G` in advance,
    we can set this flag in `G.X` to true and store the action homomorphism in `G.X`,
    and then this information should be used in the computation of the order.
    - If the flag is not known to be true then the Oscar code from
    `isomorphic_group_over_finite_field` shall be preferred.
  =#
  if GAP.Globals.HasIsHandledByNiceMonomorphism(G.X) && GAPWrap.IsHandledByNiceMonomorphism(G.X)
    # The call to `IsHandledByNiceMonomorphism` triggers an expensive
    # computation of `IsFinite` which we avoid by checking
    # `HasIsHandledByNiceMonomorphism` first.
    return fmpz(GAPWrap.Order(G.X))
  else
    order(isomorphic_group_over_finite_field(G)[1])
  end
end

function order(::Type{T}, G::MatrixGroup) where T <: IntegerUnion
   res = get_attribute!(G, :order) do
     return compute_order(G)
   end::fmpz
   return T(res)::T
end

########################################################################
#
# Constructors
#
########################################################################

function _field_from_q(q::Int)
   (n,p) = is_power(q)
   is_prime(p) || throw(ArgumentError("The field size must be a prime power"))
   return GF(p, n)
end

"""
    general_linear_group(n::Int, q::Int)
    general_linear_group(n::Int, F::FqNmodFiniteField)
    GL = general_linear_group

Return the general linear group of dimension `n` either over the field `F` or the field `GF(q)`.
"""
function general_linear_group(n::Int, F::Ring)
   G = MatrixGroup(n,F)
   G.descr = :GL
   return G
end

function general_linear_group(n::Int, q::Int)
   return general_linear_group(n, _field_from_q(q))
end

"""
    special_linear_group(n::Int, q::Int)
    special_linear_group(n::Int, F::FqNmodFiniteField)
    SL = special_linear_group

Return the special linear group of dimension `n` either over the field `F` or the field `GF(q)`.
"""
function special_linear_group(n::Int, F::Ring)
   G = MatrixGroup(n,F)
   G.descr = :SL
   return G
end

function special_linear_group(n::Int, q::Int)
   return special_linear_group(n, _field_from_q(q))
end

"""
    symplectic_group(n::Int, q::Int)
    symplectic_group(n::Int, F::FqNmodFiniteField)
    Sp = symplectic_group

Return the symplectic group of dimension `n` either over the field `F` or the
field `GF(q)`. The dimension `n` must be even.
"""
function symplectic_group(n::Int, F::Ring)
   iseven(n) || throw(ArgumentError("The dimension must be even"))
   G = MatrixGroup(n,F)
   G.descr = :Sp
   return G
end

function symplectic_group(n::Int, q::Int)
   return symplectic_group(n, _field_from_q(q))
end

"""
    orthogonal_group(e::Int, n::Int, F::Ring)
    orthogonal_group(e::Int, n::Int, q::Int)
    GO = orthogonal_group

Return the orthogonal group of dimension `n` either over the field `F` or the
field `GF(q)` of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0`
for `n` odd. If `n` is odd, `e` can be omitted.
"""
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
   return orthogonal_group(e, n, _field_from_q(q))
end

orthogonal_group(n::Int, F::Ring) = orthogonal_group(0,n,F)
orthogonal_group(n::Int, q::Int) = orthogonal_group(0,n,q)

"""
    special_orthogonal_group(e::Int, n::Int, F::Ring)
    special_orthogonal_group(e::Int, n::Int, q::Int)
    SO = special_orthogonal_group

Return the special orthogonal group of dimension `n` either over the field `F`
or the field `GF(q)` of type `e`, where `e` in {`+1`,`-1`} for `n` even and
`e`=`0` for `n` odd. If `n` is odd, `e` can be omitted.
"""
function special_orthogonal_group(e::Int, n::Int, F::Ring)
   iseven(order(F)) && return GO(e,n,F)
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
   return special_orthogonal_group(e, n, _field_from_q(q))
end

special_orthogonal_group(n::Int, F::Ring) = special_orthogonal_group(0,n,F)
special_orthogonal_group(n::Int, q::Int) = special_orthogonal_group(0,n,q)

"""
    omega_group(e::Int, n::Int, F::Ring)
    omega_group(e::Int, n::Int, q::Int)

Return the Omega group of dimension `n` over the field `GF(q)` of type `e`,
where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd,
`e` can be omitted.
"""
function omega_group(e::Int, n::Int, F::Ring)
   n==1 && return SO(e,n,F)
   if e==1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("Omega+")
   elseif e==-1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("Omega-")
   elseif e==0
      isodd(n) || throw(ArgumentError("The dimension must be odd"))
      G = MatrixGroup(n,F)
      G.descr = :Omega
   else
      throw(ArgumentError("Invalid description of orthogonal group"))
   end
   return G
end

function omega_group(e::Int, n::Int, q::Int)
   return omega_group(e, n, _field_from_q(q))
end

omega_group(n::Int, q::Int) = omega_group(0,n,q)
omega_group(n::Int, F::Ring) = omega_group(0,n,F)

"""
    unitary_group(n::Int, q::Int)
    GU = unitary_group

Return the unitary group of dimension `n` over the field `GF(q^2)`.
"""
function unitary_group(n::Int, q::Int)
   (a,b) = is_power(q)
   is_prime(b) || throw(ArgumentError("The field size must be a prime power"))
   G = MatrixGroup(n,GF(b, 2*a))
   G.descr = :GU
   return G
end

"""
    special_unitary_group(n::Int, q::Int)
    SU = special_unitary_group

Return the special unitary group of dimension `n` over the field `GF(q^2)`.
"""
function special_unitary_group(n::Int, q::Int)
   (a,b) = is_power(q)
   is_prime(b) || throw(ArgumentError("The field size must be a prime power"))
   G = MatrixGroup(n,GF(b, 2*a))
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


"""
    matrix_group(V::T...) where T<:MatrixGroup
    matrix_group(V::AbstractVector{T}) where T<:MatrixGroup

Return the matrix group generated by elements in the vector `V`.
"""
function matrix_group(V::AbstractVector{T}) where T<:Union{MatElem,MatrixGroupElem}
   if T<:MatElem      # if T <: MatrixGroupElem, we can already assume that det(V[i]) != 0 
      for v in V @assert det(v) !=0 "The matrix is not invertible" end
   end
   return MatrixGroup(nrows(V[1]),base_ring(V[1]),V)
end

matrix_group(V::T...) where T<:Union{MatElem,MatrixGroupElem} = matrix_group(collect(V))


########################################################################
#
# Subgroups
#
########################################################################

function sub(G::MatrixGroup, elements::Vector{S}) where S <: GAPGroupElem
   @assert elem_type(G) === S
   elems_in_GAP = GAP.Obj(GapObj[x.X for x in elements])
   H = GAP.Globals.Subgroup(G.X,elems_in_GAP)::GapObj
   #H is the group. I need to return the inclusion map too
   K,f = _as_subgroup(G, H)
   L = Vector{elem_type(K)}(undef, length(elements))
   for i in 1:length(L)
      L[i] = MatrixGroupElem(K, elements[i].elm, elements[i].X)
   end
   K.gens = L
   return K,f
end


########################################################################
#
# Conjugation
#
########################################################################

function Base.show(io::IO, x::GroupConjClass{T,S}) where T <: MatrixGroup where S <: MatrixGroupElem
  show(io, x.repr)
  print(io, " ^ ")
  show(io, x.X)
end

function Base.show(io::IO, x::GroupConjClass{T,S}) where T <: MatrixGroup where S <: MatrixGroup
  show(io, x.repr)
  print(io, " ^ ")
  show(io, x.X)
end

function Base.:^(H::MatrixGroup, y::MatrixGroupElem)
   if isdefined(H,:gens) && !isdefined(H,:X)
      K = MatrixGroup(H.deg, H.ring)
      K.gens = [y^-1*x*y for x in H.gens]
      for k in gens(K) k.parent = K end
   else
      K = MatrixGroup(H.deg,H.ring)
      K.X = H.X^y.X
   end

   return K
end

function Base.rand(rng::Random.AbstractRNG, C::GroupConjClass{S,T}) where S<:MatrixGroup where T<:MatrixGroup
   H = MatrixGroup(C.X.deg,C.X.ring)
   H.X = GAP.Globals.Random(GAP.wrap_rng(rng), C.CC)::GapObj
   return H
end

function Base.rand(C::GroupConjClass{S,T}) where S<:MatrixGroup where T<:MatrixGroup
   return Base.rand(Random.GLOBAL_RNG, C)
end
