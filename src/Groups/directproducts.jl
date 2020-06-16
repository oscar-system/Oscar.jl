################################################################################
#
#  Managing direct products
#
################################################################################

export
    directproduct,
    embedding,
    factorofdirectproduct,
    homomorphism_of_semidirect_product,
    isfull_direct_product,
    projection,
    semidirectproduct,
    sub,
    write_as_full


"""
    directproduct(G::S, H::T)
Return the direct product of `G` and `H`, of type ``DirectProductOfGroups{S,T}``.
"""
function directproduct(G::S, H::T) where S<:GAPGroup where T<:GAPGroup
   return DirectProductOfGroups(GAP.Globals.DirectProduct(G.X,H.X), G,H,true)
end

function factorofdirectproduct(G::DirectProductOfGroups, n::Base.Integer)
   if n==1 gr=G.G1
   elseif n==2 gr=G.G2
   else throw(ArgumentError("n must be 1 or 2"))
   end
   return gr
end

"""
    embedding(G::DirectProductOfGroups{S,T}, n::Integer)
Return the embedding of the `n`-th component of `G` into `G`, for `n` = 1,2.
"""
function embedding(G::DirectProductOfGroups{S,T}, n::Base.Integer) where S where T
   if !isfull_direct_product(G)
      throw(ArgumentError("Embedding is not a group homomorphism"))
   else
      G=write_as_full(G)
   end
   if n==1
      f=GAP.Globals.Embedding(G.X,1)
      typeG=S
      gr=G.G1
   elseif n==2
      f=GAP.Globals.Embedding(G.X,2)
      typeG=T
      gr=G.G2
   else
      throw(ArgumentError("n must be 1 or 2"))
   end
   return GAPGroupHomomorphism{typeG,DirectProductOfGroups{S,T}}(gr,G,f)
end

"""
    projection(G::DirectProductOfGroups{S,T}, n::Integer)
Return the projection of `G` into the `n`-th component of `G`, for `n` = 1,2.
"""
function projection(G::DirectProductOfGroups{S,T}, n::Base.Integer) where S where T
   if !isfull_direct_product(G)
      throw(ArgumentError("Projection is not a group homomorphism"))
   else
      G=write_as_full(G)
   end
   if n==1
      f=GAP.Globals.Projection(G.X,1)
      typeG=S
      gr=G.G1
   elseif n==2
      f=GAP.Globals.Projection(G.X,2)
      typeG=T
      gr=G.G2
   else
      throw(ArgumentError("n must be 1 or 2"))
   end
   return GAPGroupHomomorphism{DirectProductOfGroups{S,T},typeG}(G,gr,f)
end

function (G::DirectProductOfGroups{S,T})(x::GAPGroupElem{S}, y::GAPGroupElem{T}) where S where T
   return embedding(G,1)(x)*embedding(G,2)(y)
end

# start part on subgroups
function _as_subgroup_bare(G::DirectProductOfGroups{S,T}, H::GapObj) where { S , T }
  return DirectProductOfGroups(H, G.G1, G.G2, false)
end

function _as_subgroup(G::DirectProductOfGroups{S,T}, H::GapObj, ::Type{U}) where { T, S, U }
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x::U -> group_element(G, x.X))
end

function _as_subgroup(G::DirectProductOfGroups{S,T}, H::GapObj) where S <: GAPGroup where T <: GAPGroup
  return _as_subgroup(G, H, GAPGroupElem{DirectProductOfGroups{S,T}})
end

function sub(G::DirectProductOfGroups{S,T}, elms::Vector{GAPGroupElem{DirectProductOfGroups{S,T}}}) where { S, T }
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elms])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(G, H)
end

function sub(L::GAPGroupElem{DirectProductOfGroups{S,T}}...) where { S, T }
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function sub(L::Vector{GAPGroupElem{DirectProductOfGroups{S,T}}}) where { S, T }
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function Base.show(io::IO, x::DirectProductOfGroups)
   if x.IsFull
      print(io, "DirectProduct( ", GAP.gap_to_julia(GAP.Globals.StringView(x.G1.X)), " , ", GAP.gap_to_julia(GAP.Globals.StringView(x.G2.X))," )")
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.X)))
   end
end

# if a subgroup of a direct product of groups is also a direct product of groups
function write_as_full(G::DirectProductOfGroups)
   if G.IsFull
      return G
   else
      K = directproduct(G.G1,G.G2)
      K1 = image(projection(K,1),G)[1]
      K2 = image(projection(K,2),G)[1]
      H = directproduct(K1,K2)
      if index(H,G)==1
         return H
      else
         throw(ArgumentError("G is not a direct product of groups"))
      end
   end
end

function isfull_direct_product(G::DirectProductOfGroups)
   if G.IsFull
      return true
   else
      K = directproduct(G.G1,G.G2)
      K1 = image(projection(K,1),G)[1]
      K2 = image(projection(K,2),G)[1]
      H = directproduct(K1,K2)
      if index(H,G)==1
         return true
      else
         false
      end
   end
end

Base.:^(H::DirectProductOfGroups, y::GAPGroupElem) = sub([h^y for h in gens(H)])[1]


################################################################################
#
#  Semidirect products
#  
################################################################################

"""
    semidirectproduct(G::S, f::GAPGroupHomomorphism, H::T)
Return the semidirect product of `G` and `H`, of type ``SemidirectProductOfGroups{S,T}``, where `f` is a homomorphism from `G` to `Aut`(`H`).
"""
function semidirectproduct(G::S, f::GAPGroupHomomorphism{S,AutomorphismGroup{T}}, N::T) where S <: GAPGroup where T <: GAPGroup
   return SemidirectProductOfGroups(GAP.Globals.SemidirectProduct(G.X,f.map,N.X),G,N,f,true)
end

"""
    embedding(G::SemidirectProductOfGroups{S,T}, n::Integer)
Return the embedding of the `n`-th component of `G` into `G`, for `n` = 1,2.
"""
function embedding(G::SemidirectProductOfGroups{S,T}, n::Base.Integer) where S where T
   if n==1
      f=GAP.Globals.Embedding(G.X,1)
      typeG=S
      gr=G.G1
   elseif n==2
      f=GAP.Globals.Embedding(G.X,2)
      typeG=T
      gr=G.G2
   else
      throw(ArgumentError("n must be 1 or 2"))
   end
   return GAPGroupHomomorphism{typeG,SemidirectProductOfGroups{S,T}}(gr,G,f)
end



"""
    projection(G::SemidirectProductOfGroups{S,T}, n::Integer)
Return the projection of `G` into the `n`-th component of `G`, for `n` = 1,2.
"""
function projection(G::SemidirectProductOfGroups{S,T}) where S where T
 #=  if !isfull_direct_product(G)
      throw(ArgumentError("Projection is not a group homomorphism"))
   else
      G=write_as_full(G)
   end
=#
   f=GAP.Globals.Projection(G.X)
   typeG=S
   gr=G.G1
   return GAPGroupHomomorphism{SemidirectProductOfGroups{S,T},typeG}(G,gr,f)
end


# start part on subgroups
function _as_subgroup_bare(G::SemidirectProductOfGroups{S,T}, H::GapObj) where { S , T }
  return SemidirectProductOfGroups(H, G.G1, G.G2, G.f, false)
end

function _as_subgroup(G::SemidirectProductOfGroups{S,T}, H::GapObj, ::Type{U}) where { T, S, U }
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x::U -> group_element(G, x.X))
end

function _as_subgroup(G::SemidirectProductOfGroups{S,T}, H::GapObj) where S <: GAPGroup where T <: GAPGroup
  return _as_subgroup(G, H, GAPGroupElem{SemidirectProductOfGroups{S,T}})
end

function sub(G::SemidirectProductOfGroups{S,T}, elms::Vector{GAPGroupElem{SemidirectProductOfGroups{S,T}}}) where { S, T }
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elms])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(G, H)
end

function sub(L::GAPGroupElem{SemidirectProductOfGroups{S,T}}...) where { S, T }
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function sub(L::Vector{GAPGroupElem{SemidirectProductOfGroups{S,T}}}) where { S, T }
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function Base.show(io::IO, x::SemidirectProductOfGroups)
   if x.IsFull
      print(io, "SemidirectProduct( ", GAP.gap_to_julia(GAP.Globals.StringView(x.G1.X)), " , ", GAP.gap_to_julia(GAP.Globals.StringView(x.G2.X))," )")
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.X)))
   end
end

"""
    homomorphism_of_semidirect_product(G::SemidirectProductOfGroups{S,T})
If `G` is semidirect product of `C` and `N`, where `C` acts on `N`, then return the homomorphism `f` from `C` to `Aut`(`N`).
"""
homomorphism_of_semidirect_product(G::SemidirectProductOfGroups{S,T}) where {S,T} = G.f
