################################################################################
#
#  Managing direct products
#
################################################################################

export
    direct_product,
    embedding,
    factorofdirectproduct,
    homomorphism_of_semidirect_product,
    isfull_direct_product,
    isfull_semidirect_product,
    isfull_wreath_product,
    projection,
    semidirect_product,
    sub,
    wreath_product,
    write_as_full


"""
    direct_product(G::S, H::T)
Return the direct product of `G` and `H`, of type ``DirectProductOfGroups{S,T}``.
"""
function direct_product(G::S, H::T) where S<:GAPGroup where T<:GAPGroup
   X = GAP.Globals.DirectProduct(G.X,H.X)
   return DirectProductOfGroups(X,G,H,X,true)
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
   if !G.isfull
      throw(ArgumentError("Embedding is not defined for proper subgroups of direct products"))
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
function projection(G::DirectProductOfGroups{S,T}, n::Base.Integer) where {S,T}
  # @assert W.isfull "Projection not defined for proper subgroups of wreath products"
   f=GAP.Globals.Projection(G.Xfull,n)
   if n==1 H=G.G1 else H=G.G2 end
   if !G.isfull
      Gf = DirectProductOfGroups(G.Xfull,G.G1,G.G2,G.Xfull,true)
      g = embedding(Gf,G)
      p = g*GAPGroupHomomorphism{typeof(Gf),typeof(H)}(Gf,H,f)
   else
      p = GAPGroupHomomorphism{typeof(G),typeof(H)}(G,H,f)
   end
   return hom(G,H,x->p(x))
end

function (G::DirectProductOfGroups{S,T})(x::GAPGroupElem{S}, y::GAPGroupElem{T}) where S where T
   return embedding(G,1)(x)*embedding(G,2)(y)
end

# start part on subgroups
function _as_subgroup_bare(G::DirectProductOfGroups{S,T}, H::GapObj) where { S , T }
  t = H==G.X
  return DirectProductOfGroups(H, G.G1, G.G2, G.X, t)
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
   if x.isfull
      print(io, "DirectProduct( ", GAP.gap_to_julia(GAP.Globals.StringView(x.G1.X)), " , ", GAP.gap_to_julia(GAP.Globals.StringView(x.G2.X))," )")
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X)))
   end
end

# if a subgroup of a direct product of groups is also a direct product of groups
"""
    write_as_full(G::DirectProductOfGroups)
If `G` is a subgroup of the direct products `G1 x G2`, return `H1 x H2` with `H1` (resp. `H2`) subgroup of `G1` (resp. `G2`) such that `G` = `H1 x H2`. If such `H1`, `H2` do not exist, an ERROR is returned.
"""
function write_as_full(G::DirectProductOfGroups)
   if G.isfull
      return G
   else
      K = direct_product(G.G1,G.G2)
      K1 = image(projection(K,1),G)[1]
      K2 = image(projection(K,2),G)[1]
      H = direct_product(K1,K2)
      if index(H,G)==1
         return H
      else
         throw(ArgumentError("G is not a direct product of groups"))
      end
   end
end

"""
    isfull_direct_product(G::DirectProductOfGroups)
If `G` is a subgroup of the direct products `G1 x G2`, return whether `G` can be written as `H1 x H2` with `H1` (resp. `H2`) subgroup of `G1` (resp. `G2`).
"""
isfull_direct_product(G::DirectProductOfGroups) = G.isfull

Base.:^(H::DirectProductOfGroups, y::GAPGroupElem) = sub([h^y for h in gens(H)])[1]


################################################################################
#
#  Semidirect products
#  
################################################################################

"""
    semidirect_product(G::S, f::GAPGroupHomomorphism, H::T)
Return the semidirect product of `G` and `H`, of type ``SemidirectProductOfGroups{S,T}``, where `f` is a homomorphism from `G` to `Aut`(`H`).
"""
function semidirect_product(G::S, f::GAPGroupHomomorphism{S,AutomorphismGroup{T}}, N::T) where S <: GAPGroup where T <: GAPGroup
   sdp=GAP.Globals.SemidirectProduct(G.X,f.map,N.X)
   return SemidirectProductOfGroups(sdp,G,N,f,sdp,true)
end

# return the element (a,b) in G
function (G::SemidirectProductOfGroups{S,T})(a::GAPGroupElem{S},b::GAPGroupElem{T}) where {S,T}
# simply put parent(L[d+1])==W.H does not work. Example: if I want to write explicitly a permutation in H proper subgroup of Sym(n).
   Gfull = SemidirectProductOfGroups(G.Xfull,G.G1,G.G2,G.f,G.Xfull,true)
   g = embedding(Gfull,1)(a) * embedding(Gfull,2)(b)
   if !G.isfull
      @assert g in G "Element not in the group"
   end
   return g
end

"""
    isfull_semidirect_product(G::SemidirectProductOfGroups)
Return whether `G` is a semidirect product of two groups, instead of a proper subgroup.
"""
isfull_semidirect_product(G::SemidirectProductOfGroups) = G.isfull

"""
    embedding(G::SemidirectProductOfGroups{S,T}, n::Integer)
Return the embedding of the `n`-th component of `G` into `G`, for `n` = 1,2. It is not defined for proper subgroups of semidirect products.
"""
function embedding(G::SemidirectProductOfGroups{S,T}, n::Base.Integer) where S where T
   @assert G.isfull "Embedding not defined for proper subgroups of semidirect products"
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
Return the projection of `G` into the second component of `G`.
"""
function projection(G::SemidirectProductOfGroups{S,T}) where S where T
   f=GAP.Globals.Projection(G.Xfull)
   if !G.isfull
      Gf = SemidirectProductOfGroups(G.Xfull,G.G1,G.G2,G.f,G.Xfull,true)
      g = embedding(Gf,G)
      return g*GAPGroupHomomorphism{SemidirectProductOfGroups{S,T},S}(Gf,G.G1,f)
   else
      return GAPGroupHomomorphism{SemidirectProductOfGroups{S,T},S}(G,G.G1,f)
   end
end


# start part on subgroups
function _as_subgroup_bare(G::SemidirectProductOfGroups{S,T}, H::GapObj) where { S , T }
  return SemidirectProductOfGroups(H, G.G1, G.G2, G.f, G.X, false)
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
   if x.isfull
      print(io, "SemidirectProduct( ", GAP.gap_to_julia(GAP.Globals.StringViewObj(x.G1.X)), " , ", GAP.gap_to_julia(GAP.Globals.StringView(x.G2.X))," )")
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X)))
   end
end

"""
    homomorphism_of_semidirect_product(G::SemidirectProductOfGroups{S,T})
If `G` is semidirect product of `C` and `N`, where `C` acts on `N`, then return the homomorphism `f` from `C` to `Aut`(`N`).
"""
homomorphism_of_semidirect_product(G::SemidirectProductOfGroups{S,T}) where {S,T} = G.f



################################################################################
#
#  Wreath products
#  
################################################################################

"""
    wreath_product(G::T, H::PermGroup) where T<: Group
Return the wreath product of the group `G` and the permutation group `H`. The number of copies of `G` corresponds to the degree of `H` as permutation group (NOT to the number of moved points).

If `W` is a wreath product of `G` and `H`, {`g_1`, ..., `g_n`} are elements of `G` and `h` in `H`, the element `(g_1, ..., h)` of `W` can be obtained by typing
```
    W(g_1,...,g_n, h).
```
"""
function wreath_product(G::T, H::PermGroup) where T<: GAPGroup
   if Set(GAP.gap_to_julia(GAP.Globals.MovedPoints(H.X)))==Set(1:H.deg)
      Wgap=GAP.Globals.WreathProduct(G.X,H.X)
      return WreathProduct(Wgap,G,H,Wgap,true)
   else
      S=symmetric_group(H.deg)
      Wgap=GAP.Globals.WreathProduct(G.X,S.X)
      W1=GAP.Globals.PreImage(GAP.Globals.Projection(Wgap),H.X)
      return WreathProduct(W1,G,H,Wgap,true)
   end
end

# return the element (L[1],L[2],...) in W
function (W::WreathProduct)(L::Union{GAPGroupElem{T},GAPGroupElem{PermGroup}}...) where T <: GAPGroup
   d = W.H.deg
   @assert length(L)==d+1 "Wrong number of arguments"
   for i in 1:d @assert L[i] in W.G "Wrong input" end
   @assert L[d+1] in W.H "Wrong input"
# simply put parent(L[d+1])==W.H does not work. Example: if I want to write explicitly a permutation in H proper subgroup of Sym(n).
   Wfull = WreathProduct(W.Xfull,W.G,W.H,W.Xfull,true)
   g = prod([embedding(Wfull,i)(L[i]) for i in 1:d+1])
   if !W.isfull
      @assert g in W "Element not in the group"
   end
   return g
end

"""
    isfull_wreath_product(G::WreathProduct)
Return whether `G` is a wreath product of two groups, instead of a proper subgroup.
"""
isfull_wreath_product(G::WreathProduct) = G.isfull

"""
    projection(G::WreathProduct)
Return the projection of `WreathProduct(G,H)` into the permutation group `H`.
"""
function projection(W::WreathProduct)
  # @assert W.isfull "Projection not defined for proper subgroups of wreath products"
   f=GAP.Globals.Projection(W.Xfull)
   if !W.isfull
      Wf = WreathProduct(W.Xfull,W.G,W.H,W.Xfull,true)
      g = embedding(Wf,W)
      p = g*GAPGroupHomomorphism{WreathProduct,PermGroup}(Wf,W.H,f)
   else
      p = GAPGroupHomomorphism{WreathProduct,PermGroup}(W,W.H,f)
   end
   return hom(W,W.H,x->p(x))
end


"""
    embedding(G::WreathProduct, n::Integer)
Return the embedding of the `n`-th component of `G` into `G`.
"""
function embedding(W::WreathProduct, n::Base.Integer)
   if !W.isfull throw(ArgumentError("Embedding not defined for proper subgroups of wreath products")) end
   if n > W.H.deg+1
      throw(ArgumentError("n is too big"))
   else
      f = GAP.Globals.Embedding(W.Xfull,n)
      if n== W.H.deg+1 C=W.H
      else C=W.G
      end
      return GAPGroupHomomorphism{typeof(C),WreathProduct}(C,W,f)
   end
end

Base.show(io::IO, x::WreathProduct) = print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.X)))


# start part on subgroups
#TODO : to be fixed
function _as_subgroup_bare(W::WreathProduct, X::GapObj)
   if X==W.X
      return W
   else
      return WreathProduct(X, W.G, W.H, W.Xfull, false)
   end
end

function _as_subgroup(W::WreathProduct, H::GapObj, ::Type{U}) where U
  H1 = _as_subgroup_bare(W, H)
  return H1, hom(H1, W, x::U -> group_element(W, x.X))
end

function _as_subgroup(W::WreathProduct, H::GapObj)
  return _as_subgroup(W, H, GAPGroupElem{WreathProduct})
end

function sub(W::WreathProduct, elms::Vector{GAPGroupElem{WreathProduct}})
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elms])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(W, H)
end

function sub(L::GAPGroupElem{WreathProduct}...)
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function sub(L::Vector{GAPGroupElem{WreathProduct}})
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

