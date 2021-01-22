################################################################################
#
#  Managing direct products
#
################################################################################

export
    acting_subgroup,
    as_perm_group,
    as_polycyclic_group,
    cartesian_power,
    direct_product,
    embedding,
    factor_of_direct_product,
    homomorphism_of_semidirect_product,
    homomorphism_of_wreath_product,
    inner_cartesian_power,
    inner_direct_product,
    isfull_direct_product,
    isfull_semidirect_product,
    isfull_wreath_product,
    normal_subgroup,
    number_of_factors,
    projection,
    semidirect_product,
    sub,
    wreath_product,
    write_as_full


"""
    direct_product(L::AbstractVector{<:GAPGroup})
    direct_product(L::GAPGroup...)

Return the direct product of the groups in the collection `L`.
"""
function direct_product(L::AbstractVector{<:GAPGroup})
   X = GAP.Globals.DirectProduct(GAP.julia_to_gap([G.X for G in L]))
   return DirectProductGroup(X,L,X,true)
end

function direct_product(L::GAPGroup...)
   return direct_product([x for x in L])
end

"""
    inner_direct_product(L::AbstractVector{T}; morphisms)
    inner_direct_product(L::T...)

Return a direct product of groups of the same type `T` as a group of type `T`. It works for `T` of the following types:
- `PermGroup`, `PcGroup`, `FPGroup`.

The parameter `morphisms` is `false` by default. If it is set `true`, then the output is a triple (`G`, `emb`, `proj`), where `emb` and `proj` are the vectors of the embeddings (resp. projections) of the direct product `G`.
"""
function inner_direct_product(L::AbstractVector{T}; morphisms=false) where T<:Union{PcGroup,PermGroup,FPGroup}
   P = GAP.Globals.DirectProduct(GAP.julia_to_gap([G.X for G in L]))
   if T==PermGroup
      DP = T(P,GAP.Globals.NrMovedPoints(P))
   else
      DP = T(P)
   end
   if morphisms
      emb = [_hom_from_gap_map(L[i],DP,GAP.Globals.Embedding(P,i)) for i in 1:length(L)]
      proj = [_hom_from_gap_map(DP,L[i],GAP.Globals.Projection(P,i)) for i in 1:length(L)]
      return DP, emb, proj
   else
      return DP
   end
end

function inner_direct_product(L::T... ; morphisms=false) where T<:Union{PcGroup,PermGroup,FPGroup}
   return inner_direct_product([x for x in L]; morphisms=morphisms)
end

"""
    number_of_factors(G::DirectProductGroup)

Return the number of factors of `G`.
"""
number_of_factors(G::DirectProductGroup) = length(G.L)

"""
    cartesian_power(G::T, n::Int)

Return the direct product of `n` copies of `G`.
"""
function cartesian_power(G::GAPGroup, n::Base.Integer)
   L = [G for i in 1:n]
   return direct_product(L)
end

"""
    inner_cartesian_power(G::T, n::Int; morphisms)

Return the direct product of `n` copies of `G` as group of type `T`.

The parameter `morphisms` is `false` by default. If it is set `true`, then the output is a triple (`G`, `emb`, `proj`), where `emb` and `proj` are the vectors of the embeddings (resp. projections) of the direct product `G`.
"""
function inner_cartesian_power(G::T, n::Base.Integer; morphisms=false) where T <: GAPGroup
   L = [G for i in 1:n]
   return inner_direct_product(L; morphisms=morphisms)
end

"""
    factor_of_direct_product(G::DirectProductGroup, j::Int)

Return the `j`-th factor of `G`.
"""
function factor_of_direct_product(G::DirectProductGroup, j::Base.Integer)
   if j in 1:length(G.L) return G.L[j]
   else throw(ArgumentError("index not valid"))
   end
end

"""
    as_perm_group(G::DirectProductGroup)

If `G` is direct product of permutations groups, return `G` as permutation group.
"""
function as_perm_group(G::DirectProductGroup)
   if [typeof(H)==PermGroup for H in G.L]==[true for i in 1:length(G.L)]
      return PermGroup(G.X, GAP.Globals.Maximum(GAP.Globals.MovedPoints(G.X)))
   else
      throw(ArgumentError("The group is not a permutation group"))
   end
end

"""
    as_polycyclic_group(G::DirectProductGroup)

If `G` is direct product of polycyclic groups, return `G` as polycyclic group.
"""
function as_polycyclic_group(G::DirectProductGroup)
   if [typeof(H)==PcGroup for H in G.L]==[true for i in 1:length(G.L)]
      return PcGroup(G.X)
   else
      throw(ArgumentError("The group is not a polycyclic group"))
   end
end

"""
    embedding(G::DirectProductGroup, j::Integer)

Return the embedding of the `j`-th component of `G` into `G`, for `j` = 1,...,#factors of `G`.
"""
function embedding(G::DirectProductGroup, j::Base.Integer)
   j in 1:length(G.L) || throw(ArgumentError("index not valid"))
   G.isfull || throw(ArgumentError("Embedding is not defined for proper subgroups of direct products"))
   f=GAP.Globals.Embedding(G.X,j)
   gr=G.L[j]
   typeG=typeof(gr)
   return GAPGroupHomomorphism{typeG,DirectProductGroup}(gr,G,f)
end

"""
    projection(G::DirectProductGroup, j::Integer)

Return the projection of `G` into the `j`-th component of `G`, for `j` = 1,...,#factors of `G`.
"""
function projection(G::DirectProductGroup, j::Base.Integer)
   f=GAP.Globals.Projection(G.Xfull,j)
   j in 1:length(G.L) || throw(ArgumentError("index not valid"))
   H = G.L[j]
   p = GAP.Globals.GroupHomomorphismByFunction(G.X,H.X,y->GAP.Globals.Image(f,y))
   return _hom_from_gap_map(G,H,p)
end


function (G::DirectProductGroup)(V::AbstractVector{<:GAPGroupElem})
   length(V)==length(G.L) || throw(ArgumentError("Wrong number of entries"))
   arr = [GAP.Globals.Image(GAP.Globals.Embedding(G.Xfull,i),V[i].X) for i in 1:length(V)]
   xgap = prod(arr)
   xgap in G.X || throw(ArgumentError("Element not in the group"))
   return group_element(G,xgap)
end

function (G::DirectProductGroup)(V::GAPGroupElem...)
   return G([x for x in V])
end

# start part on subgroups
function _as_subgroup_bare(G::DirectProductGroup, H::GapObj)
#  t = H==G.X
  return DirectProductGroup(H, G.L, G.X, false)
end

function _as_subgroup(G::DirectProductGroup, H::GapObj, ::Type{U}) where U
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x::U -> group_element(G, x.X))
end

function _as_subgroup(G::DirectProductGroup, H::GapObj)
  return _as_subgroup(G, H, GAPGroupElem{DirectProductGroup})
end

function sub(G::DirectProductGroup, elms::Vector{<:GAPGroupElem{DirectProductGroup}})
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elms])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(G, H)
end

function sub(L::T...) where T <: GAPGroupElem{DirectProductGroup}
   length(L)>0 || throw(ArgumentError("Empty list"))
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function sub(L::Vector{<:GAPGroupElem{DirectProductGroup}})
   length(L)>0 || throw(ArgumentError("Empty list"))
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function Base.show(io::IO, G::DirectProductGroup)
   if G.isfull
      print(io, "DirectProduct of ")
      display(G.L)
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(G.X)))
   end
end

# if a subgroup of a direct product of groups is also a direct product of groups
"""
    write_as_full(G::DirectProductGroup)

If `G` is a subgroup of the direct product `G_1 x ... x G_n` such that `G = H_1 x ... x H_n` for `H_i` subgroup of `G_i`, return `G` as full direct product of the `H_i`. If such `H_i` do not exist, an ERROR is returned.
"""
function write_as_full(G::DirectProductGroup)
   if G.isfull
      return G
   else
      LK = [image(projection(G,j))[1] for j in 1:length(G.L)]
      H = direct_product(LK)
# index(H,G)==1 does not work because it does not recognize G as a subgroup of H
      order(H)==order(G) || throw(ArgumentError("G is not a direct product of groups"))
      return H
   end
end

"""
    isfull_direct_product(G::DirectProductGroup)

Return whether `G` is direct product of its factors (`false` if it is a proper subgroup).
"""
isfull_direct_product(G::DirectProductGroup) = G.isfull

Base.:^(H::DirectProductGroup, y::GAPGroupElem) = sub([h^y for h in gens(H)])[1]


################################################################################
#
#  Semidirect products
#  
################################################################################

"""
    semidirect_product(N::S, f::GAPGroupHomomorphism, H::T)

Return the semidirect product of `N` and `H`, of type ``SemidirectProductGroup{S,T}``, where `f` is a homomorphism from `H` to `Aut`(`N`).
"""
function semidirect_product(N::S, f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}, H::T) where S <: GAPGroup where T <: GAPGroup
   sdp=GAP.Globals.SemidirectProduct(H.X,f.map,N.X)
   return SemidirectProductGroup(sdp,N,H,f,sdp,true)
end

# return the element (a,b) in G
function (G::SemidirectProductGroup{S,T})(a::GAPGroupElem{S},b::GAPGroupElem{T}) where {S,T}
# simply put parent(L[d+1])==W.H does not work. Example: if I want to write explicitly a permutation in H proper subgroup of Sym(n).
   xgap = GAP.Globals.Image(GAP.Globals.Embedding(G.Xfull,1),b.X)*GAP.Globals.Image(GAP.Globals.Embedding(G.Xfull,2),a.X)
   xgap in G.X || throw(ArgumentError("Element not in the group"))
   return group_element(G,xgap)
end

"""
    normal_subgroup(G::SemidirectProductGroup)

Return `N`, where `G` is the semidirect product of the normal subgroup `N` and `H`.
"""
normal_subgroup(G::SemidirectProductGroup) = G.N

"""
    acting_subgroup(G::SemidirectProductGroup)

Return `H`, where `G` is the semidirect product of the normal subgroup `N` and `H`.
"""
acting_subgroup(G::SemidirectProductGroup) = G.H

"""
    homomorphism_of_semidirect_product(G::SemidirectProductGroup)

Return `f,` where `G` is the semidirect product of the normal subgroup `N` and the group `H` acting on `N` via the homomorphism `h`.
"""
homomorphism_of_semidirect_product(G::SemidirectProductGroup) = G.f

"""
    isfull_semidirect_product(G::SemidirectProductGroup)

Return whether `G` is a semidirect product of two groups, instead of a proper subgroup.
"""
isfull_semidirect_product(G::SemidirectProductGroup) = G.isfull

"""
    embedding(G::SemidirectProductGroup, n::Integer)

Return the embedding of the `n`-th component of `G` into `G`, for `n` = 1,2. It is not defined for proper subgroups of semidirect products.
"""
function embedding(G::SemidirectProductGroup{S,T}, n::Base.Integer) where S where T
   @assert G.isfull "Embedding not defined for proper subgroups of semidirect products"
   if n==1
      f=GAP.Globals.Embedding(G.X,2)
      typeG=S
      gr=G.N
   elseif n==2
      f=GAP.Globals.Embedding(G.X,1)
      typeG=T
      gr=G.H
   else
      throw(ArgumentError("n must be 1 or 2"))
   end
   return GAPGroupHomomorphism{typeG,SemidirectProductGroup{S,T}}(gr,G,f)
end

"""
    projection(G::SemidirectProductGroup, n::Integer)

Return the projection of `G` into the second component of `G`.
"""
function projection(G::SemidirectProductGroup)
   f=GAP.Globals.Projection(G.Xfull)
   H = G.H
   p = GAP.Globals.GroupHomomorphismByFunction(G.X,H.X,y->GAP.Globals.Image(f,y))
   return _hom_from_gap_map(G,H,p)
end

# start part on subgroups
function _as_subgroup_bare(G::SemidirectProductGroup{S,T}, H::GapObj) where { S , T }
#  t = G.X==H
  return SemidirectProductGroup(H, G.N, G.H, G.f, G.X, false)
end

function _as_subgroup(G::SemidirectProductGroup{S,T}, H::GapObj, ::Type{U}) where { T, S, U }
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x::U -> group_element(G, x.X))
end

function _as_subgroup(G::SemidirectProductGroup{S,T}, H::GapObj) where S <: GAPGroup where T <: GAPGroup
  return _as_subgroup(G, H, GAPGroupElem{SemidirectProductGroup{S,T}})
end

function sub(G::SemidirectProductGroup{S,T}, elms::Vector{<:GAPGroupElem{SemidirectProductGroup{S,T}}}) where { S, T }
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elms])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(G, H)
end

function sub(L::U...) where U <: GAPGroupElem{SemidirectProductGroup{S,T}} where { S, T }
   length(L)>0 || throw(ArgumentError("Empty list"))
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function sub(L::Vector{<:GAPGroupElem{SemidirectProductGroup{S,T}}}) where { S, T }
   length(L)>0 || throw(ArgumentError("Empty list"))
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function Base.show(io::IO, x::SemidirectProductGroup)
   if x.isfull
      print(io, "SemidirectProduct( ", GAP.gap_to_julia(GAP.Globals.StringViewObj(x.N.X)), " , ", GAP.gap_to_julia(GAP.Globals.StringView(x.H.X))," )")
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X)))
   end
end

################################################################################
#
#  Wreath products
#  
################################################################################

"""
    wreath_product(G::T, H::S, a::GAPGroupHomomorphism{S,PermGroup})
    wreath_product(G::T, H::PermGroup) where T<: Group

Return the wreath product of the group `G` and the group `H`, where `H` acts on `n` copies of `G` through the homomorphism `a` from `H` to a permutation group, and `n` is the number of moved points of `Image(a)`.

If `a` is not specified, then `H` must be a group of permutations. In this case, `n` is NOT the number of moved points, but the degree of `H`.

If `W` is a wreath product of `G` and `H`, {`g_1`, ..., `g_n`} are elements of `G` and `h` in `H`, the element `(g_1, ..., h)` of `W` can be obtained by typing
```
    W(g_1,...,g_n, h).
```
"""
function wreath_product(G::T, H::PermGroup) where T<: GAPGroup
   if Set{Int}(GAP.gap_to_julia(GAP.Globals.MovedPoints(H.X)))==Set(1:H.deg)
      Wgap=GAP.Globals.WreathProduct(G.X,H.X)
      return WreathProductGroup(Wgap,G,H,id_hom(H),Wgap,true)
   else
      S=symmetric_group(H.deg)
      Wgap=GAP.Globals.WreathProduct(G.X,S.X)
      W1=GAP.Globals.PreImage(GAP.Globals.Projection(Wgap),H.X)
   # not id_hom(H) because I need NrMovedPoints(Image(a))==degree(H), see function embedding
      return WreathProductGroup(W1,G,H,id_hom(symmetric_group(H.deg)),Wgap,true)
   end
end

function wreath_product(G::T, H::S, a::GAPGroupHomomorphism{S,PermGroup}) where S<:GAPGroup where T<:GAPGroup
   Wgap=GAP.Globals.WreathProduct(G.X,H.X,a.map)
   return WreathProductGroup(Wgap,G,H,a,Wgap,true)
end

# return the element (L[1],L[2],...) in W
function (W::WreathProductGroup)(L::Union{GAPGroupElem{T},GAPGroupElem{PermGroup}}...) where T <: GAPGroup
   d = GAP.Globals.NrMovedPoints(GAP.Globals.Image(W.a.map))
   length(L)==d+1 || throw(ArgumentError("Wrong number of arguments"))
   for i in 1:d @assert L[i] in W.G "Wrong input" end
   L[d+1] in W.H || throw(ArgumentError("Wrong input"))
# simply put parent(L[d+1])==W.H does not work. Example: if I want to write explicitly a permutation in H proper subgroup of Sym(n).
   arr = [GAP.Globals.Image(GAP.Globals.Embedding(W.Xfull,i),L[i].X) for i in 1:length(L)]
   xgap = prod(arr)
   xgap in W.X || throw(ArgumentError("Element not in the group"))
   return group_element(W,xgap)
end


"""
    normal_subgroup(W::WreathProductGroup)

Return `G`, where `W` is the wreath product of `G` and `H`.
"""
normal_subgroup(W::WreathProductGroup) = W.G

"""
    acting_subgroup(W::WreathProductGroup)

Return `H`, where `W` is the wreath product of `G` and `H`.
"""
acting_subgroup(W::WreathProductGroup) = W.H

"""
    homomorphism_of_wreath_product(G::WreathProductGroup)

If `W` is the wreath product of `G` and `H`, then return the homomorphism `f` from `H` to `Sym(n)`, where `n` is the number of copies of `G`.
"""
homomorphism_of_wreath_product(G::WreathProductGroup) = G.a

"""
    isfull_wreath_product(G::WreathProductGroup)

Return whether `G` is a wreath product of two groups, instead of a proper subgroup.
"""
isfull_wreath_product(G::WreathProductGroup) = G.isfull

"""
    projection(G::WreathProductGroup)

Return the projection of `wreath_product(G,H)` into the permutation group `H`.
"""
function projection(W::WreathProductGroup)
  # @assert W.isfull "Projection not defined for proper subgroups of wreath products"
   f=GAP.Globals.Projection(W.Xfull)
   H = W.H
   p = GAP.Globals.GroupHomomorphismByFunction(W.X,H.X,y->GAP.Globals.Image(f,y))
   return _hom_from_gap_map(W,H,p)
end

"""
    embedding(G::WreathProductGroup, n::Integer)

Return the embedding of the `n`-th component of `G` into `G`.
"""
function embedding(W::WreathProductGroup, n::Base.Integer)
   W.isfull || throw(ArgumentError("Embedding not defined for proper subgroups of wreath products"))
   n <= GAP.Globals.NrMovedPoints(GAP.Globals.Image(W.a.map))+1 || throw(ArgumentError("n is too big"))
   f = GAP.Globals.Embedding(W.Xfull,n)
   if n== GAP.Globals.NrMovedPoints(GAP.Globals.Image(W.a.map))+1 C=W.H
   else C=W.G
   end
   return GAPGroupHomomorphism{typeof(C),WreathProductGroup}(C,W,f)
end

Base.show(io::IO, x::WreathProductGroup) = print(io, GAP.gap_to_julia(GAP.Globals.StringView(x.X)))


# start part on subgroups
#TODO : to be fixed
function _as_subgroup_bare(W::WreathProductGroup, X::GapObj)
#   t = X==W.X
  return WreathProductGroup(X, W.G, W.H, W.a, W.Xfull, false)
end

function _as_subgroup(W::WreathProductGroup, H::GapObj, ::Type{U}) where U
  H1 = _as_subgroup_bare(W, H)
  return H1, hom(H1, W, x::U -> group_element(W, x.X))
end

function _as_subgroup(W::WreathProductGroup, H::GapObj)
  return _as_subgroup(W, H, GAPGroupElem{WreathProductGroup})
end

function sub(W::WreathProductGroup, elms::Vector{<:GAPGroupElem{WreathProductGroup}})
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elms])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(W, H)
end

function sub(L::T...) where T <: GAPGroupElem{WreathProductGroup}
   length(L)>0 || throw(ArgumentError("Empty list"))
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function sub(L::Vector{<:GAPGroupElem{WreathProductGroup}})
   length(L)>0 || throw(ArgumentError("Empty list"))
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end


