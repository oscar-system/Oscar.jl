################################################################################
#
#  Managing direct products
#
################################################################################

export
    as_matrix_group,
    as_perm_group,
    as_polycyclic_group,
    cartesian_power,
    direct_product,
    embedding,
    factorofdirectproduct,
    homomorphism_of_semidirect_product,
    isfull_direct_product,
    isfull_semidirect_product,
    isfull_wreath_product,
    number_of_factors,
    projection,
    semidirect_product,
    sub,
    wreath_product,
    write_as_full


"""
    direct_product(L::AbstractVector{<:GAPGroup})
    direct_product(L::GAPGroup...)
Return the direct product of the groups in the vector `L`.
"""
function direct_product(L::AbstractVector{<:GAPGroup})
   X = GAP.Globals.DirectProduct(GAP.julia_to_gap([G.X for G in L]))
   n = length(L)
   return DirectProductOfGroups(X,L,n,X,true)
end

function direct_product(L::GAPGroup...)
   return direct_product([x for x in L])
end

"""
    number_of_factors(G::DirectProductOfGroups)
Return the number of factors of `G`.
"""
number_of_factors(G::DirectProductOfGroups) = G.n

"""
    cartesian_power(G::T, n::Int)
Return the direct product of `n` copies of `G`, of type `T`.
"""
function cartesian_power(G::T, n::Base.Integer) where T <: GAPGroup
   L = [G for i in 1:n]
   return direct_product(L)
end

"""
    factorofdirectproduct(G::DirectProductOfGroups, j::Int)
Return the `j`-th factor of `G`.
"""
function factorofdirectproduct(G::DirectProductOfGroups, j::Base.Integer)
   if j in 1:G.n return G.L[j]
   else throw(ArgumentError("index not valid"))
   end
end

"""
    as_perm_group(G::DirectProductOfGroups)
If `G` is direct product of permutations groups, return `G` as permutation group.
"""
function as_perm_group(G::DirectProductOfGroups)
   if [typeof(H)==PermGroup for H in G.L]==[true for i in 1:G.n]
      return PermGroup(G.X, GAP.Globals.Maximum(GAP.Globals.MovedPoints(G.X)))
   else
      throw(ArgumentError("The group is not a permutation group"))
   end
end

"""
    as_polycyclic_group(G::DirectProductOfGroups)
If `G` is direct product of polycyclic groups, return `G` as polycyclic group.
"""
function as_polycyclic_group(G::DirectProductOfGroups)
   if [typeof(H)==PcGroup for H in G.L]==[true for i in 1:G.n]
      return PcGroup(G.X)
   else
      throw(ArgumentError("The group is not a polycyclic group"))
   end
end

"""
    as_matrix_group(G::DirectProductOfGroups)
If `G` is direct product of matrix groups over the ring `R` of dimension `n_1`, ... , `n_k` respectively, return `G` as matrix group over the ring `R` of dimension `n_1 + ... + n_k`.
"""
function as_matrix_group(G::DirectProductOfGroups)
#   if S==MatrixGroup && T==MatrixGroup && GAP.Globals.FieldOfMatrixGroup(G.G1.X)==GAP.Globals.FieldOfMatrixGroup(G.G2.X)
   if [typeof(H)==MatrixGroup for H in G.L]==[true for i in 1:G.n] && length(Set([base_ring(H) for H in G.L]))==1
      return MatrixGroup(G.X)
   else
      throw(ArgumentError("The group is not a matrix group"))
   end
end

"""
    embedding(G::DirectProductOfGroups, j::Integer)
Return the embedding of the `j`-th component of `G` into `G`, for `j` = 1,...,#factors of `G`.
"""
function embedding(G::DirectProductOfGroups, j::Base.Integer)
   if !G.isfull
      throw(ArgumentError("Embedding is not defined for proper subgroups of direct products"))
   end
   if j in 1:G.n
      f=GAP.Globals.Embedding(G.X,j)
      gr=G.L[j]
      typeG=typeof(gr)
   else
      throw(ArgumentError("index not valid"))
   end
   return GAPGroupHomomorphism{typeG,DirectProductOfGroups}(gr,G,f)
end

"""
    projection(G::DirectProductOfGroups, j::Integer)
Return the projection of `G` into the `j`-th component of `G`, for `j` = 1,...,#factors of `G`.
"""
function projection(G::DirectProductOfGroups, j::Base.Integer)
  # @assert W.isfull "Projection not defined for proper subgroups of wreath products"
   f=GAP.Globals.Projection(G.Xfull,j)
   if j in 1:G.n
      H = G.L[j]
      if !G.isfull
         Gf = DirectProductOfGroups(G.Xfull,G.L,G.n,G.Xfull,true)
         g = embedding(Gf,G)
         p = g*GAPGroupHomomorphism{DirectProductOfGroups,typeof(H)}(Gf,H,f)
      else
         p = GAPGroupHomomorphism{DirectProductOfGroups,typeof(H)}(G,H,f)
      end
      return hom(G,H,x->p(x))
   else
      throw(ArgumentError("index not valid"))
   end
end

function (G::DirectProductOfGroups)(V::Vector{GAPGroupElem})
   Gfull=DirectProductOfGroups(G.Xfull,G.L,G.n,G.Xfull,true)
   @assert length(V)==G.n "Vector of wrong length"
   x = prod([embedding(Gfull,j)(V[j]) for j in 1:length(V)])
   @assert x in G "Element not in the group"
   return group_element(G,x.X)
end

function (G::DirectProductOfGroups)(V::GAPGroupElem...)
   Gfull=DirectProductOfGroups(G.Xfull,G.L,G.n,G.Xfull,true)
   @assert length(V)==G.n "Vector of wrong length"
   x = prod([embedding(Gfull,j)(V[j]) for j in 1:length(V)])
   @assert x in G "Element not in the group"
   return group_element(G,x.X)
end

# start part on subgroups
function _as_subgroup_bare(G::DirectProductOfGroups, H::GapObj)
  t = H==G.X
  return DirectProductOfGroups(H, G.L, G.n, G.X, t)
end

function _as_subgroup(G::DirectProductOfGroups, H::GapObj, ::Type{U}) where U
  H1 = _as_subgroup_bare(G, H)
  return H1, hom(H1, G, x::U -> group_element(G, x.X))
end

function _as_subgroup(G::DirectProductOfGroups, H::GapObj)
  return _as_subgroup(G, H, GAPGroupElem{DirectProductOfGroups})
end

function sub(G::DirectProductOfGroups, elms::Vector{GAPGroupElem{DirectProductOfGroups}})
  elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elms])
  H = GAP.Globals.Group(elems_in_GAP)
  #H is the group. I need to return the inclusion map too
  return _as_subgroup(G, H)
end

function sub(L::GAPGroupElem{DirectProductOfGroups}...)
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function sub(L::Vector{GAPGroupElem{DirectProductOfGroups}})
   if length(L)==0 throw(ArgumentError("Empty list")) end
   l=collect(L)
   @assert all(x -> parent(x) == parent(l[1]), l)
   return sub(parent(l[1]),l)
end

function Base.show(io::IO, G::DirectProductOfGroups)
   if G.isfull
      print(io, "DirectProduct of ")
      display(G.L)
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(G.X)))
   end
end

#=
function Base.show(io::IO, x::GAPGroupElem{DirectProductOfGroups{S,T}}) where { S, T }
   print(io, "DirectProductElement( ", GAP.gap_to_julia(GAP.Globals.StringViewObj(projection(parent(x),1)(x).X)), " , ", GAP.gap_to_julia(GAP.Globals.StringViewObj(projection(parent(x),2)(x).X)), " )")
end
=#

# if a subgroup of a direct product of groups is also a direct product of groups
"""
    write_as_full(G::DirectProductOfGroups)
If `G` is a subgroup of the direct product `G_1 x ... x G_n` such that `G = H_1 x ... x H_n` for `H_i` subgroup of `G_i`, return `G` as full direct product of the `H_i`. If such `H_i` do not exist, an ERROR is returned.
"""
function write_as_full(G::DirectProductOfGroups)
   if G.isfull
      return G
   else
      LK = [image(projection(G,j))[1] for j in 1:G.n]
      H = direct_product(LK)
# index(H,G)==1 does not work because it does not recognize G as a subgroup of H
      if order(H)==order(G)
         return H
      else
         throw(ArgumentError("G is not a direct product of groups"))
      end
   end
end

"""
    isfull_direct_product(G::DirectProductOfGroups)
Return whether `G` is direct product of its factors (`false` if it is a proper subgroup).
"""
isfull_direct_product(G::DirectProductOfGroups) = G.isfull

Base.:^(H::DirectProductOfGroups, y::GAPGroupElem) = sub([h^y for h in gens(H)])[1]


################################################################################
#
#  Semidirect products
#  
################################################################################

"""
    semidirect_product(N::S, f::GAPGroupHomomorphism, H::T)
Return the semidirect product of `N` and `H`, of type ``SemidirectProductOfGroups{S,T}``, where `f` is a homomorphism from `H` to `Aut`(`N`).
"""
function semidirect_product(N::S, f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}, H::T) where S <: GAPGroup where T <: GAPGroup
   sdp=GAP.Globals.SemidirectProduct(H.X,f.map,N.X)
   return SemidirectProductOfGroups(sdp,N,H,f,sdp,true)
end

# return the element (a,b) in G
function (G::SemidirectProductOfGroups{S,T})(a::GAPGroupElem{S},b::GAPGroupElem{T}) where {S,T}
# simply put parent(L[d+1])==W.H does not work. Example: if I want to write explicitly a permutation in H proper subgroup of Sym(n).
   Gfull = SemidirectProductOfGroups(G.Xfull,G.N,G.H,G.f,G.Xfull,true)
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
   return GAPGroupHomomorphism{typeG,SemidirectProductOfGroups{S,T}}(gr,G,f)
end



"""
    projection(G::SemidirectProductOfGroups{S,T}, n::Integer)
Return the projection of `G` into the second component of `G`.
"""
function projection(G::SemidirectProductOfGroups{S,T}) where S where T
   f=GAP.Globals.Projection(G.Xfull)
   if !G.isfull
      Gf = SemidirectProductOfGroups(G.Xfull,G.N,G.H,G.f,G.Xfull,true)
      g = embedding(Gf,G)
      return g*GAPGroupHomomorphism{SemidirectProductOfGroups{S,T},T}(Gf,G.H,f)
   else
      return GAPGroupHomomorphism{SemidirectProductOfGroups{S,T},T}(G,G.H,f)
   end
end


# start part on subgroups
function _as_subgroup_bare(G::SemidirectProductOfGroups{S,T}, H::GapObj) where { S , T }
  return SemidirectProductOfGroups(H, G.N, G.H, G.f, G.X, false)
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
      print(io, "SemidirectProduct( ", GAP.gap_to_julia(GAP.Globals.StringViewObj(x.N.X)), " , ", GAP.gap_to_julia(GAP.Globals.StringView(x.H.X))," )")
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

