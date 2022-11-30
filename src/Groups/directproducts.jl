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
    is_full_direct_product,
    is_full_semidirect_product,
    is_full_wreath_product,
    normal_subgroup,
    number_of_factors,
    projection,
    semidirect_product,
    sub,
    wreath_product,
    write_as_full

"""
    direct_product(L::AbstractVector{<:GAPGroup}; morphisms)
    direct_product(L::GAPGroup...)

Return the direct product of the groups in the collection `L`.

The parameter `morphisms` is `false` by default. If it is set `true`, then
the output is a triple (`G`, `emb`, `proj`), where `emb` and `proj` are the
vectors of the embeddings (resp. projections) of the direct product `G`.

# Examples
```jldoctest
julia> H = symmetric_group(3)
Sym( [ 1 .. 3 ] )

julia> K = symmetric_group(2)
Sym( [ 1 .. 2 ] )

julia> G = direct_product(H,K)
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )

julia> elements(G)
12-element Vector{Oscar.BasicGAPGroupElem{DirectProductGroup}}:
 ()
 (4,5)
 (2,3)
 (2,3)(4,5)
 (1,2)
 (1,2)(4,5)
 (1,2,3)
 (1,2,3)(4,5)
 (1,3,2)
 (1,3,2)(4,5)
 (1,3)
 (1,3)(4,5)
```
"""
function direct_product(L::AbstractVector{<:GAPGroup}; morphisms=false)
   X = GAP.Globals.DirectProduct(GapObj([G.X for G in L]))
   DP = DirectProductGroup(X,L,X,true)
   if morphisms
      emb = [GAPGroupHomomorphism(L[i],DP,GAP.Globals.Embedding(X,i)) for i in 1:length(L)]
      proj = [GAPGroupHomomorphism(DP,L[i],GAP.Globals.Projection(X,i)) for i in 1:length(L)]
      return DP, emb, proj
   else
      return DP
   end
end

function direct_product(L::GAPGroup...; morphisms=false)
   return direct_product([x for x in L]; morphisms=morphisms)
end

"""
    inner_direct_product(L::AbstractVector{T}; morphisms)
    inner_direct_product(L::T...)

Return a direct product of groups of the same type `T` as a group of type
`T`. It works for `T` of the following types:
- `PermGroup`, `PcGroup`, `FPGroup`.

The parameter `morphisms` is `false` by default. If it is set `true`, then
the output is a triple (`G`, `emb`, `proj`), where `emb` and `proj` are the
vectors of the embeddings (resp. projections) of the direct product `G`.
"""
function inner_direct_product(L::AbstractVector{T}; morphisms=false) where T<:Union{PcGroup,FPGroup}
   P = GAP.Globals.DirectProduct(GapObj([G.X for G in L]))
   DP = T(P)
   if morphisms
      emb = [GAPGroupHomomorphism(L[i],DP,GAP.Globals.Embedding(P,i)) for i in 1:length(L)]
      proj = [GAPGroupHomomorphism(DP,L[i],GAP.Globals.Projection(P,i)) for i in 1:length(L)]
      return DP, emb, proj
   else
      return DP
   end
end

function inner_direct_product(L::AbstractVector{PermGroup}; morphisms=false) 
   P = GAP.Globals.DirectProductOfPermGroupsWithMovedPoints(GapObj([G.X for G in L]), GAP.Obj([collect(1:degree(G)) for G in L], recursive = true))
   DP = PermGroup(P, sum([degree(G) for G in L], init = 0))
   if morphisms
      emb = [GAPGroupHomomorphism(L[i],DP,GAP.Globals.Embedding(P,i)) for i in 1:length(L)]
      proj = [GAPGroupHomomorphism(DP,L[i],GAP.Globals.Projection(P,i)) for i in 1:length(L)]
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
function cartesian_power(G::GAPGroup, n::Int)
   L = [G for i in 1:n]
   return direct_product(L)
end

"""
    inner_cartesian_power(G::T, n::Int; morphisms)

Return the direct product of `n` copies of `G` as group of type `T`.

The parameter `morphisms` is `false` by default. If it is set `true`, then
the output is a triple (`G`, `emb`, `proj`), where `emb` and `proj` are the
vectors of the embeddings (resp. projections) of the direct product `G`.
"""
function inner_cartesian_power(G::T, n::Int; morphisms=false) where T <: GAPGroup
   L = [G for i in 1:n]
   return inner_direct_product(L; morphisms=morphisms)
end

"""
    factor_of_direct_product(G::DirectProductGroup, j::Int)

Return the `j`-th factor of `G`.
"""
function factor_of_direct_product(G::DirectProductGroup, j::Int)
   if j in 1:length(G.L) return G.L[j]
   else throw(ArgumentError("index not valid"))
   end
end

"""
    as_perm_group(G::DirectProductGroup)

If `G` is direct product of permutations groups, return `G` as permutation group.

# Examples
```jldoctest
julia> H = symmetric_group(3)
Sym( [ 1 .. 3 ] )

julia> K = symmetric_group(2)
Sym( [ 1 .. 2 ] )

julia> G = direct_product(H,K)
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )

julia> as_perm_group(G)
Group([ (1,2,3), (1,2), (4,5) ])
```
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
    embedding(G::DirectProductGroup, j::Int)

Return the embedding of the `j`-th component of `G` into `G`, for `j` = 1,...,#factors of `G`.

# Examples
```jldoctest
julia> H = symmetric_group(3)
Sym( [ 1 .. 3 ] )

julia> K = symmetric_group(2)
Sym( [ 1 .. 2 ] )

julia> G = direct_product(H,K)
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )

julia> emb1 = embedding(G,1)
Group homomorphism from
Sym( [ 1 .. 3 ] )
to
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )

julia> h = perm(H,[2,3,1])
(1,2,3)

julia> emb1(h)
(1,2,3)

julia> emb2 = embedding(G,2)
Group homomorphism from
Sym( [ 1 .. 2 ] )
to
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )

julia> k = perm(K,[2,1])
(1,2)

julia> emb2(k)
(4,5)

julia> emb1(h)*emb2(k)
(1,2,3)(4,5)
```
"""
function embedding(G::DirectProductGroup, j::Int)
   j in 1:length(G.L) || throw(ArgumentError("index not valid"))
   G.isfull || throw(ArgumentError("Embedding is not defined for proper subgroups of direct products"))
   f=GAP.Globals.Embedding(G.X,j)
   gr=G.L[j]
   typeG=typeof(gr)
   return GAPGroupHomomorphism(gr,G,f)
end

"""
    projection(G::DirectProductGroup, j::Int)

Return the projection of `G` into the `j`-th component of `G`, for `j` = 1,...,#factors of `G`.

# Examples
```jldoctest
julia> H = symmetric_group(3)
Sym( [ 1 .. 3 ] )

julia> K = symmetric_group(2)
Sym( [ 1 .. 2 ] )

julia> G = direct_product(H,K)
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )

julia> proj1 = projection(G,1)
Group homomorphism from
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )
to
Sym( [ 1 .. 3 ] )

julia> proj2 = projection(G,2)
Group homomorphism from
DirectProduct of
 Sym( [ 1 .. 3 ] )
 Sym( [ 1 .. 2 ] )
to
Sym( [ 1 .. 2 ] )

julia> g = perm([2,3,1,5,4])
(1,2,3)(4,5)

julia> proj1(g)
(1,2,3)

julia> proj2(g)
(1,2)
```
"""
function projection(G::DirectProductGroup, j::Int)
   f=GAP.Globals.Projection(G.Xfull,j)
   j in 1:length(G.L) || throw(ArgumentError("index not valid"))
   H = G.L[j]
   p = GAPWrap.GroupHomomorphismByFunction(G.X,H.X,y->GAPWrap.Image(f,y))
   return GAPGroupHomomorphism(G,H,p)
end


function (G::DirectProductGroup)(V::AbstractVector{<:GAPGroupElem})
   length(V)==length(G.L) || throw(ArgumentError("Wrong number of entries"))
   arr = [GAPWrap.Image(GAP.Globals.Embedding(G.Xfull,i),V[i].X) for i in 1:length(V)]
   xgap = prod(arr)
   xgap in G.X || throw(ArgumentError("Element not in the group"))
   return group_element(G,xgap)
end

function (G::DirectProductGroup)(V::GAPGroupElem...)
   return G([x for x in V])
end

function _as_subgroup_bare(G::DirectProductGroup, H::GapObj)
#  t = H==G.X
  return DirectProductGroup(H, G.L, G.X, false)
end

function Base.show(io::IO, G::DirectProductGroup)
   if G.isfull
      print(io, "DirectProduct of ")
      for x in G.L print(io, "\n ", x) end
   else
      print(io, String(GAPWrap.StringViewObj(G.X)))
   end
end

# if a subgroup of a direct product of groups is also a direct product of groups
@doc Markdown.doc"""
    write_as_full(G::DirectProductGroup)

If `G` is a subgroup of the direct product
$G_1 \times G_2 \times \cdots \times G_n$ such that `G` has the form
$H_1 \times H_2 \times \cdots \times H_n$, for subgroups $H_i$ of $G_i$,
return this full direct product of the $H_i$.

An exception is thrown if such $H_i$ do not exist.
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
    is_full_direct_product(G::DirectProductGroup)

Return whether `G` is direct product of its factors (`false` if it is a proper subgroup).
"""
is_full_direct_product(G::DirectProductGroup) = G.isfull

Base.:^(H::DirectProductGroup, y::GAPGroupElem) = sub([h^y for h in gens(H)]...)[1]


################################################################################
#
#  Semidirect products
#  
################################################################################

"""
    semidirect_product(N::S, f::GAPGroupHomomorphism, H::T)

Return the semidirect product of `N` and `H`, of type `SemidirectProductGroup{S,T}`,
where `f` is a group homomorphism from `H` to the automorphism group of `N`.
"""
function semidirect_product(N::S, f::GAPGroupHomomorphism{T,AutomorphismGroup{S}}, H::T) where S <: GAPGroup where T <: GAPGroup
   sdp=GAP.Globals.SemidirectProduct(H.X,f.map,N.X)
   return SemidirectProductGroup{S,T}(sdp,N,H,f,sdp,true)
end

# return the element (a,b) in G
function (G::SemidirectProductGroup{S,T})(a::GAPGroupElem{S},b::GAPGroupElem{T}) where {S,T}
# simply put parent(L[d+1])==W.H does not work. Example: if I want to write explicitly a permutation in H proper subgroup of Sym(n).
   xgap = GAPWrap.Image(GAP.Globals.Embedding(G.Xfull,1),b.X)*GAPWrap.Image(GAP.Globals.Embedding(G.Xfull,2),a.X)
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

Return `f,` where `G` is the semidirect product of the normal subgroup `N` and
the group `H` acting on `N` via the homomorphism `h`.
"""
homomorphism_of_semidirect_product(G::SemidirectProductGroup) = G.f

"""
    is_full_semidirect_product(G::SemidirectProductGroup)

Return whether `G` is a semidirect product of two groups, instead of a proper subgroup.
"""
is_full_semidirect_product(G::SemidirectProductGroup) = G.isfull

"""
    embedding(G::SemidirectProductGroup, n::Int)

Return the embedding of the `n`-th component of `G` into `G`, for `n` = 1,2.
It is not defined for proper subgroups of semidirect products.
"""
function embedding(G::SemidirectProductGroup, n::Int)
   @assert G.isfull "Embedding not defined for proper subgroups of semidirect products"
   if n==1
      f=GAP.Globals.Embedding(G.X,2)
      gr=G.N
   elseif n==2
      f=GAP.Globals.Embedding(G.X,1)
      gr=G.H
   else
      throw(ArgumentError("n must be 1 or 2"))
   end
   return GAPGroupHomomorphism(gr,G,f)
end

"""
    projection(G::SemidirectProductGroup)

Return the projection of `G` into the second component of `G`.
"""
function projection(G::SemidirectProductGroup)
   f=GAP.Globals.Projection(G.Xfull)
   H = G.H
   p = GAPWrap.GroupHomomorphismByFunction(G.X,H.X,y->GAPWrap.Image(f,y))
   return GAPGroupHomomorphism(G,H,p)
end

function _as_subgroup_bare(G::SemidirectProductGroup{S,T}, H::GapObj) where { S , T }
#  t = G.X==H
  return SemidirectProductGroup{S,T}(H, G.N, G.H, G.f, G.X, false)
end

function Base.show(io::IO, x::SemidirectProductGroup)
   if x.isfull
      print(io, "SemidirectProduct( ", String(GAPWrap.StringViewObj(x.N.X)),
                " , ", String(GAP.Globals.StringView(x.H.X))," )")
   else
      print(io, String(GAPWrap.StringViewObj(x.X)))
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

Return the wreath product of the group `G` and the group `H`, where `H` acts
on `n` copies of `G` through the homomorphism `a` from `H` to a permutation
group, and `n` is the number of moved points of `Image(a)`.

If `a` is not specified, then `H` must be a group of permutations. In this
case, `n` is NOT the number of moved points, but the degree of `H`.

If `W` is a wreath product of `G` and `H`, {`g_1`, ..., `g_n`} are elements of
`G` and `h` in `H`, the element `(g_1, ..., h)` of `W` can be obtained by
typing
```
    W(g_1,...,g_n, h).
```

# Examples
```jldoctest
julia> G = cyclic_group(3)
<pc group of size 3 with 1 generator>

julia> H = symmetric_group(2)
Sym( [ 1 .. 2 ] )

julia> W = wreath_product(G,H)
<group of size 18 with 2 generators>

julia> a = gen(W,1)
WreathProductElement(f1,<identity> of ...,())

julia> b = gen(W,2)
WreathProductElement(<identity> of ...,<identity> of ...,(1,2))

julia> a*b
WreathProductElement(f1,<identity> of ...,(1,2))
```
"""
function wreath_product(G::T, H::PermGroup) where T<: GAPGroup
   if Set{Int}(GAP.Globals.MovedPoints(H.X))==Set(1:H.deg)
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
   d = GAP.Globals.NrMovedPoints(GAPWrap.Image(W.a.map))
   length(L)==d+1 || throw(ArgumentError("Wrong number of arguments"))
   for i in 1:d @assert L[i] in W.G "Wrong input" end
   L[d+1] in W.H || throw(ArgumentError("Wrong input"))
# simply put parent(L[d+1])==W.H does not work. Example: if I want to write explicitly a permutation in H proper subgroup of Sym(n).
   arr = [GAPWrap.Image(GAP.Globals.Embedding(W.Xfull,i),L[i].X) for i in 1:length(L)]
   xgap = prod(arr)
   xgap in W.X || throw(ArgumentError("Element not in the group"))
   return group_element(W,xgap)
end


"""
    normal_subgroup(W::WreathProductGroup)

Return `G`, where `W` is the wreath product of `G` and `H`.

# Examples
```jldoctest
julia> G = cyclic_group(3)
<pc group of size 3 with 1 generator>

julia> H = symmetric_group(2)
Sym( [ 1 .. 2 ] )

julia> W = wreath_product(G,H)
<group of size 18 with 2 generators>

julia> normal_subgroup(W)
<pc group of size 3 with 1 generator>
```
"""
normal_subgroup(W::WreathProductGroup) = W.G

"""
    acting_subgroup(W::WreathProductGroup)

Return `H`, where `W` is the wreath product of `G` and `H`.

# Examples
```jldoctest
julia> G = cyclic_group(3)
<pc group of size 3 with 1 generator>

julia> H = symmetric_group(2)
Sym( [ 1 .. 2 ] )

julia> W = wreath_product(G,H)
<group of size 18 with 2 generators>

julia> acting_subgroup(W)
Sym( [ 1 .. 2 ] )
```
"""
acting_subgroup(W::WreathProductGroup) = W.H

"""
    homomorphism_of_wreath_product(G::WreathProductGroup)

If `W` is the wreath product of `G` and `H`, then return the homomorphism `f`
from `H` to `Sym(n)`, where `n` is the number of copies of `G`.
"""
homomorphism_of_wreath_product(G::WreathProductGroup) = G.a

"""
    is_full_wreath_product(G::WreathProductGroup)

Return whether `G` is a wreath product of two groups, instead of a proper subgroup.
"""
is_full_wreath_product(G::WreathProductGroup) = G.isfull

"""
    projection(G::WreathProductGroup)

Return the projection of `wreath_product(G,H)` onto the permutation group `H`.
"""
function projection(W::WreathProductGroup)
  # @assert W.isfull "Projection not defined for proper subgroups of wreath products"
   f=GAP.Globals.Projection(W.Xfull)
   H = W.H
   p = GAPWrap.GroupHomomorphismByFunction(W.X,H.X,y->GAPWrap.Image(f,y))
   return GAPGroupHomomorphism(W,H,p)
end

"""
    embedding(G::WreathProductGroup, n::Int)

Return the embedding of the `n`-th component of `G` into `G`.
"""
function embedding(W::WreathProductGroup, n::Int)
   W.isfull || throw(ArgumentError("Embedding not defined for proper subgroups of wreath products"))
   n <= GAP.Globals.NrMovedPoints(GAPWrap.Image(W.a.map))+1 || throw(ArgumentError("n is too big"))
   f = GAP.Globals.Embedding(W.Xfull,n)
   if n== GAP.Globals.NrMovedPoints(GAPWrap.Image(W.a.map))+1 C=W.H
   else C=W.G
   end
   return GAPGroupHomomorphism(C,W,f)
end

Base.show(io::IO, x::WreathProductGroup) = print(io, String(GAP.Globals.StringView(x.X)))


#TODO : to be fixed
function _as_subgroup_bare(W::WreathProductGroup, X::GapObj)
#   t = X==W.X
  return WreathProductGroup(X, W.G, W.H, W.a, W.Xfull, false)
end



