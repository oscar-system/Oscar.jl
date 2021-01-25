
# return an element in the centralizer of x in GL(n,F) with determinant d

# first: brute force way
function _elem_given_det(x,d)
   C,e = centralizer(GL(x.parent.deg, x.parent.ring),x)
   U,fa = unit_group(x.parent.ring)
   GA,ea = sub(U, [preimage(fa,det(g)) for g in gens(C)])
   l = preimage(ea,preimage(fa,d))
   return prod([C[i]^Int(l[i]) for i in 1:ngens(C)])
end



########################################################################
#
# Orders
#
########################################################################

function _GL_order(n::Int, q::fmpz)
   res = q^div(n*(n-1),2)
   for i in 1:n res *= (q^i-1) end
   return res
end

_GL_order(n::Int, F::Ring) = _GL_order(n, order(F))

function _SL_order(n::Int, q::fmpz)
   res = q^div(n*(n-1),2)
   for i in 2:n res *= (q^i-1) end
   return res
end

_SL_order(n::Int, F::Ring) = _SL_order(n, order(F))

########################################################################
#
# Centralizer in GL
#
########################################################################


# returns as matrices
function _gens_for_GL(n::Int, F::Ring)
   n !=1 || return [matrix(F,1,1,[primitive_element(F)])]
   if order(F)==2
      h1 = identity_matrix(F,n)
      h1[1,2] = 1
      h2 = zero_matrix(F,n,n)
      h2[1,n] = 1
      for i in 1:n-1 h2[i+1,i] = -1 end
      return h1,h2
   else
      h1 = identity_matrix(F,n)
      h1[1,1] = primitive_element(F)
      h2 = zero_matrix(F,n,n)
      h2[1,1] = -1
      h2[1,n] = 1
      for i in 1:n-1 h2[i+1,i] = -1 end
      return h1,h2
   end
end

# returns as matrices
# does the matrix above with F = F[x]/(f), but every entry is replaced by a diagonal join of D corresponding blocks
# ASSUMPTION: deg(f) > 1
function _gens_for_GL_matrix(f::PolyElem, n::Int, F::Ring; D=1)
   C = companion_matrix(f)
   CP = evaluate(_centralizer(f),C)            # matrix of maximal order in the centralizer of the companion matrix

   if n==1 return [diagonal_join([CP for i in 1:D])] end
   h1 = identity_matrix(F,n*degree(f)*D)
   insert_block!(h1,diagonal_join([CP for i in 1:D]),1,1)
   h2 = zero_matrix(F,n*degree(f)*D,n*degree(f)*D)
   for i in 1:(n-1)*degree(f)*D h2[i+degree(f)*D,i]=-1 end
   for i in 1:degree(f)*D
      h2[i,i]=-1
      h2[i,i+(n-1)*degree(f)*D]=1
   end
   return h1,h2      
end


# V = vector of integers of the dimensions of Jordan blocks
# return the generators for the centralizers of the unipotent element
# assumes V is sorted (e.g. [1,1,1,2,3,3])
function _centr_unipotent(F::Ring, V::AbstractVector{Int}; isSL=false) 
   n = sum(V)
   _lambda = gen(F)  # yes, gen(F) is correct; we don't need a primitive element in this case

   # L = multiset(V)
   L=[[V[1],1]]
   for i in 2:length(V)
      if V[i]==L[length(L)][1]
         L[length(L)][2] += 1
      else
         push!(L, [V[i],1])
      end
   end
   listgens = MatElem[]

   # generators for GL + internal diagonal blocks
   pos=1
   for l in L
      v_g = isSL ? _gens_for_SL(l[2],F) : _gens_for_GL(l[2],F)
      for x in v_g
         z = block_matrix(l[2],l[2],[x[i,j]*identity_matrix(F,l[1]) for i in 1:l[2] for j in 1:l[2]])
         z = insert_block(identity_matrix(F,n),z,pos,pos)
         push!(listgens,z)
      end
      if l[1]>1
         for i in 1:l[1]-1
         for j in 1:degree(F)
            z = identity_matrix(F,l[1])
            for k in 1:l[1]-i z[k,i+k]=_lambda^j end
            z = insert_block(identity_matrix(F,n),z,pos,pos)
            push!(listgens,z)
         end
         end
      end
      pos += l[1]*l[2]
   end

   # external diagonal blocks
   pos=1
   for i in 1:length(L)-1
      pos += L[i][1]*L[i][2]
      # block above diagonal
      z = identity_matrix(F,n)
      for j in 1:L[i][1] z[pos-L[i][1]+j-1,pos+L[i+1][1]-L[i][1]+j-1]=1 end
      push!(listgens,z)
      # block below diagonal
      z = identity_matrix(F,n)
      for j in 1:L[i][1] z[pos+j-1,pos-L[i][1]+j-1]=1 end
      push!(listgens,z)
   end

   # cardinality
   res = prod([_GL_order(l[2],F) for l in L])
   exp = fmpz(0)
   for i in 1:length(L)-1
   for j in i+1:length(L)
      exp += L[i][1]*L[i][2]*L[j][2]
   end
   end
   exp *= 2
   exp += sum([(L[i][1]-1)*L[i][2]^2 for i in 1:length(L)])
   res *= order(F)^exp

   return listgens, res
end


# V = vector of integers of the dimensions of Jordan blocks
# return the generators for the centralizers of the unipotent element
# assumes V is sorted (e.g. [1,1,1,2,3,3])
# does the same as above, but every entry is replaced by the corresponding block matrix
function _centr_block_unipotent(f::PolyElem, F::Ring, V::AbstractVector{Int}; isSL=false)
   d = degree(f)
   d>1 || return _centr_unipotent(F,V; isSL=isSL)
   n = sum(V)*d
   C = companion_matrix(f)

   # L = multiset(V)
   L=[[V[1],1]]
   for i in 2:length(V)
      if V[i]==L[length(L)][1]
         L[length(L)][2] += 1
      else
         push!(L, [V[i],1])
      end
   end
   listgens = MatElem[]

   # generators for GL + internal diagonal blocks
   pos=1
   for l in L
      v_g = isSL ? _gens_for_SL_matrix(f,l[2],F; D=l[1]) : _gens_for_GL_matrix(f,l[2],F; D=l[1])
      for x in v_g
         z = insert_block(identity_matrix(F,n),x,pos,pos)
         push!(listgens,z)
      end
      if l[1]>1
         for i in 1:l[1]-1
         c = identity_matrix(F,d)
         for j in 1:degree(F)*d
            z = identity_matrix(F,l[1]*d)
#            for k in 1:l[1]-i z[k,i+k]=_lambda^j end
            for k in 0:l[1]-i-1 insert_block!(z,c,d*k+1 ,d*(k+i)+1) end
            z = insert_block(identity_matrix(F,n),z,pos,pos)
            c *= C           # every time, the block C^(j-1) is inserted
            push!(listgens,z)
         end
         end
      end
      pos += l[1]*l[2]*d
   end

   # external diagonal blocks
   pos=1
   for i in 1:length(L)-1
      pos += L[i][1]*L[i][2]*d
      # block above diagonal
      z = identity_matrix(F,n)
#      for j in 1:L[i][1] z[pos-L[i][1]+j-1,pos+L[i+1][1]-L[i][1]+j-1]=1 end
      for j in 1:L[i][1]*d z[pos-L[i][1]*d+j-1,pos+(L[i+1][1]-L[i][1])*d+j-1]=1 end
      push!(listgens,z)
      # block below diagonal
      z = identity_matrix(F,n)
      for j in 1:L[i][1]*d z[pos+j-1,pos-L[i][1]*d+j-1]=1 end
      push!(listgens,z)
   end

   # cardinality
   res = prod([_GL_order(l[2],order(F)^degree(f)) for l in L])
   exp = fmpz(0)
   for i in 1:length(L)-1
   for j in i+1:length(L)
      exp += L[i][1]*L[i][2]*L[j][2]
   end
   end
   exp *= 2
   exp += sum([(L[i][1]-1)*L[i][2]^2 for i in 1:length(L)])
   exp *= degree(f)
   res *= order(F)^exp

   return listgens, res
end

# returns the list of generators
function _centralizer_GL(x::MatElem)
   _,cbm,ED = generalized_jordan_form(x; with_pol=true)    # cbm = change basis matrix
   n=nrows(x)
   listgens = MatElem[]
   res = fmpz(1)

   i=1
   pos=1
   f=ED[1][1]
   V=[ED[1][2]]
   while i <= length(ED)
      if i<length(ED) && ED[i+1][1]==f
         i+=1
         push!(V,ED[i][2])
      else
         L = _centr_block_unipotent(f,base_ring(x),V)
         for z in L[1]
            push!(listgens, insert_block(identity_matrix(base_ring(x),n),z,pos,pos))
         end
         res *= L[2]
         pos += degree(f)*sum(V)
         i+=1
         if i<=length(ED)
            f = ED[i][1]
            V = [ED[i][2]]
         end
      end
   end

   return listgens, res, cbm
end





########################################################################
#
# Generators for SL 
#
########################################################################


# returns as matrices
function _gens_for_SL(n::Int, F::Ring)
#   n !=1 || return [matrix(F,1,1,[one(F)])]
   n != 1 || return []
   if order(F)==2 || order(F)==3
      h1 = identity_matrix(F,n)
      h1[1,2] = 1
      h2 = zero_matrix(F,n,n)
      h2[1,n] = 1
      for i in 1:n-1 h2[i+1,i] = -1 end
      return h1,h2
   else
      h1 = identity_matrix(F,n)
      h1[1,1] = primitive_element(F)
      h1[2,2] = inv(h1[1,1])
      h2 = zero_matrix(F,n,n)
      h2[1,1] = -1
      h2[1,n] = 1
      for i in 1:n-1 h2[i+1,i] = -1 end
      return h1,h2
   end
end

# returns as matrices
# does the matrix above with F = F[x]/(f), but every entry is replaced by a diagonal join of D corresponding blocks
# ASSUMPTION: deg(f) > 1

function _gens_for_SL_matrix(f::PolyElem, n::Int, F::Ring; D=1)
   C = companion_matrix(f)
   CP = evaluate(_centralizer(f),C)            # matrix of maximal order in the centralizer of the companion matrix
   CPi = inv(CP)

   #if n==1 return [identity_matrix(F,n*D*degree(f))] end
   n != 1 || return []
   h1 = identity_matrix(F,n*degree(f)*D)
   insert_block!(h1,diagonal_join([CP for i in 1:D]),1,1)
   insert_block!(h1,diagonal_join([CPi for i in 1:D]),D*degree(f)+1,D*degree(f)+1)
   h2 = zero_matrix(F,n*degree(f)*D,n*degree(f)*D)
   for i in 1:(n-1)*degree(f)*D h2[i+degree(f)*D,i]=-1 end
   for i in 1:degree(f)*D
      h2[i,i]=-1
      h2[i,i+(n-1)*degree(f)*D]=1
   end
# TODO waiting for a better solution, the h3 generator is necessary because h1,h2 generate just the subgroup of SL(n*deg(f),F) isomorphic to SL(n, F^deg(f))
   h3 = identity_matrix(F,n*degree(f)*D)
   insert_block!(h3,diagonal_join([CP^(order(F)-1) for i in 1:D]),1,1)
   return h1,h2,h3      
end

# returns the list of generators
function _centralizer_SL(x::MatElem)
   _,cbm,ED = generalized_jordan_form(x; with_pol=true)    # cbm = change basis matrix
   n=nrows(x)
   listgens = MatElem[]
   _lambda = primitive_element(base_ring(x))
   res = fmpz(1)
   ind = fmpz(0)

   i=1
   pos=1
   f=ED[1][1]
   V=[ED[1][2]]
   c = evaluate(_centralizer(f), companion_matrix(f))
   c = c^(_disc_log(det(c),_lambda))  # TODO this _disc_log is bad. Don't try with large fields. 
   block_dim = [[ED[1][2],1,c]]   # list of [d,m,f]  d = dimension of the Jordan block, m = its multiplicity, c = el of max order determinant in centralizer of the companion matrix
   while i <= length(ED)
      if i<length(ED) && ED[i+1][1]==f
         i+=1
         push!(V,ED[i][2])
         if ED[i][2]==block_dim[length(block_dim)][1] block_dim[length(block_dim)][2] += 1
         else push!(block_dim, [ED[i][2], 1, c])
         end
      else
         L = _centr_block_unipotent(f,base_ring(x),V; isSL=true)
         for z in L[1]
            push!(listgens, insert_block(identity_matrix(base_ring(x),n),z,pos,pos))
         end
         ind = gcd(ind, gcd([b[1] for b in block_dim]))
         res *= L[2]
         pos += degree(f)*sum(V)
         i+=1
         if i<=length(ED)
            f = ED[i][1]
            V = [ED[i][2]]
            c = evaluate(_centralizer(f), companion_matrix(f))
            push!(block_dim, [ED[i][2],1,c])
         end
      end
   end

   # start general blocks, those which have det=1 globally, but not on every single el.div.
   Ga = abelian_group([Int(order(base_ring(x)))-1 for i in 1:length(block_dim)])
   Gb = abelian_group(Int(order(base_ring(x)))-1)
   f = hom(Ga,Gb, [Gb[1]*block_dim[i][1] for i in 1:length(block_dim)])
   K,g = kernel(f)
   for k in gens(K)
   if !iszero(g(k))               #TODO just to do not put identity matrices in the list of generators
      z = identity_matrix(base_ring(x),n)
      pos = 1
      for i in 1:length(block_dim)
         insert_block!(z,diagonal_join([block_dim[i][3]^Int(g(k)[i]) for j in 1:block_dim[i][1]]),pos,pos)
         pos += block_dim[i][1]*block_dim[i][2]*nrows(block_dim[i][3])
      end
      push!(listgens, z)
   end
   end

   ind = gcd(ind, order(base_ring(x))-1)
   res = div(res, order(base_ring(x))-1)
   res *= ind
   return listgens, res, cbm
end



#=    TODO  all of this does not work yet


# return subgroup of GL(n,F) of index d
# returns as matrices
function _gens_for_sub_GL(n::Int, F::Ring, d::Int)
   @assert mod(size(F)-1,d)==0 "Index must divide q-1"
   if n==1 return [matrix(F,1,1,[primitive_element(F)^d])] end
   if order(F)==2
      h1 = identity_matrix(F,n)
      h1[1,2] = 1
      h2 = zero_matrix(F,n,n)
      for i in 1:n-1 h2[i+1,i]=1 end
      h2[1,n] = 1
   else
      h1 = identity_matrix(F,n)
      h1[1,1] = primitive_element(F)
      h1[2,2] = h1[1,1]^(size(F)-d-2)
      h2 = zero_matrix(F,n,n)
      for i in 1:n-1 h2[i+1,i]=-1 end
      h2[1,1] = -1
      h2[1,n] = 1
   end
   return h1,h2      
end

_gens_for_SL(n::Int, F::Ring) = _gens_for_sub_GL(n, F, Int(size(F))-1)

# returns as matrices
# return subgroup of GL(n,F) of index d
# does the matrix above with F = F[x]/(f), but every entry is replaced by a diagonal join of D corresponding blocks
# ASSUMPTION: deg(f) > 1
function _gens_for_sub_GL_matrix(f::PolyElem, n::Int, F::Ring, d::Int; D=1)
   q = size(F)^degree(f)
   @assert mod(q-1,d)==0 "Index must divide q-1"
   C = companion_matrix(f)
   CP = evaluate(_centralizer(f),C)            # matrix of maximal order in the centralizer of the companion matrix

   if n==1 return [diagonal_join([CP^d for i in 1:D])] end
   h1 = identity_matrix(F,n*degree(f)*D)
   insert_block!(h1,diagonal_join([CP for i in 1:D]),1,1)
print(q-2-d)
   CP=CP^(q-2-d)
   insert_block!(h1,diagonal_join([CP for i in 1:D]),D*degree(f)+1,D*degree(f)+1)
   h2 = zero_matrix(F,n*degree(f)*D,n*degree(f)*D)
   for i in 1:(n-1)*degree(f)*D h2[i+degree(f)*D,i]=-1 end
   for i in 1:degree(f)*D
      h2[i,i]=-1
      h2[i,i+(n-1)*degree(f)*D]=1
   end
   return h1,h2      
end

=#


########################################################################
#
# User level functions
#
########################################################################

pol_elementary_divisors(x::MatrixGroupElem) = pol_elementary_divisors(x.elm)

function centralizer(G::MatrixGroup, x::MatrixGroupElem)
   if isdefined(G,:descr) && (G.descr==:GL || G.descr==:SL)
      V,card,a = G.descr==:GL ? _centralizer_GL(x.elm) : _centralizer_SL(x.elm)
      am = inv(a)
      L = [G(am*v*a) for v in V]
      H = MatrixGroup(G.deg, G.ring, L)
      H.order = card
      return H, Nothing          # do not return the embedding of the centralizer into G to do not compute G.X
   end
   C = GAP.Globals.Centralizer(G.X, x.X)
   return _as_subgroup(G, C)
end
