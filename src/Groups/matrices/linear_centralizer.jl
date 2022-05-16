# TODO: in this file are used many methods TEMPORARILY defined in files matrix_manipulation.jl and stuff_field_gen.jl
# once methods in those files will be deleted / replaced / modified, this file need to be modified too


########################################################################
#
# Orders
#
########################################################################

function _GL_order(n::Int, q::fmpz)
   res = q^div(n*(n-1),2)
   for i in 1:n
      res *= (q^i-1)
   end
   return res
end

_GL_order(n::Int, F::Ring) = _GL_order(n, order(F))

function _SL_order(n::Int, q::fmpz)
   res = q^div(n*(n-1),2)
   for i in 2:n
     res *= (q^i-1)
   end
   return res
end

_SL_order(n::Int, F::Ring) = _SL_order(n, order(F))

########################################################################
#
# Centralizer in GL
#
########################################################################

# Compute a generating set for GL(n,q) and SL(n,q)
# generators for GL(n,q) and SL(n,q) are described in [Tay87]

# returns elements of type MatElem{T}
# returns a generating set for GL(n,F) of at most two elements
function _gens_for_GL(n::Int, F::FinField)
   n !=1 || return [matrix(F,1,1,[primitive_element(F)])]
   h1 = identity_matrix(F,n)
   h2 = zero_matrix(F,n,n)
   if order(F)==2
      h1[1,2] = 1
   else
      h1[1,1] = primitive_element(F)
      h2[1,1] = -1
   end
   h2[1,n] = 1
   for i in 1:n-1
      h2[i+1,i] = -1
   end
   return [h1,h2]
end

# returns elements of type MatElem{T}
# does the same as above with the following changes:
# the field F is replaced by F = F[x]/(f);
# every entry y is replaced by a diagonal join of D copies of phi(y),
# where phi: F -> matrix_algebra(F,degree(f))
# is the ring homomorphism sending a fixed root of f into the companion matrix of f.
# (hence, elements of the final output belong to matrix_algebra(F,n*D*degree(f)))
# ASSUMPTION: deg(f) > 1
function _gens_for_GL_matrix(f::PolyElem, n::Int, F::FinField; D::Int=1)
   C = companion_matrix(f)
   CP = _centralizer(f)(C)            # matrix of maximal order in the centralizer of the companion matrix
   Df = degree(f)*D

   if n==1 return [cat([CP for i in 1:D]..., dims=(1,2))] end
   h1 = identity_matrix(F,n*Df)
   h1[1:Df,1:Df] = cat([CP for i in 1:D]..., dims=(1,2))
   h2 = zero_matrix(F,n*degree(f)*D,n*Df)
   for i in 1:(n-1)*Df
      h2[i+Df,i]=-1
   end
   for i in 1:Df
      h2[i,i]=-1
      h2[i,i+(n-1)*Df]=1
   end
   return [h1,h2]
end

# generators for centralizer of unipotent elements are described in:
# Giovanni De Franceschi, Centralizers and conjugacy classes in finite classical groups, arXiv:2008.12651

# V = vector of integers of the dimensions of Jordan blocks
# return the generators for the centralizers of the diagonal join of unipotent Jordan blocks of dimensions V
# assumes V is sorted (e.g. [1,1,1,2,3,3])
function _centr_unipotent(F::FinField, V::AbstractVector{Int}; isSL=false)
   n = sum(V)
   _lambda = gen(F)

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
      idN = identity_matrix(F,n)
      v_g = isSL ? _gens_for_SL(l[2],F) : _gens_for_GL(l[2],F)
      for x in v_g
         z = hvcat(l[2], [x[i,j]*identity_matrix(F,l[1]) for j in 1:l[2] for i in 1:l[2]]...)
         idN[pos:pos+l[1]*l[2]-1,pos:pos+l[1]*l[2]-1] = z
         push!(listgens,idN)
      end
      idN = identity_matrix(F,n)
      if l[1]>1
         for i in 1:l[1]-1, j in 1:degree(F)
            z = identity_matrix(F,l[1])
            for k in 1:l[1]-i z[k,i+k]=_lambda^j end
            idN[pos:pos+l[1]-1, pos:pos+l[1]-1] = z
            push!(listgens,idN)
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
      for j in 1:L[i][1]
         idx = pos-L[i][1]+j-1
         z[idx,idx+L[i+1][1]] = 1
      end
      push!(listgens,z)
      # block below diagonal
      z = identity_matrix(F,n)
      for j in 1:L[i][1]
         idx = pos+j-1
         z[idx,idx-L[i][1]]=1
      end
      push!(listgens,z)
   end

   # cardinality
   res = prod([_GL_order(l[2],F) for l in L])
   exp = fmpz(0)
   for i in 1:length(L)-1, j in i+1:length(L)
      exp += L[i][1]*L[i][2]*L[j][2]
   end
   exp *= 2
   exp += sum([(L[i][1]-1)*L[i][2]^2 for i in 1:length(L)])
   res *= order(F)^exp

   return listgens, res
end


# V = vector of integers of the dimensions of Jordan blocks
# assumes V is sorted (e.g. [1,1,1,2,3,3])
# does the same as above, but every entry is replaced by the corresponding block matrix in the field F/(f)
function _centr_block_unipotent(f::PolyElem, F::FinField, V::AbstractVector{Int}; isSL=false)
   d = degree(f)
   d>1 || return _centr_unipotent(F,V; isSL=isSL)
   n = sum(V)*d
   C = companion_matrix(f)
   idN = identity_matrix(F,n)
   idD = identity_matrix(F,d)

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
      idN = identity_matrix(F,n)
      for x in v_g
         idN[pos:pos+nrows(x)-1,pos:pos+ncols(x)-1] = x
         push!(listgens,idN)
      end
      if l[1]>1
         idN = identity_matrix(F,n)
         for i in 1:l[1]-1
            c = idD
            for j in 1:degree(F)*d
               z = identity_matrix(F,l[1]*d)
               for k in 0:l[1]-i-1
                  z[d*k+1:d*(k+1),d*(k+i)+1:d*(k+i+1)] = c
               end
               idN[pos:pos+l[1]*d-1,pos:pos+l[1]*d-1] = z
               c *= C           # every time, the block C^(j-1) is inserted
               push!(listgens,idN)
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
      z = idN
      for j in 1:L[i][1]*d
         idx = pos-L[i][1]*d+j-1
         z[idx,idx+L[i+1][1]*d]=1
      end
      push!(listgens,z)
      # block below diagonal
      z = idN
      for j in 1:L[i][1]*d
         idx = pos+j-1
         z[idx,idx-L[i][1]*d]=1
      end
      push!(listgens,z)
   end

   # cardinality
   res = prod([_GL_order(l[2],order(F)^degree(f)) for l in L])
   exp = fmpz(0)
   for i in 1:length(L)-1, j in i+1:length(L)
      exp += L[i][1]*L[i][2]*L[j][2]
   end
   exp *= 2
   exp += sum([(L[i][1]-1)*L[i][2]^2 for i in 1:length(L)])
   exp *= degree(f)
   res *= order(F)^exp

   return listgens, res
end

# returns the list of generators and the cardinality of the centralizer of x in GL
function _centralizer_GL(x::MatElem)
   _,a,ED = generalized_jordan_form(x; with_pol=true)    # a = change basis matrix
   am = inv(a)
   n=nrows(x)
   listgens = MatElem[]
   res = fmpz(1)
   idN = identity_matrix(base_ring(x),n)

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
            idN = identity_matrix(base_ring(x),n)
            idN[pos:pos+nrows(z)-1, pos:pos+ncols(z)-1] = z
            push!(listgens, am*idN*a)
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

   return listgens, res
end





########################################################################
#
# Generators for SL
#
########################################################################

# generators for GL(n,q) and SL(n,q) are described in [Tay87]

# returns elements of type MatElem{T}
# returns a generating set for SL(n,F) of at most two elements
function _gens_for_SL(n::Int, F::FinField)
   n != 1 || return dense_matrix_type(F)[]
   h1 = identity_matrix(F,n)
   h2 = zero_matrix(F,n,n)
   if order(F)==2 || order(F)==3
      h1[1,2] = 1
   else
      h1[1,1] = primitive_element(F)
      h1[2,2] = inv(h1[1,1])
      h2[1,1] = -1
   end
   h2[1,n] = 1
   for i in 1:n-1
      h2[i+1,i] = -1
   end
   return [h1,h2]
end

# returns elements of type MatElem{T}
# does the same as above with the following changes:
# the field F is replaced by F = F[x]/(f);
# every entry y is replaced by a diagonal join of D copies of phi(y),
# where phi: F -> matrix_algebra(F,degree(f))
# is the ring homomorphism sending a fixed root of f into the companion matrix of f.
# (hence, elements of the final output belong to matrix_algebra(F,n*D*degree(f)))
# ASSUMPTION: deg(f) > 1
function _gens_for_SL_matrix(f::PolyElem, n::Int, F::FinField; D::Int=1)
   C = companion_matrix(f)
   CP = _centralizer(f)(C)            # matrix of maximal order in the centralizer of the companion matrix
   CPi = inv(CP)
   Df = D*degree(f)

   n != 1 || return dense_matrix_type(F)[]
   h1 = identity_matrix(F,n*Df)
   h1[1:Df, 1:Df] = cat([CP for i in 1:D]..., dims=(1,2))
   h1[Df+1:2*Df, Df+1:2*Df] = cat([CPi for i in 1:D]..., dims=(1,2))
   h2 = zero_matrix(F,n*Df,n*Df)
   for i in 1:(n-1)*Df
      h2[i+Df,i]=-1
   end
   for i in 1:Df
      h2[i,i]=-1
      h2[i,i+(n-1)*Df]=1
   end
   # TODO h1,h2 generate just the subgroup K of SL(n*deg(f),F) isomorphic to SL(n, F^deg(f))
   # hence we need to add an element h3 of maximal order in SL(n*deg(f),F) \ K .
   # TODO: if in future we find out how to generate intermediate groups between GL and SL with only 2 elements,
   # then we can reduce the number of generators here from 3 to 2
   h3 = identity_matrix(F,n*Df)
   h3[1:Df, 1:Df] = cat([CP^(order(F)-1) for i in 1:D]..., dims=(1,2))
   return [h1,h2,h3]
end


# generators for centralizer of unipotent elements are described in:
# Giovanni De Franceschi, Centralizers and conjugacy classes in finite classical groups, arXiv:2008.12651

# returns the list of generators and the cardinality of the centralizer of x in SL
function _centralizer_SL(x::MatElem)
   _,a,ED = generalized_jordan_form(x; with_pol=true)    # a = change basis matrix
   am = inv(a)
   n=nrows(x)
   listgens = MatElem[]
   _lambda = primitive_element(base_ring(x))
   res = fmpz(1)
   ind = fmpz(0)
   idN = identity_matrix(base_ring(x),n)

   i=1
   pos=1
   f=ED[1][1]
   V=[ED[1][2]]
   c = _centralizer(f)(companion_matrix(f))
   # Computes d such that det(c)^d = lambda and replace c by c^d.
   # The result c^d is a primitive root for the base ring of x such that det(c^d) = lambda
   c = c^(disc_log(det(c),_lambda))
   # block_dim is a list of triples [d,m,c], where
   # d = dimension of the Jordan block, m = its multiplicity, c = el of max order determinant in centralizer of the companion matrix
   block_dim = [[ED[1][2],1,c]]
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
            temp = deepcopy(idN)
            temp[pos:pos-1+nrows(z), pos:pos-1+ncols(z)] = z
            push!(listgens, am*temp*a)
         end
         ind = gcd(ind, gcd([b[1] for b in block_dim]))
         res *= L[2]
         pos += degree(f)*sum(V)
         i+=1
         if i<=length(ED)
            f = ED[i][1]
            V = [ED[i][2]]
            c = _centralizer(f)(companion_matrix(f))
            # Computes d such that det(c)^d = lambda and replace c by c^d.
            # The result c^d is a primitive root for the base ring of x such that det(c^d) = lambda
            c = c^(disc_log(det(c),_lambda))
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
      z = deepcopy(idN)
      pos = 1
      for i in 1:length(block_dim)
         t1 = cat([block_dim[i][3]^Int(g(k)[i]) for j in 1:block_dim[i][1]]..., dims=(1,2))
         z[pos:pos-1+nrows(t1), pos:pos-1+ncols(t1)] = t1
         t1 = cat([block_dim[i][3]^Int(g(k)[i]) for j in 1:block_dim[i][1]]..., dims=(1,2))
         z[pos:pos-1+nrows(t1), pos:pos-1+ncols(t1)] = t1
         pos += block_dim[i][1]*block_dim[i][2]*nrows(block_dim[i][3])
      end
      push!(listgens, am*z*a)
   end
   end

   ind = gcd(ind, order(base_ring(x))-1)
   res = div(res, order(base_ring(x))-1)
   res *= ind
   return listgens, res
end




########################################################################
#
# User level functions
#
########################################################################

"""
    centralizer(G::MatrixGroup{T}, x::MatrixGroupElem{T})

Return (`C`,`f`), where `C` is the centralizer of `x` in `C` and `f` is the embedding of `C` into `G`.
If `G` = `GL(n,F)` or `SL(n,F)`, then `f` = `nothing`. In this case, to get the embedding homomorphism of `C` into `G`, use
> `issubgroup(G,C)[2]`
"""
function centralizer(G::MatrixGroup{T}, x::MatrixGroupElem{T}) where T <: FinFieldElem
   if isdefined(G,:descr) && (G.descr==:GL || G.descr==:SL)
      V,card = G.descr==:GL ? _centralizer_GL(x.elm) : _centralizer_SL(x.elm)
      H = MatrixGroup(G.deg, G.ring, V)
      set_attribute!(H, :order => fmpz(card))
      return H, nothing          # do not return the embedding of the centralizer into G to do not compute G.X
   end
   C = GAP.Globals.Centralizer(G.X, x.X)
   return _as_subgroup(G, C)
end
