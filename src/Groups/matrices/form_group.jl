# TODO: in this file several ways to get the preserved forms by a matrix group are implemented
# we can choose to keep all of them or only some of them

export
    invariant_alternating_forms,
    invariant_bilinear_forms,
    invariant_hermitian_forms,
    invariant_quadratic_forms,
    invariant_sesquilinear_forms,
    invariant_symmetric_forms,
    isometry_group,
    preserved_quadratic_forms,
    preserved_sesquilinear_forms

########################################################################
#
# From group to form #TODO: there are different approaches. Which is the best?
#
########################################################################

# TODO 1st approach: brute force calculation
# Algorithm furnished by Thomas Breuer, Aachen University
# extended to quadratic by Giovanni De Franceschi, TU Kaiserslautern
# WARNING: big linear system !!

# METHOD: given a group G, the condition gBg*=B is a linear condition on the coefficients of B,
# hence we write down such system of dimension n^2 (n*(n-1)/2 for quadratic forms)

"""
    invariant_bilinear_forms(G::MatrixGroup)

Return a generating set for the vector spaces of bilinear forms preserved by the group `G`.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_bilinear_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
   n = degree(G)
   M = T[]
   for mat in gens(G)
      mmat = mat^-1
      MM = zero_matrix(F,n^2,n^2)
      for i in 1:n, j in 1:n, k in 1:n
         MM[n*(i-1)+j,n*(k-1)+j] += mat[i,k]
         MM[n*(i-1)+j,n*(i-1)+k] -= mmat[j,k]
      end
      push!(M,MM)
   end

   r,K = nullspace(block_matrix(length(M),1,M))
   
   return [matrix(F,n,n,[K[i,j] for i in 1:n^2]) for j in 1:r]
end

"""
    invariant_sesquilinear_forms(G::MatrixGroup)

Return a generating set for the vector spaces of sesquilinear non-bilinear forms preserved by the group `G`.
It works only if `base_ring(G)` is a finite field with even degree on its prime subfield.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_sesquilinear_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
   @assert typeof(F)<:FinField "At the moment, only finite fields are considered"
   @assert iseven(degree(F)) "Base ring has no even degree"
   n = degree(G)
   M = T[]
   for mat in gens(G)
      mmat = map(y->frobenius(y,div(degree(F),2)),(mat.elm)^-1)
      MM = zero_matrix(F,n^2,n^2)
      for i in 1:n, j in 1:n, k in 1:n
         MM[n*(i-1)+j,n*(k-1)+j] += mat[i,k]
         MM[n*(i-1)+j,n*(i-1)+k] -= mmat[j,k]
      end
      push!(M,MM)
   end

   r,K = nullspace(block_matrix(length(M),1,M))
   
   return [matrix(F,n,n,[K[i,j] for i in 1:n^2]) for j in 1:r]
end

"""
    invariant_quadratic_forms(G::MatrixGroup)

Return a generating set for the vector spaces of quadratic forms preserved by the group `G`.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_quadratic_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
   n = degree(G)
   M = T[]
   for mat in gens(G)
      MM = zero_matrix(F,div(n*(n+1),2),div(n*(n+1),2))
      for i in 1:n
         row = div((2*n-i+2)*(i-1),2)+1
         for p in 1:n, q in p:n
            col = div((2*n-p+2)*(p-1),2)+q-p+1
            MM[row, col] = mat[i,p]*mat[i,q]
            for j in i+1:n
               MM[row+j-i, col] = mat[i,p]*mat[j,q]+mat[j,p]*mat[i,q]
            end
         end
      end
#      for i in 1:div(n*(n+1),2) MM[i,i]-=1 end
      MM -= one(MM)
      push!(M,MM)
   end

   r,K = nullspace(block_matrix(length(M),1,M))
   M = T[]
   for i in 1:r
      push!(M,upper_triangular_matrix(K[1:div(n*(n+1),2),i]))
   end
   return M
end

#TODO: do we want to keep these?

# METHOD: if B = (b_ij) is the solution matrix, then b_ji = b_ij;
# hence, the dimension of the linear system can be reduced to n(n+1)/2
"""
    invariant_symmetric_forms(G::MatrixGroup)

Return a generating set for the vector spaces of symmetric forms preserved by the group `G`.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
!!! warning "Note:"
    Work properly only in odd characteristic. In even characteristic, only alternating forms are found.
"""
invariant_symmetric_forms(G::MatrixGroup{S,T}) where {S,T} = T[x + transpose(x) for x in invariant_quadratic_forms(G)]

# METHOD: if B = (b_ij) is the solution matrix, then b_ji = -b_ij and b_ii=0;
# hence, the dimension of the linear system can be reduced to n(n-1)/2
"""
    invariant_alternating_forms(G::MatrixGroup)

Return a generating set for the vector spaces of alternating forms preserved by the group `G`.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_alternating_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
   n = degree(G)
   M = T[]

   for mat in gens(G)
      MM = zero_matrix(F,div(n*(n-1),2),div(n*(n-1),2))
      idx_r=1
      for i in 1:n, j in i+1:n
         idx_c=1
         for s in 1:n
            for t in (s+1):n
               MM[idx_r,idx_c] = mat[i,s]*mat[j,t]-mat[i,t]*mat[j,s]
               idx_c+=1
            end
         end
         idx_r+=1
      end
      MM -= one(MM)
      push!(M,MM)
   end

   r,K = nullspace(block_matrix(length(M),1,M))

   L = T[]
   for j in 1:r
      B = upper_triangular_matrix(K[1:div(n*(n-1),2),j])
      B = insert_block(zero_matrix(F,n,n),B,1,2)
#=      for i in 1:n, j in 1:i-1
         B[i,j]=-B[j,i]
      end  =#
      B -= transpose(B)
      push!(L, B)
   end
   return L
end

# METHOD: if B = (c_ij) is the solution matrix, then c_ij = (c_ji)^q (where base_ring(G) = GF(q^2))
# Let F0 be the subfield of F s.t. [F:F0] = 2 and w in F \ F0 fixed, then c_ij = x_ij +w*y_ij for x_ij, y_ij in F0
# the condition c_ij = (c_ij)^q + condition for B to be a form preserved by G
# is a F0-linear condition on the x_ij and y_ij.
# Hence, we write down the F0-linear system in the x_ij and y_ij (dimension = n*(n+1))
# NOTE: two different approaches for an appropriate w are employed for odd and even characteristic
"""
    invariant_hermitian_forms(G::MatrixGroup)

Return a generating set for the vector spaces of hermitian forms preserved by the group `G`.
It works only if `base_ring(G)` is a finite field with even degree on its prime subfield.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_hermitian_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
   @assert typeof(F)<:FinField "At the moment, only finite fields are considered"
   n = degree(G)
   M = T[]

   p = characteristic(F)
   d = degree(F)
   iseven(d) || return M

   F0 = GF(Int(p), div(d,2))[1]
   q = p^div(d,2)
   e = embed(F0,F)
   em = preimage_map(F0,F)

   if isodd(q)
      w = gen(F) - (F(2)^-1)*(gen(F)+gen(F)^q)        # w in F \ F0 s.t. w+w^q = 0
      Nw = em(w*w^q)

      # returns a,b in F0 such that x = a+bw
      coeffq(x) = F(2)^-1*(x+x^q), F(2)^-1*(x-x^q)*w^-1

      for mat in gens(G)
         A = zero_matrix(F0,n,n)
         B = zero_matrix(F0,n,n)
         for i in 1:n, j in 1:n           # mat = A + wB,  A,B in GL(n,F0)
            a,b = coeffq(mat[i,j])
            A[i,j] = em(a)
            B[i,j] = em(b)
         end
         MM = zero_matrix(F0, n^2, n^2)
         pos_r = 1
         pos_c = 1

         # the solution matrix has coefficients x_ij+w*y_ij
         # where x_ji+w*y_ji = (x_ij+w*y_ij)^q = x_ij-w*y_ij (so, we don't need i>j)
         # in the vector of solutions, the coefficients are ordered as:
         # first the x_ii, then the x_ij (i<j) and finally the y_ij (i<j)

         # coefficients for the x_ii
         for i in 1:n
            for h in 1:n
               MM[pos_r,pos_c] = A[i,h]^2+Nw*B[i,h]^2
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = 2*A[i,h]*A[i,k] + Nw*2*B[i,h]*B[i,k]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = Nw*2*(A[i,h]*B[i,k]-A[i,k]*B[i,h])
               pos_c +=1
            end
            pos_r +=1
            pos_c =1
         end

         # coefficients for the x_ij
         for i in 1:n, j in i+1:n
            for h in 1:n
               MM[pos_r,pos_c] = A[i,h]*A[j,h]+Nw*B[i,h]*B[j,h]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[i,h]*A[j,k]+A[i,k]*A[j,h] + Nw*(B[i,h]*B[j,k]+B[i,k]*B[j,h])
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = Nw*(A[i,h]*B[j,k]-A[j,k]*B[i,h]+A[j,h]*B[i,k]-A[i,k]*B[j,h])
               pos_c +=1
            end
            pos_r +=1
            pos_c =1
         end

         # coefficients for the y_ij
         for i in 1:n, j in i+1:n
            for h in 1:n
               MM[pos_r,pos_c] = A[j,h]*B[i,h]-A[i,h]*B[j,h]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[j,k]*B[i,h]-A[i,h]*B[j,k]+A[j,h]*B[i,k]-A[i,k]*B[j,h]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[i,h]*A[j,k]-A[i,k]*A[j,h]+Nw*(B[i,h]*B[j,k]-B[i,k]*B[j,h])
               pos_c +=1
            end
            pos_r +=1
            pos_c =1
         end

#=         for i in 1:n^2
            MM[i,i] -=1
         end =#
         MM -= one(MM)
         push!(M,MM)
      end

   else

      w = gen(F) * (gen(F)+gen(F)^q)^-1           # w in F \ F0 s.t. w+w^q=1
      Nw = em(w*w^q)

      # returns a,b in F0 such that x = a+bw
      coeff2(x) = x+w*(x+x^q), x+x^q

      for mat in gens(G)
         A = zero_matrix(F0,n,n)
         B = zero_matrix(F0,n,n)
         for i in 1:n, j in 1:n           # mat = A + wB,  A,B in GL(n,F0)
            a,b = coeff2(mat[i,j])
            A[i,j] = em(a)
            B[i,j] = em(b)
         end
         MM = zero_matrix(F0, n^2, n^2)
         pos_r = 1
         pos_c = 1

         # the solution matrix has coefficients x_ij+w*y_ij
         # where x_ji+w*y_ji = (x_ij+w*y_ij)^q = (x_ij+y_ij)+w*y_ij (so, we don't need i>j)
         # in the vector of solutions, the coefficients are ordered as:
         # first the x_ii, then the x_ij (i<j) and finally the y_ij (i<j)

         # coefficients for the x_ii
         for i in 1:n
            for h in 1:n
               MM[pos_r,pos_c] = A[i,h]^2+A[i,h]*B[i,h]+Nw*B[i,h]^2
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[i,k]*B[i,h]+A[i,h]*B[i,k]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[i,h]*A[i,k]+A[i,k]*B[i,h]+Nw*B[i,k]*B[i,h]
               pos_c +=1
            end
            pos_r +=1
            pos_c =1
         end

         # coefficients for the x_ij
         for i in 1:n, j in i+1:n
            for h in 1:n
               MM[pos_r,pos_c] = A[i,h]*A[j,h]+A[i,h]*B[j,h]+Nw*B[i,h]*B[j,h]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[i,h]*A[j,k]+A[i,h]*B[j,k]+A[i,k]*B[j,h]+A[i,k]*A[j,h]+Nw*(B[i,h]*B[j,k]+B[i,k]*B[j,h])
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = Nw*(A[i,h]*B[j,k]+A[j,k]*B[i,h]+B[i,k]*B[j,h]+A[i,k]*B[j,h]+A[j,h]*B[i,k])+A[i,k]*A[j,h]+A[i,k]*B[j,h]
               pos_c +=1
            end
            pos_r +=1
            pos_c =1
         end

         # coefficients for the y_ij
         for i in 1:n, j in i+1:n
            for h in 1:n
               MM[pos_r,pos_c] = A[j,h]*B[i,h]+A[i,h]*B[j,h]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[j,k]*B[i,h]+A[i,h]*B[j,k]+A[j,h]*B[i,k]+A[i,k]*B[j,h]
               pos_c +=1
            end
            for h in 1:n, k in h+1:n
               MM[pos_r,pos_c] = A[i,h]*A[j,k]+A[i,k]*A[j,h]+A[j,k]*B[i,h]+A[i,k]*B[j,h]+Nw*(B[i,h]*B[j,k]+B[i,k]*B[j,h])
               pos_c +=1
            end
            pos_r +=1
            pos_c =1
         end
 
         #= for i in 1:n^2
            MM[i,i] -=1
         end   =#
         MM -= one(MM)
         push!(M,MM)
      end
   end

   n_gens,K = nullspace(block_matrix(length(M),1,M))
   N = div(n*(n-1),2)
   L = T[]

   for r in 1:n_gens
      sol = zero_matrix(F,n,n)
      for i in 1:n
         sol[i,i] = e(K[i,r])
      end
      pos = n+1
      if isodd(q)
         for i in 1:n, j in i+1:n
            sol[i,j] = e(K[pos,r])+w*e(K[pos+N,r])       # sol[i,j] = a_ij + w*b_ij
            sol[j,i] = e(K[pos,r])-w*e(K[pos+N,r])       # sol[i,j] = sol[j,i]^q
            pos +=1
         end
      else
         for i in 1:n, j in i+1:n
            sol[i,j] = e(K[pos,r])+w*e(K[pos+N,r])        # sol[i,j] = a_ij + w*b_ij
            sol[j,i] = e(K[pos,r])+(1+w)*e(K[pos+N,r])    # sol[i,j] = sol[j,i]^q
            pos +=1
         end
      end
      push!(L,sol)
   end

   return L
end

# TODO 2nd approach: using MeatAxe GAP functionalities
# TODO: these are not exported at the moment

"""
    function invariant_bilinear_form(G::MatrixGroup)

Return an invariant bilinear form for the group `G`.
It works only if the module induced by the action of `G` is absolutely irreducible.
!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.
"""
function invariant_bilinear_form(G::MatrixGroup)
   V = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(G.X), G.mat_iso.fr.codomain)
   B = GAP.Globals.MTX.InvariantBilinearForm(V)
   return G.mat_iso(B)   
end

"""
    function invariant_sesquilinear_form(G::MatrixGroup)

Return an invariant sesquilinear (non bilinear) form for the group `G`.
It works only if the module induced by the action of `G` is absolutely irreducible.
!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.
"""
function invariant_sesquilinear_form(G::MatrixGroup)
   V = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(G.X), G.mat_iso.fr.codomain)
   B = GAP.Globals.MTX.InvariantSesquilinearForm(V)
   return G.mat_iso(B)   
end

"""
    function invariant_quadratic_form(G::MatrixGroup)

Return an invariant bilinear form for the group `G`.
It works only if the module induced by the action of `G` is absolutely irreducible.
!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.
"""
function invariant_quadratic_form(G::MatrixGroup)
   if iseven(characteristic(base_ring(G)))
      V = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(G.X), G.mat_iso.fr.codomain)
      B = GAP.Globals.MTX.InvariantQuadraticForm(V)
      return _upper_triangular_version(G.mat_iso(B))
   else
      m = invariant_bilinear_form(G)
      for i in 1:degree(G), j in i+1:degree(G)
         m[i,j]+=m[j,i]
         m[j,i]=0
      end
      return m
   end
end


# TODO 3rd approach: using GAP package "forms"
"""
    preserved_quadratic_forms(G::MatrixGroup)

Uses random methods to find all of the quadratic forms preserved by `G` up to a scalar
(i.e. such that `G` is a group of similarities for the forms). 
Since the procedure relies on a pseudo-random generator, 
the user may need to execute the operation more than once to find all invariant quadratic forms.
"""
function preserved_quadratic_forms(G::MatrixGroup{S,T}) where {S,T}
   L = GAP.Globals.PreservedQuadraticForms(G.X)
   R = SesquilinearForm{S}[]
   for f_gap in L
      f = quadratic_form(G.mat_iso(GAP.Globals.GramMatrix(f_gap)))
      f.X = f_gap
      f.mat_iso = G.mat_iso
      push!(R,f)
   end
   return R
end

"""
    preserved_sesquilinear_forms(G::MatrixGroup)

Uses random methods to find all of the sesquilinear forms preserved by `G` up to a scalar
(i.e. such that `G` is a group of similarities for the forms).
Since the procedure relies on a pseudo-random generator,
the user may need to execute the operation more than once to find all invariant sesquilinear forms.
"""
function preserved_sesquilinear_forms(G::MatrixGroup{S,T}) where {S,T}
   L = GAP.Globals.PreservedSesquilinearForms(G.X)
   R = SesquilinearForm{S}[]
   for f_gap in L
      if GAP.Globals.IsHermitianForm(f_gap)
         f = hermitian_form(G.mat_iso(GAP.Globals.GramMatrix(f_gap)))
      elseif GAP.Globals.IsSymmetricForm(f_gap)
         f = symmetric_form(G.mat_iso(GAP.Globals.GramMatrix(f_gap)))
      elseif GAP.Globals.IsAlternatingForm(f_gap)
         f = alternating_form(G.mat_iso(GAP.Globals.GramMatrix(f_gap)))
      else
         error("Invalid form")
      end
      f.X = f_gap
      f.mat_iso = G.mat_iso
      push!(R,f)
   end
   return R
end


########################################################################
#
# From form to group
#
########################################################################


# return the GAP matrix of the form preserved by the GAP standard group
function _standard_form(descr::Symbol, e::Int, n::Int, q::Int)
   if descr==:quadratic
      return GAP.Globals.InvariantQuadraticForm(GO(e,n,q).X).matrix
   elseif descr==:symmetric #|| descr==:alternating
      return GAP.Globals.InvariantBilinearForm(GO(e,n,q).X).matrix
   elseif descr==:hermitian
      return GAP.Globals.InvariantSesquilinearForm(GU(n,q).X).matrix
   elseif descr==:alternating
      return GAP.Globals.InvariantBilinearForm(Sp(n,q).X).matrix
   else
      error("unsupported description")
   end      
end

function _standard_form(descr::Symbol, e::Int, n::Int, F::Ring)
   q = order(F)
   if descr ==:hermitian q = characteristic(F)^div(degree(F),2) end
   return _standard_form(descr,e,n,Int(q))
end


"""
    isometry_group(f::SesquilinearForm{T})

Return the group of isometries for the sesquilinear form `f`.
"""
function isometry_group(f::SesquilinearForm{T}) where T
   B = gram_matrix(f)
   n = nrows(B)
   F = base_ring(B)
   r=n

   if f.descr==:quadratic
      W,phi = radical(f)
      V = VectorSpace(F,n)
      U,e = complement(V,W)
      A = zero_matrix(F,n,n)
      r = dim(U)
      for i in 1:r, j in 1:n
         A[i,j]=e(gen(U,i))[j]
      end
      for i in 1:n-r, j in 1:n
         A[i+r,j]=phi(gen(W,i))[j]
      end
      C = _upper_triangular_version(A*B*transpose(A))      
   else
      degF=0
      if f.descr==:hermitian e = div(degree(F),2) end
      C,A,r = _find_radical(B,F,n,n; e=degF, _is_symmetric=true)
   end

   if r<n
      fn = SesquilinearForm(submatrix(C,1,1,r,r),f.descr)
   else
      fn = f
   end

   e=0
   if (fn.descr==:quadratic || fn.descr==:symmetric) && iseven(r)
      if witt_index(fn)== div(r,2)
         e = 1
      else
         e = -1
      end
   end
   Xf = iscongruent(SesquilinearForm( fn.mat_iso(_standard_form(fn.descr,e,r,F)), fn.descr),fn)[2]
# if dimension is odd, fn may be congruent to a scalar multiple of the standard form
# TODO: I don't really need a primitive_element(F); I just need a non-square in F. Is there a faster way to get it?
   if Xf==nothing && isodd(r)
      Xf = iscongruent(SesquilinearForm( primitive_element(F)*fn.mat_iso(_standard_form(fn.descr,e,r,F)), fn.descr),fn)[2]
   end


   if f.descr==:hermitian G = GU(r,Int(characteristic(F)^div(degree(F),2)))
   elseif f.descr==:alternating G = Sp(r, F)
   elseif isodd(r) G = GO(0,r,F)
   elseif witt_index(f)==div(r,2) G=GO(1,r,F)
   else G=GO(-1,r,F)
   end

# 2x2 block matrix. The four blocks are: group of isometries for fn,
# everything, zero matrix, GL(n-r,F)
   if r<n
      Xfn = Xf^-1
      An=A^-1
      Idn = identity_matrix(F,n)
      L = dense_matrix_type(elem_type(F))[]
      for i in 1:ngens(G)
         push!(L, An*insert_block(Idn,Xfn*(G[i].elm)*Xf,1,1)*A)
      end
      for g in _gens_for_GL(n-r,F)
         push!(L, An*insert_block(Idn,g,r+1,r+1)*A)
      end
# TODO: not quite sure whether the last element (the one with i=r) is sufficient to generate the whole top-right block
      for i in 1:r
         y = deepcopy(Idn)
         y[i,r+1]=1
         push!(L,An*y*A)
      end
      return matrix_group(L)
   else
      Xf = GL(r,base_ring(f))(Xf)
      return G^Xf
   end
end
