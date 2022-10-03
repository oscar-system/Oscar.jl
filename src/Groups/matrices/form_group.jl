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
    preserved_sesquilinear_forms,
    automorphism_group,
    orthogonal_group,
    orthogonal_sign,
    unitary_group

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

   r,K = nullspace(vcat(M))
   
   return [matrix(F,n,n,[K[i,j] for i in 1:n^2]) for j in 1:r]
end

"""
    invariant_sesquilinear_forms(G::MatrixGroup)

Return a generating set for the vector spaces of sesquilinear non-bilinear forms preserved by the group `G`.
An exception is thrown if `base_ring(G)` is not a finite field with even degree
over its prime subfield.

!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_sesquilinear_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
   @assert F isa FinField "At the moment, only finite fields are considered"
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

   r,K = nullspace(vcat(M))
   
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

   r,K = nullspace(vcat(M))
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

   r,K = nullspace(vcat(M))

   L = T[]
   for j in 1:r
      B = zero_matrix(F,n,n)
      B[1:n-1,2:n] = upper_triangular_matrix(K[1:div(n*(n-1),2),j])
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
An exception is thrown if `base_ring(G)` is not a finite field with even degree
over its prime subfield.

!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_hermitian_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
   @assert F isa FinField "At the moment, only finite fields are considered"
   n = degree(G)
   M = T[]

   p = characteristic(F)
   d = degree(F)
   iseven(d) || return M

   F0 = GF(Int(p), div(d,2))
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

   n_gens,K = nullspace(vcat(M))
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
    invariant_bilinear_form(G::MatrixGroup)

Return an invariant bilinear form for the group `G`.
An exception is thrown if the module induced by the action of `G`
is not absolutely irreducible.

!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.
"""
function invariant_bilinear_form(G::MatrixGroup)
   V = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(G.X), codomain(G.ring_iso))
   B = GAP.Globals.MTX.InvariantBilinearForm(V)
   return preimage_matrix(G.ring_iso, B)
end

"""
    invariant_sesquilinear_form(G::MatrixGroup)

Return an invariant sesquilinear (non bilinear) form for the group `G`.
An exception is thrown if the module induced by the action of `G`
is not absolutely irreducible or if the group is defined over a finite field
of odd degree over the prime field.

!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.
"""
function invariant_sesquilinear_form(G::MatrixGroup)
   isodd(degree(base_ring(G))) && throw(ArgumentError("group is defined over a field of odd degree"))
   V = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(G.X), codomain(G.ring_iso))
   B = GAP.Globals.MTX.InvariantSesquilinearForm(V)
   return preimage_matrix(G.ring_iso, B)
end

"""
    invariant_quadratic_form(G::MatrixGroup)

Return an invariant quadratic form for the group `G`.
An exception is thrown if the module induced by the action of `G`
is not absolutely irreducible.

!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.
"""
function invariant_quadratic_form(G::MatrixGroup)
   if iseven(characteristic(base_ring(G)))
      V = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(G.X), codomain(G.ring_iso))
      B = GAP.Globals.MTX.InvariantQuadraticForm(V)
      return _upper_triangular_version(preimage_matrix(G.ring_iso, B))
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
      f = quadratic_form(preimage_matrix(G.ring_iso, GAP.Globals.GramMatrix(f_gap)))
      f.X = f_gap
      f.ring_iso = G.ring_iso
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
      if GAPWrap.IsHermitianForm(f_gap)
         f = hermitian_form(preimage_matrix(G.ring_iso, GAP.Globals.GramMatrix(f_gap)))
      elseif GAPWrap.IsSymmetricForm(f_gap)
         f = symmetric_form(preimage_matrix(G.ring_iso, GAP.Globals.GramMatrix(f_gap)))
      elseif GAPWrap.IsAlternatingForm(f_gap)
         f = alternating_form(preimage_matrix(G.ring_iso, GAP.Globals.GramMatrix(f_gap)))
      else
         error("Invalid form")
      end
      f.X = f_gap
      push!(R,f)
   end
   return R
end


"""
    orthogonal_sign(G::MatrixGroup)

For absolutely irreducible `G` of degree `n` and such that `base_ring(G)`
is a finite field, return
- `nothing` if `G` does not preserve a nonzero quadratic form,
- `0` if `n` is odd and `G` preserves a nonzero quadratic form,
- `1` if `n` is even and `G` preserves a nonzero quadratic form of `+` type,
- `-1` if `n` is even and `G` preserves a nonzero quadratic form of `-` type.
"""
function orthogonal_sign(G::MatrixGroup)
    R = base_ring(G)
    R isa FinField || error("G must be a matrix group over a finite field")
    M = GAP.Globals.GModuleByMats(GAP.Globals.GeneratorsOfGroup(G.X),
                                  codomain(iso_oscar_gap(R)))
    sign = GAP.Globals.MTX.OrthogonalSign(M)
    sign === GAP.Globals.fail && return nothing
    # If the characteristic is odd and there is an invariant
    # antisymmetric bilinear form then `GAP.Globals.MTX.OrthogonalSign`
    # does *not* return `fail`,
    # see https://github.com/gap-system/gap/issues/4936.
    if isodd(characteristic(R))
      Q = GAP.Globals.MTX.InvariantQuadraticForm(M)
      if Q === GAP.Globals.fail || Q == - GAP.Globals.TransposedMat(Q)
        return nothing
      end
    end
    return sign
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
      fn = SesquilinearForm(C[1:r, 1:r],f.descr)
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
   Xf = is_congruent(SesquilinearForm(preimage_matrix(fn.ring_iso, _standard_form(fn.descr, e, r, F)), fn.descr), fn)[2]
# if dimension is odd, fn may be congruent to a scalar multiple of the standard form
# TODO: I don't really need a primitive_element(F); I just need a non-square in F. Is there a faster way to get it?
   if Xf==nothing && isodd(r)
      Xf = is_congruent(SesquilinearForm(primitive_element(F)*preimage_matrix(fn.ring_iso, _standard_form(fn.descr, e, r, F)), fn.descr), fn)[2]
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
         temp = deepcopy(Idn)
         temp[1:r,1:r] = Xfn*(G[i].elm)*Xf
         push!(L, An*temp*A)
      end
      for g in _gens_for_GL(n-r,F)
         temp = deepcopy(Idn)
         temp[r+1:n, r+1:n] = g
         push!(L, An*temp*A)
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

"""
    isometry_group(L::AbsLat) -> MatrixGroup

Return the group of isometries of the lattice `L`.

The transformations are represented with respect to the ambient space of `L`.
"""
@attr MatrixGroup{elem_type(base_field(L)), dense_matrix_type(elem_type(base_field(L)))} function isometry_group(L::Hecke.AbsLat)
   gens = Hecke.automorphism_group_generators(L)
   G = matrix_group(gens)
   return G
end

@attr MatrixGroup{fmpq,fmpq_mat} function isometry_group(L::ZLat)
   # corner case
   if rank(L) == 0
      return matrix_group(identity_matrix(QQ,degree(L)))
   end
   @req is_definite(L) "lattice must be definite"
   G, _ = _isometry_group_via_decomposition(L)
   return G
end

"""
    _isometry_group_via_decomposition(L::ZLat) -> Tuple{MatrixGroup, Vector{fmpq_mat}}

Compute the group of isometries of the definite lattice `L` using an orthogonal decomposition.
"""
function _isometry_group_via_decomposition(L::ZLat; closed = true, direct=true)
  # TODO: adapt the direct decomposition approach for AbsLat
  # in most examples `direct=true` seems to be faster by a factor of 7
  # but in some examples it is also slower ... up to a factor of 15
  if gram_matrix(L)[1,1]<0
    L = rescale(L, -1)
  end
  # construct the sublattice M1 of L generated by the shortest vectors
  V = ambient_space(L)
  # for simplicity we work with the ambient representation
  # TODO: Swap to action on vectors once Hecke 0.15.3 is released
  _sv = shortest_vectors(L)
  if eltype(_sv) === Vector{fmpz}
    sv = fmpz_mat[matrix(ZZ, 1, length(v), v) for v in _sv]
  else
    sv = _sv
  end
  bL = basis_matrix(L)
  sv1 = fmpq_mat[v*bL for v in sv]
  h = _row_span!(sv)*bL
  M1 = lattice(V, h)
  if closed
    # basically doubles the memory usage of this function
    # a more elegant way could be to work with the corresponding projective representation
    append!(sv1, [-v for v in sv1])
  end

  # base case of the recursion
  M1primitive = primitive_closure(L, M1)

  #=
  # the following is slower than computing the automorphism_group generators of
  # M1primitive outright
  gensOM1 = Hecke.automorphism_group_generators(M1)
  OM1 = matrix_group(gensOM1)
  if M1primitive == M1
    O1 = OM1
  else
    # M1 is generated by its shortest vectors only up to finite index
    @vprint :Lattice 3 "Computing overlattice stabilizers \n"
    O1,_ = stabilizer(OM1, M1primitive, on_lattices)
  end
  =#
  O1 = matrix_group(Hecke.automorphism_group_generators(M1primitive))
  _set_nice_monomorphism!(O1, sv1, closed=closed)
  if rank(M1) == rank(L)
    @hassert :Lattice 2 M1primitive == L
    return O1, sv1
  end

  # decompose as a primitive extension: M1primitive + M2 \subseteq L
  M2 = Hecke.orthogonal_submodule(L, M1)

  @vprint :Lattice 3 "Computing orthogonal groups via an orthogonal decomposition\n"
  # recursion
  O2, sv2 = _isometry_group_via_decomposition(M2, closed=closed,direct=direct)


  # In what follows we compute the stabilizer of L in O1 x O2
  if direct
    gens12 = vcat(gens(O1),gens(O2))
    G = matrix_group(gens12)
    sv = append!(sv1, sv2)
    S,_ = stabilizer(G, L, on_lattices)
    _set_nice_monomorphism!(S, sv, closed=closed)
    return S, sv
  end

  phi, i1, i2 = glue_map(L, M1primitive, M2)
  H1 = domain(phi)
  H2 = codomain(phi)


  @vprint :Lattice 2 "glue order: $(order(H1))\n"
  # the stabilizer and kernel computations are expensive
  # alternatively we could first project to the orthogonal group of the
  # discriminant group and create an on_subgroup action
  @vprint :Lattice 3 "Computing glue stabilizers \n"
  G1, _ = stabilizer(O1, cover(H1), on_lattices)
  G2, _ = stabilizer(O2, cover(H2), on_lattices)
  # _set_nice_monomorphism!(G1, sv1, closed=closed)
  # _set_nice_monomorphism!(G2, sv2, closed=closed)

  # now we may alter sv1
  sv = append!(sv1, sv2)

  G1q =  _orthogonal_group(H1, fmpz_mat[hom(H1, H1, TorQuadModElem[H1(lift(x) * matrix(g)) for x in gens(H1)]).map_ab.map for g in gens(G1)])
  G2q =  _orthogonal_group(H2, fmpz_mat[hom(H2, H2, TorQuadModElem[H2(lift(x) * matrix(g)) for x in gens(H2)]).map_ab.map for g in gens(G2)])

  psi1 = hom(G1, G1q, gens(G1q), check=false)
  psi2 = hom(G2, G2q, gens(G2q), check=false)
  @vprint :Lattice 2 "Computing the kernel of $(psi1)\n"
  K = gens(kernel(psi1)[1])
  @vprint :Lattice 2 "Computing the kernel of $(psi2)\n"
  append!(K, gens(kernel(psi2)[1]))
  @vprint :Lattice 2 "Lifting \n"

  T = _orthogonal_group(H1, [(phi * hom(g) * inv(phi)).map_ab.map for g in gens(G2q)])
  S = _as_subgroup(G1q, GAP.Globals.Intersection(T.X, G1q.X))[1]
  append!(K, [preimage(psi1, g) * preimage(psi2, G2q(inv(phi) * hom(g) * phi)) for g in gens(S)])
  G = matrix_group(matrix.(K))
  @hassert :Lattice 2 all(on_lattices(L,g)==L for g in gens(G))
  _set_nice_monomorphism!(G, sv, closed=closed)
  @vprint :Lattice 2 "Done \n"
  return G, sv
end

function on_lattices(L::ZLat, g::MatrixGroupElem{fmpq,fmpq_mat})
  V = ambient_space(L)
  return lattice(V, basis_matrix(L) * matrix(g), check=false)
end

"""
    on_matrix(x, g::MatrixGroupElem{fmpq,fmpq_mat})

Return `x*g`.
"""
function on_matrix(x, g::MatrixGroupElem{fmpq,fmpq_mat})
  return x*matrix(g)
end

"""
    _set_nice_monomorphism!(G::MatrixGroup, short_vectors; closed=false)

Use the permutation action of `G` on `short_vectors` to represent `G` as a
finite permutation group.

Internally this sets a `NiceMonomorphism` for the underlying gap group.
No input checks whatsoever are performed.

It is assumed that the corresponding action homomorphism is injective.
Setting `closed = true` assumes that `G` actually preserves `short_vectors`.
"""
#
function _set_nice_monomorphism!(G::MatrixGroup, short_vectors; closed=false)
  phi = action_homomorphism(gset(G, on_matrix, short_vectors, closed=closed))
  GAP.Globals.SetIsInjective(phi.map, true) # fixes an infinite recursion
  GAP.Globals.SetIsHandledByNiceMonomorphism(G.X, true)
  GAP.Globals.SetNiceMonomorphism(G.X, phi.map)
end

function _row_span!(L::Vector{fmpz_mat})
  l = length(L)
  d = ncols(L[1])
  m = min(2*d,l)
  B = sparse_matrix(reduce(vcat, L[1:m]))
  h = matrix(hnf(B, truncate = true))
  for i in (m+1):l
    b = L[i]
    Hecke.reduce_mod_hnf_ur!(b, h)
    if iszero(b)
      continue
    else
      h = vcat(h, b)
      hnf!(h)
    end
  end
  return h[1:rank(h),:]
end

automorphism_group(L::Hecke.AbsLat) = isometry_group(L)

orthogonal_group(L::Hecke.ZLat) = isometry_group(L)

orthogonal_group(L::Hecke.QuadLat) = isometry_group(L)

unitary_group(L::Hecke.HermLat) = isometry_group(L)


