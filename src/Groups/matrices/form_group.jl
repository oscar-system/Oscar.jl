# TODO: in this file several ways to get the preserved forms by a matrix group are implemented
# we can choose to keep all of them or only some of them

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

   r,K = nullspace(reduce(vcat, M))
   
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
   @req F isa FinField "At the moment, only finite fields are considered"
   @req iseven(degree(F)) "Base ring has no even degree"
   n = degree(G)
   M = T[]
   for mat in gens(G)
      mmat = map(y->frobenius(y,div(degree(F),2)),matrix(mat)^-1)
      MM = zero_matrix(F,n^2,n^2)
      for i in 1:n, j in 1:n, k in 1:n
         MM[n*(i-1)+j,n*(k-1)+j] += mat[i,k]
         MM[n*(i-1)+j,n*(i-1)+k] -= mmat[j,k]
      end
      push!(M,MM)
   end

   r,K = nullspace(reduce(vcat, M))
   
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

   r,K = nullspace(reduce(vcat, M))
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

   r,K = nullspace(reduce(vcat, M))

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
   @req F isa FinField "At the moment, only finite fields are considered"
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

   n_gens,K = nullspace(reduce(vcat, M))
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

Return the Gram matrix of an invariant bilinear form for `G`.
An exception is thrown if the module induced by the action of `G`
is not absolutely irreducible.

!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> invariant_bilinear_form(Sp(4, 2))
[0   0   0   1]
[0   0   1   0]
[0   1   0   0]
[1   0   0   0]
```
"""
function invariant_bilinear_form(G::MatrixGroup)
   V = GAP.Globals.GModuleByMats(GAPWrap.GeneratorsOfGroup(GapObj(G)), codomain(_ring_iso(G)))
   B = GAP.Globals.MTX.InvariantBilinearForm(V)
   return preimage_matrix(_ring_iso(G), B)
end

"""
    invariant_sesquilinear_form(G::MatrixGroup)

Return the Gram matrix of an invariant sesquilinear (non bilinear) form for `G`.

An exception is thrown if the module induced by the action of `G`
is not absolutely irreducible or if `G` is defined over a finite field
of odd degree over the prime field.

!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> invariant_sesquilinear_form(GU(4, 2))
[0   0   0   1]
[0   0   1   0]
[0   1   0   0]
[1   0   0   0]
```
"""
function invariant_sesquilinear_form(G::MatrixGroup)
   @req iseven(degree(base_ring(G))) "group is defined over a field of odd degree"
   V = GAP.Globals.GModuleByMats(GAPWrap.GeneratorsOfGroup(GapObj(G)), codomain(_ring_iso(G)))
   B = GAP.Globals.MTX.InvariantSesquilinearForm(V)
   return preimage_matrix(_ring_iso(G), B)
end

"""
    invariant_quadratic_form(G::MatrixGroup)

Return the Gram matrix of an invariant quadratic form for `G`.
An exception is thrown if the module induced by the action of `G`
is not absolutely irreducible.

!!! warning "Note:"
    At the moment, the output is returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> invariant_quadratic_form(GO(1, 4, 2))
[0   1   0   0]
[0   0   0   0]
[0   0   0   1]
[0   0   0   0]
```
"""
function invariant_quadratic_form(G::MatrixGroup)
   if iseven(characteristic(base_ring(G)))
      V = GAP.Globals.GModuleByMats(GAPWrap.GeneratorsOfGroup(GapObj(G)), codomain(_ring_iso(G)))
      B = GAP.Globals.MTX.InvariantQuadraticForm(V)
      return _upper_triangular_version(preimage_matrix(_ring_iso(G), B))
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
   L = GAP.Globals.PreservedQuadraticForms(GapObj(G))
   R = SesquilinearForm{S}[]
   for f_gap in L
      f = quadratic_form(preimage_matrix(_ring_iso(G), GAP.Globals.GramMatrix(f_gap)))
      f.X = f_gap
      f.ring_iso = _ring_iso(G)
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
   L = GAP.Globals.PreservedSesquilinearForms(GapObj(G))
   R = SesquilinearForm{S}[]
   for f_gap in L
      if GAPWrap.IsHermitianForm(f_gap)
         f = hermitian_form(preimage_matrix(_ring_iso(G), GAP.Globals.GramMatrix(f_gap)))
      elseif GAPWrap.IsSymmetricForm(f_gap)
         f = symmetric_form(preimage_matrix(_ring_iso(G), GAP.Globals.GramMatrix(f_gap)))
      elseif GAPWrap.IsAlternatingForm(f_gap)
         f = alternating_form(preimage_matrix(_ring_iso(G), GAP.Globals.GramMatrix(f_gap)))
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

# Examples
```jldoctest
julia> orthogonal_sign(GO(1, 4, 2))
1

julia> orthogonal_sign(GO(-1, 4, 2))
-1
```
"""
function orthogonal_sign(G::MatrixGroup)
    R = base_ring(G)
    R isa FinField || error("G must be a matrix group over a finite field")
    M = GAP.Globals.GModuleByMats(GAPWrap.GeneratorsOfGroup(GapObj(G)),
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
      return GAP.Globals.InvariantQuadraticForm(GapObj(GO(e,n,q))).matrix
   elseif descr==:symmetric #|| descr==:alternating
      return GAP.Globals.InvariantBilinearForm(GapObj(GO(e,n,q))).matrix
   elseif descr==:hermitian
      return GAP.Globals.InvariantSesquilinearForm(GapObj(GU(n,q))).matrix
   elseif descr==:alternating
      return GAP.Globals.InvariantBilinearForm(GapObj(Sp(n,q))).matrix
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
      V = vector_space(F,n)
      U,e = complement(V,W)
      A = zero_matrix(F,n,n)
      r = vector_space_dim(U)
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
   Xf = is_congruent(SesquilinearForm(preimage_matrix(_ring_iso(fn), _standard_form(fn.descr, e, r, F)), fn.descr), fn)[2]
# if dimension is odd, fn may be congruent to a scalar multiple of the standard form
# TODO: I don't really need a primitive_element(F); I just need a non-square in F. Is there a faster way to get it?
   if Xf === nothing && isodd(r)
      Xf = is_congruent(SesquilinearForm(primitive_element(F)*preimage_matrix(_ring_iso(fn), _standard_form(fn.descr, e, r, F)), fn.descr), fn)[2]
   end


   if f.descr === :hermitian
      G = GU(r,Int(characteristic(F)^div(degree(F),2)))
   elseif f.descr === :alternating
      G = Sp(r, F)
   elseif isodd(r)
      G = GO(0,r,F)
   elseif 2*witt_index(f) == r
      G = GO(1,r,F)
   else
      G = GO(-1,r,F)
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
         temp[1:r,1:r] = Xfn*matrix(G[i])*Xf
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
    isometry_group(L::AbstractLat; depth::Int = -1, bacher_depth::Int = 0) -> MatrixGroup

Return the group of isometries of the lattice `L`.

The transformations are represented with respect to the ambient space of `L`.

Setting the parameters `depth` and `bacher_depth` to a positive value may improve
performance. If set to `-1` (default), the used value of `depth` is chosen
heuristically depending on the rank of `L`. By default, `bacher_depth` is set to `0`.
"""
@attr matrix_group_type(S) function isometry_group(L::Hecke.AbstractLat{S}; depth::Int=-1, bacher_depth::Int=0) where S
  gens = automorphism_group_generators(L; depth, bacher_depth)
  G = matrix_group(gens)
  return G
end

@doc raw"""
    isometry_group(
      L::ZZLat;
      algorithm::Symbol=:direct,
      depth::Int=-1,
      bacher_depth::Int=0,
      kwargs...,
     ) -> MatrixGroup

Given an integer lattice $L$ which is definite or of rank 2, return the
isometry group $O(L)$ of $L$.

One can choose which algorithm to use to compute $O(L)$. For now, we
only support the following algorithms:
- `:direct`: compute generators of $O(L)$ using an algorithm of
  Plesken--Souvignier;
- `:decomposition`: compute iteratively $O(L)$ by decomposing $L$ into
  $O(L)$-stable sublattices;
- `:dispatch`: compute $O(L)$ using heuristics which dispatch (iteratively) to
  the best algorithm between `:direct` and `:decomposition`

Setting the parameters `depth` and `bacher_depth` to a positive value may
improve performance. If set to `-1` (default), the used value of `depth`
is chosen heuristically depending on the rank of `L`. By default,
`bacher_depth` is set to `0`.
"""
@attr QQMatrixGroup function isometry_group(
  L::ZZLat;
  algorithm=:direct,
  depth::Int=-1,
  bacher_depth::Int=0,
  kwargs...,
)
  # corner cases
  if rank(L) == 0
    G = matrix_group(identity_matrix(QQ, degree(L)))
  elseif rank(L) == 1
    G = matrix_group(extend_to_ambient_space(L, -identity_matrix(QQ, rank(L)); check=false))
  elseif rank(L) == 2
    gene = automorphism_group_generators(L)
    G = matrix_group(QQMatrix[change_base_ring(QQ, m) for m in gene])
  end

  @req is_definite(L) "Lattice must be definite or of rank at most 2"

  if algorithm == :direct
    gens = automorphism_group_generators(L; depth, bacher_depth)
    G = matrix_group(gens)
  elseif algorithm == :decomposition
    G, _ = _isometry_group_via_decomposition(L; depth, bacher_depth, kwargs...)
  elseif algorithm == :dispatch
    G, _ = _isometry_group_via_heuristics(L; depth, bacher_depth, kwargs...)
  else
    error("Unknown algorithm: for the moment, we only support :direct, :decomposition and :dispatch")
  end
  return G
end

"""
    _isometry_group_via_heuristics(
      L::ZZLat;
      closed::Bool=true,
      direct::Bool=true,
      depth::Int=-1,
      bacher_depth::Int=0,
    ) -> MatrixGroup, Vector{QQMatrix}

Compute the group of isometries of the definite lattice `L` by dispatching
to Plesken--Souvignier or the decomposition algorithm, depending on the Gram
matrix of `L`.
"""
function _isometry_group_via_heuristics(
  L::ZZLat;
  closed::Bool=true,
  direct::Bool=true,
  depth::Int=-1,
  bacher_depth::Int=0,
  set_nice_mono::Bool=true,
)
  if gram_matrix(L)[1,1] < 0
    L_int = rescale(L, -1)
  else
    L_int = L
  end

  if !direct
    d = denominator(scale(L_int))
    if d > 1
      L_int = rescale(L_int, d)
    end
  end

  # for now we use very naive heuristics
  L_int = lll(L_int)
  d = diagonal(gram_matrix(L_int))
  m1 = minimum(d)
  m2 = maximum(d)

  if m2 <= 2*m1 || rank(L) <= 10
    gens = automorphism_group_generators(L_int; depth, bacher_depth)
    G = matrix_group(gens)
    sv = Vector{QQFieldElem}[]
  else
    G, sv = _isometry_group_via_decomposition(L_int; use_heuristics=true, depth, bacher_depth, closed, direct, set_nice_mono)
  end
  return G, sv
end

"""
    _isometry_group_via_decomposition(
      L::ZZLat;
      use_heuristics::Bool=false,
      closed::Bool=true,
      direct::Bool=true,
      depth::Int=-1,
      bacher_depth::Int=0,
    ) -> MatrixGroup, Vector{QQMatrix}

Compute the group of isometries of the definite lattice `L` using an orthogonal
decomposition.
"""
function _isometry_group_via_decomposition(
  L::ZZLat;
  use_heuristics::Bool=false,
  closed::Bool=true,
  direct::Bool=true,
  depth::Int=-1,
  bacher_depth::Int=0,
  set_nice_mono::Bool=true,
)
  # TODO: adapt the direct decomposition approach for AbstractLat
  # in most examples `direct=true` seems to be faster by a factor of 7
  # but in some examples it is also slower ... up to a factor of 15
  if gram_matrix(L)[1,1] < 0
    L = rescale(L, -1)
  end

  if !direct
    d = denominator(scale(L))
    if d > 1
      L = rescale(L, d)
    end
  end

  # construct the sublattice M1 of L generated by the shortest vectors
  V = ambient_space(L)
  # for simplicity we work with the ambient representation
  # TODO: Swap to action on vectors once Hecke 0.15.3 is released
  _sv1 = shortest_vectors(L)
  bL = basis_matrix(L)
  sv1 = Vector{QQFieldElem}[v*bL for v in _sv1]
  h = _row_span!(_sv1)*bL
  M1 = lattice(V, h)
  if set_nice_mono && closed
    # basically doubles the memory usage of this function
    # a more elegant way could be to work with the corresponding projective representation
    append!(sv1, [-v for v in sv1])
  end

  # base case of the recursion
  M1primitive = primitive_closure(L, M1)

 #=
  # the following is slower than computing the automorphism_group generators of
  # M1primitive outright
  gensOM1 = automorphism_group_generators(M1, depth = depth, bacher_depth = bacher_depth)
  OM1 = matrix_group(gensOM1)
  if M1primitive == M1
    O1 = OM1
  else
    # M1 is generated by its shortest vectors only up to finite index
    @vprint :Lattice 3 "Computing overlattice stabilizers \n"
    O1,_ = stabilizer(OM1, M1primitive, on_lattices)
  end
 =#

  O1 = matrix_group(automorphism_group_generators(M1primitive; depth, bacher_depth))
  if set_nice_mono
    _set_nice_monomorphism!(O1, sv1; closed)
  end
  if rank(M1) == rank(L)
    @hassert :Lattice 2 M1primitive == L
    return O1, sv1
  end

  # decompose as a primitive extension: M1primitive + M2 \subseteq L
  M2 = orthogonal_submodule(L, M1)

  @vprint :Lattice 3 "Computing orthogonal groups via an orthogonal decomposition\n"
  # recursion
  if use_heuristics
    O2, sv2 = _isometry_group_via_heuristics(M2; closed, direct, depth, bacher_depth, set_nice_mono)
  else
    O2, sv2 = _isometry_group_via_decomposition(M2; closed, direct, depth, bacher_depth, set_nice_mono)
  end

  if set_nice_mono && isempty(sv2)
    _sv2 = shortest_vectors(M2)
    sv2 = Vector{QQFieldElem}[v*basis_matrix(M2) for v in _sv2]
    if closed
      append!(sv2, [-v for v in sv2])
    end
  end
  sv = append!(sv1, sv2)

  # In what follows we compute the stabilizer of L in O1 x O2
  if direct
    gens12 = vcat(gens(O1), gens(O2))
    G = matrix_group(gens12)
    S, _ = stabilizer(G, L, on_lattices)
  else
    gensS = stabilizer_in_diagonal_action(L, M1primitive, M2, O1, O2; check=false, is_finite_known=(true, true))
    S = matrix_group(gensS)
  end
  if set_nice_mono
    _set_nice_monomorphism!(S, sv; closed)
  end
  return S, sv
end

function on_lattices(L::ZZLat, g::MatrixGroupElem{QQFieldElem,QQMatrix})
  V = ambient_space(L)
  return lattice(V, basis_matrix(L) * matrix(g); check=false)
end

"""
    on_vector(x::Vector{QQFieldElem}, g::MatrixGroupElem{QQFieldElem,QQMatrix})

Return `x*g`.
"""
function on_vector(x::Vector{QQFieldElem}, g::MatrixGroupElem{QQFieldElem,QQMatrix})
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
  phi = action_homomorphism(gset(G, on_vector, short_vectors; closed))
  GAP.Globals.SetIsInjective(phi.map, true) # fixes an infinite recursion
  GAP.Globals.SetIsHandledByNiceMonomorphism(GapObj(G), true)
  GAP.Globals.SetNiceMonomorphism(GapObj(G), phi.map)
end

function _row_span!(L::Vector{Vector{Int}})
  l = length(L)
  d = length(L[1])
  m = min(2*d,l)
  B = sparse_matrix(matrix(ZZ, L))
  h = matrix(hnf!(B))
  for i in (m+1):l
    b = matrix(ZZ, 1, d, L[i])
    Hecke.reduce_mod_hnf_ur!(b, h)
    if iszero(b)
      continue
    else
      h = vcat(h, b)
      hnf!(h)
    end
  end
  return h[1:rank(h), :]
end

automorphism_group(L::Hecke.AbstractLat; kwargs...) = isometry_group(L; kwargs...)

orthogonal_group(L::Hecke.ZZLat; kwargs...) = isometry_group(L; kwargs...)

orthogonal_group(L::Hecke.QuadLat; kwargs...) = isometry_group(L; kwargs...)

unitary_group(L::Hecke.HermLat; kwargs...) = isometry_group(L; kwargs...)

@doc raw"""
    stable_orthogonal_group(
      L::ZZLat;
      kwargs...,
    ) -> MatrixGroup, GAPGroupHomomorphism

Given an integer lattice $L$ which is definite or of rank 2, return the
subgroup $O^\#(L)$ of the orthogonal group of $L$ consisting of isometries
acting trivially on the discriminant group of $L$.

The function first computes the orthogonal group of ``L``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).

# Examples
```jldoctest
julia> A5 = root_lattice(:A, 5);

julia> H, _ = stable_orthogonal_group(A5);

julia> order(H)
720
```
"""
function stable_orthogonal_group(
    L::ZZLat;
    kwargs...,
  )
  OL = orthogonal_group(L; kwargs...)
  return stable_subgroup(L, OL; check=false)
end

@doc raw"""
    special_orthogonal_group(
      L::ZZLat;
      kwargs...,
    ) -> MatrixGroup, GAPGroupHomomorphism

Given an integer lattice $L$ which is definite or of rank 2, return the
subgroup $SO(L)$ of the orthogonal group of $L$ consisting of isometries
with determinant ``1``.

The function first computes the orthogonal group of ``L``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).

# Examples
```jldoctest
julia> D5 = root_lattice(:D, 5);

julia> H, _ = special_orthogonal_group(D5);

julia> order(H)
1920
```
"""
function special_orthogonal_group(
    L::ZZLat;
    kwargs...,
  )
  OL = orthogonal_group(L; kwargs...)
  return special_subgroup(L, OL; check=false)
end

# We do not export this one, it is just a shortcut
@doc raw"""
    _special_stable_orthogonal_group(
      L::ZZLat;
      kwargs...,
    ) -> MatrixGroup, GAPGroupHomomorphism

Given an integer lattice $L$ which is definite or of rank 2, return the
subgroup $SO^\#(L)$ of the orthogonal group of $L$ consisting of isometries
acting trivially on the discriminant group of $L$ and of determinant ``1``.

The function first computes the orthogonal group of ``L``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).

# Examples
```jldoctest
julia> A6 = root_lattice(:A, 6);

julia> H, _ = Oscar._special_stable_orthogonal_group(A6);

julia> describe(H)
"A7"
```
"""
function _special_stable_orthogonal_group(
    L::ZZLat;
    kwargs...,
  )
  OL = orthogonal_group(L; kwargs...)
  return Oscar._special_stable_subgroup(L, OL; check=false)
end
