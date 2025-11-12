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

Return a generating set for the vector space of bilinear forms preserved by `G`.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> G = sylow_subgroup(GL(3, 2), 2)[1];

julia> length(invariant_bilinear_forms(G))
2
```
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

Return a generating set for the vector space of sesquilinear forms
preserved by the group `G`.

An exception is thrown if `base_ring(G)` is not a finite field with even degree
over its prime subfield.

!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> G = sylow_subgroup(GL(3, 4), 2)[1];

julia> length(invariant_sesquilinear_forms(G))
1
```
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

Return a generating set for the vector space of quadratic forms
preserved by the group `G`.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> G = sylow_subgroup(GL(3, 2), 2)[1];

julia> length(invariant_quadratic_forms(G))
2
```
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

Return a generating set for the vector space of symmetric forms
preserved by the group `G`.

!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
!!! warning "Note:"
    Work properly only in odd characteristic. In even characteristic, only alternating forms are found.

# Examples
```jldoctest
julia> G = sylow_subgroup(GL(3, 3), 2)[1];

julia> length(invariant_symmetric_forms(G))
1
```
"""
invariant_symmetric_forms(G::MatrixGroup{S,T}) where {S,T} = T[x + transpose(x) for x in invariant_quadratic_forms(G)]

# METHOD: if B = (b_ij) is the solution matrix, then b_ji = -b_ij and b_ii=0;
# hence, the dimension of the linear system can be reduced to n(n-1)/2
"""
    invariant_alternating_forms(G::MatrixGroup)

Return a generating set for the vector space of alternating forms
preserved by the group `G`.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> G = sylow_subgroup(GL(3, 4), 2)[1];

julia> length(invariant_alternating_forms(G))
1
```
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

Return a generating set for the vector space of hermitian forms
preserved by the group `G`.
An exception is thrown if `base_ring(G)` is not a finite field with even degree
over its prime subfield.

!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.

# Examples
```jldoctest
julia> G = sylow_subgroup(GL(3, 4), 2)[1];

julia> length(invariant_hermitian_forms(G))
1
```
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
    invariant_bilinear_form(G::MatrixGroup{T}) where T <: FinFieldElem

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
function invariant_bilinear_form(G::MatrixGroup{T}) where T <: FinFieldElem
   V = GAP.Globals.GModuleByMats(GAPWrap.GeneratorsOfGroup(GapObj(G)), codomain(_ring_iso(G)))
   B = GAP.Globals.MTX.InvariantBilinearForm(V)
   return preimage_matrix(_ring_iso(G), B)
end

"""
    invariant_sesquilinear_form(G::MatrixGroup{T}) where T <: FinFieldElem

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
function invariant_sesquilinear_form(G::MatrixGroup{T}) where T <: FinFieldElem
   @req iseven(degree(base_ring(G))) "group is defined over a field of odd degree"
   V = GAP.Globals.GModuleByMats(GAPWrap.GeneratorsOfGroup(GapObj(G)), codomain(_ring_iso(G)))
   B = GAP.Globals.MTX.InvariantSesquilinearForm(V)
   return preimage_matrix(_ring_iso(G), B)
end

"""
    invariant_quadratic_form(G::MatrixGroup{T}) where T <: FinFieldElem

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
function invariant_quadratic_form(G::MatrixGroup{T}) where T <: FinFieldElem
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
    preserved_quadratic_forms(G::MatrixGroup{T}) where T <: FinFieldElem

Uses random methods to find all of the quadratic forms preserved by `G` up to a scalar
(i.e. such that `G` is a group of similarities for the forms). 
Since the procedure relies on a pseudo-random generator, 
the user may need to execute the operation more than once to find all invariant quadratic forms.
"""
function preserved_quadratic_forms(G::MatrixGroup{T}) where T <: FinFieldElem
   L = GAP.Globals.PreservedQuadraticForms(GapObj(G))
   R = QuadraticForm{T}[]
   for f_gap in L
      f = quadratic_form(preimage_matrix(_ring_iso(G), GAP.Globals.GramMatrix(f_gap)))
      f.X = f_gap
      f.ring_iso = _ring_iso(G)
      push!(R,f)
   end
   return R
end

"""
    preserved_sesquilinear_forms(G::MatrixGroup{T}) where T <: FinFieldElem

Uses random methods to find all of the sesquilinear forms preserved by `G` up to a scalar
(i.e. such that `G` is a group of similarities for the forms).
Since the procedure relies on a pseudo-random generator,
the user may need to execute the operation more than once to find all invariant sesquilinear forms.
"""
function preserved_sesquilinear_forms(G::MatrixGroup{T}) where T <: FinFieldElem
   L = GAP.Globals.PreservedSesquilinearForms(GapObj(G))
   R = SesquilinearForm{T}[]
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
    orthogonal_sign(G::MatrixGroup{T}) where T <: FinFieldElem

For absolutely irreducible `G` of degree `n`, return
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
function orthogonal_sign(G::MatrixGroup{T}) where T <: FinFieldElem
    R = base_ring(G)
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
    isometry_group(f::SesquilinearForm{T}) where T
    isometry_group(f::QuadraticForm{T}) where T

Return the group of isometries for `f`.

# Examples
```jldoctest
julia> G = symplectic_group(4, 2);

julia> f = alternating_form(invariant_alternating_forms(G)[1]);

julia> isometry_group(f) == G
true
```
"""
function isometry_group(f::Union{SesquilinearForm{T}, QuadraticForm{T}}) where T
   B = gram_matrix(f)
   n = nrows(B)
   F = base_ring(B)
   r=n

   if f isa QuadraticForm
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
      if f isa QuadraticForm
         fn = QuadraticForm(C[1:r, 1:r])
      else
         fn = SesquilinearForm(C[1:r, 1:r],f.descr)
      end
   else
      fn = f
   end

   e=0
   if (fn isa QuadraticForm || fn.descr==:symmetric) && iseven(r)
      if witt_index(fn)== div(r,2)
         e = 1
      else
         e = -1
      end
   end

   if fn isa QuadraticForm
      Xf = is_congruent(QuadraticForm(preimage_matrix(_ring_iso(fn), _standard_form(:quadratic, e, r, F))), fn)[2]
      if Xf === nothing && isodd(r)
         Xf = is_congruent(QuadraticForm(primitive_element(F)*preimage_matrix(_ring_iso(fn), _standard_form(:quadratic, e, r, F))), fn)[2]
      end
   else
      Xf = is_congruent(SesquilinearForm(preimage_matrix(_ring_iso(fn), _standard_form(fn.descr, e, r, F)), fn.descr), fn)[2]
# if dimension is odd, fn may be congruent to a scalar multiple of the standard form
# TODO: I don't really need a primitive_element(F); I just need a non-square in F. Is there a faster way to get it?
      if Xf === nothing && isodd(r)
         Xf = is_congruent(SesquilinearForm(primitive_element(F)*preimage_matrix(_ring_iso(fn), _standard_form(fn.descr, e, r, F)), fn.descr), fn)[2]
      end
   end

   if f isa SesquilinearForm && f.descr === :hermitian
      G = GU(r,Int(characteristic(F)^div(degree(F),2)))
   elseif f isa SesquilinearForm && f.descr === :alternating
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
      algorithm::Symbol=:direct
      direct::Bool=true,
      depth::Int=-1,
      bacher_depth::Int=0,
      set_nice_mono::Bool=true,
     ) -> MatrixGroup

Given an integer lattice $L$ which is definite or of rank 2, return the
isometry group $O(L)$ of $L$.

One can choose which algorithm to use to compute $O(L)$. For now, we
support the following algorithms:
- `:direct`: compute generators of $O(L)$ using an algorithm of
  Plesken--Souvignier;
- `:decomposition`: compute iteratively $O(L)$ by decomposing $L$ into
  $O(L)$-stable sublattices;
- `:default`: compute $O(L)$ use a heuristic to choose the best algorithm 
  between `:direct` and `:decomposition`

Setting the parameters `depth` and `bacher_depth` to a positive value may
improve performance. If set to `-1` (default), the used value of `depth`
is chosen heuristically depending on the rank of `L`. By default,
`bacher_depth` is set to `0`.
"""
@attr QQMatrixGroup function isometry_group(
  L::ZZLat;
  algorithm::Symbol=:default,
  direct::Bool=true,
  depth::Int=-1,
  bacher_depth::Int=0,
  set_nice_mono::Bool=true
)
  # G is represented w.r.t the basis of L
  G = Hecke._assert_has_automorphisms_ZZLat(L; algorithm, depth, bacher_depth, set_nice_mono)
  Gamb = extend_to_ambient_space(L, G; check=false)
  if set_nice_mono
    _sv = _short_vector_generators(L)
    sv = [i*basis_matrix(L) for i in _sv]
    append!(sv, [-i for i in sv])
    _set_nice_monomorphism!(Gamb, sv)
  end 
  return Gamb
end

# We overwrite the function in Hecke in order that 
# Hecke has access to the better algorithms in Oscar
function Hecke._assert_has_automorphisms_ZZLat(L::ZZLat; 
                                               algorithm=:default, 
                                               depth::Int=-1, 
                                               bacher_depth::Int=0,
                                               redo::Bool=false,  
                                               set_nice_mono::Bool=true
)
  # look in the cache
  if !redo && isdefined(L, :automorphism_group_generators)
    _gens = L.automorphism_group_generators
    G = matrix_group(_gens)
    if set_nice_mono
      Llll = lll(L)
      gram = gram_matrix(Llll)
      diag_gram = abs.(diagonal(gram))
      ma = maximum(diag_gram)
      mi = minimum(diag_gram)
      sv = first.(short_vectors(L, ma))
      append!(sv, [-i for i in sv])
      _set_nice_monomorphism!(G, sv)
    end
    return G
    end

  # corner cases
  @req rank(L) <= 2 || is_definite(L) "Lattice must be definite or of rank at most 2"
  if rank(L) <= 2
    Hecke.__assert_has_automorphisms(L; depth, bacher_depth, redo)
    _gens = L.automorphism_group_generators
    return matrix_group(_gens)
  end
  
  if algorithm == :default
    # Select an algorithm
    Llll = lll(L)
    G = gram_matrix(Llll)
    diagG = diagonal(G)
    ma = maximum(diagG)
    mi = minimum(diagG)
    r = rank(L)
    if (r < 5 && ma <10*mi) || (r < 8 && ma < 4*mi) || (r < 12 && ma < 3*mi)|| (ma < 2*mi)
      algorithm = :direct
    else
      algorithm = :decomposition
    end
  end 
  if algorithm == :direct
    Hecke.__assert_has_automorphisms(L; depth, bacher_depth, redo)
    _gens = L.automorphism_group_generators
    G = matrix_group(_gens)
    if set_nice_mono
      sv = first.(short_vectors(L, ma))
      append!(sv, [-i for i in sv])
      _set_nice_monomorphism!(G, sv)
    end
  elseif algorithm == :decomposition
    G, _ = _isometry_group_via_decomposition(L; depth=depth, bacher_depth=bacher_depth)
    # fill the cache
    L.automorphism_group_order = order(G)
    L.automorphism_group_generators = ZZMatrix[change_base_ring(ZZ,matrix(g)) for g in gens(G)]
  else
    error("Unknown algorithm: for the moment, we only support :direct, :decomposition and :default")
  end
  return G
end

"""
    _isometry_group_via_decomposition(
      L::ZZLat;
      depth::Int=-1,
      bacher_depth::Int=0,
      use_heuristics::Bool=false
      direct::Bool=true,
      set_nice_mono::Bool=true,
    ) -> MatrixGroup, Vector{QQMatrix}

Compute the group of isometries of the definite lattice `L` using an orthogonal
decomposition.
"""
function _isometry_group_via_decomposition(
  L::ZZLat;
  depth::Int=-1,
  bacher_depth::Int=0,
  use_heuristics::Bool=false,
  direct::Bool=true,
  set_nice_mono::Bool=true
)
  # TODO: adapt the decomposition approach for AbstractLat
  # in most examples `direct=true` seems to be faster by a factor of 7
  # but in some examples it is also slower ... up to a factor of 15
  
  L = lattice(rational_span(L))
  if gram_matrix(L)[1,1] < 0
    L = rescale(L, -1)
  end

  if !direct
    # need an integral lattice to work with discriminant groups
    d = denominator(scale(L))
    if d > 1
      L = rescale(L, d)
    end
  end

  # construct the sublattice M1 of L generated by the shortest vectors
  V = ambient_space(L)
  # for simplicity we work with the ambient representation
  M1, M1primitive, sv1 = _shortest_vector_primitive_sublattice(L)
  basisM1prim = ZZ.(basis_matrix(M1primitive))
  # basically doubles the memory usage of this function
  # a more elegant way could be to work with the corresponding projective representation
  append!(sv1, [-v for v in sv1]) # given in the coordinates of L
    
  # prepare to select algorithm
  GM = gram_matrix(M1primitive)
  diagGM = diagonal(GM)
  ma = maximum(diagGM)
  mi = minimum(diagGM)
  if ma < 2*mi
    # compute O(M1primitive) directly in Hecke
    @vprintln :Isometry 3 "Computing orthogonal group in Hecke"
    Hecke.__assert_has_automorphisms(M1primitive; depth, bacher_depth) # avoid an infinite recursion
    O1 = matrix_group(M1primitive.automorphism_group_generators)
  else
    # first compute O(M1) then the stabiliser of M1primitive
    @vprintln :Isometry 3 "Computing orthogonal group of shortest sublattice in Hecke"
    Hecke.__assert_has_automorphisms(M1; depth, bacher_depth) # avoid an infinite recursion
    OM1 = matrix_group(M1.automorphism_group_generators)
    if M1primitive == M1
      _O1 = OM1
    else
      # M1 is generated by its shortest vectors only up to finite index
      @vprintln :Isometry 3 "Computing overlattice stabilizers \n"
      #O1,_ = stabilizer(OM1, M1primitive, on_lattices) # needs ambient representation
      @time _O1, _ = _overlattice_stabilizer(OM1, M1, M1primitive)
    end
    # transform to the basis of M1prim
    T = coordinates(basis_matrix(M1primitive), M1)
    invT = inv(T)
    O1 = matrix_group([ZZ.(T*matrix(g)*invT) for g in gens(_O1)])
    @hassert :Isometry 2 is_isometry_group(M1primitive, O1, false)
  end

  if rank(M1) == rank(L)
    B = basisM1prim 
    Binv = inv(B) 
    O1 = matrix_group([Binv*matrix(i)*B for i in gens(O1)])
    if set_nice_mono
      _set_nice_monomorphism!(O1, sv1)
    end
    @hassert :Isometry 2 is_isometry_group(L, O1, false)
    return O1, sv1
  end
  
  # decompose as a primitive extension: M1primitive + M2 \subseteq L
  M2 = orthogonal_submodule(L, M1)

  # recursion
  @vprintln :Isometry 3 "Recursion\n"
  O2, sv2 = _isometry_group_via_decomposition(M2; direct, depth, bacher_depth, set_nice_mono)

  # go to to the basis of L
  basisM2 = ZZ.(basis_matrix(M2))
  sv = append!(sv1, [i*basisM2 for i in sv2])
  BB = vcat(basisM1prim,basisM2)
  
  # In what follows we compute the stabilizer of L in O1 x O2
  if direct
    # Cook up O1 x O2
    r1 = rank(M1primitive)
    r2 = rank(M2) 
    I1 = identity_matrix(ZZ,r1)
    I2 = identity_matrix(ZZ,r2)
    gens12 = [block_diagonal_matrix(ZZMatrix[matrix(i),I2]) for i in gens(O1)] 
    append!(gens12, [block_diagonal_matrix(ZZMatrix[I1,matrix(i)]) for i in gens(O2)]) 
    O1timesO2 = matrix_group(gens12)
    B = vcat(basis_matrix(M1primitive),basis_matrix(M2))
    SL = lattice_in_same_ambient_space(L, B)
    _S,_ = _overlattice_stabilizer(O1timesO2, SL, L)
    BBs = solve_init(BB)
    # transform _S to the basis of L 
    _gens = [solve(BBs,matrix(i)*BB;side=:right) for i in gens(_S)]
    S = matrix_group(_gens)
    
    #=
    O1amb = extend_to_ambient_space(M1primitive, O1;check=false)
    O2amb = extend_to_ambient_space(M2 ,O2;check=false)
    gens12 = vcat(gens(O1amb),gens(O2amb))
    G = matrix_group(gens12)
    (_S, _) = stabilizer(G, L, on_lattices)
    S = matrix_group([ZZ.(matrix(i)) for i in gens(_S)])
    =#
  else
    gensS = stabilizer_in_diagonal_action(L, M1primitive, M2, O1, O2; check=false, is_finite_known=(true, true))
    S = matrix_group(gensS)
  end

  if set_nice_mono
    _set_nice_monomorphism!(S, sv)
  end
  @hassert :Isometry 2 is_isometry_group(L, S; ambient_representation=false)
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
    _set_nice_monomorphism!(G::MatrixGroup, short_vectors)

Use the permutation action of `G` on `short_vectors` to represent `G` as a
finite permutation group.

Internally this sets a `NiceMonomorphism` for the underlying gap group.
No input checks whatsoever are performed.

It is assumed that the corresponding action homomorphism is injective.
"""
#
function _set_nice_monomorphism!(G::MatrixGroup, short_vectors; julia::Bool=true)
  short_vectors = ZZMatrix[matrix(ZZ,1, degree(G), i) for i in short_vectors]
  if julia 
    phi = _nice_hom(G, short_vectors)
  else 
    phi = action_homomorphism(gset(G, on_vector, short_vectors; closed=true))
  end
  GAP.Globals.SetIsInjective(phi.map, true) # fixes an infinite recursion
  GAP.Globals.SetIsHandledByNiceMonomorphism(GapObj(G), true)
  GAP.Globals.SetNiceMonomorphism(GapObj(G), phi.map)
end

_set_nice_monomorphism(G, _short_vectors) = _set_nice_monomorphism(G, [matrix(ZZ,1,d,i) for i in short_vectors])
# return X as permutation on X
# assumes that X is sorted
function _as_perm(w, g::ZZMatrix, X::Vector{ZZMatrix})
  n = length(X)
  per = Vector{Int}(undef, n)
  i = 0
  for x in X
    i +=1
    mul!(w, x, g)
    j = searchsortedfirst(X, w,lt=Hecke._isless)
    per[i] = j
  end 
  return per
end 

# inefficient conversion...
_as_perm(w, g::QQMatrix, X::Vector{ZZMatrix}) = _as_perm(w, ZZ.(g), X)


function _nice_hom(G, _short_vectors::Vector{ZZMatrix})
  _short_vectors = copy(_short_vectors)
  sort!(_short_vectors, lt=Hecke._isless);
  w = similar(_short_vectors[1])
  n = length(_short_vectors)
  Sn = symmetric_group(n)
  act_func(g) = perm(Sn, _as_perm(w, matrix(g), _short_vectors))
  return hom(G, Sn, act_func)
end 

# stabilizer of L in G < O(L)
function _overlattice_stabilizer(G::MatrixGroup{ZZRingElem,ZZMatrix}, S::ZZLat, L::ZZLat)
  _BL = coordinates(basis_matrix(L),S)
  n = denominator(_BL) 
  BL = ZZ.(n*_BL)
  R,iR = residue_ring(ZZ, Int(n))
  BLmod = change_base_ring(R, BL)
  howell_form!(BLmod)
  tmp = similar(BLmod)
  return stabilizer(G, BLmod, on_howell_form)
end 

function on_howell_form(M::zzModMatrix, g::MatrixGroupElem{ZZRingElem,ZZMatrix})
  Mg = M*matrix(g)
  howell_form!(Mg)
  return Mg
end 

function _row_span!(L::Vector{Vector{ZZRingElem}})
  l = length(L)
  d = length(L[1])
  m = min(2*d,l)
  B = sparse_matrix(matrix(ZZ, L))
  h = matrix(hnf(B; truncate=true))
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

function _shortest_vector_primitive_sublattice(L::ZZLat)
  if rank(L) == 0 
    return L
  end
  V = ambient_space(L)
  sv = shortest_vectors(L)
  h = _row_span!(sv)*basis_matrix(L)
  M1 = lattice(V, h)
  M1primitive = primitive_closure(L, M1)
  return lll(M1), lll(M1primitive), sv
end 

function _short_vector_generators(L::ZZLat)
  sv = shortest_vectors(L)
  B = _row_span!(sv)*basis_matrix(L) 
  if nrows(B) == rank(L)
    return sv
  end
  M = orthogonal_submodule(L, B)
  svM = _short_vector_generators(M)
  T = ZZ.(coordinates(basis_matrix(M), L))
  append!(sv, [i*T for i in svM])
  return sv
end 

function Hecke._is_isometric_with_isometry_definite(L1::ZZLat, 
                                                    L2::ZZLat; 
                                                    kwargs...)
   return _is_isometric_with_isometry_definite_via_decomposition(L1, L2; kwargs...)
end

function _is_isometric_with_isometry_definite_via_decomposition(L1::ZZLat, L2::ZZLat; depth::Int = -1, bacher_depth::Int = 0)
  L1 = lattice(rational_span(L1))
  L2 = lattice(rational_span(L2))
  
  # algorithm selection
  
  L1lll = lll(L1)
  G = gram_matrix(L1lll)
  diagG = diagonal(G)
  ma = maximum(diagG)
  mi = minimum(diagG)
  r = rank(L1)
  if r < 4 || (r < 5 && ma <50*mi) || (r < 8 && ma < 4*mi) || (r < 12 && ma < 3*mi)|| (ma < 2*mi)
    b,fL = Hecke.__is_isometric_with_isometry_definite(L1, L2; depth, bacher_depth)
    @hassert :Isometry 3 fL*gram_matrix(L2)*transpose(fL) == gram_matrix(L1)
    return b, fL
  end
  
  # todo: if degree > rank go to lattice(rational_span) and transform back?
  rank(L1) != rank(L2) && return false, zero_matrix(QQ, 0, 0)
  # TODO: compute short vectors only once
  # TODO: start with the non-primitive one
  M1s, M1 = _shortest_vector_primitive_sublattice(L1)
  M2s, M2 = _shortest_vector_primitive_sublattice(L2)
  
  # algorithm selection
  G = gram_matrix(M1)
  diagG = diagonal(G)
  mi = minimum(diagG)
  ma = maximum(diagG)
  if ma >= 2*mi
    b, fMs = Hecke.__is_isometric_with_isometry_definite(M1s, M2s;  depth, bacher_depth)
    b || return false, zero_matrix(QQ, 0, 0)
    OM1s,_ = _isometry_group_via_decomposition(M1s)
    # make it work for now. Go through hecke later?
    @vprint :Isometry 2 "computing orbit of an overlattice..."
    DM1s = discriminant_group(M1s)
    dM1s = discriminant_representation(M1s, OM1s; check=false, full=false, ambient_representation=false)
    OM1sq = image(dM1s)[1]
    to_gapM1s = get_attribute(OM1sq,:to_gap)
    B1 = basis_matrix(M1)
    J1gap = sub(codomain(to_gapM1s), [to_gapM1s(DM1s(B1[i,:])) for i in 1:nrows(B1)])[1]
    fMsinvB = inv(fMs)*basis_matrix(M1s)
    B2 = basis_matrix(M2)
    J2gap = sub(codomain(to_gapM1s), [to_gapM1s(DM1s(coordinates(B2[i,:],M2s)*fMsinvB)) for i in 1:nrows(B2)])[1]
    XM = orbit(OM1sq, on_subgroups, J1gap)
    b, tmp = is_conjugate_with_data(XM, J1gap, J2gap)
    b || return false, zero_matrix(QQ, 0, 0)
    @vprintln :Isometry 2 "done"
    # transform back to the bases of M1 and M2
    T1 = solve(basis_matrix(M1), basis_matrix(M1s); side=:left)
    T2 = solve(basis_matrix(M2), basis_matrix(M2s); side=:left)
    fM = inv(T1)*matrix(dM1s\tmp)*fMs*T2
  else 
    b, fM = Hecke.__is_isometric_with_isometry_definite(M1, M2;  depth, bacher_depth)
    b || return false, zero_matrix(QQ, 0, 0)
  end
  @hassert :Isometry 3 fM*gram_matrix(M2)*transpose(fM) == gram_matrix(M1)
  
  # base case of the recursion 
  if rank(M1) == rank(L1)
    @assert degree(M1)==rank(M1)
    @assert degree(M2)==rank(M2)
    f = inv(basis_matrix(M1))*fM*basis_matrix(M2)
    @hassert :Isometry 1 lattice_in_same_ambient_space(L2, basis_matrix(L1)*f) == L2
    @hassert :Isometry 1 gram_matrix(ambient_space(L2), f) == gram_matrix(ambient_space(L1))
    return b, f
  end
  # recurse 
  N1 = orthogonal_submodule(L1, M1)
  N2 = orthogonal_submodule(L2, M2)
  b, fN = is_isometric_with_isometry(N1, N2;  depth, bacher_depth)
  b || return false, zero_matrix(QQ, 0, 0)
  @hassert :Isometry 3 fN*gram_matrix(N2)*transpose(fN) == gram_matrix(N1)
  #b, fN = _is_isometric_via_decomposition(N1, N2)
  

  
  # modify (fM,fN) so that it extends to L
  @vtime :Isometry 4 OM1,_ = _isometry_group_via_decomposition(M1)
  @vtime :Isometry 4 ON1,_ = _isometry_group_via_decomposition(N1)
  
  DM1 = discriminant_group(M1)
  dM1 = discriminant_representation(M1, OM1; check=false, full=false, ambient_representation=false)
  OM1q = image(dM1)[1]
  DN1 = discriminant_group(N1)
  dN1 = discriminant_representation(N1, ON1; check=false, full=false, ambient_representation=false)
  ON1q = image(dN1)[1]
  
  phi1, iHM1, iHN1 = glue_map(L1, M1,N1;_snf=true, check=false)
  phi2, iHM2, iHN2 = glue_map(L2, M2,N2;_snf=true, check=false)
  
  # TODO: early out if direct sum. 
  
  # Modify (fM,fN) so that fM(HM1) = HM2 and fN(HN1) = HN2
  HM1 = domain(phi1)
  HN1 = codomain(phi1)
  HM2 = domain(phi2)
  HN2 = codomain(phi2)
  
  to_gapM1 = get_attribute(OM1q,:to_gap)
  #to_oscarM = get_attribute(OM1q,:to_oscar)
  HM1gap = sub(codomain(to_gapM1), [to_gapM1(iHM1(i)) for i in gens(HM1)])[1]
  fMinvB = inv(fM)*basis_matrix(M1)
  HM2gap = sub(codomain(to_gapM1), [to_gapM1(DM1(coordinates(lift(i),M2)*fMinvB)) for i in gens(HM2)])[1]
  @vtime :Isometry 4 XM = orbit(OM1q, on_subgroups, HM1gap)
  @vtime :Isometry 4 (b, tmp) = is_conjugate_with_data(XM, HM1gap, HM2gap)
  b || return false, zero_matrix(QQ, 0, 0)
  
  fM = matrix(dM1\tmp)*fM
  fMB = fM*basis_matrix(M2)
  # confirm result
  @hassert :Isometry 3 all(coordinates(lift(i),M1)*fMB in cover(HM2) for i in gens(HM1))
  @hassert :Isometry 3 fM*gram_matrix(M2)*transpose(fM) == gram_matrix(M1)

  to_gapN1 = get_attribute(ON1q,:to_gap)
  #to_oscarN = get_attribute(ON1q,:to_oscar)
  HN1gap = sub(codomain(to_gapN1), [to_gapN1(iHN1(i)) for i in gens(HN1)])[1]
  fNinvB = inv(fN)*basis_matrix(N1)
  HN2gap = sub(codomain(to_gapN1), [to_gapN1(DN1(coordinates(lift(i),N2)*fNinvB)) for i in gens(HN2)])[1]
  XN = orbit(ON1q, on_subgroups, HN1gap)
  
  b, tmp = is_conjugate_with_data(XN, HN1gap, HN2gap)
  b || return false, zero_matrix(QQ, 0, 0)
  # fN is represented w.r.t. the bases of N1 and N2
  fN = matrix(dN1\tmp)*fN
  @hassert :Isometry 3 fN*gram_matrix(N2)*transpose(fN) == gram_matrix(N1)
  fNB = fN*basis_matrix(N2)
  # confirm result
  @hassert :Isometry 3 all(coordinates(lift(i),N1)*fNB in cover(HN2) for i in gens(HN1))
  
  fHM1_HM2 = hom(HM1, HM2, [HM2(coordinates(lift(i), M1)*fMB) for i in gens(HM1)])
  fHN1_HN2 = hom(HN1, HN2, [HN2(coordinates(lift(i), N1)*fNB) for i in gens(HN1)])
  # go in a square and obtain an isometry of HN1 
  # if the square commutes, g is trivial 
  # thus we check if we can make g trivial by modifying f
  _g = fHM1_HM2*phi2*inv(fHN1_HN2)*inv(phi1)
  OHM1 = orthogonal_group(HM1) 
  g = OHM1(_g)
  
  # setup the induced orthogonal groups
  stabHM1,inc_stabHM1 = stabilizer(OM1q, HM1gap, on_subgroups)
  stabHN1,inc_stabHN1 = stabilizer(ON1q, HN1gap, on_subgroups)
  stabHM1_on_HM1,res_stabM = restrict_automorphism_group(stabHM1, iHM1; check=false)
  stabHN1_on_HN1,res_stabN= restrict_automorphism_group(stabHN1, iHN1;check=false)
  phi_stabHN1_on_HN1_phiinv,_ = sub(OHM1, elem_type(OHM1)[OHM1(phi1*hom(i)*inv(phi1)) for i in gens(stabHN1_on_HN1)])
  
  
  # notation D:= UgV \subseteq G
  # find u,v with 1 = u g v
  G = OHM1
  U = stabHM1_on_HM1
  V = phi_stabHN1_on_HN1_phiinv
  
  # permutation groups are faster
  GtoGperm = isomorphism(PermGroup, G) # there are other ways to get a permutation group
  Gperm = codomain(GtoGperm)
  Uperm = GtoGperm(U)[1]
  Vperm = GtoGperm(V)[1]
  g1perm = one(Gperm)
  gperm = GtoGperm(g)
  D = double_coset(Uperm, gperm, Vperm)
  b, uperm, vperm = _decompose(D, g1perm)
  b || return false, zero_matrix(QQ, 0, 0)
  v = GtoGperm\vperm 
  u = GtoGperm\uperm
  @hassert :Isometry 3 u*v == g
  
  gMbar = inv(u) 
  gNbar = inv(v)
  @hassert :Isometry 3 one(G) == gMbar * g * gNbar 
  gNbar = stabHN1_on_HN1(inv(phi1)*hom(v)*phi1)
  f1 = inc_stabHM1\(res_stabM\gMbar)
  f2 = inc_stabHN1\(res_stabN\gNbar)
    
  fM = (dM1\f1)*fM
  fN = (dN1\f2)*fN
  @hassert :Isometry 3 fM*gram_matrix(M2)*transpose(fM) == gram_matrix(M1)
  @hassert :Isometry 3 fN*gram_matrix(N2)*transpose(fN) == gram_matrix(N1)
  # assemble the isometry 
  
  B1 = vcat(basis_matrix(M1),basis_matrix(N1))
  B2 = vcat(basis_matrix(M2),basis_matrix(N2))
  f = solve(B1, diagonal_matrix([fM,fN]); side=:right)*B2
  @hassert :Isometry 3 lattice_in_same_ambient_space(L2, basis_matrix(M1)*f) == M2
  @hassert :Isometry 3 lattice_in_same_ambient_space(L2, basis_matrix(N1)*f) == N2

  @hassert :Isometry 1 lattice_in_same_ambient_space(L2, basis_matrix(L1)*f) == L2
  @hassert :Isometry 1 gram_matrix(ambient_space(L2), f) == gram_matrix(ambient_space(L1))
  return b, f
end
