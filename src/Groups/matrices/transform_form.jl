# TODO: in this file are used many methods TEMPORARILY defined in files matrix_manipulation.jl and stuff_field_gen.jl
# once methods in those files will be deleted / replaced / modified, this file need to be modified too



import AbstractAlgebra: Field

export
    iscongruent

###############################################################################################################

# Algorithm

###############################################################################################################

# The functions _find_radical, _block_anisotropic_elim and _block_herm_elim are based on the article
# James B. Wilson, Optimal Algorithms of Gram-Schmidt type, Linear Algebra and its Applications 438 (2013) 4573-4583

# returns x,y such that ax^2+by^2 = c
# at the moment, it simply search for x such that (c-ax^2)/b is a square
# TODO: is there a faster way to find x and y?
# TODO: it would be better if this is deterministic. This depends on gen(F) and issquare(F).

function _solve_eqn(a::T, b::T, c::T) where T <: FinFieldElem
   F = parent(a)  
   for x in F
      s = (c - a*x^2)*b^-1
      if issquare(s) return x,square_root(s) end
   end
   return nothing, nothing
end



# if _is_symmetric, returns C,A,d where A*B*transpose(frobenius(A,e)) = C, C = block_matrix(2,2,[C,0,0,0]) and d = rank(C)
# else returns C,A,d where B*A = C, C = [C 0] and d = rank(C)
# Assumption: if _is_symmetric==true, then nr=nc always
# Assumption: e = deg(F)/2 in the hermitian case, e=0 otherwise
function _find_radical(B::MatElem{T}, F::Field, nr::Int, nc::Int; e::Int=0, _is_symmetric::Bool=false) where T <: FinFieldElem

   @assert !_is_symmetric || nr==nc
   V1 = VectorSpace(F,nc)
   V2 = VectorSpace(F,nr)
   K, embK = kernel(ModuleHomomorphism(V1, V2, _is_symmetric ? B : transpose(B)))
   U, embU = complement(V1,K)
   d = dim(U)
   elemT = elem_type(V1)
   A = matrix(elemT[embU.(gens(U)) ; embK.(gens(K))])
#   type_vector = elem_type(V1)
#   A = matrix(vcat(type_vector[embU(v) for v in gens(U)], type_vector[embK(v) for v in gens(K)] ))

   if _is_symmetric
      return A*B*transpose(map(y -> frobenius(y,e),A)), A, d
   else
      A = transpose(A)
      return B*A, A, d
   end

end





# returns D, A such that A*B*transpose(frobenius(A)) = D and 
# D is diagonal matrix (or with blocks [0 1 s 0])
# f = dimension of the zero block in B in the isotropic case
function _block_anisotropic_elim(B::MatElem{T}, _type::Symbol; isotr=false, f=0)  where T <: FinFieldElem

   d = nrows(B)
   F = base_ring(B)
   if d <= 1
      return B, identity_matrix(F,d)
   end

   if _type==:symmetric
      degF=0
      s=1
   elseif _type==:alternating
      degF=0
      s=-1
   elseif _type==:hermitian
      degF=div(degree(F),2)
      s=1
   end

   # conjugate transpose in hermitian case
   # transpose in the other cases
   star(X) = transpose(map(y -> frobenius(y,degF),X))

   if isotr
      q = characteristic(F)^degF
      g = d-f
      U = submatrix(B,1,f+1,f,g)
      V = submatrix(B,f+1,f+1,g,g)
      C,A,e = _find_radical(U,F,f,g)
      # I expect C always to be of rank f
      C = submatrix(C,1,1,f,f)                     
      Vprime = star(A)*V*A
      Z = submatrix(Vprime, 1,1,f,f)
      Y = submatrix(Vprime,f+1,1,g-f,f)
      Bprime = submatrix(Vprime,f+1,f+1,g-f,g-f)
      TR = zero_matrix(F,f,f)
      D = zero_matrix(F,f,f)
      for i in 1:f
         D[i,i] = Z[i,i]
         for j in i+1:f  TR[i,j] = Z[i,j] end
      end
      Aarray = MatElem{T}[]  #TODO type?
      Barray = MatElem{T}[]
      for i in 1:f
         alpha = D[i,i]
         if alpha != 0
            push!(Barray, matrix(F,2,2,[-alpha^-1,0,0,alpha]))
            push!(Aarray, matrix(F,2,2,[1,-alpha^-1,0,1]))
         else
            push!(Barray, matrix(F,2,2,[0,1,s,0]))
            push!(Aarray, matrix(F,2,2,[1,0,0,1]))
         end
      end
      B0,A0 = _block_anisotropic_elim(Bprime,_type)
      B1 = diagonal_join(Barray)
      B1 = diagonal_join(B1,B0)
      C = C^-1
      Temp = vcat(C,-TR*C)
      Temp = vcat(Temp,-Y*C)
      Temp = hcat(Temp, vcat(zero_matrix(F,f,g),star(A)))
      P = zero_matrix(F,2*f,2*f)
      for i in 1:f
         P[2*i-1,i] = 1
         P[2*i,i+f] = 1
      end
      A1 = diagonal_join(Aarray)*P
      A1 = diagonal_join(A1,A0)
      return B1, A1*Temp
   else
      c,f = Int(ceil(d/2)), Int(floor(d/2))
      B0 = submatrix(B,1,1,c,c)
      U = submatrix(B,c+1,1,f,c)
      V = submatrix(B,c+1,c+1,f,f)
      B1,A0,e = _find_radical(B0,F,c,c; e=degF, _is_symmetric=true)
      B1 = submatrix(B1,1,1,e,e)
      U = U*star(A0)
      U1 = submatrix(U,1,1,f,e)
      U2 = submatrix(U,1,e+1,f,c-e)
      Z = V-s*U1*B1^-1*star(U1)
      D1,A1 = _block_anisotropic_elim(B1,_type)
      Temp = zero_matrix(F,d-e,d-e)
      insert_block!(Temp,s*star(U2),1,c-e+1)
      insert_block!(Temp,U2,c-e+1,1)
      insert_block!(Temp,Z,c-e+1,c-e+1)
      if c-e==0
         D2,A2 = _block_anisotropic_elim(Temp,_type)
      else
         D2,A2 = _block_anisotropic_elim(Temp, _type; isotr=true, f=c-e)
      end
      Temp = hcat(-U1*B1^-1, zero_matrix(F,f,c-e))*A0
      Temp = vcat(A0,Temp)
      Temp = insert_block(identity_matrix(F,d),Temp,1,1)

      return diagonal_join(D1,D2), diagonal_join(A1,A2)*Temp
   end
end



# assume B is nondegenerate


# returns D, A such that A*B*transpose(frobenius(A)) = D and 
# D is diagonal matrix (or with blocks [0 1 s 0])
# f = dimension of the zero block in B in the isotropic case
function _block_herm_elim(B::MatElem{T}, _type) where T <: FinFieldElem
   d = nrows(B)
   F = base_ring(B)

   if d==1
      return B, identity_matrix(F,1)
   end

   c = Int(ceil(d/2))
   B2 = submatrix(B,1,1,c,c)
   if B2==0
      D,A = _block_anisotropic_elim(B,_type; isotr=true, f=c)
   else
      D,A = _block_anisotropic_elim(B,_type)
   end

   return D,A
end



# returns D such that D*B*conjugatetranspose(D) is the standard basis
# it modifies the basis_change_matrix of the function _block_herm_elim
# TODO: not done for orthogonal

function _to_standard_form(B::MatElem{T}, _type::Symbol)  where T <: FinFieldElem
   F = base_ring(B)
   n = nrows(B)
   A,D = _block_herm_elim(B, _type)

   if _type==:alternating
      our_perm = vcat(1:2:n, reverse(2:2:n))
      D = permutation_matrix(F,our_perm)*D
   elseif _type==:hermitian
      w = primitive_element(F)
      q = Int(sqrt(order(F)))
      Z = identity_matrix(F,n)
      # turn the elements on the main diagonal into 1
      for i in 1:n
         if A[i,i]!=0
            lambda = _disc_log(w^(q+1),A[i,i])
            Z[i,i] = w^-lambda
            A[i,i] = 1
         end
      end
      D = Z*D
      # moving all hyperbolic lines at the end
      Z = identity_matrix(F,n)
      our_permut = Array(1:n)
      NOZ = 0      # Number Of Zeros on the diagonal before A[i,i]
      for i in 1:n
         if A[i,i]==0
            NOZ += 1
         else
            j = i
            while j>i-NOZ
               swap_cols!(A,j,j-1)
               j -= 1
            end
            j = i
            while j > i-NOZ
               swap_rows!(A,j,j-1)
               j -= 1
            end
            for j in 1:NOZ
               our_permut[i+1-j] = our_permut[i-j]
            end
            our_permut[i-NOZ] = i
         end
      end
      Z = permutation_matrix(F, our_permut)
      D = Z*D
      # turn 2x2 identities into 2x2 anti-diagonal blocks
      Z = identity_matrix(F,n)
      if isodd(q)
         b = (1+w^(div((q-1)^2,2)))^-1
         a = b*w^(div((1-q),2))
         d = 1
         c = w^(div((q-1),2))
      else
         b = (1+w^(q-1))^-1
         a = b*w^(q-1)
         d = 1
         c = 1
      end
      S = matrix(F,2,2,[a,b,c,d])
      if div(n-NOZ,2)==0
         S = zero_matrix(F,0,0)
      else
         S = diagonal_join([S for i in 1:div(n-NOZ,2)])
      end
      # turn into standard GAP form
      sec_perm = Int[]
      for i in 1:div(n,2)
         sec_perm = vcat([i,n+1-i],sec_perm)
      end
      if isodd(n)
         insert_block!(Z,S,2,2)
         sec_perm = vcat([div(n+1,2)],sec_perm)
      else
         insert_block!(Z,S,1,1)
      end
      D = transpose(permutation_matrix(F,sec_perm))*Z*D
   end

   return D
end



###############################################################################################################

# Change of basis between two matrices

###############################################################################################################

# modifies A by eliminating all hyperbolic lines and turning A into a diagonal matrix
# return the matrix Z such that Z*A*transpose(Z) is diagonal;
# works only in odd characteristic
function _elim_hyp_lines(A::MatElem{T}) where T <: FinFieldElem
   F = base_ring(A)
   n = nrows(A)
   b = matrix(F,2,2,[1,1,1,-1])  # change of basis from matrix([0,1,1,0]) to matrix([2,0,0,-2])
   Z = identity_matrix(F,n)

   i = 1
   while i <= n
      if A[i,i]==0
         A[i,i]=2
         A[i,i+1]=0
         A[i+1,i]=0
         A[i+1,i+1]=-2
         insert_block!(Z,b,i,i)
         i+=2
      else
         i+=1
      end
   end

   return Z
end


# return true, D such that D*B2*conjugatetranspose(D)=B1
# return false, nothing if D does not exist
# TODO: orthogonal only in odd char, at the moment
function _change_basis_forms(B1::MatElem{T}, B2::MatElem{T}, _type::Symbol)  where T <: FinFieldElem

   if _type==:alternating || _type==:hermitian
      D1 = _to_standard_form(B1,_type)
      D2 = _to_standard_form(B2,_type)
      return true, D1^-1*D2
   elseif _type==:symmetric
      F = base_ring(B1)
      isodd(characteristic(F)) || error("Even characteristic not supported")
      n = nrows(B1)
      A1,D1 = _block_herm_elim(B1, _type)
      A2,D2 = _block_herm_elim(B2, _type)
      q = order(F)
      # eliminate all hyperbolic lines and turn A1,A2 into diagonal matrices
      # TODO: assure that the function _elim_hyp_lines actually modifies A1 and A2
      D1 = _elim_hyp_lines(A1)*D1
      D2 = _elim_hyp_lines(A2)*D2
      issquare( prod(diagonal(A1))*prod(diagonal(A2)) )[1] || return false, nothing
      # move all the squares on the diagonal at the begin
      _squares = [i for i in 1:n if issquare(A1[i,i])[1]]
      our_perm = vcat(_squares, setdiff(1:n, _squares))
      P = permutation_matrix(F,our_perm)
      s1 = length(_squares)
      D1 = P*D1
      A1 = P*A1*transpose(P)
      _squares = [i for i in 1:n if issquare(A2[i,i])[1]]
      our_perm = vcat(_squares, setdiff(1:n, _squares))
      P = permutation_matrix(F,our_perm)
      s2 = length(_squares)
      D2 = P*D2
      A2 = P*A2*transpose(P)
      # get same number of squares on the two diagonals of A1 and A2 by modifying A1
      if s1!=s2
         s = min(s1,s2)+1
         w = A1[s,s]*A2[s,s]     # I'm sure this is not a square
         a,b = _solve_eqn(F(1),F(1),w)
         L = identity_matrix(F,n)
         for i in 0:div(abs(s1-s2),2)-1
            k = s+2*i
            r = square_root(A1[k,k]*A1[k+1,k+1]^-1)
            L[k:k+1,k:k+1] = [a b*r ; b -a*r]
         end
         D1 = L*D1
         A1 = L*A1*transpose(L)
      end
      # change matrix from A1 to A2
      S = diagonal_matrix([square_root(A1[i,i]*A2[i,i]^-1) for i in 1:n])
      return true, D1^-1*S*D2
   end

end

"""
    iscongruent(f::SesquilinearForm{T}, g::SesquilinearForm{T}) where T <: RingElem

If `f` and `g` are sesquilinear forms, return (`true`, `C`) if there exists a matrix `C` such that `f^C = g`, or equivalently, `CBC* = A`, where `A` and `B` are the Gram matrices of `f` and `g` respectively, and `C*` is the transpose-conjugate matrix of `C`. If such `C` does not exist, then return (`false`, `nothing`).
If `f` and `g` are quadratic forms, return (`true`, `C`) if there exists a matrix `C` such that `f^A = ag` for some scalar `a`. If such `C` does not exist, then return (`false`, `nothing`).
"""
function iscongruent(f::SesquilinearForm{T}, g::SesquilinearForm{T}) where T <: RingElem

   @assert base_ring(f)==base_ring(g) "The forms have not the same base ring"
   @assert nrows(gram_matrix(f))==nrows(gram_matrix(g)) "The forms act on vector spaces of different dimensions"
   f.descr==g.descr || return false, nothing
   n = nrows(gram_matrix(f))
   F = base_ring(f)
   
   if f.descr==:quadratic
      if iseven(characteristic(F))            # in this case we use the GAP algorithms
         Bg = g.mat_iso(GAP.Globals.BaseChangeToCanonical(g.X))
         Bf = f.mat_iso(GAP.Globals.BaseChangeToCanonical(f.X))
         UTf = _upper_triangular_version(Bf*gram_matrix(f)*transpose(Bf))
         UTg = _upper_triangular_version(Bg*gram_matrix(g)*transpose(Bg))
         if _is_scalar_multiple_mat(UTf, UTg)[1]
            return true, Bf^-1*Bg
         else
            return false, nothing
         end
      else
         return iscongruent(corresponding_bilinear_form(f), corresponding_bilinear_form(g))
      end
   else
      rank_f = rank(gram_matrix(f))
      rank_f==rank(gram_matrix(g)) || return false, nothing
      if rank_f<n
         degF=0
         if f.descr==:hermitian degF=div(degree(F),2) end
         Cf,Af,d = _find_radical(gram_matrix(f),F,n,n; e=degF, _is_symmetric=true)
         Cg,Ag,_ = _find_radical(gram_matrix(g),F,n,n; e=degF, _is_symmetric=true)
         _is_true, Z = _change_basis_forms( submatrix(Cf,1,1,d,d), submatrix(Cg,1,1,d,d), f.descr)
         _is_true || return false, nothing
         Z = insert_block(identity_matrix(F,n), Z, 1,1)
         return true, Af^-1*Z*Ag
      else
         return _change_basis_forms(gram_matrix(f), gram_matrix(g), f.descr)
      end
   end

   return false, nothing
end
