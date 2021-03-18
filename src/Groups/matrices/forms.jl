# TODO: in this file are used many methods TEMPORARILY defined in files matrix_manipulation.jl and stuff_field_gen.jl
# once methods in those files will be deleted / replaced / modified, this file need to be modified too



import AbstractAlgebra: FieldElem, PolyElem, MPolyElem
import Hecke: base_ring, defining_polynomial, gram_matrix, radical

export
    alternating_form,
    corresponding_bilinear_form,
    corresponding_quadratic_form,
    hermitian_form,
    isalternating_form,
    isdegenerate,
    ishermitian_form,
    isometry_group,
    isquadratic_form,
    issingular,
    issymmetric_form,
    preserved_quadratic_forms,
    preserved_sesquilinear_forms,
    quadratic_form,
    SesquilinearForm,
    symmetric_form,
    witt_index



# descr is always defined
# matrix is always defined except when descr="quadratic"; in such a case, at least one of matrix and pol is defined
# NOTE: the fields ring_iso and mat_iso are always defined if the field X is
"""
    SesquilinearForm{T<:RingElem}

Type of groups `G` of `n x n` matrices over the ring `R`, where `n = degree(G)` and `R = base_ring(G)`.
At the moment, only rings of type `FqNmodFiniteField` are supported.
"""
mutable struct SesquilinearForm{T<:RingElem}
   matrix::MatElem{T}
   descr::Symbol       # quadratic, symmetric, alternating or hermitian
   pol::MPolyElem{T}     # only for quadratic forms
   X::GapObj
   mat_iso::GenMatIso

   function SesquilinearForm{T}(B::MatElem{T},sym) where T
      if sym==:hermitian
         @assert ishermitian_matrix(B) "The matrix is not hermitian"
      elseif sym==:symmetric
         @assert issymmetric(B) "The matrix is not symmetric"
      elseif sym==:alternating
         @assert isskewsymmetric_matrix(B) "The matrix is not skew-symmetric"
      elseif sym != :quadratic
         error("Unsupported description")
      end

      if sym==:quadratic
         return new{T}(_upper_triangular_version(B),sym)
      else
         return new{T}(B,sym)
      end
   end

   function SesquilinearForm{T}(f::MPolyElem{T},sym) where T
      @assert sym==:quadratic "Only quadratic forms are described by polynomials"
      @assert Set([total_degree(x) for x in monomials(f)])==Set(2) "The polynomials is not homoneneous of degree 2"
      r = new{T}()
      r.pol = f
      r.descr = :quadratic
      return r
   end
end

SesquilinearForm(B::MatElem{T}, sym) where T = SesquilinearForm{T}(B,sym)
SesquilinearForm(f::MPolyElem{T},sym) where T = SesquilinearForm{T}(f,sym)



########################################################################
#
# Properties
#
########################################################################


"""
    isalternating_form(f::SesquilinearForm)

Return whether the form `f` is an alternating form.
"""
isalternating_form(f::SesquilinearForm) = f.descr==:alternating

"""
    ishermitian_form(f::SesquilinearForm)

Return whether the form `f` is a hermitian form.
"""
ishermitian_form(f::SesquilinearForm) = f.descr==:hermitian

"""
    isquadratic_form(f::SesquilinearForm)

Return whether the form `f` is a quadratic form.
"""
isquadratic_form(f::SesquilinearForm) = f.descr==:quadratic

"""
    issymmetric_form(f::SesquilinearForm)

Return whether the form `f` is a symmetric form.
"""
issymmetric_form(f::SesquilinearForm) = f.descr==:symmetric

"""
    preserved_quadratic_forms(G::MatrixGroup)

Return a generating set for the vector space of quadratic forms preserved by `G`.
!!! warning "Note:"
    The process involves random procedures, so the function may return different outputs every time.
"""
function preserved_quadratic_forms(G::MatrixGroup)
   L = GAP.Globals.PreservedQuadraticForms(G.X)
   return [G.mat_iso(GAP.Globals.GramMatrix(L[i])) for i in 1:length(L)]
end

"""
    preserved_sesquilinear_forms(G::MatrixGroup)

Return a generating set for the vector space of sesquilinear forms preserved by `G`.
!!! warning "Note:"
    The process involves random procedures, so the function may return different outputs every time.
"""
function preserved_sesquilinear_forms(G::MatrixGroup)
   L = GAP.Globals.PreservedSesquilinearForms(G.X)
   return [G.mat_iso(GAP.Globals.GramMatrix(L[i])) for i in 1:length(L)]
end




########################################################################
#
# Constructors
#
########################################################################


"""
    alternating_form(B::MatElem{T})

Return the alternating form with Gram matrix `B`.
"""
alternating_form(B::MatElem{T}) where T <: FieldElem = SesquilinearForm(B, :alternating)

"""
    symmetric_form(B::MatElem{T})

Return the symmetric form with Gram matrix `B`.
"""
symmetric_form(B::MatElem{T}) where T <: FieldElem = SesquilinearForm(B, :symmetric)

"""
    hermitian_form(B::MatElem{T})

Return the hermitian form with Gram matrix `B`.
"""
hermitian_form(B::MatElem{T}) where T <: FieldElem = SesquilinearForm(B, :hermitian)

# turns the matrix of a quadratic form into an upper triangular matrix of the same form
# (two matrices A,B represent the same quadratic form iff A-B is skew-symmetric)
function _upper_triangular_version(C::MatElem)
   B = deepcopy(C)
   for i in 1:nrows(B), j in i+1:nrows(B)
      B[i,j]+=B[j,i]
      B[j,i]=0
   end
   return B
end

"""
    quadratic_form(B::MatElem{T})

Return the quadratic form with Gram matrix `B`.
"""
quadratic_form(B::MatElem{T}) where T <: FieldElem = SesquilinearForm(B, :quadratic)

"""
    quadratic_form(f::MPolyElem{T}; check=true)

Return the quadratic form described by the polynomial `f`. Here, `f` must be a homogeneous polynomial of degree 2. If `check` is set as `false`, it does not check whether the polynomial is homogeneous of degree 2.
To define quadratic forms of dimension 1, `f` can also have type `PolyElem{T}`.
"""
quadratic_form(f::MPolyElem{T}) where T <: FieldElem = SesquilinearForm(f, :quadratic)
# TODO : neither ishomogeneous or ishomogenous works for variables of type MPolyElem{T}

# just to allow quadratic forms over vector fields of dimension 1, so defined over polynomials in 1 variable
function quadratic_form(f::PolyElem{T}) where T <: FieldElem
   @assert degree(f)==2 && coefficients(f)[0]==0 && coefficients(f)[1]==0 "The polynomials is not homoneneous of degree 2"
   R1 = PolynomialRing(base_ring(f), [string(parent(f).S)])[1]

   return SesquilinearForm(R1[1]^2*coefficients(f)[2], :quadratic)
end

########################################################################
#
# Show
#
########################################################################


function _assign_description(sym::Symbol)
   if sym== :alternating print("Alternating")
   elseif sym== :hermitian print("Hermitian")
   elseif sym== :symmetric print("Symmetric")
   elseif sym== :quadratic print("Quadratic")
   else error("unsupported description")
   end
end


function Base.show(io::IO, f::SesquilinearForm)
   _assign_description(f.descr)
   println(" form with Gram matrix ")
   show(io, "text/plain", gram_matrix(f))
end




########################################################################
#
# Basic
#
########################################################################

#TODO: checking whether two quadratic forms coincide by checking their polynomials is not possible yet.
==(B::SesquilinearForm, C::SesquilinearForm) = gram_matrix(B)==gram_matrix(C) && B.descr==C.descr

function base_ring(B::SesquilinearForm)
   if isdefined(B,:matrix) return base_ring(gram_matrix(B))
   else return base_ring(B.pol)
   end
end

"""
    corresponding_bilinear_form(Q::SesquilinearForm)

Given a quadratic form `Q`, return the bilinear form `B` defined by `B(u,v) = Q(u+v)-Q(u)-Q(v)`.
"""
function corresponding_bilinear_form(B::SesquilinearForm)
   B.descr==:quadratic || throw(ArgumentError("The form must be a quadratic form"))
   M = gram_matrix(B)+transpose(gram_matrix(B))
   if characteristic(base_ring(B))==2 return alternating_form(M)
   else return symmetric_form(M)
   end
end

"""
    corresponding_quadratic_form(Q::SesquilinearForm)

Given a symmetric form `f`, returns the quadratic form `Q` defined by `Q(v) = f(v,v)/2`. It is defined only in odd characteristic.
"""
function corresponding_quadratic_form(B::SesquilinearForm)
   B.descr==:symmetric || throw(ArgumentError("The form must be a symmetric form"))
   characteristic(base_ring(B))!=2 || throw(ArgumentError("Corresponding quadratic form not uniquely determined"))
   M = deepcopy(gram_matrix(B))
   l = inv(base_ring(B)(2))
   for i in 1:nrows(M)
      for j in i+1:nrows(M)
         M[j,i]=0
      end
      M[i,i]*=l
   end
   return quadratic_form(M)
end



########################################################################
#
# Fields of the variable
#
########################################################################

"""
    gram_matrix(B::SesquilinearForm)

Return the Gram matrix of a sesquilinear or quadratic form `B`.
"""
function gram_matrix(f::SesquilinearForm)
   isdefined(f,:matrix) && return f.matrix

   @assert f.descr==:quadratic && isdefined(f,:pol) "Cannot determine Gram matrix"
   d = nvars(parent(f.pol))
   B = zero_matrix( base_ring(f.pol), d, d )
   V = collect(exponent_vectors(f.pol))
   C = collect(coeffs(f.pol))
   for i in 1:length(V)
      x = y = 0
      for j in 1:d
         if V[i][j] !=0
            x = j
            break
         end
      end
      for j in 1:d
         if V[i][d+1-j] !=0
            y = d+1-j
            break
         end
      end
      B[x,y] = C[i]
   end
   f.matrix = B
   return B
end

"""
    defining_polynomial(f::SesquilinearForm)

Return the polynomial that defines the quadratic form `f`.
"""
function defining_polynomial(f::SesquilinearForm)
   isdefined(f,:pol) && return f.pol

   @assert f.descr == :quadratic "Polynomial defined only for quadratic forms"
   R = PolynomialRing(base_ring(f.matrix), nrows(f.matrix) )[1]
   p = zero(R)
   for i in 1:nrows(f.matrix), j in i:nrows(f.matrix)
      p += f.matrix[i,j] * R[i]*R[j]
   end
   f.pol = p
   return p
end


function assign_from_description(f::SesquilinearForm)
   if f.descr==:quadratic f.X=GAP.Globals.QuadraticFormByMatrix(f.mat_iso(gram_matrix(f)),f.mat_iso.fr.codomain)
   elseif f.descr==:symmetric || f.descr==:alternating f.X=GAP.Globals.BilinearFormByMatrix(f.mat_iso(gram_matrix(f)),f.mat_iso.fr.codomain)
   elseif f.descr==:hermitian f.X=GAP.Globals.HermitianFormByMatrix(f.mat_iso(gram_matrix(f)),f.mat_iso.fr.codomain)
   else error("unsupported description")
   end
end


function Base.getproperty(f::SesquilinearForm, sym::Symbol)

   if isdefined(f,sym) return getfield(f,sym) end

   if sym === :mat_iso
      f.mat_iso = gen_mat_iso(nrows(gram_matrix(f)), base_ring(f))

   elseif sym == :X
      if !isdefined(f, :X)
         if !isdefined(f,:mat_iso)
            f.mat_iso = gen_mat_iso(nrows(gram_matrix(f)), base_ring(f))
         end
         assign_from_description(f)
      end

   end

   return getfield(f, sym)

end



########################################################################
#
# Operations
#
########################################################################

function Base.:*(f::SesquilinearForm, l::FieldElem)
   l !=0 || throw(ArgumentError("Zero is not admitted"))
   parent(l)==base_ring(f) || throw(ArgumentError("The scalar does not belong to the base ring of the form"))
   if !isdefined(f,:matrix)
      return SesquilinearForm(l*f.pol, f.descr)
   else
      g = SesquilinearForm(l*gram_matrix(f), f.descr)
      if isdefined(f,:pol) g.pol=f.pol end
      return g
   end
end

Base.:*(l::FieldElem, f::SesquilinearForm) = f*l

function Base.:^(f::SesquilinearForm{T}, x::MatElem{T}; check=false) where T <: RingElem
   @assert base_ring(f)==base_ring(x) "Incompatible base rings"
   @assert nrows(gram_matrix(f))==nrows(x) "Incompatible dimensions"

   if check @assert rank(x)==nrows(x) "Matrix not invertible" end
   if f.descr==:hermitian
      m = x^-1*gram_matrix(f)*conjugate_transpose(x^-1)
   else
      m = x^-1*gram_matrix(f)*transpose(x^-1)
   end
   return SesquilinearForm(m, f.descr)
end

Base.:^(f::SesquilinearForm{T}, x::MatrixGroupElem{T}; check=false) where T <: RingElem = f^x.elm

function (f::SesquilinearForm{T})(v::AbstractAlgebra.Generic.FreeModuleElem{T},w::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: RingElem
   @assert f.descr!=:quadratic "Quadratic forms requires only one argument"

   if f.descr==:hermitian
      return v*gram_matrix(f)*map( y->frobenius(y,div(degree(base_ring(w)),2)),w)
   else
      return v*gram_matrix(f)*w
   end
end

function (f::SesquilinearForm{T})(v::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: RingElem
   @assert f.descr==:quadratic "Sesquilinear forms requires two arguments"
   return v*gram_matrix(f)*v
end




########################################################################
#
# Functionalities
#
########################################################################

"""
    radical(f::SesquilinearForm{T})

Return the radical of the sesquilinear form `f`, i.e. the subspace of all `v` such that `f(u,v)=0` for all `u`. The radical of a quadratic form `Q` is the set of vectors `v` such that `Q(v)=0` and `v` lies in the radical of the corresponding bilinear form.
"""
function radical(f::SesquilinearForm{T}) where T
   V = VectorSpace(base_ring(f), nrows(gram_matrix(f)) )
   R = GAP.Globals.RadicalOfForm(f.X)
   GAP.Globals.Dimension(R)==0 && return sub(V,[])
   L = AbstractAlgebra.Generic.FreeModuleElem{T}[]
   for l in GAP.Globals.GeneratorsOfVectorSpace(R)
      v = V([f.mat_iso.fr(t) for t in l])
      push!(L,v)
   end
   return sub(V,L)
end

"""
    witt_index(f::SesquilinearForm{T})

Return the Witt index of the form induced by `f` on `V/Rad(f)`. The Witt Index is the dimension of a maximal totally isotropic (singular for quadratic forms) subspace.
"""
witt_index(f::SesquilinearForm{T}) where T = GAP.Globals.WittIndex(f.X)

"""
    isdegenerate(f::SesquilinearForm{T})

Return whether `f` is degenerate, i.e. `f` has nonzero radical. A quadratic form is degenerate if the corresponding bilinear form is.
"""
function isdegenerate(f::SesquilinearForm{T}) where T
   f.descr != :quadratic && return det(gram_matrix(f))==0
   return det(gram_matrix(f)+transpose(gram_matrix(f)))==0
end

"""
    issingular(Q::SesquilinearForm{T})

For a quadratic form `Q`, return whether `Q` is singular, i.e. `Q` has nonzero radical.
"""
function issingular(f::SesquilinearForm{T}) where T
   f.descr != :quadratic && throw(ArgumentError("The form is not quadratic"))
   return GAP.Globals.IsSingularForm(f.X)
end




########################################################################
#
# From group to form #TODO: there are different approaches. Which is the best?
#
########################################################################

# TODO 1st approach: brute force calculation
# Algorithm furnished by Thomas Breuer, Aachen University
# extended to quadratic by Giovanni De Franceschi, TU Kaiserslautern
# WARNING: huge linear system !!

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
    invariant_hermitian_forms(G::MatrixGroup)
Return a generating set for the vector spaces of sesquilinear non-bilinear forms preserved by the group `G`. It works only if `base_ring(G)` has even degree.
!!! warning "Note:"
    At the moment, elements of the generating set are returned of type `mat_elem_type(G)`.
"""
function invariant_hermitian_forms(G::MatrixGroup{S,T}) where {S,T}
   F = base_ring(G)
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
    invariant_sesquilinear_forms(G::MatrixGroup)
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
      for i in 1:n, j in i:n
      for p in 1:n, q in p:n
         MM[div((2*n-i+2)*(i-1),2)+1, div((2*n-p+2)*(p-1),2)+q-p+1] = mat[i,p]*mat[i,q]
      end
      for j in i+1:n, p in 1:n, q in p:n
         MM[div((2*n-i+2)*(i-1),2)+j-i+1, div((2*n-p+2)*(p-1),2)+q-p+1] = mat[i,p]*mat[j,q]+mat[j,p]*mat[i,q]
      end
      end
      for i in 1:div(n*(n+1),2) MM[i,i]-=1 end
      push!(M,MM)
   end

   r,K = nullspace(block_matrix(length(M),1,M))
   M = T[]
   for i in 1:r
      push!(M,upper_triangular_matrix([K[j,i] for j in 1:div(n*(n+1),2)]))
   end
   return M
end


# TODO 2nd approach: using MeatAxe GAP functionalities

"""
    function invariant_bilinear_form(G::MatrixGroup)
Return an invariant bilinear form for the group `G`. It works only if the module induced by the action of `G` is absolutely irreducible.
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
Return an invariant sesquilinear (non bilinear) form for the group `G`. It works only if the module induced by the action of `G` is absolutely irreducible.
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
Return an invariant bilinear form for the group `G`. It works only if the module induced by the action of `G` is absolutely irreducible.
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
Uses random methods to find all of the quadratic forms preserved by `G` up to a scalar (i.e. such that `G` is a group of similarities for the forms). Since the procedure relies on a pseudo-random generator, the user may need to execute the operation more than once to find all invariant quadratic forms.
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
Uses random methods to find all of the sesquilinear forms preserved by `G` up to a scalar (i.e. such that `G` is a group of similarities for the forms). Since the procedure relies on a pseudo-random generator, the user may need to execute the operation more than once to find all invariant sesquilinear forms.
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
#      @assert rank(B+transpose(B))==n "At the moment, only nondegenerate quadratic forms are considered"
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
      C,A,r = find_radical(B,F,n,n; e=degF, Symmetric=true)
   end

   if r<n
      fn = SesquilinearForm(submatrix(C,1,1,r,r),f.descr)
   else
      fn = f
   end

   e=0
   if (fn.descr==:quadratic || fn.descr==:symmetric) && iseven(r)
      if witt_index(fn)== div(r,2) e=1
      else e=-1
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
# TODO: not quite sure whether the last element is sufficient to generate the whole top-right block
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
