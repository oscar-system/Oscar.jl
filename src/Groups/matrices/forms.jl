# TODO: in this file are used many methods TEMPORARILY defined in files matrix_manipulation.jl and stuff_field_gen.jl
# once methods in those files will be deleted / replaced / modified, this file need to be modified too




# descr is always defined
# matrix is always defined
# NOTE: the field ring_iso is always defined if the field X is
"""
    SesquilinearForm{T<:RingElem}

Type of alternating and symmetric bilinear forms, hermitian forms,
defined by a Gram matrix.
"""
mutable struct SesquilinearForm{T<:RingElem}
   matrix::MatElem{T}
   descr::Symbol       # symmetric, alternating or hermitian
   X::GapObj
   ring_iso::MapFromFunc

   function SesquilinearForm{T}(B::MatElem{T},sym) where T
      if sym==:hermitian
         @assert is_hermitian(B) "The matrix is not hermitian"
      elseif sym==:symmetric
         @assert is_symmetric(B) "The matrix is not symmetric"
      elseif sym==:alternating
         @assert is_alternating(B) "The matrix does not correspond to an alternating form"
      else
         error("Unsupported description")
      end

      return new{T}(B,sym)
   end
end

SesquilinearForm(B::MatElem{T}, sym) where T = SesquilinearForm{T}(B,sym)


# descr is always defined
# at least one of matrix and pol is defined
# NOTE: the field ring_iso is always defined if the field X is
"""
    QuadraticForm{T<:RingElem}

Type of quadratic forms, defined by a Gram matrix or by a polynomial.
"""
mutable struct QuadraticForm{T<:RingElem}
   matrix::MatElem{T}
   descr::Symbol  # perhaps unnecessary
   pol::MPolyRingElem{T}
   X::GapObj
   ring_iso::MapFromFunc

   function QuadraticForm{T}(B::MatElem{T}) where T
      return new{T}(_upper_triangular_version(B))
   end

   function QuadraticForm{T}(f::MPolyRingElem{T}) where T
      @assert Set([total_degree(x) for x in AbstractAlgebra.monomials(f)])==Set(2) "The polynomials is not homogeneous of degree 2"
      r = new{T}()
      r.pol = f
      return r
   end
end

QuadraticForm(B::MatElem{T}) where T = QuadraticForm{T}(B)
QuadraticForm(f::MPolyRingElem{T}) where T = QuadraticForm{T}(f)


########################################################################
#
# Properties
#
########################################################################


"""
    is_alternating(f::SesquilinearForm)

Return whether `f` has been created as an alternating form.

The Gram matrix `M` of an alternating form satisfies `M = -transpose(M)`
and `M` has zeros on the diagonal,
see [`is_alternating(M::MatrixElem)`](@ref).

# Examples
```jldoctest
julia> M = matrix(GF(2), [0 1; 1 0])
[0   1]
[1   0]

julia> is_alternating(alternating_form(M))
true

julia> is_alternating(symmetric_form(M))
false
```
"""
is_alternating(f::SesquilinearForm) = f.descr==:alternating


"""
    is_hermitian(f::SesquilinearForm)

Return whether `f` has been created as a hermitian form.

The Gram matrix `M` of a hermitian form satisfies
`M = conjugate_transpose(M)`,
see [`is_hermitian(B::MatElem{T}) where T <: FinFieldElem`](@ref).

# Examples
```jldoctest
julia> F, z = finite_field(2, 2); M = matrix(F, [0 z; z+1 0])
[    0   o]
[o + 1   0]

julia> is_hermitian(hermitian_form(M))
true
```
"""
is_hermitian(f::SesquilinearForm) = f.descr==:hermitian

"""
    is_quadratic(f::SesquilinearForm)
    is_quadratic(f::QuadraticForm)

Return whether `f` has been created as a quadratic form.

# Examples
```jldoctest
julia> M = matrix(GF(2), [0 1; 0 0])
[0   1]
[0   0]

julia> is_quadratic(quadratic_form(M))
true
```
"""
is_quadratic(f::SesquilinearForm) = false
is_quadratic(f::QuadraticForm) = true
#T needed?

"""
    is_symmetric(f::SesquilinearForm)

Return whether `f` has been created as a symmetric form.
The Gram matrix `M` of a symmetric form satisfies `M = transpose(M)`,
see [`is_symmetric(M::MatrixElem)`](@ref).

# Examples
```jldoctest
julia> M = matrix(GF(2), [0 1; 1 0])
[0   1]
[1   0]

julia> is_symmetric(symmetric_form(M))
true

julia> is_symmetric(alternating_form(M))
false
```
"""
is_symmetric(f::SesquilinearForm) = f.descr==:symmetric


########################################################################
#
# Constructors
#
########################################################################


"""
    alternating_form(B::MatElem{T})

Return the alternating form with Gram matrix `B`.
An exception is thrown if `B` is not square or does not have zeros on the
diagonal or does not satisfy `B = -transpose(B)`.

# Examples
```jldoctest
julia> f = alternating_form(matrix(GF(3), 2, 2, [0, 1, -1, 0]))
alternating form with Gram matrix
[0   1]
[2   0]

julia> describe(isometry_group(f))
"SL(2,3)"
```
"""
alternating_form(B::MatElem{T}) where T <: FieldElem = SesquilinearForm(B, :alternating)

"""
    symmetric_form(B::MatElem{T})

Return the symmetric form with Gram matrix `B`.
An exception is thrown if `B` is not square or not symmetric.

# Examples
```jldoctest
julia> f = symmetric_form(matrix(GF(3), 2, 2, [0, 1, 1, 0]))
symmetric form with Gram matrix
[0   1]
[1   0]

julia> describe(isometry_group(f))
"C2 x C2"
```
"""
symmetric_form(B::MatElem{T}) where T <: FieldElem = SesquilinearForm(B, :symmetric)

"""
    hermitian_form(B::MatElem{T})

Return the hermitian form with Gram matrix `B`.
An exception is thrown if `B` is not square or does not satisfy
`B = conjugate_transpose(B)`, see [`conjugate_transpose`](@ref).

# Examples
```jldoctest
julia> F = GF(4);  z = gen(F);

julia> f = hermitian_form(matrix(F, 2, 2, [0, z, z^2, 0]))
hermitian form with Gram matrix
[    0   o]
[o + 1   0]

julia> describe(isometry_group(f))
"C3 x S3"
```
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

# Examples
```jldoctest
julia> f = quadratic_form(matrix(GF(3), 2, 2, [2, 2, 0, 1]))
quadratic form with Gram matrix
[2   2]
[0   1]

julia> describe(isometry_group(f))
"D8"
```
"""
quadratic_form(B::MatElem{T}) where T <: FieldElem = QuadraticForm(B)

"""
    quadratic_form(f::MPolyRingElem{T}; check=true)

Return the quadratic form described by the polynomial `f`.
Here, `f` must be a homogeneous polynomial of degree 2.
If `check` is set to `false`, it is not checked whether `f` is homogeneous of degree 2.
To define quadratic forms of dimension 1, `f` can also have type `PolyRingElem{T}`.

# Examples
```jldoctest
julia> _, (x1, x2) = polynomial_ring(GF(3), [:x1, :x2]);

julia> f = quadratic_form(2*x1^2 + 2*x1*x2 + x2^2)
quadratic form with Gram matrix
[2   2]
[0   1]

julia> describe(isometry_group(f))
"D8"
```
"""
quadratic_form(f::MPolyRingElem{T}) where T <: FieldElem = QuadraticForm(f)
# TODO : neither is_homogeneous or is_homogeneous works for variables of type MPolyRingElem{T}

# just to allow quadratic forms over vector spaces of dimension 1, so defined over polynomials in 1 variable
function quadratic_form(f::PolyRingElem{T}) where T <: FieldElem
   @assert degree(f)==2 && coefficients(f)[0]==0 && coefficients(f)[1]==0 "The polynomials is not homogeneous of degree 2"
   R1 = polynomial_ring(base_ring(f), [string(parent(f).S)]; cached=false)[1]

   return QuadraticForm(R1[1]^2*coefficients(f)[2])
end

########################################################################
#
# Show
#
########################################################################

function Base.show(io::IO, f::SesquilinearForm)
   println(io, "$(f.descr) form with Gram matrix")
   show(io, "text/plain", gram_matrix(f))
end

function Base.show(io::IO, f::QuadraticForm)
   println(io, "quadratic form with Gram matrix")
   show(io, "text/plain", gram_matrix(f))
end

########################################################################
#
# Basic
#
########################################################################

#TODO: checking whether two quadratic forms coincide by checking their polynomials is not possible yet.
==(B::SesquilinearForm, C::SesquilinearForm) = B.descr == C.descr && gram_matrix(B) == gram_matrix(C)
==(B::QuadraticForm, C::QuadraticForm) = gram_matrix(B) == gram_matrix(C)

function Base.hash(f::Union{SesquilinearForm, QuadraticForm}, h::UInt)
   b = 0xf64440baac005f8c % UInt
   h = hash(f.descr, h)
   h = hash(gram_matrix(f), h)
   return xor(h, b)
end

base_ring(B::SesquilinearForm) = base_ring(gram_matrix(B))

function base_ring(B::QuadraticForm)
   isdefined(B,:matrix) && return base_ring(gram_matrix(B))
   return base_ring(B.pol)
end

"""
    corresponding_bilinear_form(Q::QuadraticForm)

Return the symmetric form `B` defined by `B(u,v) = Q(u+v)-Q(u)-Q(v)`.

# Examples
```jldoctest; filter = r"[0-9]"
julia> Q = quadratic_form(invariant_quadratic_form(GO(3,3)))
quadratic form with Gram matrix
[0   1   0]
[0   0   0]
[0   0   1]

julia> corresponding_bilinear_form(Q)
symmetric form with Gram matrix
[0   1   0]
[1   0   0]
[0   0   2]
```
"""
function corresponding_bilinear_form(Q::QuadraticForm)
   gram = gram_matrix(Q)
   M = gram + transpose(gram)
   characteristic(base_ring(Q)) == 2 && return alternating_form(M)
   return symmetric_form(M)
end

"""
    corresponding_quadratic_form(f::SesquilinearForm)

Given a symmetric form `f`, return the quadratic form `Q`
defined by `Q(v) = f(v,v)/2`.
It is defined only in odd characteristic.

# Examples
```jldoctest; filter = r"[0-9]"
julia> f = symmetric_form(invariant_bilinear_form(GO(3, 3)))
symmetric form with Gram matrix
[0   2   0]
[2   0   0]
[0   0   1]

julia> corresponding_quadratic_form(f)
quadratic form with Gram matrix
[0   2   0]
[0   0   0]
[0   0   2]
```
"""
function corresponding_quadratic_form(B::SesquilinearForm)
   @req B.descr==:symmetric "The form must be a symmetric form"
   @req characteristic(base_ring(B))!=2 "Corresponding quadratic form not uniquely determined"
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
    gram_matrix(f::SesquilinearForm)
    gram_matrix(f::QuadraticForm)

Return the Gram matrix that defines `f`.

# Examples
```jldoctest
julia> M = invariant_bilinear_forms(GO(3, 5))[1];

julia> M == gram_matrix(symmetric_form(M))
true

julia> M = invariant_quadratic_forms(GO(3, 5))[1];

julia> M == gram_matrix(quadratic_form(M))
true
```
"""
gram_matrix(f::SesquilinearForm) = f.matrix

function gram_matrix(f::QuadraticForm)
   isdefined(f,:matrix) && return f.matrix
   @req isdefined(f,:pol) "Cannot determine Gram matrix"
   d = nvars(parent(f.pol))
   B = zero_matrix( base_ring(f.pol), d, d )
   V = collect(AbstractAlgebra.exponent_vectors(f.pol))
   C = collect(AbstractAlgebra.coefficients(f.pol))
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
    defining_polynomial(f::QuadraticForm)

Return the polynomial that defines `f`.

# Examples
```jldoctest
julia> g = GO(5, 2);

julia> f = quadratic_form(invariant_quadratic_forms(g)[1]);

julia> defining_polynomial(f)
x1^2 + x2*x4 + x3*x5
```
"""
function defining_polynomial(f::QuadraticForm)
   isdefined(f,:pol) && return f.pol
   R = polynomial_ring(base_ring(f.matrix), nrows(f.matrix); cached=false)[1]
   p = zero(R)
   for i in 1:nrows(f.matrix), j in i:nrows(f.matrix)
      p += f.matrix[i,j] * R[i]*R[j]
   end
   f.pol = p
   return p
end


function assign_from_description(f::SesquilinearForm)
  if f.descr == :symmetric || f.descr == :alternating
    f.X = GAP.Globals.BilinearFormByMatrix(map_entries(_ring_iso(f), gram_matrix(f)), codomain(_ring_iso(f)))
  elseif f.descr == :hermitian
    f.X = GAP.Globals.HermitianFormByMatrix(map_entries(_ring_iso(f), gram_matrix(f)), codomain(_ring_iso(f)))
  else
    error("unsupported description")
  end
end

function assign_from_description(f::QuadraticForm)
  f.X = GAP.Globals.QuadraticFormByMatrix(map_entries(_ring_iso(f), gram_matrix(f)), codomain(_ring_iso(f)))
end

function _ring_iso(f::Union{SesquilinearForm, QuadraticForm})
  if !isdefined(f, :ring_iso)
      f.ring_iso = iso_oscar_gap(base_ring(f))
  end
  return f.ring_iso
end

GAP.@install function GapObj(f::Union{SesquilinearForm, QuadraticForm})
  !isdefined(f, :X) && assign_from_description(f)
  return f.X
end


########################################################################
#
# Operations
#
########################################################################

function Base.:*(f::SesquilinearForm, l::FieldElem)
   @req l != 0 "Zero is not admitted"
   @req parent(l)==base_ring(f) "The scalar does not belong to the base ring of the form"
   return SesquilinearForm(l*gram_matrix(f), f.descr)
end

function Base.:*(f::QuadraticForm, l::FieldElem)
   @req l != 0 "Zero is not admitted"
   @req parent(l)==base_ring(f) "The scalar does not belong to the base ring of the form"
   if !isdefined(f,:matrix)
      return SesquilinearForm(l*f.pol, f.descr)
   else
      g = QuadraticForm(l*gram_matrix(f), f.descr)
      if isdefined(f,:pol) g.pol=f.pol end
#T really???
      return g
   end
end

Base.:*(l::FieldElem, f::Union{SesquilinearForm, QuadraticForm}) = f*l

function Base.:^(f::SesquilinearForm{T}, x::MatElem{T}; check=false) where T <: RingElem
   @req base_ring(f)==base_ring(x) "Incompatible base rings"
   @req nrows(gram_matrix(f))==nrows(x) "Incompatible dimensions"

   if check @assert rank(x)==nrows(x) "Matrix not invertible" end
   if f.descr==:hermitian
      m = x^-1*gram_matrix(f)*conjugate_transpose(x^-1)
   else
      m = x^-1*gram_matrix(f)*transpose(x^-1)
   end
   return SesquilinearForm(m, f.descr)
end

function Base.:^(f::QuadraticForm{T}, x::MatElem{T}; check=false) where T <: RingElem
   @req base_ring(f)==base_ring(x) "Incompatible base rings"
   @req nrows(gram_matrix(f))==nrows(x) "Incompatible dimensions"

   if check @assert rank(x)==nrows(x) "Matrix not invertible" end
    m = x^-1*gram_matrix(f)*transpose(x^-1)
   return QuadraticForm(m)
end

Base.:^(f::Union{SesquilinearForm{T}, QuadraticForm{T}}, x::MatrixGroupElem{T}; check=false) where T <: RingElem = f^matrix(x)

function (f::SesquilinearForm{T})(v::AbstractAlgebra.Generic.FreeModuleElem{T},w::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: RingElem
   if f.descr==:hermitian
      w = conjugate_transpose(w.v)
   else
      w = transpose(w.v)
   end
   return v.v*gram_matrix(f)*w
end

function (f::QuadraticForm{T})(v::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: RingElem
   return v.v*gram_matrix(f)*transpose(v.v)
end



########################################################################
#
# Functionalities
#
########################################################################

"""
    radical(f::SesquilinearForm{T})
    radical(f::QuadraticForm{T})

Return `(U, emb)` where `U` is the radical of `f` and `emb` is the
embedding of `U` into the full vector space on which `f` is defined.

For a quadratic form `f`, the radical is defined as the set of vectors `v`
such that `f(v) = 0` and `v` lies in the radical of the corresponding
bilinear form.

Otherwise the radical of `f` is defined as the subspace of all `v`
such that `f(u, v) = 0` for all `u`.

# Examples
```jldoctest
julia> g = GO(7, 2);

julia> f = symmetric_form(invariant_symmetric_forms(g)[1]);

julia> U, emb = radical(f);

julia> vector_space_dim(U)
1
```
"""
function radical(f::Union{SesquilinearForm{T}, QuadraticForm{T}}) where T
   V = vector_space(base_ring(f), nrows(gram_matrix(f)) )
   R = GAP.Globals.RadicalOfForm(GapObj(f))
   GAPWrap.Dimension(R) == 0 && return sub(V, [])
   L = AbstractAlgebra.Generic.FreeModuleElem{T}[]
   for l in GAP.Globals.GeneratorsOfVectorSpace(R)
      v = V([preimage(_ring_iso(f), t) for t in l])
      push!(L,v)
   end
   return sub(V,L)
end

"""
    witt_index(f::SesquilinearForm)
    witt_index(f::QuadraticForm)

Return the Witt index of the form induced by `f` on `V/Rad(f)`.
The Witt index is the dimension of a maximal totally isotropic
(singular for quadratic forms) subspace.

# Examples
```jldoctest
julia> g = GO(1, 6, 2);

julia> witt_index(quadratic_form(invariant_quadratic_forms(g)[1]))
3

julia> g = GO(-1, 6, 2);

julia> witt_index(quadratic_form(invariant_quadratic_forms(g)[1]))
2
```
"""
witt_index(f::Union{SesquilinearForm, QuadraticForm}) = GAP.Globals.WittIndex(GapObj(f))

"""
    is_degenerate(f::SesquilinearForm)
    is_degenerate(f::QuadraticForm)

Return whether `f` is degenerate, i.e. `f` has nonzero radical.
A quadratic form is degenerate if the corresponding bilinear form is.

# Examples
```jldoctest
julia> g = GO(5, 2);

julia> f = symmetric_form(invariant_symmetric_forms(g)[1]);

julia> is_degenerate(f)
true
```
"""
is_degenerate(f::SesquilinearForm) = det(gram_matrix(f)) == 0

function is_degenerate(f::QuadraticForm)
   gram = gram_matrix(f)
   return det(gram + transpose(gram)) == 0
end

"""
    is_singular(Q::QuadraticForm)

Return whether `Q` is singular, i.e. `Q` has nonzero radical.

# Examples
```jldoctest
julia> G = GO(5, 2);

julia> f = quadratic_form(invariant_quadratic_forms(G)[1]);

julia> is_singular(f)
false
```
"""
is_singular(f::QuadraticForm) = GAPWrap.IsSingularForm(GapObj(f))
