################################################################################
## Some more getter functions
################################################################################

@doc raw"""
    defining_matrix(X::DeterminantalGerm) -> MatElem

Return the defining matrix `A` of the derterminantal germ `X` over the ring of the ambient germ of `X`. 
!!! note
    The returned matrix `A` is not a matrix over a polynomial ring, but a matrix over a localization of a polynomial ring at the complement of a maximal ideal. (Hence each entry of `A` has a numerator and a denominator.)

# Examples:
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> A = R[x 0 z;  0 y  z]
[x   0   z]
[0   y   z]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0])
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x, y, z
        over rational field
      by ideal (x*y, x*z, -y*z)
    at complement of maximal ideal of point (0, 0, 0)

julia> defining_matrix(X_A)
[x   0   z]
[0   y   z]

julia> base_ring(A) == base_ring(defining_matrix(X_A))
false

julia> base_ring(A)
Multivariate polynomial ring in 3 variables x, y, z
  over rational field

julia> base_ring(defining_matrix(X_A))
Localization
  of multivariate polynomial ring in 3 variables x, y, z
    over rational field
  at complement of maximal ideal of point (0, 0, 0)
```
"""
defining_matrix(X::DeterminantalGerm) = X.A

@doc raw"""
    determinantal_type(X::DeterminantalGerm) -> Int, Int, Int

Return the determinantal type `(n,m,t)` of the derterminantal germ `X`, where `n x m` is the size of the defining matrix and `t` the size of the minors. 


# Examples:
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> A = R[x 0 z;  0 y  z]
[x   0   z]
[0   y   z]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0])
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 3 variables x, y, z
        over rational field
      by ideal (x*y, x*z, -y*z)
    at complement of maximal ideal of point (0, 0, 0)

julia> determinantal_type(X_A)
(2, 3, 2)
```
"""
function determinantal_type(X::DeterminantalGerm)
  n, m = size(X.A)
  return n, m, X.t
end

_matrix_type(X::DeterminantalGerm{<:Ring, <:Ring, <:AffineScheme, T}) where {T} = T

################################################################################
## More constructors
################################################################################

#TODO: To add more
@doc raw"""
    DeterminantalGerm(A::MatElem{<:MPolyRingElem}, t::Int, p::Vector{T}; mat_type::Symbol = :generic, check::Bool=true)

Return the DeterminantalGerm $(X_A^t, p)$. 

# Examples:
```jldoctest
julia> R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z];

julia> A = R[v w x y;  w x y z]
[v   w   x   y]
[w   x   y   z]

julia> B = R[v w x;  w x y;  x y z]
[v   w   x]
[w   x   y]
[x   y   z]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0,0,0])
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 5 variables v, w, x, y, z
        over rational field
      by ideal (v*x - w^2, v*y - w*x, w*y - x^2, v*z - w*y, w*z - x*y, x*z - y^2)
    at complement of maximal ideal of point (0, 0, 0, 0, 0)

julia> X_B = DeterminantalGerm(B, 2, [0,0,0,0,0], mat_type=:symmetric)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 5 variables v, w, x, y, z
        over rational field
      by ideal (v*x - w^2, v*y - w*x, w*y - x^2, v*y - w*x, v*z - x^2, w*z - x*y, w*y - x^2, w*z - x*y, x*z - y^2)
    at complement of maximal ideal of point (0, 0, 0, 0, 0)

julia> X_A == X_B
false

julia> SpaceGerm(X_A) == SpaceGerm(X_B)
true
```
"""
function DeterminantalGerm(A::MatElem{<:MPolyRingElem}, t::Int, p::Vector{T};
                           mat_type::Symbol = :generic, check::Bool=true
                          ) where T<:Union{Integer, FieldElem}
  R = base_ring(A)
  kk = coefficient_ring(R)
  point = [kk.(v) for v in p]  ## throws an error, if vector entries are not compatible
  L, _ = localization(R, complement_of_point_ideal(R, point))
  return DeterminantalGerm(L.(A), t, mat_type=mat_type, check=check)
end

################################################################################
## basic functionality for determinantal germs, which needs to be overwriten
################################################################################

function ==(X::DeterminantalGerm, Y::DeterminantalGerm)
  X === Y && return true
  _matrix_type(X) == _matrix_type(Y) || return false
  determinantal_type(X) == determinantal_type(Y) || return false
  ambient_coordinate_ring(X) === ambient_coordinate_ring(Y) || return false
  defining_matrix(X) == defining_matrix(Y) && return true
  return underlying_scheme(X) == underlying_scheme(Y)
end

################################################################################
## T1_GL module
################################################################################

function _R_ij(A::MatElem, i::Integer, j::Integer)
  R_ij = zero(A)
  R_ij[i, :] = A[j, :]
  return R_ij
end

function _C_ij(A::MatElem, i::Integer, j::Integer)
  C_ij = zero(A)
  C_ij[:, i] = A[:, j]
  return C_ij
end

# TODO: move to a more fitting place (maybe AbstractAlgebra)
function _sym_mat_gens(A::MatElem)
  n, m = size(A)
  @req n == m "matrix 'A' must be a quadratic."
  R = base_ring(A)
  gens = typeof(A)[]
  sizehint!(gens, div(n*(n+1), 2))
  for i in 1:n
    for j in i:n
      tmp = zero(A)
      tmp[i,j] = one(R)
      tmp[j,i] = one(R)
      push!(gens, tmp)
    end
  end
  return gens
end

# TODO: move to a more fitting place (maybe AbstractAlgebra)
function _skew_sym_mat_gens(A::MatElem)
  n, m = size(A)
  @req n == m "matrix 'A' must be a quadratic."
  R = base_ring(A)
  characteristic(R) == 2 && return _sym_mat_gens(A)
  gens = typeof(A)[]
  sizehint!(gens, div(n*(n-1), 2))
  for i in 1:n
    for j in i+1:n
      tmp = zero(A)
      tmp[i,j] = -one(R)
      tmp[j,i] = one(R)
      push!(gens, tmp)
    end
  end
  return gens
end

@doc raw"""
    T1_GL_module(X::DeterminantalGerm) -> SubquoModule

Return the $T^1_{GL}$-module of the defining matrix `A` of the determinantal germ of `X`. 
!!! note
    Different determinantal structures for the same underlying space germ may yield different $T^1_{GL}$-modules.

# Examples:
```jldoctest
julia> R, (x,y) = QQ[:x, :y];

julia> A = R[x 0; 0 y^2+x^2]
[x           0]
[0   x^2 + y^2]

julia> X_A = DeterminantalGerm(A, 2, [0,0]);

julia> X_A_sym = DeterminantalGerm(A, 2, [0,0], mat_type=:symmetric);

julia> X_A == X_A_sym
false

julia> underlying_space_germ(X_A) == underlying_space_germ(X_A_sym)
true

julia> T1_A = T1_GL_module(X_A)
Subquotient of submodule with 4 generators
  1: E[1,1]
  2: E[1,2]
  3: E[2,1]
  4: E[2,2]
by submodule with 10 generators
  1: E[1,1] + 2*x*E[2,2]
  2: 2*y*E[2,2]
  3: x*E[1,1]
  4: (x^2 + y^2)*E[1,2]
  5: x*E[2,1]
  6: (x^2 + y^2)*E[2,2]
  7: x*E[1,1]
  8: (x^2 + y^2)*E[2,1]
  9: x*E[1,2]
  10: (x^2 + y^2)*E[2,2]

julia> vector_space_dim(T1_A)
6

julia> T1_A_sym = T1_GL_module(X_A_sym)
Subquotient of submodule with 3 generators
  1: E[1,1]
  2: E[1,2] + E[2,1]
  3: E[2,2]
by submodule with 6 generators
  1: E[1,1] + 2*x*E[2,2]
  2: 2*y*E[2,2]
  3: 2*x*E[1,1]
  4: (x^2 + y^2)*E[1,2] + (x^2 + y^2)*E[2,1]
  5: x*E[1,2] + x*E[2,1]
  6: (2*x^2 + 2*y^2)*E[2,2]

julia> vector_space_dim(T1_A_sym)
4
```
"""
@attr SubquoModule function T1_GL_module(X::DeterminantalGerm)
  # transposing, since '_vec' vcats the columms of A and we would rather read rowwise
  A = transpose(defining_matrix(X))
  m, n = size(A)
  L = base_ring(parent(A))
  N = ngens(L)
  # defining_matrix has size n x m
  F = FreeMod(L, [Symbol("E[$i,$j]") for i in 1:n for j in 1:m])

  if _matrix_type(X) === Val{:generic}
    rels = vcat([derivative.(A, i) for i in 1:N],
                [_C_ij(A, i, j) for i in 1:n for j in 1:n],
                [_R_ij(A, i, j) for i in 1:m for j in 1:m]
               )
    return SubquoModule(F, gens(F), F.(_vec.(rels)))
  end

  # Case: _matrix_type(X) === Val{:symmetric} or # _matrix_type(X) === Val{:skew_symmetric}
  erz = _matrix_type(X) === Val{:symmetric} ? _sym_mat_gens(A) : _skew_sym_mat_gens(A)
  rels = vcat([derivative.(A, i) for i in 1:N],
              [_R_ij(A, i, j) + _C_ij(A, i, j) for i in 1:n for j in 1:n],
             )
  return SubquoModule(F, F.(_vec.(erz)), F.(_vec.(rels)))
end

@doc raw"""
    tjurina_GL_number(X::DeterminantalGerm)

Return the tjurina_GL_number of the determinantal germ $(X_A^t, p)$. 

!!! note
    Different determinantal structures for the same underlying space germ may yield different `tjurina_GL_number`s.

# Examples:
```jldoctest
julia> R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z];

julia> A = R[v w x y;  w x y z]
[v   w   x   y]
[w   x   y   z]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0,0,0]);

julia> B = R[v w x;  w x y;  x y z]
[v   w   x]
[w   x   y]
[x   y   z]

julia> X_B = DeterminantalGerm(B, 2, [0,0,0,0,0], mat_type=:symmetric);

julia> X_A == X_B
false

julia> underlying_space_germ(X_A) == underlying_space_germ(X_B)
true

julia> tjurina_GL_number(X_A)
3

julia> tjurina_GL_number(X_B)
1
```
"""
tjurina_GL_number(X::DeterminantalGerm) = vector_space_dim(T1_GL_module(X))

is_determinantally_rigid(X::DeterminantalGerm) = is_zero(T1_GL_module(X))

@attr Bool is_EIDS(X::DeterminantalGerm{<:Any, <:Any, <:Any, Val{:generic}}) = krull_dim(T1_GL_module(X)) <= 0

@doc raw"""
    basis_versal_det_unfolding(X::DeterminantalGerm) -> Vector{SubquoModuleElem}

Return a basis of a versal determinantal unfolding of the determinantal germ `X`. 

!!! note
    Different determinantal structures for the same underlying space germ may yield different versal determinantal unfoldings.

# Examples:
```jldoctest
julia> R, (v,w,x,y,z) = QQ[:v,:w,:x,:y,:z];

julia> A = R[v w x y;  w x y z]
[v   w   x   y]
[w   x   y   z]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0,0,0]);

julia> B = R[v w x;  w x y;  x y z]
[v   w   x]
[w   x   y]
[x   y   z]

julia> X_B = DeterminantalGerm(B, 2, [0,0,0,0,0], mat_type=:symmetric);

julia> X_A == X_B
false

julia> underlying_space_germ(X_A) == underlying_space_germ(X_B)
true

julia> basis_versal_det_unfolding(X_A)
3-element Vector{SubquoModuleElem{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}}:
 E[1,2]
 E[1,3]
 E[1,4]

julia> basis_versal_det_unfolding(X_B)
1-element Vector{SubquoModuleElem{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}}:
 E[1,3] + E[3,1]
```
"""
basis_versal_det_unfolding(X::DeterminantalGerm) = vector_space_basis(T1_GL_module(X))

