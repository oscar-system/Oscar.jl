################################################################################
## generators of ambient_module of T1_GL_module and T1_SL_module
################################################################################

# TODO: move some of them to a more fitting place (maybe AbstractAlgebra)

function _mat_space_gens(M::MatSpace)
  R = base_ring(M)
  gens = sizehint!(elem_type(M)[], nrows(M)*ncols(M))
  for j in 1:ncols(M) 
    for i in 1:nrows(M)
      tmp = zero(M)
      tmp[i,j] = one(R)
      push!(gens, tmp)
    end
  end
  return gens
end

function _sym_mat_gens(M::MatSpace)
  AbstractAlgebra.check_square(M)
  n = nrows(M)
  R = base_ring(M)
  gens = sizehint!(elem_type(M)[], div(n*(n+1), 2))
  for i in 1:n
    for j in i:n
      tmp = zero(M)
      tmp[i,j] = one(R)
      tmp[j,i] = one(R)
      push!(gens, tmp)
    end
  end
  return gens
end

function _skew_sym_mat_gens(M::MatSpace)
  AbstractAlgebra.check_square(M)
  R = base_ring(M)
  characteristic(R) == 2 && return _sym_mat_gens(M)
  n = nrows(M)
  gens = sizehint!(elem_type(M)[], div(n*(n-1), 2))
  for i in 1:n
    for j in i+1:n
      tmp = zero(M)
      tmp[i,j] = -one(R)
      tmp[j,i] = one(R)
      push!(gens, tmp)
    end
  end
  return gens
end


## dispatch generators of 'ambient_module' based on the Val type

_T1_mat_gens(A::MatElem, ::Type{Val{:generic}}) = _mat_space_gens(parent(A))
_T1_mat_gens(A::MatElem, ::Type{Val{:symmetric}}) = _sym_mat_gens(parent(A))
_T1_mat_gens(A::MatElem, ::Type{Val{:skew_symmetric}}) = _skew_sym_mat_gens(parent(A))


################################################################################
## relations of T1_GL_module and T1_SL_module
################################################################################

## internal helper functions for relations of T1_GL_module and T1_SL_module

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

_J(A::MatElem) = [derivative.(A, i) for i in 1:ngens(base_ring(A))]



## dispatch relations based on the val type

function _T1_GL_rels(A::MatElem, ::Type{Val{:generic}}) 
  return vcat(_J(A), [_C_ij(A, i, j) for i in 1:ncols(A) for j in 1:ncols(A)], 
                     [_R_ij(A, i, j) for i in 1:nrows(A) for j in 1:nrows(A)])
end

function _T1_GL_rels(A::MatElem, ::Union{Type{Val{:symmetric}}, Type{Val{:skew_symmetric}}})
  AbstractAlgebra.check_square(A)
  return vcat(_J(A), [_R_ij(A, i, j) + _C_ij(A, i, j) for i in 1:nrows(A) for j in 1:ncols(A)])
end



function _T1_SL_rels(A::MatElem, ::Type{Val{:generic}})
  rels = _J(A)    # partial derivatives of matrix
  n, m = size(A)
  sizehint!(rels, length(rels) + n^2-1 + m^2-1)
  for i in 1:m
    for j in 1:m
      if i == j   # entry on diagonal - use bottom right entry to make trace zero
        i == m && continue    # skip bottom right entry
        push!(rels, _C_ij(A, i, i) - _C_ij(A, m, m))
      else
        push!(rels, _C_ij(A, i, j))
      end
    end
  end
  for j in 1:n
    for i in 1:n
      if i == j   # entry on diagonal - use bottom right entry to make trace zero
        i == n && continue    # skip bottom right entry
        push!(rels, _R_ij(A, i, i) - _R_ij(A, n, n))
      else
        push!(rels, _R_ij(A, i, j))
      end
    end
  end
  return rels
end

function _T1_SL_rels(A::MatElem, ::Union{Type{Val{:symmetric}}, Type{Val{:skew_symmetric}}})
  AbstractAlgebra.check_square(A)
  rels = _J(A)    # partial derivatives of matrix
  n = nrows(A)
  sizehint!(rels, length(rels) + n^2-1)
  for i in 1:n
    for j in 1:n
      if i == j   # entry on diagonal - use bottom right entry to make trace zero
        i == n && continue    # skip bottom right entry
        push!(rels, _R_ij(A, i, i) + _C_ij(A, i, i) - _R_ij(A, n, n) - _C_ij(A, n, n))
      else
        push!(rels, _R_ij(A, i, j) + _C_ij(A, i, j))
      end
    end
  end
  return rels
end



################################################################################
## internal workhorse for 'T1_GL module' and 'T1_SL_module'
################################################################################

function _T1_GL_module(A::MatElem, val_type::Type{<:Val})
  @req val_type <: Union{Val{:generic}, Val{:symmetric}, Val{:skew_symmetric}} "$val_type is not supported"
  M = parent(A)
  n, m = size(A)
  # transposing, since '_vec' vcats the columms of A and we would rather read rowwise
  A = transpose(A)
  # transitive function calls by '_T1_gens' throws an errow for a non-square matrix in the (skew-)symmetric case, therefore call it already here
  erz = _T1_mat_gens(A, val_type)
  L = base_ring(A)
  function symb_fun()
    return [Symbol("E[$i,$j]") for i in 1:n for j in 1:m]
  end
  F = FreeMod{elem_type(L)}(n*m, L, symb_fun) 
  T1_GL = SubquoModule(F, F.(_vec.(erz)), F.(_vec.(_T1_GL_rels(A, val_type))))
  # interpretation map
  im(v::FreeModElem) = M(transpose(reshape(Vector(coordinates(v), n*m), (m,n))))
  pre(B::MatElem) = F(_vec(transpose(B)))
  interp = MapFromFunc(ambient_free_module(T1_GL), M, im, pre)
  return T1_GL, interp
end



function _T1_SL_module(A::MatElem, val_type::Type{<:Val})
  @req val_type <: Union{Val{:generic}, Val{:symmetric}, Val{:skew_symmetric}} "$val_type is not supported"
  M = parent(A)
  n, m = size(A)
  # transposing, since '_vec' vcats the columms of A and we would rather read rowwise
  A = transpose(A)
  # transitive function calls by '_T1_gens' throws an errow for a non-square matrix in the (skew-)symmetric case, therefore call it already here
  erz = _T1_mat_gens(A, val_type)
  L = base_ring(A)
  function symb_fun()
    return [Symbol("E[$i,$j]") for i in 1:n for j in 1:m]
  end
  F = FreeMod{elem_type(L)}(n*m, L, symb_fun) 
  T1_SL = SubquoModule(F, F.(_vec.(erz)), F.(_vec.(_T1_SL_rels(A, val_type))))
  # interpretation map
  im(v::FreeModElem) = M(transpose(reshape(Vector(coordinates(v), n*m), (m,n))))
  pre(B::MatElem) = F(_vec(transpose(B)))
  interp = MapFromFunc(ambient_free_module(T1_SL), M, im, pre)
  return T1_SL, interp
end



################################################################################
## T1_GL module
################################################################################

@doc raw"""
    T1_GL_module(A::MatElem; mat_type::Symbol = :generic) -> SubquoModule, MapFromFunc

Return for a matrix `A` with entries in either a `MPolyRing` or a `MPolyLocRing` the $T^1_{GL}(A)$-module of `A` and an interpretation map, which is an isomorphism from the `ambient_free_module` of the $T1_GL(A)$-module to the matrix space of `A`.

# Examples:
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> A = R[x 0 z; 0 y z]
[x   0   z]
[0   y   z]

julia> T1_GL, intr = T1_GL_module(A, mat_type=:generic);

julia> B = vector_space_basis(T1_GL)
3-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 E[1,2]
 E[1,3]
 E[2,1]

julia> intr.(repres.(B))
3-element Vector{AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}}:
 [0 1 0; 0 0 0]
 [0 0 1; 0 0 0]
 [0 0 0; 1 0 0]
```
"""
function T1_GL_module(A::MatElem; mat_type::Symbol = :generic)
  @req mat_type in (:generic, :symmetric, :skew_symmetric) "'mat_type' must be either ':generic', ':symmetric' or 'skew_symmetric'"
  return _T1_GL_module(A, Val{mat_type})
end



@doc raw"""
    T1_GL_module(X::DeterminantalGerm) -> SubquoModule, MapFromFunc

Return for a determinantal germ of `X` the $T^1_{GL}(A)$-module of the defining matrix `A` and an interpretation map, which is an isomorphism from the `ambient_free_module` of the $T1_GL(A)$-module to the matrix space of `A`. 
!!! note
    Different determinantal structures for the same underlying space germ may yield different $T^1_{GL}$-modules.

# Examples:
```jldoctest
julia> R, (x,y) = QQ[:x, :y];

julia> A = R[x 0;  0 y^2+x^2]
[x           0]
[0   x^2 + y^2]

julia> X_A = DeterminantalGerm(A, 2, [0,0]);

julia> X_A_sym = DeterminantalGerm(A, 2, [0,0], mat_type = :symmetric);

julia> X_A == X_A_sym
false

julia> SpaceGerm(representative(X_A), point(X_A)) == SpaceGerm(representative(X_A_sym), point(X_A_sym))
true

julia> T1_A, _ = T1_GL_module(X_A);

julia> T1_A
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

julia> T1_A_sym, _ = T1_GL_module(X_A_sym);

julia> T1_A_sym
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



julia> P, (t,) = polynomial_ring(QQ, [:t])
(Multivariate polynomial ring in 1 variable over QQ, QQMPolyRingElem[t])

julia> B = P[0 t^3;  -t^3 0]
[   0   t^3]
[-t^3     0]

julia> X_B_skew = DeterminantalGerm(B, 1, [0], mat_type = :skew_symmetric)
Spectrum
  of localization
    of quotient
      of multivariate polynomial ring in 1 variable t
        over rational field
      by ideal (t^3)
    at complement of maximal ideal of point (0)

julia> T1_B_skew, _  = T1_GL_module(X_B_skew);

julia> T1_B_skew
Subquotient of submodule with 1 generator
  1: E[1,2] - E[2,1]
by submodule with 5 generators
  1: 3*t^2*E[1,2] - 3*t^2*E[2,1]
  2: t^3*E[1,2] - t^3*E[2,1]
  3: 0
  4: 0
  5: t^3*E[1,2] - t^3*E[2,1]
```
"""
@attr Tuple{SubquoModule, MapFromFunc} T1_GL_module(X::DeterminantalGerm) = _T1_GL_module(defining_matrix(X), _mat_type(X))



@doc raw"""
    tjurina_GL_number(X::DeterminantalGerm) -> Union{Integer, PosInf}

Return the `tjurina_GL_number` of the determinantal germ `X`. 

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

julia> X_B = DeterminantalGerm(B, 2, [0,0,0,0,0], mat_type = :symmetric);

julia> X_A == X_B
false

julia> SpaceGerm(representative(X_A), point(X_A)) == SpaceGerm(representative(X_B), point(X_B))
true

julia> tjurina_GL_number(X_A)
3

julia> tjurina_GL_number(X_B)
1



julia> P, (a, b) = QQ[:a,:b];

julia> C = P[0 0 a b^2;  0 0 b^2 a^2;  -a -b^2 0 0;  -b^2 -a^2 0 0]
[   0      0     a   b^2]
[   0      0   b^2   a^2]
[  -a   -b^2     0     0]
[-b^2   -a^2     0     0]

julia> X_C_skew = DeterminantalGerm(C, 2, [0,0], mat_type = :skew_symmetric);

julia> tjurina_GL_number(X_C_skew)
12
```
"""
tjurina_GL_number(X::DeterminantalGerm{<:Field, <:Ring, <:AffineScheme, <:Val}) = vector_space_dim(T1_GL_module(X)[1])



@doc raw"""
    is_determinantally_rigid(X::DeterminantalGerm)

Return whether the determinantal germ `X` is determinantally rigid. 

!!! note
    Different determinantal structures for the same underlying space germ may yield different results.

# Examples:
```jldoctest
julia> R, (w,x,y,z) = QQ[:w,:x,:y,:z];

julia> A = R[x z;  w y]
[x   z]
[w   y]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0,0]);

julia> B = matrix([x*y-w*z])
[-w*z + x*y]

julia> X_B = DeterminantalGerm(B, 1, [0,0,0,0]);

julia> X_A == X_B
false

julia> SpaceGerm(representative(X_A), point(X_A)) == SpaceGerm(representative(X_B), point(X_B))
true

julia> is_determinantally_rigid(X_A)
true

julia> is_determinantally_rigid(X_B)
false
```
"""
is_determinantally_GL_rigid(X::DeterminantalGerm) = is_zero(T1_GL_module(X)[1])



@doc raw"""
    is_EIDS(X::DeterminantalGerm{<:Ring, <:Ring, <:AffineScheme, Val{:generic}}) -> Bool

Return whether the determinantal germ `X` is an essentially isolated determinantal singularity (EIDS). 

# Examples:
```jldoctest
julia> R, (w,x,y,z) = QQ[:w,:x,:y,:z];

julia> A = R[x y z;  w x y]
[x   y   z]
[w   x   y]

julia> X_A = DeterminantalGerm(A, 2, [0,0,0,0]);

julia> is_EIDS(X_A)
true

julia> B = R[x y z;  w x 0]
[x   y   z]
[w   x   0]

julia> X_B = DeterminantalGerm(B, 2, [0,0,0,0]);

julia> is_EIDS(X_B)
false
```
"""
@attr Bool is_EIDS(X::DeterminantalGerm{<:Ring, <:Ring, <:AffineScheme, Val{:generic}}) = krull_dim(T1_GL_module(X)[1]) <= 0



@doc raw"""
    basis_GL_versal_determinantal_unfolding(X::DeterminantalGerm) -> Vector{SubquoModuleElem}

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

julia> X_B_sym = DeterminantalGerm(B, 2, [0,0,0,0,0], mat_type = :symmetric);

julia> X_A == X_B_sym
false

julia> SpaceGerm(representative(X_A), point(X_A)) == SpaceGerm(representative(X_B_sym), point(X_B_sym))
true

julia> basis_GL_versal_determinantal_unfolding(X_A)
3-element Vector{SubquoModuleElem{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}}:
 E[1,2]
 E[1,3]
 E[1,4]

julia> basis_GL_versal_determinantal_unfolding(X_B_sym)
1-element Vector{SubquoModuleElem{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}}:
 E[1,3] + E[3,1]



julia> P, (a, b) = QQ[:a,:b];

julia> C = P[0 0 a b;  0 0 b a;  -a -b 0 0; -b -a 0 0]
[ 0    0   a   b]
[ 0    0   b   a]
[-a   -b   0   0]
[-b   -a   0   0]

julia> X_C_skew = DeterminantalGerm(C, 2, [0,0], mat_type = :skew_symmetric);

julia> basis_GL_versal_determinantal_unfolding(X_C_skew)
4-element Vector{SubquoModuleElem{MPolyLocRingElem{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem, MPolyComplementOfKPointIdeal{QQField, QQFieldElem, QQMPolyRing, QQMPolyRingElem}}}}:
 E[1,2] - E[2,1]
 E[1,3] - E[3,1]
 E[1,4] - E[4,1]
 E[3,4] - E[4,3]
```
"""
basis_GL_versal_determinantal_unfolding(X::DeterminantalGerm{<:Field, <:Ring, <:AffineScheme, <:Val}) = vector_space_basis(T1_GL_module(X)[1])



################################################################################
## T1_SL module
################################################################################

@doc raw"""
    T1_SL_module(A::MatElem; mat_type::Symbol = :generic) -> SubquoModule, MapFromFunc

Return for a matrix `A` with entries in either a `MPolyRing` or a `MPolyLocRing` the $T^1_{GL}(A)$-module of `A` and an interpretation map, which is an isomorphism from the `ambient_free_module` of the $T1_GL(A)$-module to the matrix space of `A`.

# Examples:
```jldoctest
julia> R, (x,y,z) = QQ[:x, :y, :z];

julia> A = R[z y x^2; 0 x y]
[z   y   x^2]
[0   x     y]

julia> T1_SL, intr = T1_SL_module(A, mat_type=:generic);

julia> B = vector_space_basis(T1_SL)
4-element Vector{SubquoModuleElem{QQMPolyRingElem}}:
 E[1,2]
 E[1,3]
 E[2,1]
 E[2,2]

julia> intr.(repres.(B))
4-element Vector{AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}}:
 [0 1 0; 0 0 0]
 [0 0 1; 0 0 0]
 [0 0 0; 1 0 0]
 [0 0 0; 0 1 0]
```
"""
function T1_SL_module(A::MatElem; mat_type::Symbol = :generic)
  @req mat_type in (:generic, :symmetric, :skew_symmetric) "'mat_type' must be either ':generic', ':symmetric' or 'skew_symmetric'"
  return _T1_SL_module(A, Val{mat_type})
end



@doc raw"""
    T1_SL_module(X::DeterminantalGerm) -> SubquoModule, MapFromFunc

Return the $T^1_{SL}(A)$-module of the defining matrix `A` of the determinantal germ of `X` and an interpretation map, which maps the representative of elements of the $T1_{SL}(A)$-module to their corresponding matrices in the matrix space of `A`. 
!!! note
    Different determinantal structures for the same underlying space germ may yield different $T^1_{SL}$-modules.

# Examples:
```jldoctest
julia> R, (x,y) = QQ[:x, :y];

julia> A = R[x 0;  0 x*y+y^3]
[x           0]
[0   x*y + y^3]

julia> X_A = DeterminantalGerm(A, 2, [0,0]);

julia> X_A_sym = DeterminantalGerm(A, 2, [0,0], mat_type = :symmetric);

julia> X_A == X_A_sym
false

julia> SpaceGerm(representative(X_A), point(X_A)) == SpaceGerm(representative(X_A_sym), point(X_A_sym))
true

julia> T1_A, _ = T1_GL_module(X_A);

julia> T1_A
Subquotient of submodule with 4 generators
  1: E[1,1]
  2: E[1,2]
  3: E[2,1]
  4: E[2,2]
by submodule with 10 generators
  1: E[1,1] + y*E[2,2]
  2: (x + 3*y^2)*E[2,2]
  3: x*E[1,1]
  4: (x*y + y^3)*E[1,2]
  5: x*E[2,1]
  6: (x*y + y^3)*E[2,2]
  7: x*E[1,1]
  8: (x*y + y^3)*E[2,1]
  9: x*E[1,2]
  10: (x*y + y^3)*E[2,2]

julia> T1_A_sym, _ = T1_GL_module(X_A_sym);

julia> T1_A_sym
Subquotient of submodule with 3 generators
  1: E[1,1]
  2: E[1,2] + E[2,1]
  3: E[2,2]
by submodule with 6 generators
  1: E[1,1] + y*E[2,2]
  2: (x + 3*y^2)*E[2,2]
  3: 2*x*E[1,1]
  4: (x*y + y^3)*E[1,2] + (x*y + y^3)*E[2,1]
  5: x*E[1,2] + x*E[2,1]
  6: (2*x*y + 2*y^3)*E[2,2]

julia> vector_space_dim(T1_A)
9

julia> vector_space_dim(T1_A_sym)
6
```
"""
@attr Tuple{SubquoModule, MapFromFunc} T1_SL_module(X::DeterminantalGerm) = _T1_SL_module(defining_matrix(X), _mat_type(X))



@doc raw"""
    tjurina_SL_number(X::DeterminantalGerm) -> Union{Integer, PosInf}

Return the `tjurina_SL_number` of the determinantal germ `X`. 

!!! note
    Different determinantal structures for the same underlying space germ may yield different `tjurina_SL_number`s.

# Examples:
```jldoctest
julia> R, (x,y) = QQ[:x, :y];

julia> A = R[0 x y;  x y 0; y 0 x^2]
[0   x     y]
[x   y     0]
[y   0   x^2]

julia> X_A = DeterminantalGerm(A, 3, [0,0]);

julia> X_A_sym = DeterminantalGerm(A, 3, [0,0], mat_type = :symmetric);

julia> X_A == X_A_sym
false

julia> SpaceGerm(representative(X_A), point(X_A)) == SpaceGerm(representative(X_A_sym), point(X_A_sym))
true

julia> tjurina_SL_number(X_A)
9

julia> tjurina_SL_number(X_A_sym)
6
```
"""
tjurina_SL_number(X::DeterminantalGerm{<:Field, <:Ring, <:AffineScheme, <:Val}) = vector_space_dim(T1_SL_module(X)[1])


# TODO: doctest
@doc raw"""
    is_determinantally_SL_rigid(X::DeterminantalGerm)

Return whether the determinantal germ `X` is determinantally rigid. 

!!! note
    Different determinantal structures for the same underlying space germ may yield different results.

"""
is_determinantally_SL_rigid(X::DeterminantalGerm) = is_zero(T1_SL_module(X)[1])

# TODO: doctest
@doc raw"""
    basis_SL_versal_determinantal_unfolding(X::DeterminantalGerm) -> Vector{SubquoModuleElem}

Return a basis of a versal determinantal unfolding of the determinantal germ `X`. 

!!! note
    Different determinantal structures for the same underlying space germ may yield different versal determinantal unfoldings.

"""
basis_SL_versal_determinantal_unfolding(X::DeterminantalGerm{<:Field, <:Ring, <:AffineScheme, <:Val}) = vector_space_basis(T1_SL_module(X)[1])
