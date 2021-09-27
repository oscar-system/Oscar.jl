export invariant_ring, primary_invariants, secondary_invariants, irreducible_secondary_invariants, fundamental_invariants, affine_algebra
export coefficient_ring, polynomial_ring, action, group
export ismodular
export reynolds_operator, molien_series

################################################################################
#
#  Field access
#
################################################################################

coefficient_ring(I::InvRing) = I.field

polynomial_ring(I::InvRing) = I.poly_ring

action(I::InvRing) = I.action

_action_singular(I::InvRing) = I.action_singular

group(I::InvRing) = I.group

ismodular(I::InvRing) = I.modular

function invariant_ring(M::Vector{<: MatrixElem})
  return invariant_ring(base_ring(M[1]), M)
end

invariant_ring(matrices::MatrixElem{T}...) where {T} = invariant_ring(collect(matrices))

function invariant_ring(K::Field, M::Vector{<: MatrixElem})
  return invariant_ring(matrix_group([change_base_ring(K, g) for g in M]))
end

#######################################################

@doc Markdown.doc"""
    invariant_ring(G::MatrixGroup)

Return the invariant ring of the finite matrix group `G`.
    
CAVEAT: The creation of invariant rings is lazy in the sense that no explicit computations are done until specifically invoked (for example by the `primary_invariants` function).

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [M1, M2]);

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]
```
"""
function invariant_ring(G::MatrixGroup)
  n = degree(G)
  action = mat_elem_type(typeof(G))[g.elm for g in gens(G)]
  return InvRing(base_ring(G), G, action)
end

#######################################################

function Base.show(io::IO, IR::InvRing)
  print(io, "Invariant ring of\n")
  print(io, group(IR), "\n")
  print(io, "with generators\n")
  print(io, action(IR))
end

# Returns a map performing the right action of M on the ring R
function right_action(R::MPolyRing{T}, M::MatrixElem{T}) where T
  @assert nvars(R) == ncols(M)
  @assert nrows(M) == ncols(M)
  n = nvars(R)

  vars = zeros(R, n)
  x = gens(R)
  for i = 1:n
    for j = 1:n
      if iszero(M[i, j])
        continue
      end
      vars[i] = addeq!(vars[i], M[i, j]*x[j])
    end
  end

  right_action_by_M = (f::MPolyElem{T}) -> evaluate(f, vars)

  return MapFromFunc(right_action_by_M, R, R)
end

right_action(R::MPolyRing{T}, M::MatrixGroupElem{T}) where T = right_action(R, M.elm)
right_action(f::MPolyElem{T}, M::MatrixElem{T}) where T = right_action(parent(f), M)(f)
right_action(f::MPolyElem{T}, M::MatrixGroupElem{T}) where T = right_action(f, M.elm)

function reynolds_molien_via_singular(IR::InvRing{T}) where {T <: Union{FlintRationalField, AnticNumberField}}
  if !isdefined(IR, :reynolds_singular) || !isdefined(IR, :molien_singular)
    singular_matrices = _action_singular(IR)

    rey, mol = Singular.LibFinvar.reynolds_molien(singular_matrices...)
    IR.reynolds_singular = rey
    IR.molien_singular = mol
  end
  return IR.reynolds_singular, IR.molien_singular
end

function reynolds_molien_via_singular(IR::InvRing{T}) where {T <: Union{Nemo.GaloisField, Nemo.GaloisFmpzField}}
  @assert !ismodular(IR)
  if !isdefined(IR, :reynolds_singular) || !isdefined(IR, :molien_singular)
    singular_matrices = _action_singular(IR)

    rey = Singular.LibFinvar.reynolds_molien(singular_matrices..., "")
    mol = Singular.lookup_library_symbol("Finvar", "newring")[2][:M]
    IR.reynolds_singular = rey
    IR.molien_singular = mol
  end
  return IR.reynolds_singular, IR.molien_singular
end

function reynolds_via_singular(IR::InvRing{T}) where {T <: Union{FlintRationalField, AnticNumberField, Nemo.GaloisField, Nemo.GaloisFmpzField}}
  return reynolds_molien_via_singular(IR)[1]
end

# Singular.LibFinvar.reynolds_molien does not work for finite fields which are
# not prime fields.
function reynolds_via_singular(IR::InvRing{T}) where {T <: Union{FqNmodFiniteField, FqFiniteField}}
  @assert !ismodular(IR)
  if !isdefined(IR, :reynolds_singular)
    singular_matrices = _action_singular(IR)

    rey = Singular.LibFinvar.group_reynolds(singular_matrices...)[1]
    IR.reynolds_singular = rey
  end
  return IR.reynolds_singular
end

# Need a special function for char 0 at the moment
# TODO: Remove this function once we can iterate over and call `order` on char 0
# matrix groups
function _prepare_reynolds_operator(IR::InvRing{FldT, GrpT, T}) where {FldT <: Union{FlintRationalField, AnticNumberField}, GrpT <: MatrixGroup, T}
  if isdefined(IR, :reynolds_operator)
    return nothing
  end

  H, GtoH = Oscar.isomorphic_group_over_finite_field(group(IR))
  elements_of_G = [ preimage(GtoH, h) for h in H ]
  actions = [ right_action(polynomial_ring(IR), g) for g in elements_of_G ]
  function reynolds(f::T)
    g = parent(f)()
    for action in actions
      g = addeq!(g, action(f))
    end
    return g*base_ring(f)(1//order(H))
  end

  IR.reynolds_operator = MapFromFunc(reynolds, polynomial_ring(IR), polynomial_ring(IR))
  return nothing
end

function _prepare_reynolds_operator(IR::InvRing{FldT, GrpT, PolyElemT}) where {FldT, GrpT, PolyElemT}
  @assert !ismodular(IR)

  if isdefined(IR, :reynolds_operator)
    return nothing
  end

  actions = [ right_action(polynomial_ring(IR), g) for g in group(IR) ]
  function reynolds(f::PolyElemT)
    g = parent(f)()
    for action in actions
      g = addeq!(g, action(f))
    end
    return g*base_ring(f)(1//order(group(IR)))
  end

  IR.reynolds_operator = MapFromFunc(reynolds, polynomial_ring(IR), polynomial_ring(IR))
  return nothing
end

function reynolds_operator_via_oscar(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem}
  @assert !ismodular(IR)

  if !isdefined(IR, :reynolds_operator)
    _prepare_reynolds_operator(IR)
  end
  return IR.reynolds_operator(f)
end

function reynolds_operator_via_singular(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem}
   @assert parent(f) === polynomial_ring(IR)

   rey = reynolds_via_singular(IR)
   fSing = singular_ring(polynomial_ring(IR))(f)
   fReySing = Singular.LibFinvar.evaluate_reynolds(rey, fSing)
   # fReySing is an ideal...
   @assert length(gens(fReySing)) == 1
   return polynomial_ring(IR)(gens(fReySing)[1])
end

@doc Markdown.doc"""
     reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem}

In the non-modular case, return the image of `f` under the Reynolds operator projecting onto `IR`.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> R = polynomial_ring(IR)
Multivariate Polynomial Ring in x[1], x[2], x[3] over Cyclotomic field of order 3 graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

julia> x = gens(R)
3-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]
 x[2]
 x[3]

julia> f = x[1]^3
x[1]^3

julia> reynolds_operator(IR, f)
1//3*x[1]^3 + 1//3*x[2]^3 + 1//3*x[3]^3
```
```jldoctest
julia> M = matrix(GF(3), [0 1 0; -1 0 0; 0 0 -1])
[0   1   0]
[2   0   0]
[0   0   2]

julia> G = MatrixGroup(3, GF(3), [M])
Matrix group of degree 3 over Galois field with characteristic 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Galois field with characteristic 3
with generators
gfp_mat[[0 1 0; 2 0 0; 0 0 2]]

julia> R = polynomial_ring(IR)
Multivariate Polynomial Ring in x[1], x[2], x[3] over Galois field with characteristic 3 graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

julia> x = gens(R)
3-element Vector{MPolyElem_dec{gfp_elem, gfp_mpoly}}:
 x[1]
 x[2]
 x[3]

julia> f = x[1]^2
x[1]^2

julia> reynolds_operator(IR, f)
2*x[1]^2 + 2*x[2]^2

julia> f = x[1]^3
x[1]^3

julia> reynolds_operator(IR, f)
0
```
"""
reynolds_operator(IR::InvRing{FldT, GrpT, T}, f::T) where {FldT, GrpT, T <: MPolyElem} = reynolds_operator_via_oscar(IR, f)

function reynolds_operator(IR::InvRing, f::MPolyElem)
  @assert parent(f) === polynomial_ring(IR).R
  return reynolds_operator(IR, polynomial_ring(IR)(f))
end

function basis_via_reynolds(IR::InvRing, d::Int)
  @assert d >= 0 "Degree must be non-negative"
  @assert !ismodular(IR)
  R = polynomial_ring(IR)
  if d == 0
    return elem_type(R)[ one(R) ]
  end

  rey = reynolds_via_singular(IR)
  basisSing = Singular.LibFinvar.invariant_basis_reynolds(rey, d)
  res = Vector{elem_type(R)}()
  # [ 0 ] is not a basis, let's return [ ]
  if length(gens(basisSing)) == 1 && iszero(gens(basisSing)[1])
    return res
  end
  for f in gens(basisSing)
    push!(res, R(f))
  end
  return res
end

function basis_via_linear_algebra(IR::InvRing, d::Int)
  @assert d >= 0 "Degree must be non-negative"
  R = polynomial_ring(IR)
  if d == 0
    return elem_type(R)[ one(R) ]
  end

  basisSing = Singular.LibFinvar.invariant_basis(d, _action_singular(IR)...)
  res = Vector{elem_type(R)}()
  # [ 0 ] is not a basis, let's return [ ]
  if length(gens(basisSing)) == 1 && iszero(gens(basisSing)[1])
    return res
  end
  for f in gens(basisSing)
    push!(res, R(f))
  end
  return res
end

@doc Markdown.doc"""
     basis(IR::InvRing, d::Int)

Given an invariant ring `IR` and an integer `d`, return a basis for the invariants in degree `d`.

# Examples
```
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> basis(IR, 6)
4-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]^6 + x[2]^6 + x[3]^6
 x[1]^4*x[2]*x[3] + x[1]*x[2]^4*x[3] + x[1]*x[2]*x[3]^4
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
 x[1]^2*x[2]^2*x[3]^2
```
```
julia> M = matrix(GF(3), [0 1 0; -1 0 0; 0 0 -1])
[0   1   0]
[2   0   0]
[0   0   2]

julia> G = MatrixGroup(3, GF(3), [M])
Matrix group of degree 3 over Galois field with characteristic 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Galois field with characteristic 3
with generators
gfp_mat[[0 1 0; 2 0 0; 0 0 2]]

julia> basis(IR, 2)
2-element Vector{MPolyElem_dec{gfp_elem, gfp_mpoly}}:
 x[1]^2 + x[2]^2
 x[3]^2

julia> basis(IR, 3)
2-element Vector{MPolyElem_dec{gfp_elem, gfp_mpoly}}:
 x[1]^2*x[3] + 2*x[2]^2*x[3]
 x[1]*x[2]*x[3]
```
"""
basis(IR::InvRing, d::Int, algo = :default) = collect(iterate_basis(IR, d, algo))

function primary_invariants_via_singular(IR::InvRing)
  if !isdefined(IR, :primary_singular)
    IR.primary_singular = Singular.LibFinvar.primary_invariants(_action_singular(IR)...)
    P = IR.primary_singular[1]
    R = polynomial_ring(IR)
    p = Vector{elem_type(R)}()
    for i = 1:ncols(P)
      push!(p, R(P[1, i]))
    end
    IR.primary = p
  end
  return IR.primary
end

#######################################################

@doc Markdown.doc"""
    primary_invariants(IR::InvRing)

Return a system of primary invariants for `IR`.

If a system of primary invariants for `IR` is already cached, return the cached system. 
Otherwise, compute and cache such a system first.

NOTE: The primary invariants are sorted by increasing degree.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [M1, M2]);

julia> IR = invariant_ring(G);

julia> primary_invariants(IR)
3-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]*x[2]*x[3]
 x[1]^3 + x[2]^3 + x[3]^3
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
```
"""    
function primary_invariants(IR::InvRing)
  if !isdefined(IR, :primary)
    primary_invariants_via_singular(IR)
  end
  return copy(IR.primary)
end

function secondary_invariants_via_singular(IR::InvRing)
  if !isdefined(IR, :secondary)
    rey, mol = reynolds_molien_via_singular(IR)
    primary_invariants_via_singular(IR)
    P = IR.primary_singular
    if iszero(characteristic(coefficient_ring(IR)))
      S, IS = Singular.LibFinvar.secondary_char0(P[1], rey, mol)
    else
      S, IS = Singular.LibFinvar.secondary_charp(P...)
    end
    R = polynomial_ring(IR)
    s = Vector{elem_type(R)}()
    for i = 1:ncols(S)
      push!(s, R(S[1, i]))
    end
    is = Vector{elem_type(R)}()
    for i = 1:ncols(IS)
      push!(is, R(IS[1, i]))
    end
    IR.secondary = s
    IR.irreducible_secondary = is
  end
  return IR.secondary
end

#######################################################

@doc Markdown.doc"""
    secondary_invariants(IR::InvRing)

Return a system of secondary invariants for `IR` with respect to the currently cached system of primary invariants for `IR`
(if no system of primary invariants for `IR` is cached, compute and cache such a system first).

If a system of secondary invariants is already cached, return the cached system. 
Otherwise, compute and cache such a system first.

NOTE: The secondary invariants are sorted by increasing degree.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [M1, M2]);

julia> IR = invariant_ring(G);

julia> secondary_invariants(IR)
2-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 1
 x[1]^6*x[3]^3 + x[1]^3*x[2]^6 + x[2]^3*x[3]^6
```
"""    
function secondary_invariants(IR::InvRing)
  if !isdefined(IR, :secondary)
    secondary_invariants_via_singular(IR)
  end
  return copy(IR.secondary)
end

@doc Markdown.doc"""
    irreducible_secondary_invariants(IR::InvRing)

From among a system of secondary invariants for `IR` (with respect to the currently cached system of primary invariants for `IR`), return the irrreducible secondary invariants.

If a system of secondary invariants is already cached, return the irreducible ones from that system. 
Otherwise, compute and cache a system of secondary invariants first.

NOTE: A secondary invariant is *irreducible* if it cannot be written as a polynomial expession in the primary invariants and the other secondary invariants. The multiplicative unit 1 is not irreducible: It is considered to be the empty power product.

# Examples
```
julia> M = matrix(QQ, [0 -1 0 0 0; 1 -1 0 0 0; 0 0 0 0 1; 0 0 1 0 0; 0 0 0 1 0]);

julia> G = MatrixGroup(5, QQ, [M]);

julia> IR = invariant_ring(G);

julia> secondary_invariants(IR)
12-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 1
 x[1]*x[3] - x[1]*x[5] - x[2]*x[3] + x[2]*x[4]
 x[1]^2 - x[1]*x[2] + x[2]^2
 x[3]^2*x[5] + x[3]*x[4]^2 + x[4]*x[5]^2
 x[3]^3 + x[4]^3 + x[5]^3
 x[1]*x[3]*x[4] - x[1]*x[3]*x[5] - x[2]*x[3]*x[4] + x[2]*x[4]*x[5]
 x[1]*x[3]^2 - x[1]*x[4]^2 + x[2]*x[4]^2 - x[2]*x[5]^2
 x[1]^2*x[3] + x[1]^2*x[5] - 2*x[1]*x[2]*x[3] + x[2]^2*x[3] + x[2]^2*x[4]
 x[1]^2*x[3] - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[1]*x[2]*x[5] + x[2]^2*x[4]
 x[1]^3*x[3] - x[1]^3*x[5] - 2*x[1]^2*x[2]*x[3] + x[1]^2*x[2]*x[4] + x[1]^2*x[2]*x[5] + 2*x[1]*x[2]^2*x[3] - x[1]*x[2]^2*x[4] - x[1]*x[2]^2*x[5] - x[2]^3*x[3] + x[2]^3*x[4]
 x[1]^4 - 2*x[1]^3*x[2] + 3*x[1]^2*x[2]^2 - 2*x[1]*x[2]^3 + x[2]^4
 x[1]^5*x[3] - x[1]^5*x[5] - 3*x[1]^4*x[2]*x[3] + x[1]^4*x[2]*x[4] + 2*x[1]^4*x[2]*x[5] + 5*x[1]^3*x[2]^2*x[3] - 2*x[1]^3*x[2]^2*x[4] - 3*x[1]^3*x[2]^2*x[5] - 5*x[1]^2*x[2]^3*x[3] + 3*x[1]^2*x[2]^3*x[4] + 2*x[1]^2*x[2]^3*x[5] + 3*x[1]*x[2]^4*x[3] - 2*x[1]*x[2]^4*x[4] - x[1]*x[2]^4*x[5] - x[2]^5*x[3] + x[2]^5*x[4]

julia> irreducible_secondary_invariants(IR)
8-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x[1]*x[3] - x[1]*x[5] - x[2]*x[3] + x[2]*x[4]
 x[1]^2 - x[1]*x[2] + x[2]^2
 x[3]^2*x[5] + x[3]*x[4]^2 + x[4]*x[5]^2
 x[3]^3 + x[4]^3 + x[5]^3
 x[1]*x[3]*x[4] - x[1]*x[3]*x[5] - x[2]*x[3]*x[4] + x[2]*x[4]*x[5]
 x[1]*x[3]^2 - x[1]*x[4]^2 + x[2]*x[4]^2 - x[2]*x[5]^2
 x[1]^2*x[3] + x[1]^2*x[5] - 2*x[1]*x[2]*x[3] + x[2]^2*x[3] + x[2]^2*x[4]
 x[1]^2*x[3] - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[1]*x[2]*x[5] + x[2]^2*x[4]
```
"""
function irreducible_secondary_invariants(IR::InvRing)
  if !isdefined(IR, :irreducible_secondary)
    secondary_invariants_via_singular(IR)
  end
  return copy(IR.irreducible_secondary)
end

# Doesn't belong here...
# Matrices act from the left here!
function heisenberg_group(n::Int)
  K, a = CyclotomicField(n, "a")
  M1 = zero_matrix(K, n, n)
  M1[1, n] = one(K)
  for i = 2:n
    M1[i, i - 1] = one(K)
  end

  M2 = zero_matrix(K, n, n)
  M2[1, 1] = one(K)
  for i = 2:n
    M2[i, i] = M2[i - 1, i - 1]*a
  end
  return MatrixGroup(n, K, [ M1, M2 ])
end

################################################################################
#
#  Molien series
#
################################################################################

# Some functionality via Singular
function __molien_series_via_singular(IR::InvRing{T}) where {T <: Union{FlintRationalField, AnticNumberField, Nemo.GaloisField, Nemo.GaloisFmpzField}}
  return reynolds_molien_via_singular(IR)[2]
end

function _molien_series_via_singular(S::PolyRing, IR::InvRing{T}) where {T <: Union{FlintRationalField, AnticNumberField}}
  mol = __molien_series_via_singular(IR)
  # Singular does not build a new polynomial ring for the univariate Hilbert series
  # (how could it after all), but uses the first variable of the given ring.

  R = polynomial_ring(IR)
  K = coefficient_ring(IR)
  Kx, _ = PolynomialRing(K, "x", cached = false)
  # Need an extra coercion here while waiting for https://github.com/Nemocas/AbstractAlgebra.jl/pull/1009
  #return to_univariate(S, R(mol[1, 1]))//to_univariate(S, R(mol[1, 2]))
  _num = Kx(to_univariate(Kx, R(mol[1, 1])))
  _den = Kx(to_univariate(Kx, R(mol[1, 2])))
  num = change_coefficient_ring(coefficient_ring(S), _num, parent = S)
  den = change_coefficient_ring(coefficient_ring(S), _den, parent = S)
  return num//den
end

function _molien_series_via_singular(S::PolyRing, IR::InvRing{T}) where {T <: Union{Nemo.GaloisField, Nemo.GaloisFmpzField}}
  @assert !ismodular(IR)
  mol = __molien_series_via_singular(IR)

  # TODO: Write a conversion from Singular number fields to Nemo ones
  Qx, x = PolynomialRing(FlintQQ, "x", cached = false)
  K, a = number_field(Qx(Singular.n_transExt_to_spoly(Singular.modulus(coefficient_ring(parent(mol[1, 1]))))), "a", cached = false)
  R, y = PolynomialRing(K, ["y"], cached = false)
  Kx, _ = PolynomialRing(K, "x", cached = false)
  _num = to_univariate(Kx, R(mol[1, 1]))
  _den = to_univariate(Kx, R(mol[1, 2]))
  num = change_coefficient_ring(coefficient_ring(S), _num, parent = S)
  den = change_coefficient_ring(coefficient_ring(S), _den, parent = S)
  return num//den
end

function _molien_series_char0(S::PolyRing, I::InvRing)
  G = group(I)
  n = degree(G)
  Gp, GtoGp = isomorphic_group_over_finite_field(G)
  K = coefficient_ring(I)
  Kt, _ = PolynomialRing(K, "t", cached = false)
  C = conjugacy_classes(Gp)
  res = zero(FractionField(Kt))
  for c in C
    g = (GtoGp\(representative(c)))::elem_type(G)
    f = charpoly(Kt, g.elm)
    res = res + length(c)::fmpz * 1//reverse(f)
  end
  res = divexact(res, order(Gp)::fmpz)
  num = change_coefficient_ring(coefficient_ring(S),
                                numerator(res), parent = S)
  den = change_coefficient_ring(coefficient_ring(S),
                                denominator(res), parent = S)
  return num//den
end

function _molien_series_charp_nonmodular_via_gap(S::PolyRing, I::InvRing)
  G = group(I)
  @assert G isa MatrixGroup
  t = GAP.Globals.CharacterTable(G.X)
  chi = [GAP.Globals.BrauerCharacterValue(GAP.Globals.Representative(c))
         for c in GAP.Globals.ConjugacyClasses(t)]
  info = GAP.Globals.MolienSeriesInfo(GAP.Globals.MolienSeries(t,
                                                               GAP.GapObj(chi)))
  num = S(Vector{fmpz}(GAP.Globals.CoefficientsOfUnivariatePolynomial(info.numer))::Vector{fmpz})
  den = S(Vector{fmpz}(GAP.Globals.CoefficientsOfUnivariatePolynomial(info.denom))::Vector{fmpz})
  return num//den
end

function molien_series(S::PolyRing, I::InvRing)
  if characteristic(coefficient_ring(I)) == 0
    return _molien_series_char0(S, I)
  else
    if !ismodular(I)
      return _molien_series_charp_nonmodular_via_gap(S, I)
    else
      throw(NotImplemented())
    end
  end
end

@doc doc"""
    molien_series([S::PolyRing], I::InvRing)

Return the Molien series of `I` as a rational function. The invariant ring must
be non-modular.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a");

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0]);

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1]);

julia> G = MatrixGroup(3, K, [M1, M2]);

julia> IR = invariant_ring(G);

julia> molien_series(IR)
(-t^6 + t^3 - 1)//(t^9 - 3*t^6 + 3*t^3 - 1)
```
"""
function molien_series(I::InvRing)
  if !isdefined(I, :molien_series)
    S, t = PolynomialRing(QQ, "t", cached = false)
    I.molien_series = molien_series(S, I)
  end
  return I.molien_series
end


@doc Markdown.doc"""
    fundamental_invariants(IR::InvRing)

Return a system of fundamental invariants for `IR`.

If a system of fundamental invariants is already cached, return the cached system. 
Otherwise, compute and cache such a system first.

NOTE: The fundamental invariants are sorted by increasing degree.

# Examples
```jldoctest
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> fundamental_invariants(IR)
4-element Vector{MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}}:
 x[1]^3 + x[2]^3 + x[3]^3
 x[1]*x[2]*x[3]
 x[1]^6 + x[2]^6 + x[3]^6
 x[1]^6*x[3]^3 + x[1]^3*x[2]^6 + x[2]^3*x[3]^6
```
"""
function fundamental_invariants(IR::InvRing, algo::Symbol = :king)
  if !isdefined(IR, :fundamental)
    if algo == :king
      IR.fundamental = fundamental_invariants_via_king(IR)
    elseif algo == :minimal_subalgebra
      IR.fundamental = fundamental_invariants_via_minimal_subalgebra(IR)
    else
      error("Unsupported argument :$(algo) for algo")
    end
  end
  return IR.fundamental
end

function fundamental_invariants_via_minimal_subalgebra(IR::InvRing)
  V = primary_invariants(IR)
  append!(V, irreducible_secondary_invariants(IR))
  return minimal_subalgebra_generators(V)
end

function fundamental_invariants_via_king(IR::InvRing)
  rey = reynolds_via_singular(IR)
  F = Singular.LibFinvar.invariant_algebra_reynolds(rey)::Singular.smatrix{<: Singular.spoly}
  R = polynomial_ring(IR)
  f = Vector{elem_type(R)}()
  for i = 1:ncols(F)
    push!(f, R(F[1, i]))
  end
  return f
end

@doc Markdown.doc"""
    affine_algebra(IR::InvRing)

Given an invariant ring `IR` with underlying graded polynomial ring, say `R`,
return a graded affine algebra, say `A`, together with a graded algebra homomomorphism 
$A\rightarrow R$ which maps `A` isomorphically onto `IR`.

NOTE: If a system of fundamental invariants for `IR` is already cached, the function makes use of that system. 
Otherwise, such a system is computed and cached first. The algebra `A` is graded according to the 
degrees of the fundamental invariants, the modulus of `A` is generated by the algebra relations
on these invariants, and the algebra homomomorphism $A\rightarrow R$ is defined by sending the
$i$th generator of `A` to the $i$th fundamental invariant.

# Examples 
```jldoctest 
julia> K, a = CyclotomicField(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = MatrixGroup(3, K, [ M1, M2 ])
Matrix group of degree 3 over Cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring of
Matrix group of degree 3 over Cyclotomic field of order 3
with generators
AbstractAlgebra.Generic.MatSpaceElem{nf_elem}[[0 0 1; 1 0 0; 0 1 0], [1 0 0; 0 a 0; 0 0 -a-1]]

julia> affine_algebra(IR)
(Quotient of Multivariate Polynomial Ring in y[1], y[2], y[3], y[4] over Cyclotomic field of order 3 graded by
  y[1] -> [3]
  y[2] -> [3]
  y[3] -> [6]
  y[4] -> [9] by ideal(y[1]^6 - 3*y[1]^4*y[3] - 16*y[1]^3*y[2]^3 - 4*y[1]^3*y[4] + 3*y[1]^2*y[3]^2 + 24*y[1]*y[2]^3*y[3] + 4*y[1]*y[3]*y[4] + 72*y[2]^6 + 24*y[2]^3*y[4] - y[3]^3 + 8*y[4]^2), Algebra homomorphism with

domain: Quotient of Multivariate Polynomial Ring in y[1], y[2], y[3], y[4] over Cyclotomic field of order 3 graded by
  y[1] -> [3]
  y[2] -> [3]
  y[3] -> [6]
  y[4] -> [9] by ideal(y[1]^6 - 3*y[1]^4*y[3] - 16*y[1]^3*y[2]^3 - 4*y[1]^3*y[4] + 3*y[1]^2*y[3]^2 + 24*y[1]*y[2]^3*y[3] + 4*y[1]*y[3]*y[4] + 72*y[2]^6 + 24*y[2]^3*y[4] - y[3]^3 + 8*y[4]^2)

codomain: Multivariate Polynomial Ring in x[1], x[2], x[3] over Cyclotomic field of order 3 graded by
  x[1] -> [1]
  x[2] -> [1]
  x[3] -> [1]

defining images of generators: MPolyElem_dec{nf_elem, AbstractAlgebra.Generic.MPoly{nf_elem}}[x[1]^3 + x[2]^3 + x[3]^3, x[1]*x[2]*x[3], x[1]^6 + x[2]^6 + x[3]^6, x[1]^6*x[3]^3 + x[1]^3*x[2]^6 + x[2]^3*x[3]^6]
)
```
"""
function affine_algebra(IR::InvRing)
    C = polynomial_ring(IR)
    V = fundamental_invariants(IR)
    d = Int[]
    for i = 1:length(V)
        push!(d, degree(V[i])[1])
    end
    D,  = PolynomialRing(coefficient_ring(IR), "y" => (1:length(d)))
    D, = grade(D, d)
    PHI = hom(D, C, V)
    I = kernel(PHI)
    A, = quo(D, I)
    PHIBAR = hom(A, C, V)
    return (A,  PHIBAR)
end
