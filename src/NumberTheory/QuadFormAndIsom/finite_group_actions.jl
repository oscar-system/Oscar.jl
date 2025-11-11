###############################################################################
#
#  Testing whether a collection of matrices define isometries
#
###############################################################################

@doc raw"""
    is_isometry(
      Q::QuadSpace,
      f::QQMatrix;
      is_special::Bool=false,
    ) -> Bool

Return whether ``f`` defines an isometry of the quadratic space ``Q``.

# Arguments
- If `is_special` is set to `true`, return whether in addition ``f`` has
  determinant 1.
"""
function is_isometry(
    Q::Hecke.QuadSpace,
    f::QQMatrix;
    is_special::Bool=false,
  )
  !(rank(Q) == nrows(f) == ncols(f)) && return false
  if gram_matrix(Q, f) != gram_matrix(Q)
    return false
  end
  return !is_special || isone(det(f))
end

@doc raw"""
    is_isometry_list(
      Q::QuadSpace,
      V::Vector{QQMatrix};
      is_special::Bool=false,
    ) -> Bool

Return whether the matrices in ``V`` define isometries of the quadratic space
``Q``.

# Arguments
- If `is_special` is set to `true`, return whether in addition all isometries
  in ``V`` have determinant 1.
"""
function is_isometry_list(
    Q::Hecke.QuadSpace,
    V::Vector{QQMatrix};
    is_special::Bool = false,
  )
  return all(f -> is_isometry(Q, f; is_special), V)
end

@doc raw"""
    is_isometry_group(
      Q::QuadSpace
      G::MatrixGroup;
      is_special::Bool=false,
    ) -> Bool

Return whether the matrix ``G`` defines a group of isometries of the quadratic
space ``Q``.

# Arguments
- If `is_special` is set to `true`, return whether in addition all isometries
  in ``G`` have determinant 1.
"""
function is_isometry_group(
    Q::Hecke.QuadSpace,
    G::MatrixGroup;
    is_special::Bool=false,
  )
  return all(f -> is_isometry(Q, matrix(f); is_special), gens(G))
end

@doc raw"""
    is_isometry(
      L::ZZLat,
      f::QQMatrix,
      ambient_representation::Bool = false;
      is_special::Bool=false,
      is_stable::Bool=false,
    ) -> Bool

Return whether ``f`` defines an isometry of the lattice ``L``.

# Arguments
- If `is_special` is set to `true`, return whether in addition ``f``
  has determinant 1.

- If `is_stable` is set to `true`, return whether in addition ``f``
  acts trivially on the discriminant group of ``L`` (note: this
  requires ``L`` being integral).

- If `ambient_representation` is set to `true`, ``f`` is considered as a
  $\mathbb{Q}$-linear map on the ambient space of ``L``.
"""
function is_isometry(
    L::ZZLat,
    f::QQMatrix,
    ambient_representation::Bool = false;
    is_special::Bool=false,
    is_stable::Bool=false,
  )
  if is_special && !isone(det(f))
    return false
  end

  if ambient_representation
    !is_isometry(ambient_space(L), f) && return false
    ok, fL =  can_solve_with_solution(basis_matrix(L), basis_matrix(L)*f)
  else
    ok = true
    fL = f
    if fL*gram_matrix(L)*transpose(fL) != gram_matrix(L)
      return false
    end
  end
  if !ok || !isone(denominator(fL))
    return false
  elseif !(rank(L) == nrows(fL) == ncols(fL))
    return false
  elseif !is_stable
    return true
  end
  q = discriminant_group(L)
  fq = hom(q, q, elem_type(q)[q(solve(basis_matrix(L), lift(t))*fL*basis_matrix(L)) for t in gens(q)])
  return isone(matrix(fq))
end

@doc raw"""
    is_isometry_list(
      L::ZZLat,
      V::Vector{QQMatrix},
      ambient_representation::Bool = false;
      is_special::Bool=false,
      is_stable::Bool = false,
    ) -> Bool

Return whether the matrices in ``V`` define isometries of the lattice ``L``.

# Arguments
- For the optional arguments, see
  [`is_isometry(::ZZLat, ::QQMatrix)`].
"""
function is_isometry_list(
    L::ZZLat,
    V::Vector{QQMatrix},
    ambient_representation::Bool = false;
    is_special::Bool=false,
    is_stable::Bool=false,
  )
  return all(f -> is_isometry(L, f, ambient_representation; is_special, is_stable), V)
end

@doc raw"""
    is_isometry_group(
      L::ZZLat,
      G::MatrixGroup,
      ambient_representation::Bool = false;
      is_special::Bool=false,
      is_stable::Bool=false,
    ) -> Bool

Return whether the matrix ``G`` defines a group of isometries of the lattice
``L``.

# Arguments
- For the optional arguments, see
  [`is_isometry(::ZZLat, ::QQMatrix)`].
"""
function is_isometry_group(
    L::ZZLat,
    G::MatrixGroup,
    ambient_representation::Bool = false;
    is_special::Bool=false,
    is_stable::Bool=false,
  )
  return all(f -> is_isometry(L, matrix(f), ambient_representation; is_special, is_stable), gens(G))
end

###############################################################################
#
#  Change of representation
#
###############################################################################

###############################################################################
# Lattice to ambient

@doc raw"""
    extend_to_ambient_space(
      L::ZZLat,
      F::T;
      check::Bool=false,
    ) where T <: Union{QQMatrix, Vector{QQMatrix}, MatrixGroup} -> T

Given ``F`` being either:
  * a matrix with rational entries;
  * a list of matrices with rational entries;
  * a matrix group over the rationals,
defining a collection of isometries of ``L``, represented in the fixed basis of
``L``, return the same collection but represented in the standard basis of the
ambient space of ``L``.

If ``L`` is not of full rank in its ambient quadratic space ``V``, each
isometry is extended to be the identity on the quadratic space orthogonal to
$L\otimes\mathbb{Q}$ inside ``V``.

# Arguments
- If `check` is set to `true`, the function checks whether ``F`` indeed defines
  a collection of isometries of ``L``.
"""
extend_to_ambient_space

function extend_to_ambient_space(
    L::ZZLat,
    f::QQMatrix;
    check::Bool=true,
  )
  @req !check || is_isometry(L, f, false) "Matrix does not define an isometry of the lattice"
  V = ambient_space(L)
  B = basis_matrix(L)
  B2 = orthogonal_complement(V, B)
  C = vcat(B, B2)
  f_ambient = block_diagonal_matrix(QQMatrix[f, identity_matrix(QQ, nrows(B2))])
  f_ambient = inv(C)*f_ambient*C
  return f_ambient
end

function extend_to_ambient_space(
    L::ZZLat,
    F::Vector{QQMatrix};
    check::Bool=true,
  )
  @req !check || is_isometry_list(L, F, false) "Matrices do not define isometries of the lattice"
  V = ambient_space(L)
  B = basis_matrix(L)
  B2 = orthogonal_complement(V, B)
  C = vcat(B, B2)
  iC = inv(C)
  I = identity_matrix(QQ, nrows(B2))
  F_ambient = eltype(F)[]
  for f in F
    f_ambient = iC*block_diagonal_matrix(QQMatrix[f, I])*C
    push!(F_ambient, f_ambient)
  end
  return F_ambient
end

function extend_to_ambient_space(
    L::ZZLat,
    F::MatrixGroup{QQFieldElem, QQMatrix};
    check::Bool=true,
  )
  F_ambient_gens = extend_to_ambient_space(L, matrix.(gens(F)); check)
  return matrix_group(F_ambient_gens)
end

###############################################################################
# Ambient to lattice

@doc raw"""
    restrict_to_lattice(
      L::ZZLat,
      F::T;
      check::Bool=false,
    ) where T <: Union{QQMatrix, Vector{QQMatrix}, MatrixGroup} -> T

Given ``F`` being either:
  * a matrix with rational entries;
  * a list of matrices with rational entries;
  * a matrix group over the rationals,
defining a collection of isometries of the ambient space ``V`` of ``L``
preserving ``L``, represented in the standard basis of ``V``, return the
associated collection of isometries of ``L``, represented in the fixed basis of
``L``, obtained by restricting all the isometries in ``F`` to ``L``.

# Arguments
- If `check` is set to `true`, the function checks whether ``F`` indeed defines
  a collection of isometries of the ambient space of ``L`` which preserve
  ``L``.
"""
restrict_to_lattice

function restrict_to_lattice(
    L::ZZLat,
    f::QQMatrix;
    check::Bool=true,
  )
  @req !check || is_isometry(L, f, true) "Matrix does not define an isometry of the lattice"
  return solve(basis_matrix(L), basis_matrix(L)*f)
end

function restrict_to_lattice(
    L::ZZLat,
    F::Vector{QQMatrix};
    check::Bool=true,
  )
  @req !check || is_isometry_list(L, F, true) "Matrices do not define isometries of the lattice"
  F_lattice = QQMatrix[]
  for f in F
    f_lattice = solve(basis_matrix(L), basis_matrix(L)*f)
    push!(F_lattice, f_lattice)
  end
  return F_lattice
end

function restrict_to_lattice(
    L::ZZLat,
    F::MatrixGroup{QQFieldElem, QQMatrix};
    check::Bool=true,
  )
  F_lattice_gens = restrict_to_lattice(L, matrix.(gens(F)); check)
  return matrix_group(F_lattice_gens)
end

###############################################################################
#
#  Stabilizers in orthogonal group
#
###############################################################################

@doc raw"""
    stabilizer_discriminant_subgroup(
      L::ZZLat,
      G::MatrixGroup,
      H::TorQuadModule;
      pointwise::Bool=false,
      ambient_representation::Bool=true,
      check::Bool=true,
   ) -> MatrixGroup, GAPGroupHomomorphism

Given an integral lattice ``L``, a group ``G`` of isometries of ``L`` and a
submodule ``H`` of the discriminant group $D_L$ of ``L``, return the largest
subgroup of ``G`` preserving ``H`` under its action on $D_L$. This is also
the largest subgroup of ``G`` acting on the overlattice lattice of ``L``
determined by ``H``.

# Arguments
- If `pointwise` is set to `true`, return the largest subgroup of ``G``
  fixing ``H`` pointwise;

- If `ambient_representation` is set to `true`, the group ``G`` is considered
  as a group of $\mathbb{Q}$-linear automorphisms on the ambient space of
  ``L``.

- If `check` is set to `true`, the functions tests whether ``G`` defines a
  group of isometries of ``L`` and whether ``H`` is a submodule of ``D_L``.
"""
function stabilizer_discriminant_subgroup(
  L::ZZLat,
  G::MatrixGroup,
  H::TorQuadModule;
  pointwise::Bool=false,
  ambient_representation::Bool=true,
  check::Bool=true,
)
  if check
    @req is_isometry_group(L, G, ambient_representation) "Group does not define isometries of the lattice"
    @req relations(H) == L "Torsion module does not define a subgroup of the discriminant of the lattice"
    @req is_sublattice(dual(L), cover(H)) "Torsion module does not define a subgroup of the discriminant of the lattice"
  end

  discL = discriminant_representation(L, G; ambient_representation, check=false, full=false)
  GL, _ = image(discL)
  qL = domain(GL)
  j = hom(H, qL, TorQuadModuleElem[qL(lift(h)) for h in gens(H)])
  stabL, _ = stabilizer(GL, j)
  if pointwise
    _, resL = restrict_automorphism_group(stabL, j; check=false)
    stabL, _ = kernel(resL)
  end
  return preimage(discL, stabL)
end

@doc raw"""
    stabilizer_in_diagonal_action(
      L::ZZLat,
      K::ZZLat,
      N::ZZLat,
      OK::MatrixGroup,
      ON::MatrixGroup;
      check::Bool=true,
    ) -> Vector{QQMatrix}

Given a primitive extension ``K\oplus N \subseteq L`` of integral lattices,
and given two groups of isometries ``OK`` and ``ON`` of ``K`` and ``N``,
respectively, return generators for the setwise stabilizer of ``L`` in
``OK\times ON``, seen as a group of isometries of the ambient space of ``L``.

# Arguments
- If `check` is set to `true`, the function tests whether ``L`` is integral,
  ``K`` and ``N`` are orthogonal and both primitive in ``L``, and
  ``K\oplus N`` and ``L`` are rationally equal.
"""
function stabilizer_in_diagonal_action(
    L::ZZLat,
    K::ZZLat,
    N::ZZLat,
    OK::MatrixGroup,
    ON::MatrixGroup;
    check::Bool=true,
    is_finite_known=(false, false),
  )
  if check
    @req is_integral(L) "Only available for integral lattices"
    @req is_primitive(L, N) && is_primitive(L, K) "Sublattices are not both primitive in L"
    @req iszero(basis_matrix(N) * gram_matrix(ambient_space(L)) * transpose(basis_matrix(K))) "Sublattices are not orthogonal"
    @req rank(N) + rank(K) == rank(L) "Incompatible ranks"
  end

  # Can speed up kernel computations
  if is_finite_known[1]
    set_is_finite(OK, true)
  else
    is_finite(OK)
  end

  if is_finite_known[2]
    set_is_finite(ON, true)
  else
    is_finite(ON)
  end

  # Need the glue map to determine isometries of OK\times ON preserving L in
  # K^\vee\oplus N^\vee
  phi, HKinqK, HNinqN = glue_map(L, K, N)
  iphi = inv(phi)

  gen = QQMatrix[]
  discN = discriminant_representation(N, ON; check=false, full=false)
  GN, _ = image(discN)
  stabN, _ = stabilizer(GN, HNinqN)
  SN, resN = restrict_automorphism_group(stabN, HNinqN; check=false)

  # Largest subgroup of ON satisfying that {id_K}\times kerN preserves L
  kerN, _ = preimage(discN, first(kernel(resN)))
  append!(gen, matrix.(gens(kerN)))

  discK = discriminant_representation(K, OK; check=false, full=false)
  GK, _ = image(discK)
  stabK, _ = stabilizer(GK, HKinqK)
  SK, resK = restrict_automorphism_group(stabK, HKinqK; check=false)

  # Largest subgroup of OK satisfying that kerK\times {id_N} preserves L
  kerK, _ = preimage(discK, first(kernel(resK)))
  append!(gen, matrix.(gens(kerK)))

  SNphi = Oscar._orthogonal_group(domain(SK), TorQuadModuleMap[phi * hom(g) * iphi for g in gens(SN)])
  # C is isomorphic to kerK\P/kerN where P is the group we aim to construct
  C, _ = intersect(SK, SNphi)

  for g in gens(C)
    f = matrix(preimage(discK, preimage(resK, g)))
    g = matrix(preimage(discN, preimage(resN, SN(iphi * hom(g) * phi; check=false))))
    push!(gen, f*g)
  end
  unique!(gen)
  return gen
end

@doc raw"""
    maximal_extension(
      L::ZZLat,
      S::ZZLat,
      OS::MatrixGroup;
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Given a primitive sublattice ``S`` of an integral lattice ``L`` whose
orthogonal complement ``K`` is definite or of rank 2, and given a group``OS``
of isometries ``S``, return the largest group ``P`` of isometries of ``L``
preserving ``S`` and whose restriction to ``S`` is contained in ``OS``.

By definition, ``P`` is the setwise stabilizer of ``L`` in ``O(K)\times OS``
where ``O(K)`` is the full orthogonal group of ``K``.

# Arguments
-  If `check` is set to `true`, the function tests whether ``L`` is integral,
  ``S`` is primitive in ``L`` and ``K`` is definite or of rank 2.

- The keyword arguments in `kwargs` are optional arguments for the computation
  of ``O(K)``; see [`isometry_group(::ZZLat)`](@ref).
"""
function maximal_extension(
    L::ZZLat,
    S::ZZLat,
    OS::MatrixGroup;
    check::Bool=true,
    kwargs...,
  )
  T = orthogonal_submodule(L, S)
  if check
    @req is_integral(L) "Only available for integral lattices"
    @req is_primitive(L, S) "Lattice S is not primitive in L"
    @req is_definite(T) || rank(T) == 2 "Orthogonal complement is not definite or of rank 2"
  end

  if rank(T) == 0
    P = OS
  else
    OT = isometry_group(T; kwargs...)
    gen = stabilizer_in_diagonal_action(L, S, T, OS, OT)
    P = matrix_group(gen)
  end
  return P
end

@doc raw"""
    stabilizer_in_orthogonal_group(
      L::ZZLat,
      B::QQMatrix;
      stable::Bool=false,
      special::Bool=false,
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Given ``B`` a matrix with rational entries whose set of rows represent a finite
collection ``F`` of vectors in the rational span ``L\otimes \mathbb{Q}`` of
the integral lattice ``L``, return the subgroup

```math
S_F := \{ f \in O(L) | \forall v\in F : f(v) = v\} \subseteq O(L)
```
defining the joint stabilizers, in the orthogonal group of ``L``, of the
vectors in ``F``.

The implementation requires that the largest saturated submodule ``K`` of ``L``
orthogonal to ``B`` is definite or of rank 2.

# Arguments
- If `stable` is set to `true`, return the intersection of ``S_F`` with the
  stable orthogonal group $O^\#(L)$ of ``L``.

- If `special` is set to `true`, return the intersection of ``S_F`` with the
  special orthogonal group $SO(L)$ of ``L``.

- If both of the previous are `true`, then return the intersection
  $S_F\cap O^\#(L)\cap SO(L)$.

- If `check` is set to true, the function tests whether the lattice ``L`` is
  integral, the matrix ``B`` defines a set of vectors in $L\otimes \mathbb{Q}$
  and whether the lattice ``K`` is definite or of rank 2.

- The function first computes the orthogonal group of ``K``: the extra keyword
  arguments in `kwargs` are optional arguments in the computations of such a
  group (see [`isometry_group(::ZZLat)`](@ref)).

# Examples
```jldoctest
julia> A2 = root_lattice(:A, 2)
Integer lattice of rank 2 and degree 2
with gram matrix
[ 2   -1]
[-1    2]

julia> v = QQ[1 1;]
[1   1]

julia> H = stabilizer_in_orthogonal_group(A2, v)
Matrix group of degree 2
  over rational field

julia> order(H)
2
```
"""
function stabilizer_in_orthogonal_group(
    L::ZZLat,
    B::QQMatrix;
    stable::Bool=false,
    special::Bool=false,
    check::Bool=true,
    kwargs...,
  )
  K = orthogonal_submodule(L, B)
  if check
    @req is_integral(L) "Only available for integral lattices"
    @req can_solve(basis_matrix(L), B) "B does not define a subspace of the rational span of L"
    @req is_definite(K) || rank(K) == 2 "The orthogonal complement of B in L is not definite or of rank 2"
  end

  if rank(K) == 0
    return matrix_group(identity_matrix(QQ, degree(K)))
  end

  # The different orthogonal group functions for K produce isometries of
  # K which are the identity on the orthogonal quadratic space
  if stable
    # Stable isometries of K are in bijection with stable isometries of L
    # acting trivially on the orthogonal complement K^\perp_L (K is saturated)
    if special
      P, _ = Oscar._special_stable_orthogonal_group(K; kwargs...)
    else
      P, _ = stable_orthogonal_group(K; kwargs...)
    end
  else
    N = orthogonal_submodule(L, K)
    _, j, _ = glue_map(L, K, N)
    if special
      OK, _ = special_orthogonal_group(K; kwargs...)
    else
      OK = orthogonal_group(K; kwargs...)
    end
    P, _ = stabilizer_discriminant_subgroup(K, OK, domain(j); check=false, pointwise=true)
  end
  return P
end

@doc raw"""
    pointwise_stabilizer_in_orthogonal_group(
      L::ZZLat,
      S::ZZLat;
      kwargs...,
    ) -> MatrixGroup

Given a sublattice ``S`` of an integral lattice ``L``, return the subgroup

```math
\{ f \in O(L) | \forall v\in S : f(v) = v\} \subseteq O(L)
```
defining the pointwise stabilizer, in the orthogonal group of ``L``, of the
sublattice ``S`` of ``L``.

# Arguments
- The function requires that the largest saturated submodule ``K`` of ``L``
  orthogonal to ``S`` is definite or of rank 2.

- For the other keyword arguments, see
  [`stabilizer_in_orthogonal_group(::ZZLat, ::QQMatrix)`](@ref).

# Examples
```jldoctest
julia> A4 = root_lattice(:A, 4)
Integer lattice of rank 4 and degree 4
with gram matrix
[ 2   -1    0    0]
[-1    2   -1    0]
[ 0   -1    2   -1]
[ 0    0   -1    2]

julia> A2 = lattice_in_same_ambient_space(A4, QQ[1 0 0 0; 0 1 0 0])
Integer lattice of rank 2 and degree 4
with gram matrix
[ 2   -1]
[-1    2]

julia> H = pointwise_stabilizer_in_orthogonal_group(A4, A2)
Matrix group of degree 4
  over rational field

julia> order(H)
2
```
"""
function pointwise_stabilizer_in_orthogonal_group(
    L::ZZLat,
    S::ZZLat;
    kwargs...
  )
  @req ambient_space(L) === ambient_space(S) "Lattices are not contained in the same quadratic space"
  return stabilizer_in_orthogonal_group(L, basis_matrix(S); kwargs...)
end

@doc raw"""
    setwise_stabilizer_in_orthogonal_group(
      L::ZZLat,
      S::Union{QQMatrix, ZZLat};
      stable::Bool=false,
      special::Bool=false,
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Given a sublattice ``S`` of an integral lattice ``L``, or a generating set of
vectors given as rows in a matrix with rational entries, return the subgroup

```math
\{ f \in O(L) | f(S) \subseteq S\} \subseteq O(L)
```
defining the setwise stabilizer, in the orthogonal group of ``L``, of the
sublattice ``S``.

# Arguments
- The function requires that ``S`` and its orthogonal complement ``T`` in ``L``
  are both definite or of rank 2. However ``S`` need not be primitive in ``L``.

- If `check` is set to `true`, the function tests whether ``L`` is integral,
  whether ``S`` is a sublattice of ``L``, and whether ``S`` and ``T`` are
  definite or of rank 2.

- For the other keyword arguments, see
  [`stabilizer_in_orthogonal_group(::ZZLat, ::QQMatrix)`](@ref).

# Examples
```jldoctest
julia> A2 = root_lattice(:A, 2)
Integer lattice of rank 2 and degree 2
with gram matrix
[ 2   -1]
[-1    2]

julia> L = lattice_in_same_ambient_space(A2, QQ[1 1; 1 -1])
Integer lattice of rank 2 and degree 2
with gram matrix
[2   0]
[0   6]

julia> H1 = setwise_stabilizer_in_orthogonal_group(A2, L)
Matrix group of degree 2
  over rational field

julia> describe(H1)
"C2 x C2"

julia> H2 = setwise_stabilizer_in_orthogonal_group(A2, L; stable=true, special=true)
Matrix group of degree 2
  over rational field

julia> order(H2)
1
```
"""
function setwise_stabilizer_in_orthogonal_group(
    L::ZZLat,
    S::ZZLat;
    stable::Bool=false,
    special::Bool=false,
    check::Bool=true,
    kwargs...,
  )
  if check
    @req is_integral(L) "Only available for integral lattices"
    @req ambient_space(L) === ambient_space(S) "Lattices are not contained in the same quadratic space"
    @req is_sublattice(L, S) "Matrix does not define a sublattice"
    @req is_definite(S) || rank(S) == 2 "The lattice S is not definite or of rank 2"
  end

  # Compute the isometries preserving the rational span of S
  # If S is not saturated, we then filter the ones actually preserving S 
  OS = orthogonal_group(S; kwargs...)
  if !is_primitive(L, S)
    is_sat = false
    Ssat = primitive_closure(L, S)
    qS = discriminant_group(S)
    H, _ = sub(qS, elem_type(qS)[qS(vec(collect(basis_matrix(Ssat)[i:i, :]))) for i in 1:rank(S)])
    OS, _ = stabilizer_discriminant_subgroup(S, OS, H; check=false)
  else
    is_sat = true
    Ssat = S
  end

  P = maximal_extension(L, Ssat, OS; kwargs...)

  if stable
    if special
      P, _ = Oscar._special_stable_subgroup(L, P; check=false)
    else
      P, _ = stable_subgroup(L, P; check=false)
    end
  elseif special
    P, _ = special_subgroup(L, P; check=false)
  end
  return P
end

function setwise_stabilizer_in_orthogonal_group(
    L::ZZLat,
    B::QQMatrix;
    kwargs...,
  )
  S = lattice(ambient_space(L), B; isbasis=(rank(B)==nrows(B)))
  return setwise_stabilizer_in_orthogonal_group(L, S; kwargs...)
end

@doc raw"""
    pointwise_stabilizer_orthogonal_complement_in_orthogonal_group(
      L::ZZLat,
      S::Union{QQMatrix, ZZLat};
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Given a sublattice ``S`` of an integral lattice ``L``, or a generating set of
vectors given as rows in a matrix with rational entries, return the
pointwise stabilizer, in the orthogonal group of ``L``, of the orthogonal
complement of ``S`` in ``L``.

# Arguments
- The function requires that ``S`` is definite or of rank 2. However
  ``S`` need not be primitive in ``L``.

- If `check` is set to `true`, the function tests whether ``L`` is integral,
  whether ``S`` is a sublattice of ``L``, and whether ``S`` is definite or
  of rank 2.

- For the other keyword arguments, see
  [`stabilizer_in_orthogonal_group(::ZZLat, ::QQMatrix)`](@ref).

# Examples
```jldoctest
julia> A2 = root_lattice(:A, 2)
Integer lattice of rank 2 and degree 2
with gram matrix
[ 2   -1]
[-1    2]

julia> v = QQ[1 1;]
[1   1]

julia> H1 = pointwise_stabilizer_orthogonal_complement_in_orthogonal_group(A2, v)
Matrix group of degree 2
  over rational field

julia> H2 = pointwise_stabilizer_orthogonal_complement_in_orthogonal_group(A2, v; special=true)
Matrix group of degree 2
  over rational field

julia> index(H1, H2)
2
```
"""
function pointwise_stabilizer_orthogonal_complement_in_orthogonal_group(
    L::ZZLat,
    B::QQMatrix;
    kwargs...,
  )
  K = orthogonal_submodule(L, B)
  return stabilizer_in_orthogonal_group(L, basis_matrix(K); kwargs...)
end

function pointwise_stabilizer_orthogonal_complement_in_orthogonal_group(
    L::ZZLat,
    S::ZZLat;
    kwargs...,
  )
  @req (ambient_space(L) === ambient_space(S)) "Lattices are not contained in the same quadratic space"
  K = orthogonal_submodule(L, S)
  return stabilizer_in_orthogonal_group(L, basis_matrix(K); kwargs...)
end

###############################################################################
#
#  Saturation
#
###############################################################################

@doc raw"""
    saturation(
      L::ZZLat,
      G::MatrixGroup;
      ambient_representation::Bool=true,
      check::Bool=true,
      special::Bool=false,
      stable::Bool=false,
      kwargs...,
    ) -> MatrixGroup

Given a lattice ``L`` and a group of isometries ``G`` of ``L`` whose associated
coinvariant sublattice ``S`` is definite or of rank 2, return the saturation of
``G`` in the orthogonal group ``O(L)`` of ``L``, that is the pointwise
stabilizer in ``O(L)`` of its invariant sublattice.

See [`pointwise_stabilizer_in_orthogonal_group`](@ref).

# Arguments
- If `ambient_representation` is set to `true`, the group ``G`` is considered as
  a group of $\mathbb{Q}$-linear automorphisms on the ambient space of ``L``.

- If `check` is set to `true`, the function tests whether ``G`` defines a group
  of isometries of the lattice ``L`` and whether the associated coinvariant
  sublattice ``S`` is definite or of rank 2.

- If `special` is set to `true`, the function returns the saturation of ``G``
  in the special orthogonal group ``SO(L)`` of ``L``. Note that this requires
  ``G`` to consist also of special isometries. If `check` is set to `true`, the
  function tests in addition whether all isometries in ``G`` have determinant
  1.

- If `stable` is set to `true`, the function returns the saturation of ``G``
  in the stable orthogonal group ``O^\#(L)`` of ``L``. Note that this requires
  ``G`` to consist also of stable isometries. If `check` is set to `true`, the
  function tests in addition whether all isometries in ``G`` act trivially on
  the discriminant group of ``L``.

- The keyword arguments in `kwargs` are optional arguments for the computation
  of isometry group of lattices (see [`isometry_group(::ZZLat)`](@ref)).
"""
function saturation(
  L::ZZLat,
  G::MatrixGroup;
  ambient_representation::Bool=true,
  check::Bool=true,
  special::Bool=false,
  stable::Bool=false,
  kwargs...,
)
  if check
    @req is_isometry_group(L, G, ambient_representation; is_special=special, is_stable=stable) "Group does not define a group of isometries of the lattice L with the given properties"
  end

  d = denominator(scale(L))
  if d > 1
    L = rescale(L, d)
  end
  F = invariant_lattice(L, G; ambient_representation, check=false)
  Gsat = stabilizer_in_orthogonal_group(L, basis_matrix(F); check, special, stable, kwargs...)

  if !ambient_representation
    return restrict_to_lattice(L, Gsat; check=false)
  else
    return Gsat
  end
end

@doc raw"""
    saturation(
      L::ZZLat,
      G::MatrixGroup,
      H::MatrixGroup;
      ambient_representation::Bool=true,
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Given a lattice ``L`` and two finite group of isometries $H \leq G$
of ``L`` such that the coinvariant sublattice ``S`` of ``H`` is definite,
return the saturation of ``H`` in the group ``G``, that is the pointwise
stabilizer in ``G`` of its invariant sublattice associated to ``H``.

See [`saturation(::ZZLat, ::MatrixGroup)`](@ref).

# Arguments
- If `ambient_representation` is set to `true`, the groups ``G`` and ``H`` are
  considered as a group of $\mathbb{Q}$-linear automorphisms on the ambient
  space of ``L``.

- If `check` is set to `true`, the function tests whether ``G`` is a finite
  group of isometries of ``L``, whether ``H`` is contained in ``G`` and
  whether the associated coinvariant sublattice ``S`` is definite.

- The keyword arguments in `kwargs` are optional arguments for the computation
  of isometry group of lattices (see [`isometry_group(::ZZLat)`](@ref)).
"""
function saturation(
  L::ZZLat,
  G::MatrixGroup,
  H::MatrixGroup;
  ambient_representation::Bool=true,
  check::Bool=true,
  kwargs...,
)
  if check
    @req is_definite(first(coinvariant_lattice(L, H))) "Coinvariant sublattice of the second group is not definite"
    @req is_isometry_group(L, G, ambient_representation) "first group does not define a group of isometries of the lattice"
    @req is_subgroup(G, H) "Second group is not a subgroup of the first group"
    @req is_finite(G) "First group is not finite"
  end
  Hsat = saturation(L, H; ambient_representation)
  Gsat, _ = intersect(G, Hsat)
  return Gsat
end

@doc raw"""
    is_saturated_with_saturation(
      L::ZZLat,
      [G::MatrixGroup],
      H::MatrixGroup;
      kwargs...,
    ) -> Bool, MatrixGroup

Given a lattice ``L`` and a group of isometries ``H`` of ``L`` with
definite coinvariant sublattice``S`` , return whether the group ``H`` is
saturated in the orthogonal group ``O(L)`` of ``L``, meaning ``H`` is the
pointwise stabilizer in ``O(L)`` of its invariant sublattice. The second
returned output is the saturation of ``H`` in ``O(L)``.

Alternatively, one can ask whether ``H`` is saturated in a given finite group
``G`` of isometries of ``L`` containing ``H``.

# Arguments
- See [`saturation(::ZZLat, ::MatrixGroup)`](@ref) and
  [`saturation(::ZZLat, ::MatrixGroup, ::MatrixGroup)`](@ref) for the optional
  keyword arguments.
"""
function is_saturated_with_saturation(
  L::ZZLat,
  G::MatrixGroup,
  H::MatrixGroup;
  kwargs...,
)
  Hsat = saturation(L, G, H; kwargs...)
  return H == Hsat, Hsat
end

function is_saturated_with_saturation(
  L::ZZLat,
  H::MatrixGroup;
  ambient_representation::Bool=true,
  kwargs...,
)
  @req is_definite(first(coinvariant_lattice(L, H; ambient_representation))) "Associated coinvariant sublattice is not definite"
  Hsat = saturation(L, H; ambient_representation, kwargs...)
  return H == Hsat, Hsat
end
