###############################################################################
#
#  Testing whether a collection of matrices define isometries
#
###############################################################################

@doc raw"""
    is_isometry(
      Q::QuadSpace,
      f::QQMatrix
    ) -> Bool

Return whether ``f`` defines an isometry of the quadratic space ``Q``.
"""
function is_isometry(
    Q::Hecke.QuadSpace,
    f::QQMatrix,
  )
  !(degree(Q) == nrows(f) == ncols(f)) && return false
  return f*gram_matrix(Q)*transpose(f) == gram_matrix(Q)
end

@doc raw"""
    is_isometry_list(
      Q::QuadSpace,
      V::Vector{QQMatrix},
    ) -> Bool

Return whether the matrices in ``V`` define isometries of the quadratic space
``Q``.
"""
function is_isometry_list(
    Q::Hecke.QuadSpace,
    V::Vector{QQMatrix},
    ambient_representation::Bool = false,
  )
  return all(f -> is_isometry(Q, f), V)
end

@doc raw"""
    is_isometry_group(
      Q::QuadSpace
      G::MatrixGroup,
    ) -> Bool

Return whether the matrix ``G`` defines a group of isometries of the quadratic
space ``Q``.
"""
function is_isometry_group(
    Q::Hecke.QuadSpace,
    G::MatrixGroup,
    ambient_representation::Bool = false,
  )
  return all(f -> is_isometry(Q, matrix(f)), gens(G))
end

@doc raw"""
    is_isometry(
      L::ZZLat,
      f::QQMatrix,
      ambient_representation::Bool = false,
    ) -> Bool

Return whether ``f`` defines an isometry of the lattice ``L``.

If `ambient_representation` is set to `true`, ``f`` is considered as a
$\mathbb{Q}$-linear map on the ambient space of ``L``.
"""
function is_isometry(
    L::ZZLat,
    f::QQMatrix,
    ambient_representation::Bool = false,
  )
  if ambient_representation
    !is_isometry(ambient_space(L), f) && return false
    return can_solve(basis_matrix(L), basis_matrix(L)*f)
  end
  !(rank(L) == nrows(f) == ncols(f)) && return false
  return f*gram_matrix(L)*transpose(f) == gram_matrix(L)
end

@doc raw"""
    is_isometry_list(
      L::ZZLat,
      V::Vector{QQMatrix},
      ambient_representation::Bool = false,
    ) -> Bool

Return whether the matrices in ``V`` define isometries of the lattice ``L``.

If `ambient_representation` is set to `true`, the matrices in ``V`` are
considered as $\mathbb{Q}$-linear maps on the ambient space of ``L``.
"""
function is_isometry_list(
    L::ZZLat,
    V::Vector{QQMatrix},
    ambient_representation::Bool = false,
  )
  return all(f -> is_isometry(L, f, ambient_representation), V)
end

@doc raw"""
    is_isometry_group(
      L::ZZLat,
      G::MatrixGroup,
      ambient_representation::Bool = false,
    ) -> Bool

Return whether the matrix ``G`` defines a group of isometries of the lattice
``L``.

If `ambient_representation` is set to `true`, the group ``G`` is considered as
a group of $\mathbb{Q}$-linear automorphisms on the ambient space of ``L``.
"""
function is_isometry_group(
    L::ZZLat,
    G::MatrixGroup,
    ambient_representation::Bool = false,
  )
  return all(f -> is_isometry(L, matrix(f), ambient_representation), gens(G))
end

###############################################################################
#
#  Change of representation
#
###############################################################################

###############################################################################
# Lattice to ambient

@doc raw"""
    representation_in_ambient_coordinates(
      L::ZZLat,
      F::T;
      check::Bool=false,
    ) where T <: Union{QQMatrix, Vector{QQMatrix}, MatrixGroup} -> T

Given ``F`` being either:
  * a matrix with rational entries;
  * a list of matrices of rational entries;
  * a matrix group over the rationals,
representing a collection of isometries of ``L`` in the fixed basis of ``L``,
return the same collection but represented in the standard basis of the ambient
space of ``L``.

If ``L`` is not of full rank in its ambient quadratic space ``V``, each
isometry is extended to be the identity on the quadratic space orthogonal to
$L\otimes\mathbb{Q}$ inside ``V``.

If `check` is set to `true`, the function checks whether ``F`` indeed defines
a collection of isometries of ``L``.
"""
representation_in_ambient_coordinates

function representation_in_ambient_coordinates(
    L::ZZLat,
    f::QQMatrix;
    check::Bool=true,
  )
  if check
    @req is_isometry(L, f, false) "Matrix must define isometry of the lattice"
  end
  V = ambient_space(L)
  B = basis_matrix(L)
  B2 = orthogonal_complement(V, B)
  C = vcat(B, B2)
  f_ambient = block_diagonal_matrix(QQMatrix[f, identity_matrix(QQ, nrows(B2))])
  f_ambient = inv(C)*f_ambient*C
  return f_ambient
end

function representation_in_ambient_coordinates(
    L::ZZLat,
    F::Vector{QQMatrix};
    check::Bool=true,
  )
  if check
    @req is_isometry_list(L, F, false) "Matrices must define isometries of the lattice"
  end
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

function representation_in_ambient_coordinates(
    L::ZZLat,
    F::MatrixGroup{QQFieldElem, QQMatrix};
    check::Bool=true,
  )
  F_ambient_gens = representation_in_ambient_coordinates(L, matrix.gens(F); check)
  return matrix_group(F_ambient_gens)
end

###############################################################################
# Ambient to lattice

@doc raw"""
    representation_in_lattice_coordinates(
      L::ZZLat,
      F::T;
      check::Bool=false,
    ) where T <: Union{QQMatrix, Vector{QQMatrix}, MatrixGroup} -> T

Given ``F`` being either:
  * a matrix with rational entries;
  * a list of matrices of rational entries;
  * a matrix group over the rationals,
representing a collection of isometries of ``L`` in the standard basis of the
ambient space of ``L``, return the same collection but represented in the fixed
basis of ``L``.

If `check` is set to `true`, the function checks whether ``F`` indeed defines
a collection of isometries of the ambient space of ``L`` which preserve ``L``.
"""
representation_in_lattice_coordinates

function representation_in_lattice_coordinates(
    L::ZZLat,
    f::QQMatrix;
    check::Bool=true,
  )
  if check
    @req is_isometry(L, f, true) "Matrix must define isometry ambient space of the lattice"
  end
  return solve(basis_matrix(L), basis_matrix(L)*f)
end

function representation_in_lattice_coordinates(
    L::ZZLat,
    F::Vector{QQMatrix};
    check::Bool=true,
  )
  if check
    @req is_isometry_list(L, F, true) "Matrices must define isometries ambient space of the lattice"
  end
  F_lattice = QQMatrix[]
  for f in F
    f_lattice = solve(basis_matrix(L), basis_matrix(L)*f)
    push!(F_lattice, f_lattice)
  end
  return F_lattice
end

function representation_in_lattice_coordinates(
    L::ZZLat,
    F::MatrixGroup{QQFieldElem, QQMatrix};
    check::Bool=true,
  )
  F_lattice_gens = representation_in_lattice_coordinates(L, matrix.(gens(F)); check)
  return matrix_group(F_lattice_gens)
end

###############################################################################
#
#  Stabilizer in diagonal actions
#
###############################################################################

@doc raw"""
    stabilizer_in_orthogonal_group(
      L::ZZLat,
      B::QQMatrix,
      stable::Bool = false;
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Return the joint stabilizer in the orthogonal group of ``L`` of the vectors in
the rational span of ``L`` given by the rows of the matrix ``B``.

The implementation requires that the largest saturated submodule ``K`` of ``L``
orthogonal to ``B`` is definite.

If `stable` is set to true, compute the joint stabilizer in the stable subgroup
$O^\#(L)$, consisting of isometries acting trivially on the discriminant group
of ``L``.

If `check` is set to true, the function tests whether ``B`` defines a set of
vectors in $L\otimes \mathbb{Q}$ and whether the lattice ``K`` is definite.

The function first computes the orthogonal group of ``K``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).
"""
function stabilizer_in_orthogonal_group(
    L::ZZLat,
    B::QQMatrix,
    stable::Bool = false;
    check::Bool=true,
    kwargs...,
  )
  K = orthogonal_submodule(L, B)
  if check
    @req can_solve(basis_matrix(L), B) "B does not define a lattice in the rational span of L"
    @req is_definite(K) "The orthogonal of B in L must be definite"
  end
  if stable
    # Stable isometries of K are in bijection with stable isometries of L
    # acting trivially on the orthogonal complement K^\perp_L (K is saturated)
    P = stable_orthogonal_group(K; kwargs...)
  else
    OK = orthogonal_group(K; kwargs...) # These are identity on K^\perp_L
    P, _ = stabilizer(OK, L, on_lattices) # There come from isometries of L
  end
  return P
end

@doc raw"""
    pointwize_stabilizer_in_orthogonal_group(
      L::ZZLat,
      B::QQMatrix,
      stable::Bool = false;
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Return the pointwize stabilizer in the orthogonal group of ``L`` of the lattice
``S``, contained in the rational span of ``L``.

The implementation requires that the largest saturated submodule ``K`` of ``L``
orthogonal to ``S`` is definite.

If `stable` is set to true, compute the pointwize stabilizer in the stable
subgroup $O^\#(L)$, consisting of isometries acting trivially on the
discriminant group of ``L``.

If `check` is set to true, the function tests whether ``S`` is contained in
$L\otimes \mathbb{Q}$ and whether the lattice ``K`` is definite.

The function first computes the orthogonal group of ``K``: the extra keyword
arguments in `kwargs` are optional arguments in the computations of such a
group (see [`isometry_group(::ZZLat)`](@ref)).
"""
function pointwize_stabilizer_in_orthogonal_group(
    L::ZZLat,
    S::ZZLat,
    stable::Bool = false;
    check::Bool=true,
    kwargs...
  )
  @req ambient_space(L) === ambient_space(S) "Lattices are not contained in the same quadratic space"
  return stabilizer_in_orthogonal_group(L, basis_matrix(S), stable; kwargs...)
end

@doc raw"""
    _stabilizer_in_diagonal_action(
      L::ZZLat,
      K::ZZLat,
      N::ZZLat,
      OK::MatrixGroup,
      ON::MatrixGroup,
    ) -> Vector{QQMatrix}

Given a primitive extension ``K\oplus N \subseteq L`` of even lattices,
and given two finite groups of isometries ``OK`` and ``ON`` of ``K`` and ``N``,
respectively, return generators for the stabilizer of ``L`` in ``OK\times ON``,
seen as a group of isometries of the ambient space of ``L``.
"""
function _stabilizer_in_diagonal_action(
    L::ZZLat,
    K::ZZLat,
    N::ZZLat,
    OK::MatrixGroup,
    ON::MatrixGroup,
  )
  # Necessary tests and can speed up kernel computations
  @assert is_finite(OK)
  @assert is_finite(ON)

  # Need the glue map to determine isometries of OK\times ON preserving L in
  # K^\vee\oplus N^\vee
  phi, HKinqK, HNinqN = glue_map(L, K, N)
  iphi = inv(phi)

  gen = QQMatrix[]
  qN = codomain(HNinqN)
  genN = ZZMatrix[matrix(hom(qN, qN, elem_type(qN)[qN(lift(a)*matrix(g)) for a in gens(qN)])) for g in gens(ON)]
  GN = Oscar._orthogonal_group(qN, genN; check=false)
  discN = hom(ON, GN, gens(GN); check=false)
  stabN, _ = stabilizer(GN, HNinqN)
  SN, resN = restrict_automorphism_group(stabN, HNinqN; check=false)
  # Largest subgroup of ON satisfying that {id_K}\times kerN preserves L
  kerN, _ = preimage(discN, first(kernel(resN)))
  append!(gen, matrix.(small_generating_set(kerN))) # Might be overkill, but we want to avoid having gazillions of generators

  qK = codomain(HKinqK)
  genK = ZZMatrix[matrix(hom(qK, qK, elem_type(qK)[qK(lift(a)*matrix(g)) for a in gens(qK)])) for g in gens(OK)]
  GK = Oscar._orthogonal_group(qK, genK; check=false)
  discK = hom(OK, GK, gens(GK); check=false)
  stabK, _ = stabilizer(GK, HKinqK)
  SK, resK = restrict_automorphism_group(stabK, HKinqK; check=false)
  # Largest subgroup of OK satisfying that kerK\times {id_N} preserves L
  kerK, _ = preimage(discK, first(kernel(resK)))
  append!(gen, matrix.(small_generating_set(kerK)))

  SNphi = Oscar._orthogonal_group(domain(SK), ZZMatrix[matrix(phi * hom(g) * iphi) for g in gens(SN)])
  # C is isomorphic to kerK\P/kerN where P is the group we aim to construct
  C, _ = _as_subgroup(SK, GAP.Globals.Intersection(GapObj(SK), GapObj(SNphi)))

  for g in small_generating_set(C)
    f = matrix(preimage(discK, preimage(resK, g)))
    g = matrix(preimage(discN, preimage(resN, SN(iphi * hom(g) * phi; check=false))))
    push!(gen, f*g)
  end
  unique!(gen)
  return gen
end

@doc raw"""
    _maximal_extension(
      L::ZZLat,
      S::ZZLat,
      OS::MatrixGroup;
      kwargs...,
    ) -> MatrixGroup, Vector{Vector{QQFieldElem}}

Given a primitive sublattice ``S`` of an even lattice ``L`` with definite
orthogonal complement ``K``, and given a finite group of isometries ``OS``
of ``S``, return the largest (finite) group ``P`` of isometries of ``L``
preserving ``S`` and whose restriction to ``S`` is contained in ``OS``.

By definition, ``P`` is the stabilizer of ``L`` in ``O(K)\times OS``
where ``O(K)`` is the full orthogonal group of ``K``.

The second output is the list of vectors of shortest length in ``K``.

The keyword arguments in `kwargs` are optional arguments for calling the
function `_isometry_group_via_decomposition(K; kwargs...)`.
"""
function _maximal_extension(
    L::ZZLat,
    S::ZZLat,
    OS::MatrixGroup;
    kwargs...,
  )
  T = orthogonal_submodule(L, T)
  @req is_definite(T) "Orthogonal complement must be definite"
  OT, svT = _isometry_group_via_decomposition(T; kwargs...)

  gen = _stabilizer_in_diagonal_action(L, S, T, OS, OT)
  GL = matrix_group(gen)
  return GL, svT
end

@doc raw"""
    stabilizer_sublattice_in_orthogonal_group(
      L::ZZLat,
      S::Union{QQMatrix, ZZLat},
      pointwize::Bool = false;
      check::Bool=true,
      kwargs...,
    ) -> MatrixGroup

Given a sublattice ``S`` of an even lattice ``L``, or a generating set of
vectors given as rows in a matrix with rational entries, return the stabilizer
of ``S`` in the orthogonal group of ``L``.

Note that the function requires the orthogonal complement of ``S`` in ``L``
is definite. However ``S`` need not be primitive in ``L``.

If `pointwize` is set to `true`, compute the pointwize stabilizer of ``S`` in
``O(L)``, i.e. the isometries of ``L`` acting trivially on ``S``. It it is set
to `false`, then the lattice ``S`` must be definite too.

If `check` is set to `true`, the function tests whether ``S`` is a sublattice
of ``L``, and whether ``S`` is definite whenever `pointwize = false`.

The keyword arguments in `kwargs` are optional arguments for the computation
of isometry group of definite lattices (see [`isometry_group(::ZZLat)`](@ref)).
"""
function stabilizer_sublattice_in_orthogonal_group(
    L::ZZLat,
    S::ZZLat,
    pointwize::Bool = false;
    check::Bool=true,
    kwargs...,
  )
  if check
    @req is_sublattice(L, S) "Matrix does not define a sublattice"
    if !pointwize
      @req is_definite(S) "The lattice defined by B must be definite for non-pointwize stabilizer"
    end
  end
  # Compute the isometries preserving the rational span of S
  # If S is not saturated, we then filter the ones actually preserving S 
  if !is_primitive(L, S)
    is_sat = false
    Ssat = primitive_closure(L, S)
  else
    is_sat = true
    Ssat = S
  end
  if pointwize
    OSsat = matrix_group(QQMatrix[identity_matrix(QQ, degree(S))])
  else
    OSsat = orthogonal_group(Ssat; kwargs...)
  end
  P, _ = _maximal_extension(L, Ssat, OSsat; kwargs...)
  if !is_sat
    P, _ = stabilizer(P, S, on_lattices)
  end
  return P
end

function stabilizer_sublattice_in_orthogonal_group(
    L::ZZLat,
    B::QQMatrix,
    pointwize::Bool = false;
    check::Bool=true,
    kwargs...,
  )
  S = lattice(ambient_space(L), B; isbasis=(rank(B)==nrows(B)))
  return stabilizer_in_orthogonal_group(L, S, pointwize; check, kwargs...)
end
