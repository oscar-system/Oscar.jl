################################################################################
#
#  HomBasisBuilder (used in cox_ring(::LinearQuotient))
#
################################################################################

base_ring(HBB::HomBasisBuilder) = HBB.R

base_ring_type(
  ::Type{HomBasisBuilder{RingType,RingElemType}}
) where {RingType,RingElemType} = RingType

power_product_cache(HBB::HomBasisBuilder) = HBB.C

group(HBB::HomBasisBuilder) = HBB.G

group_action(HBB::HomBasisBuilder) = HBB.action

eigenvalues_to_group(HBB::HomBasisBuilder) = HBB.eigenvalues_to_group

preimage_ring(HBB::HomBasisBuilder) = HBB.T

# not a great name, but should only be used internally
hom_gens(HBB::HomBasisBuilder) = HBB.hom_gens
inhom_gens(HBB::HomBasisBuilder) = power_base(power_product_cache(HBB))

degrees_of_hom_gens(HBB::HomBasisBuilder) = HBB.degrees

power_base_to_hom_gens(HBB::HomBasisBuilder) = HBB.power_base_to_hom_gens

# Computes for any element g in gens(G) a matrix N such that
# (e_i*N)*polys == G_action(polys[i], g)
function action_on_basis(
  G::FinGenAbGroup, G_action::Function, polys::Vector{<:MPolyRingElem{T}}
) where {T<:FieldElem}
  @req !isempty(polys) "Input must not be empty"

  R = parent(polys[1])
  K = coefficient_ring(R)

  monomial_to_column = enumerate_monomials(polys)
  M = polys_to_smat(polys, monomial_to_column)

  res = SMat{elem_type(K)}[]
  for g in gens(G)
    polys_g = [G_action(f, g) for f in polys]
    N = polys_to_smat(polys_g, monomial_to_column)

    # Little bit of a waste to recompute the rref of M all the time.
    # But I don't see how to do it better and mats should not contain many
    # elements anyways.
    sol = solve(M, N; side=:left)

    push!(res, sol)
  end
  return res
end

function process_eigenvector!(
  HBB::HomBasisBuilder, d::Int, e::Vector{<:FieldElem}, v::MatElem, is_gen::BitVector
)
  C = power_product_cache(HBB)
  T = preimage_ring(HBB)
  basisd = all_power_products_of_degree!(C, d, true)

  res = elem_type(T)[]
  for i in 1:nrows(v)
    # We need to distinguish whether the eigenvector only involves elements
    # of basisd, which do NOT involve any element of inhom_gens(HBB) of degree d.
    if any(j -> !is_zero(v[i, j]) && is_gen[j], 1:ncols(v))
      push!(HBB.hom_gens, dot(v[i, :], basisd))

      # The G-degree of the new element is given by the corresponding eigenvalues
      push!(HBB.degrees, eigenvalues_to_group(HBB)(e))

      # We added a new generator, so the corresponding element of T is just
      # the next variable.
      push!(res, gen(T, length(hom_gens(HBB))))
    else
      # In this case we do not need to add an element to hom_gens, since
      # the eigenvector corresponds to powers of lower degree elements (for
      # which we already have a homogeneous generating system).
      g = T()
      for j in 1:ncols(v)
        is_zero(v[i, j]) && continue

        m = T(v[i, j])
        expj = get_exponent_vector(C, basisd[j])
        # inhom_gens[1]^expj[1] \cdots inhom_gens[end]^expj[end] == basisd[j]
        for k in 1:length(inhom_gens(HBB))
          is_zero(expj[k]) && continue

          # We know that inhom_gens[k] must be of degree lower than d, so this
          # key must exist
          m *= power_base_to_hom_gens(HBB)[inhom_gens(HBB)[k]]^expj[k]
        end
        g += m
      end
      push!(res, g)
    end
  end
  return res
end

# Compute a group(HBB)-homogeneous basis for the degree d component
function fill_degree!(HBB::HomBasisBuilder, d::Int)
  T = preimage_ring(HBB)
  C = power_product_cache(HBB)
  K = coefficient_ring(T)

  basisd = all_power_products_of_degree!(C, d, true)

  # Whether basisd[i] is an element of inhom_gens(HBB) or not
  is_gen = BitVector(is_one(sum(get_exponent_vector(C, f))) for f in basisd)

  gensGd = [matrix(M) for M in action_on_basis(group(HBB), group_action(HBB), basisd)]

  ceig = common_eigenspaces(gensGd; side=:left)

  # Now translate every eigenvector into a polynomial in base_ring(HBB) and build
  # a basis of the degree d component by homogeneous polynomials (via their
  # preimages in T).
  hom_basisd = elem_type(T)[]
  V = zero_matrix(K, 0, length(basisd)) # Collect all eigenvectors in a matrix
  for (e, v) in ceig
    append!(hom_basisd, process_eigenvector!(HBB, d, e, v, is_gen))
    V = vcat(V, v)
  end

  # Now translate basisd into hom_basisd to be able to translate inhom_gens in
  # hom_gens.

  # We only need to do this for elements of basisd which are elements of inhom_gens.
  row_to_gen = Int[]
  M = zero_matrix(K, 0, length(basisd))
  for i in 1:length(basisd)
    if !is_gen[i]
      continue
    end
    # Find k with inhom_gens[k] == basisd[i]
    k = findfirst(isequal(1), get_exponent_vector(C, basisd[i]))
    push!(row_to_gen, k)

    N = zero_matrix(K, 1, length(basisd))
    N[1, i] = one(K)
    M = vcat(M, N)
  end

  sol = solve(V, M; side=:left)

  # The i-th row of sol gives the coefficient of inhom_gens[row_to_gen[i]] in
  # the basis hom_basisd.
  for i in 1:nrows(sol)
    HBB.power_base_to_hom_gens[inhom_gens(HBB)[row_to_gen[i]]] = dot(sol[i, :], hom_basisd)
  end
  return nothing
end

# Let R = K[x_1, ..., x_n] be a polynomial ring and inhom_gens(HBB) polynomials in R
# generating a standard graded K-subalgebra A of R.
# Despite the name, the generators inhom_gens are assumed to be homogeneous with
# respect to the standard grading!
# Let action(HBB) be a linear action by an abelian group G on the subalgebra A.
# This action induces a grading of the group of (linear) characters Hom(G, C^\times)
# of G on A. As G is isomorphic to Hom(G, C^\times), we interpret this as a grading
# by G.
# Return generators of A which are homogeneous with respect to this grading, their
# degrees and an automorphism phi of a polynomial ring T over K in #inhom_gens
# variables such that the diagram
#   T ---> A
#   |      ||
#   v      ||
#   T ---> A
# commutes, where the upper horizontal map is the evaluation at inhom_gens, the
# lower horizontal map the evaluation at the homogeneous generators and the left
# vertical map is phi.
function homogeneous_generators(HBB::HomBasisBuilder)
  T = preimage_ring(HBB)

  # The degrees in which A (the algebra generated by inhom_gens) is generated
  # as a K-algebra
  degree_needed = falses(maximum(total_degree, inhom_gens(HBB)))
  for f in inhom_gens(HBB)
    @assert !is_zero(total_degree(f))
    degree_needed[total_degree(f)] = true
  end

  for d in 1:length(degree_needed)
    !degree_needed[d] && continue
    fill_degree!(HBB, d)
  end

  phi = hom(T, T, [power_base_to_hom_gens(HBB)[f] for f in inhom_gens(HBB)])
  return hom_gens(HBB), degrees_of_hom_gens(HBB), phi
end

################################################################################
#
#  Cox ring of a linear quotient
#
################################################################################

# Let G < GL(V) be a finite group and GtoA a surjective map G -> A to an abelian
# group A. Let H be the kernel of GtoA. The invariant ring K[V]^H has a
# well-defined action by A and is hence graded by the characters of A. The group
# of irreducible characters of the abelian group A is isomorphic to A, so we
# assume that K[V]^H is graded by A. To choose a canonical isomorphism, we
# require a root of unity of order at least exponent(A). This is the input zeta:
# it is the root of unity and its order (this is in the context of
# fixed_root_of_unity(::LinearQuotient)).
# This is intended to be used in the situation, where f is an element of the Cox
# ring of the linear quotient V/G and GtoA is the map from G to the class group
# of V/G.
# Return the degree of the polynomial f as an element of codomain(GtoA),
# assuming that f is homogeneous.
function ab_g_degree(GtoA::Map, f::MPolyRingElem, zeta::Tuple{<:FieldElem,Int})
  A = codomain(GtoA)
  K = coefficient_ring(parent(f))
  @assert K === parent(zeta[1])

  powers_of_zeta = _powers_of_root_of_unity(zeta...)
  l = zeta[2]

  c = zeros(ZZRingElem, ngens(A))
  eldivs = elementary_divisors(A)
  for i in 1:ngens(A)
    fi = right_action(f, GtoA \ A[i])
    q, r = divrem(fi, f)
    @assert is_zero(r) "Polynomial is not homogeneous"
    z = first(AbstractAlgebra.coefficients(q))
    @assert parent(q)(z) == q "Polynomial is not homogeneous"
    k = (powers_of_zeta[z] * eldivs[i])//l
    @assert is_integral(k)
    c[i] = numerator(k)
  end
  return A(c)
end

@doc raw"""
    cox_ring(L::LinearQuotient)

Return the Cox ring of the linear quotient `L` in a presentation as a graded affine
algebra (`MPolyQuoRing`) and an injective map from this ring into a polynomial ring.

Let `G = group(L)` and let `H` be the subgroup generated by the pseudo-reflections
contained in `G`. By a theorem of Arzhantsev--Gaifullin [AG10](@cite), the Cox
ring is graded isomorphic to the invariant ring of the group `H[G,G]`, where
`[G,G]` is the derived subgroup of `G`.
We use ideas from [DK17](@cite) to find homogeneous generators of the invariant ring.
To get a map from `group(G)` to the grading group of the returned ring, use
[`class_group`](@ref).
"""
function cox_ring(
  L::LinearQuotient; algo_gens::Symbol=:default, algo_rels::Symbol=:groebner_basis
)
  G = group(L)
  Grefl, GrefltoG = subgroup_of_pseudo_reflections(G)
  H, HtoG = derived_subgroup(G)
  # Compute H*Grefl
  H = matrix_group(
    base_ring(G), degree(G), vcat(map(HtoG, gens(H)), map(GrefltoG, gens(Grefl)))
  )
  A, GtoA = class_group(L)

  RH = invariant_ring(H)
  Q, QtoR = affine_algebra(RH; algo_gens=algo_gens, algo_rels=algo_rels)
  S = base_ring(Q)
  StoR = hom(S, codomain(QtoR), [QtoR(f) for f in gens(Q)])

  Rgraded = codomain(StoR)
  R = forget_grading(Rgraded)

  invars = elem_type(R)[forget_grading(StoR(x)) for x in gens(S)]

  if is_trivial(A)
    T = grade(forget_grading(S), [zero(A) for i in 1:ngens(S)])[1]

    relsT = elem_type(T)[T(forget_grading(r)) for r in gens(modulus(Q))]

    Q, TtoQ = quo(T, ideal(T, relsT))
    QtoR = hom(Q, Rgraded, [Rgraded(f) for f in invars])

    return Q, QtoR
  end

  # We now have generators of C[x]^H (invars) and a presentation as affine
  # algebra via Q = S/<relations>.
  # We have to transform this into generators which are homogeneous with respect
  # to the grading induced by the action of A.

  # The action of A on RH
  A_action = (f, g) -> right_action(f, matrix(preimage(GtoA, g)))

  # We compute eigenvectors of the linear action of A on the vector space of
  # polynomials of a fixed degree of RH. The corresponding eigenvalues are the
  # A-degrees of these eigenvectors; to translate them to group elements, we need
  # to fix a root of unity.
  zeta = fixed_root_of_unity(L)
  powers_of_zeta = _powers_of_root_of_unity(zeta...)
  l = zeta[2]
  eldiv = elementary_divisors(A)
  @assert all(is_zero, [mod(l, e) for e in eldiv])
  l_div_eldiv = ZZRingElem[div(l, e) for e in eldiv]
  function eig_to_group(e::Vector{<:FieldElem})
    return A(ZZRingElem[div(powers_of_zeta[e[k]], l_div_eldiv[k]) for k in 1:length(e)])
  end

  HBB = HomBasisBuilder(PowerProductCache(R, invars), A, A_action, eig_to_group)
  hom_invars, degrees, phi = homogeneous_generators(HBB)

  T = domain(phi)

  # Move the relations from S to T
  StoT = hom(S, T, [phi(gen(T, i)) for i in 1:ngens(T)])
  relsT = elem_type(T)[StoT(r) for r in gens(modulus(Q))]
  T = grade(T, degrees)[1]

  # The relations are in general not homogeneous
  relsT = reduce(
    vcat,
    [collect(values(homogeneous_components(T(f)))) for f in relsT];
    init=elem_type(T)[],
  )
  Q, TtoQ = quo(T, ideal(T, relsT))
  QtoR = hom(Q, Rgraded, [Rgraded(f) for f in hom_invars])

  return Q, QtoR
end

################################################################################
#
#  Cox ring of a QQ-factorial terminalization (aka minimal model) of a
#  linear quotient
#
################################################################################

@doc raw"""
    cox_ring_of_qq_factorial_terminalization(L::LinearQuotient)

Return the Cox ring of a QQ-factorial terminalization of the linear quotient `L`
in a presentation as a graded affine algebra (`MPolyQuoRing`) and an injective
map from this ring into a Laurent polynomial ring using the algorithm from
[Yam18](@cite).
"""
function cox_ring_of_qq_factorial_terminalization(L::LinearQuotient)
  # TODO: can we handle this without error (in a type-stable way)?
  @assert !has_terminal_singularities(L) "Variety is already terminal"

  G = group(L)
  is_subgroup_of_sl(G) || throw(Hecke.NotImplemented())

  # Compute the Cox ring of L itself
  RVG, RVGtoR = cox_ring(L)
  @vprint :LinearQuotients "Computed K[V]^[G, G]\n"

  K = base_ring(L)
  R = codomain(RVGtoR)
  SAbG = base_ring(RVG) # RVG is a quotient ring graded by Ab(G)
  S = forget_grading(SAbG)
  StoR = hom(S, R, [RVGtoR(x) for x in gens(RVG)])
  I = forget_grading(modulus(RVG))

  Sdeg, _ = grade(S, [degree(StoR(x)) for x in gens(S)])

  juniors = representatives_of_junior_elements(G, fixed_root_of_unity(L))
  # TODO: make sure the case where G is not generated by junior elements
  # is handled correctly in what follows

  vals = Dict{elem_type(G),Tuple}()
  for g in juniors
    vals[g] = (
      weights_of_action(R, g, fixed_root_of_unity(L)),
      monomial_valuation(R, g, fixed_root_of_unity(L)),
    )
  end

  # Phase 1: Assure (*{i}) for all i (see [Yam18], [Sch23])
  for i in 1:length(juniors)
    @vprint :LinearQuotients "Starting search for junior element #$i\n"
    @vprint :LinearQuotients "Number of variables: $(ngens(S))\n"
    new_gens = new_generators_phase_1(StoR, I, Sdeg, SAbG, juniors[i], vals)
    while !isempty(new_gens)
      @vprint :LinearQuotients "Found $(length(new_gens)) new generators\n"
      S, Sdeg, SAbG, StoR, I = add_generators(StoR, I, new_gens, SAbG, L)
      @vprint :LinearQuotients "Number of variables: $(ngens(S))\n"
      new_gens = new_generators_phase_1(StoR, I, Sdeg, SAbG, juniors[i], vals)
    end
  end

  # Phase 2: Assure (*{1, 2}), (*{1, 3}), ..., (*{1, ..., length(juniors)})
  for i in 2:length(juniors)
    for j in 1:(i - 1)
      @vprint :LinearQuotients "Running phase 2 for i = $i and i' = $j\n"
      @vprint :LinearQuotients "Number of variables: $(ngens(S))\n"
      new_gens = new_generators_phase_2(StoR, I, Sdeg, SAbG, juniors, vals, i, j)
      while !isempty(new_gens)
        @vprint :LinearQuotients "Found $(length(new_gens)) new generators\n"
        S, Sdeg, SAbG, StoR, I = add_generators(StoR, I, new_gens, SAbG, L)
        @vprint :LinearQuotients "Number of variables: $(ngens(S))\n"
        new_gens = new_generators_phase_2(StoR, I, Sdeg, SAbG, juniors, vals, i, j)
      end
    end
  end

  # We now have elements of R that give rise to generators of the Cox ring as a
  # subring of a Laurent polynomial ring over R, see [Sch23, Section 6.1].

  # Transform the collected generators into Laurent polynomials
  Rt, t = laurent_polynomial_ring(
    forget_grading(codomain(StoR)), ["t$i" for i in 1:length(juniors)]
  )
  gensRt = Vector{elem_type(Rt)}()
  degsRt = Vector{Vector{ZZRingElem}}()
  for x in gens(S)
    f = StoR(x)
    ft = Rt(forget_grading(f))
    d = Vector{ZZRingElem}()
    for i in 1:length(juniors)
      v = vals[juniors[i]][2](f)
      ft *= t[i]^v
      push!(d, v)
    end
    push!(gensRt, ft)
    push!(degsRt, d)
  end

  # Add the additional generators t_i^{-r_i} corresponding to the exceptional
  # divisors
  r = [order(s) for s in juniors]
  for i in 1:length(juniors)
    push!(gensRt, t[i]^-r[i])
    d = zeros(ZZRingElem, length(juniors))
    d[i] = -r[i]
    push!(degsRt, d)
  end

  # Relations
  rels = relations(StoR, I, juniors, vals, degsRt)
  Q, _ = quo(base_ring(rels), rels)
  QtoRt = hom(Q, Rt, gensRt)

  return Q, QtoRt
end

# Let base_ring(I) be graded with weights w and let min(h) be the homogeneous
# component of a polynomial h of minimal degree.
# Return the ideal \langle min(h) | h in I \rangle.
# See [Sch23, Algorithm 6.3.1]
function minimal_parts(I::MPolyIdeal, w::Vector{ZZRingElem})
  if all(iszero, w)
    return I
  end

  R = base_ring(I)
  K = coefficient_ring(R)

  @assert nvars(R) == length(w)

  w1 = push!(copy(w), ZZRingElem(-1))
  S, t = graded_polynomial_ring(K, "t#" => 1:(nvars(R) + 1), w1)

  Ihom = homogenize_at_last_variable(I, S)

  StoR = hom(S, R, push!(copy(gens(R)), R()))
  return StoR(Ihom)
end

# Let I be a non-homogeneous ideal in forget_grading(SAbG).
# Return the ideal \langle h | h in I is homogeneous in SAbG \rangle.
# See [Sch23, Algorithm 6.3.2]
function group_homogeneous_ideal(I::MPolyIdeal, SAbG::MPolyDecRing)
  A = grading_group(SAbG)
  @assert is_finite(A)

  eldivs = elementary_divisors(A)
  for i in 1:ngens(A)
    I = g_homogeneous_ideal(I, ZZRingElem[degree(x)[i] for x in gens(SAbG)], eldivs[i])
    if is_zero(I)
      break
    end
  end
  return I
end

# Helper function for group_homogeneous_ideal, see [Sch23, Algorithm 6.3.2]
function g_homogeneous_ideal(I::MPolyIdeal, weights::Vector{ZZRingElem}, order::ZZRingElem)
  R = base_ring(I)

  w = push!(copy(weights), ZZRingElem(1))
  S, t = graded_polynomial_ring(coefficient_ring(R), "t#" => 1:(nvars(R) + 1), w)

  Ihom = homogenize_at_last_variable(I, S)

  # Have to dance around the fact that one cannot build an inhomogeneous ideal
  # in a graded ring...
  T = forget_grading(S)
  J =
    ideal(T, [forget_grading(x) for x in gens(Ihom)]) + ideal(T, gen(T, nvars(T))^order - 1)
  phi = hom(R, T, gens(T)[1:nvars(R)])
  return preimage(phi, J)
end

# Assumes I\subseteq J and returns h_1,\dots, h_l \in J with
# J = I + (h_1,\dots, h_l).
# If both ideals are homogeneous (with respect to any weights), then the h_i
# are homogeneous too.
function as_subideal(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where {T<:MPolyDecRingElem}
  R = base_ring(I)
  Q, RtoQ = quo(R, I)
  gensJQ = minimal_generating_set(ideal(Q, [RtoQ(f) for f in gens(J)]))
  return elem_type(R)[RtoQ \ f for f in gensJQ]
end

# Ensures (after iterative calls) (*{k}) where k is the index of the element
# junior in the array juniors
# See [Sch23, Algorithm 6.2.2]
function new_generators_phase_1(
  StoR::MPolyAnyMap,
  I::MPolyIdeal,
  Sdeg::MPolyDecRing,
  SAbG::MPolyDecRing,
  junior::MatrixGroupElem,
  vals::Dict{<:MatrixGroupElem,<:Tuple},
)
  S = domain(StoR)
  R = codomain(StoR)

  val = vals[junior][2]
  RtoReig = vals[junior][1][2]
  weights = ZZRingElem[]
  min_parts = elem_type(R)[]
  for x in gens(S)
    v = val(StoR(x))
    push!(weights, v)
    push!(min_parts, forget_grading(homogeneous_component(RtoReig(StoR(x)), v))(gens(R)...))
  end

  @vprint :LinearQuotients "Phase 1: Computing minI\n"
  minI = ideal(Sdeg, minimal_parts(I, weights))

  @vprint :LinearQuotients "Phase 1: Computing kernel of beta\n"
  beta = hom(Sdeg, R, min_parts)
  J = forget_grading(kernel(beta))

  @vprint :LinearQuotients "Phase 1: Computing minJ\n"
  minJ = ideal(Sdeg, group_homogeneous_ideal(J, SAbG))

  @vprint :LinearQuotients "Phase 1: Representing as subideal\n"
  return [forget_grading(f) for f in as_subideal(minI, minJ)]
end

# Ensures (after iterative calls) (*A\cup {l,k}) assuming (*A\cup {l}) and (*A\cup {k}), where A is any set of indices
# See [Sch23, Algorithm 6.2.3]
function new_generators_phase_2(
  StoR::MPolyAnyMap,
  I::MPolyIdeal,
  Sdeg::MPolyDecRing,
  SAbG::MPolyDecRing,
  juniors::Vector{<:MatrixGroupElem},
  vals::Dict{<:MatrixGroupElem,<:Tuple},
  k::Int,
  l::Int,
)
  @req k > l "first index $k is not larger than the second index $l"

  S = domain(StoR)
  R = codomain(StoR)

  T = S
  Ihom = I
  inds = [l, k]
  for i in 1:2
    @vprint :LinearQuotients "Phase 2: Homogenizing at $(inds[i])\n"
    weights = ZZRingElem[vals[juniors[inds[i]]][2](StoR(x)) for x in gens(S)]
    append!(weights, zeros(ZZRingElem, i))
    weights[end] = -1
    T, _ = graded_polynomial_ring(
      coefficient_ring(S), ["t$j" for j in 1:(nvars(S) + i)], weights
    )
    Ihom = homogenize_at_last_variable(Ihom, T)
  end

  @vprint :LinearQuotients "Phase 2: Setting up quotient ring\n"
  Ilk = Ihom * ideal(T, gen(T, ngens(S) + 1)) + Ihom * ideal(T, gen(T, ngens(S) + 2))
  Q, TtoQ = quo(T, Ilk)
  IhomQ = ideal(Q, [TtoQ(x) for x in gens(Ihom)])
  tlkQ = ideal(Q, [TtoQ(gen(T, ngens(S) + 1)), TtoQ(gen(T, ngens(S) + 2))])
  @vprint :LinearQuotients "Phase 2: Computing intersection\n"
  J = intersect(IhomQ, tlkQ)
  @vprint :LinearQuotients "Phase 2: Simplifying\n"
  J = simplify(J)

  @vprint :LinearQuotients "Phase 2: Moving generators around\n"
  TtoS = hom(T, S, append!([x for x in gens(S)], [one(S) for i in 1:2]))
  gensJ = elem_type(Q)[]
  for f in gens(J)
    if !iszero(f)
      push!(gensJ, f)
    end
  end
  new_gens = elem_type(S)[TtoS(TtoQ \ x) for x in gensJ]
  Sk, _ = grade(S, ZZRingElem[vals[juniors[k]][2](StoR(x)) for x in gens(S)])
  for i in 1:length(new_gens)
    f = new_gens[i]
    hc = homogeneous_components(Sk(f))
    j = argmin(x -> x[1], keys(hc))
    new_gens[i] = forget_grading(hc[j])
  end
  return new_gens
end

# Update the data structures.
# TODO: Should there be a struct for all this stuff one moves around?
function add_generators(
  StoRold::MPolyAnyMap,
  Iold::MPolyIdeal,
  new_gens::Vector{<:MPolyRingElem},
  SAbGold::MPolyDecRing,
  L::LinearQuotient,
)
  R = codomain(StoRold)
  K = coefficient_ring(R)
  Sold = domain(StoRold)

  _, GtoAbG = class_group(L)
  zeta = fixed_root_of_unity(L)

  new_gens_imgs = elem_type(R)[StoRold(x) for x in new_gens]
  imgs = append!(elem_type(R)[StoRold(x) for x in gens(Sold)], new_gens_imgs)
  S, _ = polynomial_ring(K, :t => 1:(ngens(Sold) + length(new_gens)))
  Sdeg, _ = grade(S, [degree(f) for f in imgs])

  AbG_degrees = append!(
    [degree(x) for x in gens(SAbGold)],
    [ab_g_degree(GtoAbG, f, zeta) for f in new_gens_imgs],
  )
  SAbG, _ = grade(S, AbG_degrees)

  StoR = hom(S, R, imgs)

  SoldtoS = hom(Sold, S, gens(S)[1:ngens(Sold)])

  gensI = gens(SoldtoS(Iold))
  for i in 1:length(new_gens)
    push!(gensI, SoldtoS(new_gens[i]) - gen(S, ngens(Sold) + i))
  end
  I = ideal(S, gensI)

  return S, Sdeg, SAbG, StoR, I
end

# Relations of the Cox ring, see [Sch23, Algorithm 6.3.3].
function relations(
  StoR::MPolyAnyMap,
  I::MPolyIdeal,
  juniors::Vector{<:MatrixGroupElem},
  vals::Dict{<:MatrixGroupElem,<:Tuple},
  degsRt::Vector{Vector{ZZRingElem}},
)
  R = codomain(StoR)
  S = domain(StoR)

  T = S
  Ihom = I
  for i in 1:length(juniors)
    @vprint :LinearQuotients "Relations: Homogenizing at $i\n"
    weights = ZZRingElem[vals[juniors[i]][2](StoR(x)) for x in gens(S)]
    append!(weights, zeros(ZZRingElem, i))
    weights[end] = -1
    T, _ = graded_polynomial_ring(
      coefficient_ring(S), ["t$j" for j in 1:(nvars(S) + i)], weights
    )
    Ihom = homogenize_at_last_variable(Ihom, T)
  end

  @vprint :LinearQuotients "Relations: Computing Gröbner basis\n"
  # We need homogeneous generators (with respect to all gradings) of Ihom
  gb = groebner_basis(Ihom; complete_reduction=true)
  T, _ = graded_polynomial_ring(
    coefficient_ring(S),
    append!(["X$i" for i in 1:ngens(S)], ["Y$i" for i in 1:length(juniors)]),
    degsRt,
  )
  @assert ngens(T) == ngens(base_ring(Ihom))

  @vprint :LinearQuotients "Relations: Finishing up\n"
  r = [order(s) for s in juniors]
  rels = Vector{elem_type(T)}()
  for f in gb
    F = MPolyBuildCtx(T)
    for (c, e) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
      ee = deepcopy(e)
      for i in (ngens(S) + 1):ngens(T)
        @assert is_zero(mod(ee[i], r[i - ngens(S)]))
        ee[i] = div(ee[i], r[i - ngens(S)])
      end
      push_term!(F, c, ee)
    end
    push!(rels, finish(F))
  end

  return ideal(T, rels)
end
