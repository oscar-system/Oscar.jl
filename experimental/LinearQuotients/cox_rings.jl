################################################################################
#
#  Cox ring of a linear quotient
#
################################################################################

# Computes for any matrix M in mats a matrix N such that
# (e_i*N)*polys == polys[i]^M
function matrices_in_basis(mats::Vector{<: MatElem{T}}, polys::Vector{<: MPolyRingElem{T}}) where {T <: FieldElem}
  @assert !isempty(mats) && !isempty(polys)

  R = parent(polys[1])
  K = coefficient_ring(R)

  monomial_to_column = enumerate_monomials(polys)
  M = polys_to_smat(polys, monomial_to_column)

  res = SMat{elem_type(K)}[]
  for A in mats
    actionA = right_action(R, A)
    polysA = [ actionA(f) for f in polys ]
    N = polys_to_smat(polysA, monomial_to_column)

    # Little bit of a waste to recompute the rref of M all the time.
    # But I don't see how to do it better and mats should not contain many
    # elements anyways.
    fl, sol = can_solve_with_solution(M, N, side = :left)
    @assert fl

    push!(res, sol)
  end
  return res
end

# Let G < GL(V) be a finite group and GtoAbG the map G -> Ab(G) = G/[G, G].
# The invariant ring K[V]^{[G,G]} has a well-defined action by Ab(G) and is hence
# graded by the characters of Ab(G). The group of irreducible characters of the
# abelian group Ab(G) is isomorphic to Ab(G), so we assume that K[V]^{[G, G]} is
# graded by Ab(G). To choose a canonical isomorphism, we require a root of unity
# of order at least exponent(Ab(G)). This is the input zeta: it is the root of
# unity and its order (this is in the context of fixed_root_of_unity(::LinearQuotient)).
# Return the degree of the polynomial f as an element of codomain(GtoAbG),
# assuming that f is homogeneous.
function ab_g_degree(GtoAbG::Map, f::MPolyRingElem, zeta::Tuple{<:FieldElem, Int})
  AbG = codomain(GtoAbG)
  K = coefficient_ring(parent(f))
  @assert K === parent(zeta[1])

  powers_of_zeta = _powers_of_root_of_unity(zeta...)
  l = zeta[2]

  c = zeros(ZZRingElem, ngens(AbG))
  eldivs = elementary_divisors(AbG)
  for i = 1:ngens(AbG)
    fi = right_action(f, GtoAbG\AbG[i])
    q, r = divrem(fi, f)
    @assert is_zero(r) "Polynomial is not homogeneous"
    z = first(AbstractAlgebra.coefficients(q))
    @assert parent(q)(z) == q "Polynomial is not homogeneous"
    k = (powers_of_zeta[z]*eldivs[i])//l
    @assert is_integral(k)
    c[i] = numerator(k)
  end
  return AbG(c)
end

@doc raw"""
    cox_ring(L::LinearQuotient)

Return the Cox ring of the linear quotient `L` in a presentation as a graded affine
algebra (`MPolyQuoRing`) and an injective map from this ring into a polynomial ring.

By a theorem of Arzhantsev--Gaifullin [AG10](@cite) the Cox ring is graded isomorphic
to the invariant ring of the derived subgroup of `group(L)`.
We use ideas from [DK17](@cite) to find homogeneous generators of the invariant ring.
To get a map from `group(G)` to the grading group of the returned ring, use
[`class_group`](@ref).
"""
function cox_ring(L::LinearQuotient; algo_gens::Symbol = :default, algo_rels::Symbol = :groebner_basis)
  G = group(L)

  is_small_group(G) || error("Only implemented for groups not containing reflections")

  H, HtoG = derived_subgroup(G)
  A, GtoA = class_group(L)

  K = base_ring(L)
  zeta, l = fixed_root_of_unity(L)

  eldiv = elementary_divisors(A)
  @assert all(is_zero, [ mod(l, e) for e in eldiv ])
  l_div_eldiv = ZZRingElem[ div(l, e) for e in eldiv ]

  powers_of_zeta = Dict{elem_type(K), Int}()
  t = one(K)
  for i = 0:l - 1
    powers_of_zeta[t] = i
    t *= zeta
  end

  RH = invariant_ring(H)
  Q, QtoR = affine_algebra(RH, algo_gens = algo_gens, algo_rels = algo_rels)
  S = base_ring(Q)
  StoR = hom(S, codomain(QtoR), [ QtoR(f) for f in gens(Q) ])

  Rgraded = codomain(StoR)
  R = Rgraded.R
  K = coefficient_ring(R)

  invars = elem_type(R)[ StoR(x).f for x in gens(S) ]

  if istrivial(A)
    T = grade(forget_grading(S), [ zero(A) for i = 1:ngens(S) ])[1]

    relsT = elem_type(T)[ T(r.f) for r in gens(modulus(Q)) ]

    Q, TtoQ = quo(T, ideal(T, relsT))
    QtoR = hom(Q, Rgraded, [ Rgraded(f) for f in invars ])

    return Q, QtoR
  end

  gensA = dense_matrix_type(K)[]
  for g in gens(A)
    M = preimage(GtoA, g).elm
    push!(gensA, M)
  end

  # We now have generators of C[x]^H (invars) and a presentation as affine algebra
  # via Q = S/<relations>.
  # We have to transform this into generators which are homogeneous with respect to
  # the grading induced by the action of gensA.
  # For the presentation we already set up the polynomial ring T which gets a
  # grading by A in the end (we don't know the degrees yet).
  # We now compute homogeneous generators (elements of hom_invars) and remember
  # how an element of invars is written in hom_invars, that is, the image of
  # gens(S)[i] under the isomorphism S -> T that makes the obvious diagram commute
  #  S --> R, gens(S)[i] \mapsto invars[i]
  #  |
  #  |     ||
  #  v
  #  T --> R, gens(T)[i] \mapsto hom_invars[i]
  # NOTE: Whenever it says homogeneous in the comments, it is understood to mean
  # homogeneous w.r.t. the grading by the action of gensA!

  T, t = polynomial_ring(K, "t" => 1:length(invars))

  # The degrees in which C[x]^H is generated as a C-algebra
  degree_needed = falses(maximum(total_degree, invars))
  for f in invars
    @assert !iszero(total_degree(f))
    degree_needed[total_degree(f)] = true
  end

  hom_invars = elem_type(R)[]
  degrees = elem_type(A)[] # the A-degrees of hom_invars
  invars_in_hom_invars = Dict{elem_type(R), elem_type(T)}()
  C = PowerProductCache(R, invars)
  for d = 1:length(degree_needed)
    !degree_needed[d] && continue

    # Find a basis of common eigenvectors for the action of gensA on the degree d
    # component of C[x]^H. This is going to give us the homogeneous invariants.

    # A basis of the degree d component of C[x]^H
    basisd = all_power_products_of_degree!(C, d, true)

    # Whether basisd[i] is an element of invars or not
    is_gen = BitVector( isone(sum(C.exponent_vectors[f])) for f in basisd )

    gensAd = [ matrix(M) for M in matrices_in_basis(gensA, basisd) ]

    ceig = common_eigenspaces(gensAd, side = :left)

    # Now translate every eigenvector into a polynomial in R and build a basis
    # of the degree d component by homogeneous polynomials (via their preimages
    # in T).
    hom_basisd = elem_type(T)[]
    V = zero_matrix(K, 0, length(basisd)) # Collect all eigenvectors in a matrix
    for (e, v) in ceig
      for i = 1:nrows(v)

        # Whether the eigenvector only involves elements of basisd, which do NOT
        # involve any element of invars of degree d.
        from_lower_degree = true
        for j = 1:ncols(v)
          if !iszero(v[i, j]) && is_gen[j]
            from_lower_degree = false
            break
          end
        end

        if from_lower_degree
          # In this case we do not need to add an element to hom_invars, since
          # the eigenvector corresponds to powers of lower degree elements (for
          # which we already have a homogeneous generating system).
          # We still need to add it to hom_basisd though.
          g = T()
          for j = 1:ncols(v)
            if iszero(v[i, j])
              continue
            end

            m = T(v[i, j])
            expj = C.exponent_vectors[basisd[j]]
            # invars[1]^expj[1] \cdots invars[end]^expj[end] == basisd[j]
            for k = 1:length(invars)
              if iszero(expj[k])
                continue
              end
              # We know that invars[k] must be of degree lower than d, so this
              # key must exist
              m *= invars_in_hom_invars[invars[k]]^expj[k]
            end
            g += m
          end
          push!(hom_basisd, g)
        else
          push!(hom_invars, dot(v[i, :], basisd))

          # The A-degree of the new element is given by the corresponding eigenvalues
          push!(degrees, A(ZZRingElem[ div(powers_of_zeta[e[k]], l_div_eldiv[k]) for k = 1:length(e) ]))

          # We added a new generator, so the corresponding element of T is just
          # the next variable.
          push!(hom_basisd, gens(T)[length(hom_invars)])
        end
      end
      V = vcat(V, v)
    end

    # Now translate basisd into hom_basisd to be able to translate invars in
    # hom_invars.

    # We only need to do this for elements of basisd which are elements of invars.
    row_to_invar = Int[]
    M = zero_matrix(K, 0, length(basisd))
    for i = 1:length(basisd)
      if !is_gen[i]
        continue
      end
      # Find k with invars[k] == basid[i]
      k = findfirst(isequal(1), C.exponent_vectors[basisd[i]])
      push!(row_to_invar, k)

      N = zero_matrix(K, 1, length(basisd))
      N[1, i] = one(K)
      M = vcat(M, N)
    end

    fl, sol = can_solve_with_solution(V, M, side = :left)
    @assert fl

    # The i-th row of sol gives the coefficient of invars[row_to_invar[i]] in
    # the basis hom_basisd.
    for i = 1:nrows(sol)
      invars_in_hom_invars[invars[row_to_invar[i]]] = dot(sol[i, :], hom_basisd)
    end
  end

  # Now grade T, build the maps, and move the relations from S to T

  T = grade(T, degrees)[1]

  # Wants the images to be homogeneous, so have to turn off the check.
  StoT = hom(S, T, [ invars_in_hom_invars[f] for f in invars ], check = false)
  relsT = elem_type(T)[ StoT(r) for r in gens(modulus(Q)) ]

  # TODO: Funnily enough ideal(relsT) would raise an error about the generators
  # not being homogeneous (which is correct, I still want to build this ideal)
  Q, TtoQ = quo(T, ideal(T, relsT))
  QtoR = hom(Q, Rgraded, [ Rgraded(f) for f in hom_invars ])

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
function cox_ring_of_qq_factorial_terminalization(L::LinearQuotient; verbose::Bool = false)
  # TODO: can we handle this without error (in a type-stable way)?
  @assert !has_terminal_singularities(L) "Variety is already terminal"

  G = group(L)
  is_subgroup_of_sl(G) || throw(Hecke.NotImplemented())

  # Compute the Cox ring of L itself
  RVG, RVGtoR = cox_ring(L)
  verbose && @info "Computed K[V]^[G, G]"

  K = base_ring(L)
  R = codomain(RVGtoR)
  SAbG = base_ring(RVG) # RVG is a quotient ring graded by Ab(G)
  S = forget_grading(SAbG)
  StoR = hom(S, R, [ RVGtoR(gens(RVG)[i]) for i = 1:ngens(RVG) ])
  I = forget_grading(modulus(RVG))

  Sdeg, _ = grade(S, [ degree(StoR(x)) for x in gens(S) ])

  juniors = representatives_of_junior_elements(G, fixed_root_of_unity(L))
  # TODO: make sure the case where G is not generated by junior elements
  # is handled correctly in what follows

  vals = Dict{elem_type(G), Tuple}()
  for g in juniors
    vals[g] = (weights_of_action(R, g, fixed_root_of_unity(L)), monomial_valuation(R, g, fixed_root_of_unity(L)))
  end

  # Phase 1: Assure (*{i}) for all i (see [Yam18], [Sch23])
  for i = 1:length(juniors)
    verbose && @info "Starting search for junior element #$i"
    verbose && @info "Number of variables: $(ngens(S))"
    new_gens = new_generators_phase_1(StoR, I, Sdeg, SAbG, juniors[i], vals, verbose = verbose)
    while !isempty(new_gens)
      verbose && @info "Found $(length(new_gens)) new generators"
      S, Sdeg, SAbG, StoR, I = add_generators(StoR, I, new_gens, SAbG, L)
      verbose && @info "Number of variables: $(ngens(S))"
      new_gens = new_generators_phase_1(StoR, I, Sdeg, SAbG, juniors[i], vals, verbose = verbose)
    end
  end

  # Phase 2: Assure (*{1, 2}), (*{1, 3}), ..., (*{1, ..., length(juniors)})
  for i = 2:length(juniors)
    for j = 1:i - 1
      verbose && @info "Running phase 2 for i = $i and i' = $j"
      verbose && @info "Number of variables: $(ngens(S))"
      new_gens = new_generators_phase_2(StoR, I, Sdeg, SAbG, juniors, vals, i, j, verbose = verbose)
      while !isempty(new_gens)
        verbose && @info "Found $(length(new_gens)) new generators"
        S, Sdeg, SAbG, StoR, I = add_generators(StoR, I, new_gens, SAbG, L)
        verbose && @info "Number of variables: $(ngens(S))"
        new_gens = new_generators_phase_2(StoR, I, Sdeg, SAbG, juniors, vals, i, j, verbose = verbose)
      end
    end
  end

  # We now have elements of R that give rise to generators of the Cox ring as a
  # subring of a Laurent polynomial ring over R, see [Sch23, Section 6.1].

  # Transform the collected generators into Laurent polynomials
  Rt, t = LaurentPolynomialRing(forget_grading(codomain(StoR)), [ "t$i" for i = 1:length(juniors) ])
  gensRt = Vector{elem_type(Rt)}()
  degsRt = Vector{Vector{ZZRingElem}}()
  for x in gens(S)
    f = StoR(x)
    ft = Rt(forget_grading(f))
    d = Vector{ZZRingElem}()
    for i = 1:length(juniors)
      v = vals[juniors[i]][2](f)
      ft *= t[i]^v
      push!(d, v)
    end
    push!(gensRt, ft)
    push!(degsRt, d)
  end

  # Add the additional generators t_i^{-r_i} corresponding to the exceptional
  # divisors
  r = [ order(s) for s in juniors ]
  for i = 1:length(juniors)
    push!(gensRt, t[i]^-r[i])
    d = zeros(ZZRingElem, length(juniors))
    d[i] = -r[i]
    push!(degsRt, d)
  end

  # Relations
  rels = relations(StoR, I, juniors, vals, degsRt, verbose = verbose)
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
  S, t = graded_polynomial_ring(K, [ "t$i" for i = 1:nvars(R) + 1 ], w1)

  Ihom = homogenize_at_last_variable(I, S)

  StoR = hom(S, R, push!(copy(gens(R)), R()))
  return StoR(Ihom)
end

# Let I be a non-homogeneous ideal in forget_grading(SAbG).
# Return the ideal \langle h | h in I is homogeneous in SAbG \rangle.
# See [Sch23, Algorithm 6.3.2]
function group_homogeneous_ideal(I::MPolyIdeal, SAbG::MPolyDecRing)
  A = grading_group(SAbG)
  @assert isfinite(A)

  eldivs = elementary_divisors(A)
  for i = 1:ngens(A)
    I = g_homogeneous_ideal(I, ZZRingElem[ degree(x)[i] for x in gens(SAbG) ], eldivs[i])
    if iszero(I)
      break
    end
  end
  return I
end

# Helper function for group_homogeneous_ideal, see [Sch23, Algorithm 6.3.2]
function g_homogeneous_ideal(I::MPolyIdeal, weights::Vector{ZZRingElem}, order::ZZRingElem)
  R = base_ring(I)

  w = push!(copy(weights), ZZRingElem(1))
  S, t = graded_polynomial_ring(coefficient_ring(R), [ "t$i" for i = 1:nvars(R) + 1 ], w)

  Ihom = homogenize_at_last_variable(I, S)

  # Have to dance around the fact that one cannot build an inhomogeneous ideal
  # in a graded ring...
  T = forget_grading(S)
  J = ideal(T, [ forget_grading(x) for x in gens(Ihom) ]) + ideal(T, gen(T, nvars(T))^order - 1)
  phi = hom(R, T, gens(T)[1:nvars(R)])
  return preimage(phi, J)
end

# Assumes I\subseteq J and returns h_1,\dots, h_l \in J with
# J = I + (h_1,\dots, h_l).
# If both ideals are homogeneous (with respect to any weights), then the h_i
# are homogeneous too.
function as_subideal(I::MPolyIdeal{T}, J::MPolyIdeal{T}) where {T <: MPolyDecRingElem}
  R = base_ring(I)
  Q, RtoQ = quo(R, I)
  gensJQ = minimal_generating_set(ideal(Q, [ RtoQ(f) for f in gens(J) ]))
  return elem_type(R)[ RtoQ\f for f in gensJQ ]
end

# Ensures (after iterative calls) (*{k}) where k is the index of the element
# junior in the array juniors
# See [Sch23, Algorithm 6.2.2]
function new_generators_phase_1(StoR::MPolyAnyMap, I::MPolyIdeal, Sdeg::MPolyDecRing, SAbG::MPolyDecRing, junior::MatrixGroupElem, vals::Dict{<:MatrixGroupElem, <:Tuple}; verbose::Bool = false)

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

  verbose && @info "Phase 1: Computing minI"
  minI = ideal(Sdeg, minimal_parts(I, weights))

  verbose && @info "Phase 1: Computing kernel of beta"
  beta = hom(Sdeg, R, min_parts)
  J = forget_grading(kernel(beta))

  verbose && @info "Phase 1: Computing minJ"
  minJ = ideal(Sdeg, group_homogeneous_ideal(J, SAbG))

  verbose && @info "Phase 1: Representing as subideal"
  return [ forget_grading(f) for f in as_subideal(minI, minJ) ]
end

# Ensures (after iterative calls) (*A\cup {l,k}) assuming (*A\cup {l}) and (*A\cup {k}), where A is any set of indices
# See [Sch23, Algorithm 6.2.3]
function new_generators_phase_2(StoR::MPolyAnyMap, I::MPolyIdeal, Sdeg::MPolyDecRing, SAbG::MPolyDecRing, juniors::Vector{<:MatrixGroupElem}, vals::Dict{<:MatrixGroupElem, <:Tuple}, k::Int, l::Int; verbose::Bool = false)
  @assert k > l

  S = domain(StoR)
  R = codomain(StoR)

  T = S
  Ihom = I
  inds = [ l, k ]
  for i = 1:2
    verbose && @info "Phase 2: Homogenizing at $(inds[i])"
    weights = ZZRingElem[ vals[juniors[inds[i]]][2](StoR(x)) for x in gens(S) ]
    append!(weights, zeros(ZZRingElem, i))
    weights[end] = -1
    T, _ = graded_polynomial_ring(coefficient_ring(S), [ "t$j" for j = 1:nvars(S) + i ], weights)
    Ihom = homogenize_at_last_variable(Ihom, T)
  end

  verbose && @info "Phase 2: Setting up quotient ring"
  Ilk = Ihom*ideal(T, gens(T)[ngens(S) + 1]) + Ihom*ideal(T, gens(T)[ngens(S) + 2])
  Q, TtoQ = quo(T, Ilk)
  IhomQ = ideal(Q, [ TtoQ(x) for x in gens(Ihom) ])
  tlkQ = ideal(Q, [ TtoQ(gens(T)[ngens(S) + 1]), TtoQ(gens(T)[ngens(S) + 2]) ])
  verbose && @info "Phase 2: Computing intersection"
  J = intersect(IhomQ, tlkQ)
  verbose && @info "Phase 2: Simplifying"
  J = simplify(J)

  verbose && @info "Phase 2: Moving generators around"
  TtoS = hom(T, S, append!([ x for x in gens(S) ], [ one(S) for i = 1:2 ]))
  gensJ = elem_type(Q)[]
  for f in gens(J)
    if !iszero(f)
      push!(gensJ, f)
    end
  end
  new_gens = elem_type(S)[ TtoS(TtoQ\x) for x in gensJ ]
  Sk, _ = grade(S, ZZRingElem[ vals[juniors[k]][2](StoR(x)) for x in gens(S) ])
  for i = 1:length(new_gens)
    f = new_gens[i]
    hc = homogeneous_components(Sk(f))
    j = argmin(x -> x[1], keys(hc))
    new_gens[i] = forget_grading(hc[j])
  end
  return new_gens
end

# Update the data structures.
# TODO: Should there be a struct for all this stuff one moves around?
function add_generators(StoRold::MPolyAnyMap, Iold::MPolyIdeal, new_gens::Vector{<:MPolyRingElem}, SAbGold::MPolyDecRing, L::LinearQuotient)
  R = codomain(StoRold)
  K = coefficient_ring(R)
  Sold = domain(StoRold)

  _, GtoAbG = class_group(L)
  zeta = fixed_root_of_unity(L)

  new_gens_imgs = elem_type(R)[ StoRold(x) for x in new_gens]
  imgs = append!(elem_type(R)[ StoRold(x) for x in gens(Sold) ], new_gens_imgs)
  S, _ = polynomial_ring(K, "t" => 1:(ngens(Sold) + length(new_gens)))
  Sdeg, _ = grade(S, [ degree(f) for f in imgs ])

  AbG_degrees = append!([ degree(x) for x in gens(SAbGold) ], [ ab_g_degree(GtoAbG, f, zeta) for f in new_gens_imgs ])
  SAbG, _ = grade(S, AbG_degrees)

  StoR = hom(S, R, imgs)

  SoldtoS = hom(Sold, S, gens(S)[1:ngens(Sold)])

  gensI = gens(SoldtoS(Iold))
  for i = 1:length(new_gens)
    push!(gensI, SoldtoS(new_gens[i]) - gens(S)[ngens(Sold) + i])
  end
  I = ideal(S, gensI)

  return S, Sdeg, SAbG, StoR, I
end

# Relations of the Cox ring, see [Sch23, Algorithm 6.3.3].
function relations(StoR::MPolyAnyMap, I::MPolyIdeal, juniors::Vector{<:MatrixGroupElem}, vals::Dict{<:MatrixGroupElem, <:Tuple}, degsRt::Vector{Vector{ZZRingElem}}; verbose::Bool = false)
  R = codomain(StoR)
  S = domain(StoR)

  T = S
  Ihom = I
  for i = 1:length(juniors)
    verbose && @info "Relations: Homogenizing at $i"
    weights = ZZRingElem[ vals[juniors[i]][2](StoR(x)) for x in gens(S) ]
    append!(weights, zeros(ZZRingElem, i))
    weights[end] = -1
    T, _ = graded_polynomial_ring(coefficient_ring(S), [ "t$j" for j = 1:nvars(S) + i ], weights)
    Ihom = homogenize_at_last_variable(Ihom, T)
  end

  verbose && @info "Relations: Computing Gröbner basis"
  # We need homogeneous generators (with respect to all gradings) of Ihom
  gb = groebner_basis(Ihom, complete_reduction = true)
  T, _ = graded_polynomial_ring(coefficient_ring(S), append!([ "X$i" for i = 1:ngens(S) ], [ "Y$i" for i = 1:length(juniors) ]), degsRt)
  @assert ngens(T) == ngens(base_ring(Ihom))

  verbose && @info "Relations: Finishing up"
  r = [ order(s) for s in juniors ]
  rels = Vector{elem_type(T)}()
  for f in gb
    F = MPolyBuildCtx(T)
    for (c, e) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_vectors(f))
      ee = deepcopy(e)
      for i = ngens(S) + 1:ngens(T)
        @assert is_zero(mod(ee[i], r[i - ngens(S)]))
        ee[i] = div(ee[i], r[i - ngens(S)])
      end
      push_term!(F, c, ee)
    end
    push!(rels, finish(F))
  end

  return ideal(T, rels)
end
