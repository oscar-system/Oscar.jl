# Computes for any matrix M in mats a matrix N such that
# (e_i*N)*polys == polys[i]^M
function matrices_in_basis(mats::Vector{<: MatElem{T}}, polys::Vector{<: MPolyElem{T}}) where {T <: FieldElem}
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
function ab_g_degree(GtoAbG::Map, f::MPolyElem, zeta::Tuple{<:FieldElem, Int})
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

@doc Markdown.doc"""
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
