
################################################################################
#
#  King's algorithm
#
################################################################################

# Computes a d-truncated Gröbner basis of I
# TODO: This should not be in this file and integrated in the general groebner_basis
# functionality
function _groebner_basis(
  I::MPolyIdeal, d::Int; ordering::MonomialOrdering=default_ordering(base_ring(I))
)
  sI = singular_generators(I.gens, ordering)
  R = base_ring(sI)
  J = Singular.Ideal(R, gens(sI)...)
  G = Singular.with_degBound(d) do
    return Singular.std(J)
  end
  BA = IdealGens(base_ring(I), G)
  return BA
end

# [Kin13, p. 5] See also [DK15, Algorithm 3.8.2]
function fundamental_invariants_via_king(RG::FinGroupInvarRing, beta::Int=0)
  @assert !is_modular(RG)

  Rgraded = _internal_polynomial_ring(RG)
  R = forget_grading(Rgraded)
  ordR = degrevlex(gens(R))

  S = elem_type(R)[]
  G = IdealGens(R, elem_type(R)[], ordR)
  GO = elem_type(R)[]

  g = order(Int, group(RG))
  if is_cyclic(group(RG))
    dmax = g
  else
    # We get a somewhat better bound if the group is not cyclic, see [DK15, Theorem 3.2.8]
    if iseven(g)
      dmax = floor(Int, 3//4 * g)
    else
      dmax = floor(Int, 5//8 * g)
    end
  end
  if beta > 0 && beta < dmax
    dmax = beta
  end
  d = 1
  gb_is_full = true # whether we have a full or truncated Gröbner basis
  while d <= dmax
    if length(S) >= ngens(R) && total_degree(S[end]) == d - 2
      # We haven't added any invariants in the last round, so there is a chance
      # that we are done.
      # DK15 never compute a full Gröbner basis, but experience shows that it
      # can save many truncated Gröbner bases if we don't add any new invariants
      # in the following rounds.
      I = ideal(R, GO)
      GO = gens(groebner_basis(I; ordering=ordR))
      if is_zero(dim(I))
        mons = gens(ideal(R, Singular.kbase(I.gb[ordR].gensBiPolyArray.S)))
        dmax = maximum(total_degree(f) for f in mons)
        d > dmax ? break : nothing
      end
      G = I.gb[ordR]
      gb_is_full = true
    end

    if is_zero(dimension_via_molien_series(Int, RG, d))
      d += 1
      continue
    end

    if !gb_is_full
      # We don't have a full Gröbner basis, so we compute a degree truncated one.
      G = _groebner_basis(ideal(R, GO), d; ordering=ordR)
      GO = collect(G)
    end

    # There are two possible strategies to find new candidates in degree d
    # 1) apply the Reynolds operator to monomials of R of degree d which are
    #    not divisible by any leading monomial in G.S (this is what [Kin13]
    #    proposes)
    # 2) iterate the basis of the degree d component of RG in the usual fashion
    #
    # We now estimate the complexity of both approaches using the formulas in
    # [KS99, Section 17.2].
    # However, experience shows that we have to add a generous 1/2*|G| to the
    # Reynolds operator runtime to get closer to reality.

    # Compute the monomials which would be input for the Reynolds operator
    # TODO: Properly wrap kbase (or reimplement it; an iterator would be lovely)
    mons = gens(ideal(R, Singular.kbase(singular_generators(G), d)))
    if isempty(mons) || (length(mons) == 1 && is_zero(mons[1]))
      break
    end

    # Runtime estimates, see [KS99, Section17.2]
    time_rey = length(mons) * d * order(group(RG))
    time_lin_alg = ngens(group(RG)) * length(monomials_of_degree(R, d))^2
    X = 1 / 2 * order(Int, group(RG)) # magical extra factor (see above)

    if X * time_rey < time_lin_alg
      # Reynolds approach
      invs = (
        _cast_in_internal_poly_ring(
          RG, reynolds_operator(RG, _cast_in_external_poly_ring(RG, Rgraded(m)))
        ) for m in mons
      )
    else
      # Linear algebra approach
      invs = (
        _cast_in_internal_poly_ring(RG, f) for f in iterate_basis(RG, d, :linear_algebra)
      )
    end

    for m in invs
      f = forget_grading(m)
      if is_zero(f)
        continue
      end

      _, g = divrem(f, GO)
      if is_zero(g)
        continue
      end

      push!(S, f)
      push!(GO, g)
      gb_is_full = false
    end

    d += 1
  end

  invars_cache = FundamentalInvarsCache{elem_type(Rgraded),typeof(Rgraded)}()
  polys_ext = [_cast_in_external_poly_ring(RG, Rgraded(f)) for f in S]
  # Cancelling the leading coefficient is not mathematically necessary and
  # should be done with the ordering that is used for the printing
  invars_cache.invars = [inv(AbstractAlgebra.leading_coefficient(f)) * f for f in polys_ext]
  invars_cache.via_primary_and_secondary = false
  invars_cache.S = graded_polynomial_ring(
    coefficient_ring(R), ["y$i" for i in 1:length(S)], [total_degree(f) for f in S];
    cached=false,
  )[1]
  return invars_cache
end

################################################################################
#
#  Via primary and secondary invariants
#
################################################################################

# By definition, an element of the irreducible secondary invariants cannot be
# written as a polynomial expression in the primary invariants and the other
# irreducible secondary invariants (note that this is also true in the modular
# case in OSCAR!).
# Hence if we take the primary and irreducible secondary invariants, we only have
# to make sure, that none of the *primary* invariants is in the algebra
# generated by the others.
function fundamental_invariants_via_primary_and_secondary(IR::FinGroupInvarRing)
  R = polynomial_ring(IR)
  K = coefficient_ring(R)

  invars_cache = FundamentalInvarsCache{elem_type(R),typeof(R)}()
  invars_cache.via_primary_and_secondary = true

  if isempty(irreducible_secondary_invariants(IR))
    # The easy case: there are no irreducible secondary invariants, so `IR` is
    # generated by the primary invariants

    invars_cache.invars = primary_invariants(IR)
    invars_cache.S = graded_polynomial_ring(
      K,
      ["y$i" for i in 1:length(invars_cache.invars)],
      [total_degree(forget_grading(f)) for f in invars_cache.invars];
      cached=false,
    )[1]
    invars_cache.toS = Dict{elem_type(R),elem_type(invars_cache.S)}(
      invars_cache.invars[i] => gen(invars_cache.S, i) for
      i in 1:length(invars_cache.invars)
    )

    return invars_cache
  end

  invars = append!(irreducible_secondary_invariants(IR), primary_invariants(IR))

  res, rels = _minimal_subalgebra_generators_with_relations(
    invars,
    ideal(R, [zero(R)]);
    check=false,
    start=length(irreducible_secondary_invariants(IR)),
  )

  # Sort the result by degree
  sp = sortperm(
    res; lt=(x, y) -> total_degree(forget_grading(x)) < total_degree(forget_grading(y))
  )
  res = res[sp]

  # Bookkeeping: we need to transform the relations in rels to the new ordering
  # (and potentially less variables)
  T, _ = graded_polynomial_ring(
    K, ["y$i" for i in 1:length(res)], [total_degree(forget_grading(x)) for x in res];
    cached=false,
  )

  invars_cache.invars = res
  invars_cache.S = T

  t = gens(T)[invperm(sp)]
  invars_cache.toS = Dict{elem_type(R),elem_type(T)}()
  for (f, g) in zip(invars, rels)
    invars_cache.toS[f] = g(t...)
  end
  return invars_cache
end

################################################################################
#
#  User functions
#
################################################################################

@doc raw"""
    fundamental_invariants(IR::FinGroupInvarRing, algorithm::Symbol = :default; beta::Int = 0)

Return a system of fundamental invariants for `IR`.

The result is cached, so calling this function again with argument `IR` 
will be fast and give the same result.

# Implemented Algorithms

In the non-modular case the function relies on King's algorithm [Kin13](@cite) which
finds a system of fundamental invariants directly, without computing primary and
secondary invariants.
If an upper bound for the degrees of fundamental invariants is known, this can be
supplied by the keyword argument `beta` and might result in an earlier termination
of the algorithm. By default, the algorithm uses the bounds from [DH00](@cite)
and [Sez02](@cite).

Alternatively, if specified by `algorithm = :primary_and_secondary`, the function computes
fundamental invariants from a collection of primary and irreducible secondary
invariants.
The optional keyword argument `beta` is ignored for this algorithm.

In the modular case, only the second method is available for theoretical reasons.

# Examples
```jldoctest
julia> K, a = cyclotomic_field(3, "a")
(Cyclotomic field of order 3, a)

julia> M1 = matrix(K, [0 0 1; 1 0 0; 0 1 0])
[0   0   1]
[1   0   0]
[0   1   0]

julia> M2 = matrix(K, [1 0 0; 0 a 0; 0 0 -a-1])
[1   0        0]
[0   a        0]
[0   0   -a - 1]

julia> G = matrix_group(M1, M2)
Matrix group of degree 3
  over cyclotomic field of order 3

julia> IR = invariant_ring(G)
Invariant ring
  of matrix group of degree 3 over K

julia> fundamental_invariants(IR)
4-element Vector{MPolyDecRingElem{AbsSimpleNumFieldElem, AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}}:
 x[1]*x[2]*x[3]
 x[1]^3 + x[2]^3 + x[3]^3
 x[1]^3*x[2]^3 + x[1]^3*x[3]^3 + x[2]^3*x[3]^3
 x[1]^3*x[2]^6 + x[1]^6*x[3]^3 + x[2]^3*x[3]^6
```
"""
function fundamental_invariants(
  IR::FinGroupInvarRing, algorithm::Symbol=:default; beta::Int=0
)
  if !isdefined(IR, :fundamental)
    if algorithm == :default
      algorithm = is_modular(IR) ? :primary_and_secondary : :king
    end

    if algorithm == :king
      IR.fundamental = fundamental_invariants_via_king(IR, beta)
    elseif algorithm == :primary_and_secondary
      IR.fundamental = fundamental_invariants_via_primary_and_secondary(IR)
    else
      error("Unsupported argument :$(algorithm) for algorithm")
    end
  end
  return copy(IR.fundamental.invars)
end
