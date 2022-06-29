export fundamental_invariants

################################################################################
#
#  King's algorithm
#
################################################################################

# Computes a d-truncated Gröbner basis of I
# TODO: This should not be in this file and integrated in the general groebner_basis
# functionality
function _groebner_basis(I::MPolyIdeal, d::Int; ordering::MonomialOrdering = default_ordering(base_ring(I)))
  singular_assure(I, ordering)
  R = I.gens.Sx
  J = Singular.Ideal(R, gens(I.gens.S)...)
  G = Singular.with_degBound(d) do
        return Singular.std(J)
      end
  BA = BiPolyArray(base_ring(I), G)
  return BA
end

# [Kin13, p. 5] See also [DK15, Algorithm 3.8.2]
function fundamental_invariants_via_king(RG::InvRing, beta::Int = 0)
  @assert !is_modular(RG)

  Rgraded = polynomial_ring(RG)
  R = Rgraded.R
  # R needs to have the correct ordering for application of divrem
  @assert ordering(R) == :degrevlex
  ordR = degrevlex(gens(R))

  S = elem_type(R)[]
  G = BiPolyArray(R, elem_type(R)[])
  singular_assure(G, ordR)
  GO = elem_type(R)[]

  # TODO: Should we use [DK15, Theorem 3.2.8] to get a somewhat better bound
  # for non-cyclic groups?
  dmax = order(group(RG))
  if beta > 0 && beta < dmax
    dmax = beta
  end
  d = 1
  while d <= dmax
    if !isempty(S)
      I = ideal(R, GO)

      if total_degree(S[end]) == d - 2
        GO = groebner_basis(I, ordering = ordR)
        if is_zero(dim(I))
          mons = gens(ideal(R, Singular.kbase(I.gb[ordR].S)))
          dmax = maximum( total_degree(f) for f in mons )
          d > dmax ? break : nothing
        end
        G = I.gb[ordR]
      elseif total_degree(S[end]) == d - 1
        G = _groebner_basis(I, d, ordering = ordR)
        GO = collect(G)
      end
    end

    # TODO: Properly wrap kbase (or reimplement it; an iterator would be lovely)
    mons = gens(ideal(R, Singular.kbase(G.S, d)))
    if isempty(mons)
      break
    end

    for m in mons
      f = reynolds_operator(RG, Rgraded(m)).f
      if is_zero(f)
        continue
      end

      # Workaround while waiting for https://github.com/Nemocas/AbstractAlgebra.jl/pull/1192
      if isempty(GO)
        g = f
      else
        _, g = divrem(f, GO)
      end
      if is_zero(g)
        continue
      end

      push!(S, inv(leading_coefficient(f))*f)
      push!(GO, g)
    end

    d += 1
  end
  return [ Rgraded(f) for f in S ]
end

################################################################################
#
#  Via primary and secondary invariants
#
################################################################################

function fundamental_invariants_via_minimal_subalgebra(IR::InvRing)
  V = primary_invariants(IR)
  append!(V, irreducible_secondary_invariants(IR))
  return minimal_subalgebra_generators(V)
end

################################################################################
#
#  Singular
#
################################################################################

function fundamental_invariants_via_singular(IR::InvRing)
  @assert !is_modular(IR)
  rey = reynolds_via_singular(IR)
  F = Singular.LibFinvar.invariant_algebra_reynolds(rey)::Singular.smatrix{<: Singular.spoly}
  R = polynomial_ring(IR)
  f = Vector{elem_type(R)}()
  for i = 1:ncols(F)
    push!(f, R(F[1, i]))
  end
  return f
end

################################################################################
#
#  User functions
#
################################################################################

@doc Markdown.doc"""
    fundamental_invariants(IR::InvRing, algo::Symbol = :default; beta::Int = 0)

Return a system of fundamental invariants for `IR`.

The result is cached, so calling this function again with argument `IR` 
will be fast and give the same result.

# Implemented Algorithms

In the non-modular case the function relies on King's algorithm [Kin13](@cite) which
finds a system of fundamental invariants directly, without computing primary and
secondary invariants.
If an upper bound for the degrees of fundamental invariants better than the order
of the group is known, this can be supplied by the keyword argument `beta` and might
result in an earlier termination of the algorithm.

Alternatively, if specified by `algo = :minimal_subalgebra`, the function computes
fundamental invariants from a collection of primary and irreducible secondary
invariants using the function `minimal_subalgebra_generators`.
The optional keyword argument `beta` is ignored for this algorithm.

In the modular case, only the second method is available for theoretical reasons.

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
 x[1]^3*x[2]^6 + x[1]^6*x[3]^3 + x[2]^3*x[3]^6
```
"""
function fundamental_invariants(IR::InvRing, algo::Symbol = :default; beta::Int = 0)
  if !isdefined(IR, :fundamental)
    if algo == :default
      algo = is_modular(IR) ? :minimal_subalgebra : :king
    end

    if algo == :king
      IR.fundamental = fundamental_invariants_via_king(IR, beta = beta)
    elseif algo == :minimal_subalgebra
      IR.fundamental = fundamental_invariants_via_minimal_subalgebra(IR)
    else
      error("Unsupported argument :$(algo) for algo")
    end
  end
  return IR.fundamental
end
