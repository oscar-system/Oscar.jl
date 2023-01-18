export secondary_invariants, irreducible_secondary_invariants,
semi_invariants, relative_invariants

function add_invariant!(C::SecondaryInvarsCache{T}, f::T, isirred::Bool, exps::Vector{Int}) where T
  push!(C.invars, f)
  push!(C.is_irreducible, isirred)
  if isirred
    for exp in C.sec_in_irred
      push!(exp, 0)
    end
  end
  push!(C.sec_in_irred, exps)
  return nothing
end

################################################################################
#
#  Modular case
#
################################################################################

# DK15, Algorithm 3.7.5
function secondary_invariants_modular(RG::InvRing)
  Rgraded = polynomial_ring(RG)
  # We have to compute a lot with these polynomials and the grading only
  # gets in the way (one cannot ask for total_degree and even if one could
  # the answer would be a GrpAbFinGenElem where it could be an Int)
  R = Rgraded.R
  K = coefficient_ring(R)

  p_invars = elem_type(R)[ f.f for f in primary_invariants(RG) ]

  s_invars = elem_type(R)[ one(R) ]
  s_invars_cache = SecondaryInvarsCache{elem_type(Rgraded)}()
  add_invariant!(s_invars_cache, one(Rgraded), false, Int[])

  # The maximal degree in which we have to look for secondary invariants by [Sym11].
  maxdeg = sum( total_degree(f) - 1 for f in p_invars )

  # Store the secondary invariants sorted by their total degree.
  # We only store the indices in s_invars_cache.invars.
  s_invars_sorted = Vector{Vector{Int}}(undef, maxdeg)
  for d = 1:maxdeg
    s_invars_sorted[d] = Int[]
  end

  # Store the indices of the irreducible secondary invariants in
  # s_invars_cache.invars .
  is_invars = Vector{Int}()

  C = PowerProductCache(R, p_invars)
  for d = 1:maxdeg
    Md = generators_for_given_degree!(C, s_invars, d, false)[1]

    # We have to find invariants of degree d which are not in the linear span of Md.
    # We could call iterate_basis(RG, d) and reduce any element by the basis matrix
    # of Md, but we can do a little bit better.

    # Build the basis matrix of Md compatible with Bd
    # We need to reverse the columns of this matrix, see below.
    Bd = iterate_basis_linear_algebra(RG, d)
    ncB1 = length(Bd.monomials_collected) + 1
    mons_to_cols = Dict{Vector{Int}, Int}(first(AbstractAlgebra.exponent_vectors(Bd.monomials_collected[i].f)) => ncB1 - i for i = 1:length(Bd.monomials_collected))
    B = BasisOfPolynomials(R, Md, mons_to_cols)

    # Do a slight detour and first try to build invariants as products of ones
    # of smaller degree. This is no speed-up (it's in fact more work), but this
    # way we have a set of irreducible secondary invariants of a much smaller
    # cardinality.
    # Use Kin07, Lemma 2: We do not need to build all power products of lower
    # degree secondary invariants, but only products of the form i*s, where
    # i is an irreducible secondary invariant of degree < d and s a secondary
    # invariant of degree d - deg(i).
    products = Set{elem_type(R)}()
    for i = 1:length(is_invars)
      f = s_invars_cache.invars[is_invars[i]].f
      @assert total_degree(f) < d
      dd = d - total_degree(f)
      for j in s_invars_sorted[dd]
        g = s_invars_cache.invars[j].f
        lp = length(products)
        fg = f*g
        push!(products, fg)
        if lp == length(products)
          # We already constructed this product from other factors.
          continue
        end

        if add_to_basis!(B, fg)
          # fg is a product of monic polynomials, so monic itself
          exp = copy(s_invars_cache.sec_in_irred[j])
          exp[i] += 1
          add_invariant!(s_invars_cache, Rgraded(fg), false, exp)
          push!(s_invars_sorted[total_degree(fg)], length(s_invars_cache.invars))
          push!(s_invars, fg)
        end
      end
    end

    # Now look for "new" invariants in this degree (the irreducible ones)
    B = B.M
    N = Bd.kernel

    # Now find all columns of N which are not in the span of B.
    # (The rows of B must span a subspace of the columns of N.)
    # B is in upper-right reduced row echelon form, N in upper-right reduced
    # column echelon form and column i of B corresponds to row nrows(N) - i + 1
    # of N in terms of basis elements since we reversed the columns of B.
    # We iterate N from the bottom right corner and B from the top left.

    c = ncols(N)
    b = 1
    for r = nrows(N):-1:1
      if c < 1
        break
      end
      if iszero(N[r, c])
        continue
      end
      # N[r, c] is a pivot of N

      if b <= nrows(B)
        k = searchsortedfirst(B.rows[b].pos, ncB1 - r)
        if k <= length(B.rows[b].pos) && B.rows[b].pos[k] == ncB1 - r
          # B has also has a pivot at the corresponding position, so we can skip
          # this column of N.
          c -= 1
          b += 1
          continue
        end
      end

      f = R()
      for j = 1:nrows(N)
        if iszero(N[j, c])
          continue
        end
        f += N[j, c]*Bd.monomials_collected[j].f
      end
      f = inv(AbstractAlgebra.leading_coefficient(f))*f
      push!(s_invars, f)
      add_invariant!(s_invars_cache, Rgraded(f), true, push!(zeros(Int, length(is_invars)), 1))
      push!(s_invars_sorted[total_degree(f)], length(s_invars_cache.invars))
      push!(is_invars, length(s_invars_cache.invars))
      c -= 1
    end
  end
  return s_invars_cache
end

################################################################################
#
#  Non modular case
#
################################################################################

# DK15, Algorithm 3.7.2 and Kin07, Section 4 "Improved new algorithm"
function secondary_invariants_nonmodular(RG::InvRing)
  @assert !is_modular(RG)
  p_invars = primary_invariants(RG)
  I = ideal_of_primary_invariants(RG)
  LI = leading_ideal(I, ordering = default_ordering(base_ring(I)))

  h = reduce_hilbert_series_by_primary_degrees(RG)

  Rgraded = polynomial_ring(RG)
  R = Rgraded.R
  # R needs to have the correct ordering for application of divrem
  @assert ordering(R) == :degrevlex

  K = coefficient_ring(R)
  s_invars_cache = SecondaryInvarsCache{elem_type(Rgraded)}()
  add_invariant!(s_invars_cache, one(Rgraded), false, Int[])

  # Store the secondary invariants sorted by their total degree.
  # We only store the indices in s_invars_cache.invars.
  s_invars_sorted = Vector{Vector{Int}}(undef, degree(h))
  for d = 1:degree(h)
    s_invars_sorted[d] = Int[]
  end

  # Store the indices of the irreducible secondary invariants in
  # s_invars_cache.invars .
  is_invars = Vector{Int}()

  # The Groebner basis should already be cached
  gbI = [ f.f for f in groebner_basis(I, ordering = degrevlex(gens(base_ring(I)))) ]

  for d = 1:degree(h)
    k = coeff(h, d) # number of invariants we need in degree d
    if iszero(k)
      continue
    end
    invars_found = 0

    gb = copy(gbI)

    # Try to build as many invariants as possible as power products of lower
    # degree ones.
    # Use Kin07, Lemma 2: We do not need to build all power products of lower
    # degree secondary invariants, but only products of the form i*s, where
    # i is an irreducible secondary invariant of degree < d and s a secondary
    # invariant of degree d - deg(i).
    products = Set{elem_type(R)}()
    for i = 1:length(is_invars)
      f = s_invars_cache.invars[is_invars[i]].f
      @assert total_degree(f) < d
      dd = d - total_degree(f)
      for j in s_invars_sorted[dd]
        g = s_invars_cache.invars[j].f
        lp = length(products)
        fg = f*g
        push!(products, fg)
        if lp == length(products)
          # We already constructed this product from other factors.
          continue
        end

        # DK15 propose to check containment via linear algebra; this approach
        # from Kin07 using d-truncated Groebner bases appears to be faster.
        _, r = divrem(fg, gb) # degrevlex from assert ordering(R) == :degrevlex
        if !is_zero(r)
          # fg is a product of monic polynomials, so monic itself
          exp = copy(s_invars_cache.sec_in_irred[j])
          exp[i] += 1
          add_invariant!(s_invars_cache, Rgraded(fg), false, exp)
          invars_found += 1
          push!(s_invars_sorted[total_degree(fg)], length(s_invars_cache.invars))
          push!(gb, r)
          invars_found == k && break
        end
      end
      invars_found == k && break
    end

    invars_found == k && continue

    # We have to find more invariants of this degree not coming from power products
    # of smaller degree ones.
    mons = all_monomials(Rgraded, d)
    for m in mons

      # Can exclude some monomials, see DK15, Remark 3.7.3 (b)
      skip = false
      for g in gens(LI)
        if mod(m.f, g.f) == 0
          skip = true
          break
        end
      end
      skip && continue

      f = reynolds_operator(RG, m).f
      if iszero(f)
        continue
      end
      _, r = divrem(f, gb)  # degrevlex from assert ordering(R) == :degrevlex
      if !is_zero(r)
        f = inv(AbstractAlgebra.leading_coefficient(f))*f
        add_invariant!(s_invars_cache, Rgraded(f), true, push!(zeros(Int, length(is_invars)), 1))
        push!(s_invars_sorted[total_degree(f)], length(s_invars_cache.invars))
        push!(is_invars, length(s_invars_cache.invars))
        push!(gb, r)
        invars_found += 1
        invars_found == k && break
      end
    end
  end

  return s_invars_cache
end

################################################################################
#
#  User functions
#
################################################################################

function _secondary_invariants(IR::InvRing)
  if isdefined(IR, :secondary)
    return nothing
  end
  if is_modular(IR)
    IR.secondary = secondary_invariants_modular(IR)
  else
    IR.secondary = secondary_invariants_nonmodular(IR)
  end
  return nothing
end

@doc Markdown.doc"""
    secondary_invariants(IR::InvRing)

Return a system of secondary invariants for `IR` as a `Vector` sorted by
increasing degree. The result is cached, so calling this function again 
with argument `IR` will be fast and give the same result.
Note that the secondary invariants are defined with respect to the currently
cached system of primary invariants for `IR` (if no system of primary invariants
for `IR` is cached, such a system is computed and cached first).

# Implemented Algorithms

For the non-modular case, the function relies on Algorithm 3.7.2 in [DK15](@cite), 
enhanced by ideas from [Kin07](@cite). In the modular case, Algorithm 3.7.5 in 
[DK15](@cite) is used.

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
 x[1]^3*x[2]^6 + x[1]^6*x[3]^3 + x[2]^3*x[3]^6
```
"""
function secondary_invariants(IR::InvRing)
  _secondary_invariants(IR)
  return copy(IR.secondary.invars)
end

@doc Markdown.doc"""
    irreducible_secondary_invariants(IR::InvRing)

Return a system of irreducible secondary invariants for `IR` as a `Vector` sorted
by increasing degree. The result is cached, so calling this function again will
be fast and give the same result.
Here, a secondary invariant is called irreducible, if it cannot be written as a
polynomial expression in the primary invariants and the other secondary
invariants.

Note that the secondary invariants and hence the irreducible secondary invariants
are defined with respect to the currently cached system of primary invariants for
`IR` (if no system of primary invariants for `IR` is cached, such a system is 
computed and cached first).

# Examples
```jldoctest
julia> M = matrix(QQ, [0 -1 0 0 0; 1 -1 0 0 0; 0 0 0 0 1; 0 0 1 0 0; 0 0 0 1 0]);

julia> G = MatrixGroup(5, QQ, [M]);

julia> IR = invariant_ring(G);

julia> secondary_invariants(IR)
12-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 1
 x[1]*x[3] - x[2]*x[3] + x[2]*x[4] - x[1]*x[5]
 x[3]^2 + x[4]^2 + x[5]^2
 x[1]^3 - 3*x[1]*x[2]^2 + x[2]^3
 x[1]^2*x[3] - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[2]^2*x[4] + x[1]*x[2]*x[5]
 x[1]*x[3]^2 - x[2]*x[3]^2 + x[2]*x[4]^2 - x[1]*x[5]^2
 x[1]^2*x[3] + x[1]^2*x[4] - 2*x[1]*x[2]*x[4] + x[2]^2*x[4] + x[2]^2*x[5]
 x[1]*x[3]*x[4] - x[2]*x[3]*x[4] - x[1]*x[3]*x[5] + x[2]*x[4]*x[5]
 x[3]*x[4]^2 + x[3]^2*x[5] + x[4]*x[5]^2
 x[1]*x[3]^3 - x[2]*x[3]^3 + x[2]*x[3]^2*x[4] + x[1]*x[3]*x[4]^2 - x[2]*x[3]*x[4]^2 + x[2]*x[4]^3 - x[1]*x[3]^2*x[5] - x[1]*x[4]^2*x[5] + x[1]*x[3]*x[5]^2 - x[2]*x[3]*x[5]^2 + x[2]*x[4]*x[5]^2 - x[1]*x[5]^3
 x[3]^4 + 2*x[3]^2*x[4]^2 + x[4]^4 + 2*x[3]^2*x[5]^2 + 2*x[4]^2*x[5]^2 + x[5]^4
 x[1]*x[3]^5 - x[2]*x[3]^5 + x[2]*x[3]^4*x[4] + 2*x[1]*x[3]^3*x[4]^2 - 2*x[2]*x[3]^3*x[4]^2 + 2*x[2]*x[3]^2*x[4]^3 + x[1]*x[3]*x[4]^4 - x[2]*x[3]*x[4]^4 + x[2]*x[4]^5 - x[1]*x[3]^4*x[5] - 2*x[1]*x[3]^2*x[4]^2*x[5] - x[1]*x[4]^4*x[5] + 2*x[1]*x[3]^3*x[5]^2 - 2*x[2]*x[3]^3*x[5]^2 + 2*x[2]*x[3]^2*x[4]*x[5]^2 + 2*x[1]*x[3]*x[4]^2*x[5]^2 - 2*x[2]*x[3]*x[4]^2*x[5]^2 + 2*x[2]*x[4]^3*x[5]^2 - 2*x[1]*x[3]^2*x[5]^3 - 2*x[1]*x[4]^2*x[5]^3 + x[1]*x[3]*x[5]^4 - x[2]*x[3]*x[5]^4 + x[2]*x[4]*x[5]^4 - x[1]*x[5]^5

julia> irreducible_secondary_invariants(IR)
8-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x[1]*x[3] - x[2]*x[3] + x[2]*x[4] - x[1]*x[5]
 x[3]^2 + x[4]^2 + x[5]^2
 x[1]^3 - 3*x[1]*x[2]^2 + x[2]^3
 x[1]^2*x[3] - x[1]*x[2]*x[3] - x[1]*x[2]*x[4] + x[2]^2*x[4] + x[1]*x[2]*x[5]
 x[1]*x[3]^2 - x[2]*x[3]^2 + x[2]*x[4]^2 - x[1]*x[5]^2
 x[1]^2*x[3] + x[1]^2*x[4] - 2*x[1]*x[2]*x[4] + x[2]^2*x[4] + x[2]^2*x[5]
 x[1]*x[3]*x[4] - x[2]*x[3]*x[4] - x[1]*x[3]*x[5] + x[2]*x[4]*x[5]
 x[3]*x[4]^2 + x[3]^2*x[5] + x[4]*x[5]^2
```
"""
function irreducible_secondary_invariants(IR::InvRing)
  _secondary_invariants(IR)
  is_invars = elem_type(polynomial_ring(IR))[]
  for i = 1:length(IR.secondary.invars)
    IR.secondary.is_irreducible[i] ? push!(is_invars, IR.secondary.invars[i]) : nothing
  end
  return is_invars
end

################################################################################
#
#  Semi-invariants / relative invariants
#
################################################################################

# Gat96, Algorithim 3.16 and DK15, Algorithm 3.7.2
@doc Markdown.doc"""
    semi_invariants(IR::InvRing, chi::GAPGroupClassFunction)
    relative_invariants(IR::InvRing, chi::GAPGroupClassFunction)

Given an irreducible character `chi` of the underlying group, return a system of
semi-invariants (or relative invariants) with respect to `chi`.
By this, we mean a set of free generators of the isotypic component of the of the
polynomial ring with respect to `chi` as a module over the algebra generated
by primary invariants for `IR`.
See also [Gat96](@cite) and [Sta79](@cite).

!!! note
    If `coefficient_ring(IR)` does not contain all character values of `chi`, an error is raised.

This function is so far only implemented in the case of characteristic zero.

# Examples
```jldoctest
julia> S2 = symmetric_group(2);

julia> RS2 = invariant_ring(S2);

julia> F = abelian_closure(QQ)[1];

julia> chi = Oscar.group_class_function(S2, [ F(sign(representative(c))) for c in conjugacy_classes(S2) ])
group_class_function(character_table(Sym( [ 1 .. 2 ] )), QQAbElem{nf_elem}[1, -1])

julia> semi_invariants(RS2, chi)
1-element Vector{MPolyElem_dec{fmpq, fmpq_mpoly}}:
 x[1] - x[2]

```
"""
function semi_invariants(RG::InvRing, chi::GAPGroupClassFunction)
  @assert is_zero(characteristic(coefficient_ring(RG)))
  @assert is_irreducible(chi)

  p_invars = primary_invariants(RG)
  I = ideal_of_primary_invariants(RG)
  LI = leading_ideal(I, ordering = default_ordering(base_ring(I)))

  h = reduce_hilbert_series_by_primary_degrees(RG, chi)

  Rgraded = polynomial_ring(RG)
  R = Rgraded.R

  B = BasisOfPolynomials(R, elem_type(R)[])

  semi_invars = elem_type(Rgraded)[]

  rey_op = reynolds_operator(RG, chi)

  for d = 0:degree(h)
    k = coeff(h, d) # number of invariants we need in degree d
    if iszero(k)
      continue
    end
    invars_found = 0

    mons = all_monomials(Rgraded, d)
    for m in mons

      # Can exclude some monomials, see DK15, Remark 3.7.3 (b)
      skip = false
      for g in gens(LI)
        if mod(m.f, g.f) == 0
          skip = true
          break
        end
      end
      skip && continue

      f = rey_op(m).f
      if iszero(f)
        continue
      end
      nf = normal_form(f, I).f
      if add_to_basis!(B, nf)
        f = inv(AbstractAlgebra.leading_coefficient(f))*f
        push!(semi_invars, Rgraded(f))
        invars_found += 1
        invars_found == k && break
      end
    end
  end

  return semi_invars
end

relative_invariants(RG::InvRing, chi::GAPGroupClassFunction) = semi_invariants(RG, chi)
