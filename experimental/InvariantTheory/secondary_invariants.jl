################################################################################
#
#  Helpers
#
################################################################################

# Add a base element to the cache.
function add_base_element!(C::PowerProductCache{S, T}, f::T, remember_exponents::Bool = true) where {S, T}
  @assert !iszero(f)

  push!(C.base, f)
  if remember_exponents
    for (g, exps) in C.exponent_vectors
      push!(exps, 0)
    end
  end
  for d in sort!(collect(keys(C.power_products)))
    extend!(C, d, length(C.base), remember_exponents)
  end
  return nothing
end

# Extend the already computed power products by powers of C.base[i].
# Only used internally.
function extend!(C::PowerProductCache{S, T}, d::Int, i::Int, remember_exponents::Bool = true) where {S, T <: MPolyElem}
  @assert d >= 0
  @assert 1 <= i && i <= length(C.base)

  if d <= 0
    return nothing
  end

  exps = zeros(Int, length(C.base))

  di = total_degree(C.base[i])

  # Larger degree, ignore
  if di > d
    return nothing
  end

  # Exactly matching degree
  if di == d
    push!(C.power_products[d], C.base[i])
    if remember_exponents
      exps[i] = 1
      C.exponent_vectors[C.base[i]] = copy(exps)
      exps[i] = 0
    end
    C.last_factor[C.base[i]] = i
    return nothing
  end

  # Build all products with elements of degree d - total_degree(C.base[i])
  dd = d - di
  all_power_products_of_degree!(C, dd, remember_exponents)
  for j = 1:length(C.power_products[dd])
    # We only need the products for which the second factor (of degree dd) does
    # not involve any factor of index > i. Otherwise we will get duplicates.
    if C.last_factor[C.power_products[dd][j]] > i
      continue
    end

    f = C.power_products[dd][j]*C.base[i]
    push!(C.power_products[d], f)
    if remember_exponents
      C.exponent_vectors[f] = copy(C.exponent_vectors[C.power_products[dd][j]])
      C.exponent_vectors[f][i] += 1
    end
    C.last_factor[f] = i
  end
  return nothing
end

# Computes all power products of elements of C.base of total degree exactly d.
# This function is meant to be used in the case where one needs the result for
# different degrees d as results for lower degrees are reused (and cached in C).
# If one is only interested in all monomials of a certain degree, one should use
# all_monomials. Note that however also for a single degree d a naive evaluation
# of the result of all_monomials at the elements of C.base is in general slower.
# If remember_exponents is true, the exponent vectors of the computed power
# products will be stored in the dictionary C.exponent_vectors, that is, if
# C.base = [ f_1, ..., f_k ] and we compute f = f_1^e_1 \cdots f_k^e_k, then
# C.exponent_vectors[f] = [ e_1, ..., e_k ].
function all_power_products_of_degree!(C::PowerProductCache{S, T}, d::Int, remember_exponents::Bool = true) where {S, T <: MPolyElem}
  @assert d >= 0

  # Degree d has already been computed?
  if haskey(C.power_products, d) && (C.exponent_vectors_known[d] || !remember_exponents)
    return C.power_products[d]
  end

  C.power_products[d] = Vector{T}()

  if d == 0
    push!(C.power_products[d], one(C.ring))
    remember_exponents ? C.exponent_vectors[one(C.ring)] = zeros(Int, length(C.base)) : nothing
    C.last_factor[one(C.ring)] = 0
    return C.power_products[d]
  end

  for i = 1:length(C.base)
    extend!(C, d, i, remember_exponents)
  end
  C.exponent_vectors_known[d] = remember_exponents
  return C.power_products[d]
end

# Returns a generating set for the degree d component of a module
# \oplus_k K[f_1,\dots, f_n]*g_k where the f_i are given by the entries of C.base
# and module_gens is the list of the g_k.
# If remember_exponents is true, the second return value is a dictionary D
# containing the exponent vectors of the computed elements, that is, if
# f = f_1^e_1 \cdots f_n^e_k \cdot g_1^e_{k + 1} \cdots g_m^{e_{k + m}} was
# computed, then D[f] = [ e_1, ..., e_{k + m} ] where e_l in { 0, 1 } for l > k.
# (If remember_exponents is false, the dictionary is empty.)
function generators_for_given_degree!(C::PowerProductCache{S, T}, module_gens::Vector{T}, d::Int, remember_exponents) where {S, T <: MPolyElem}
  @assert !isempty(module_gens)
  @assert all(!iszero, module_gens)

  R = C.ring
  generators = elem_type(R)[]
  exps = Dict{T, Vector{Int}}()
  coords = zeros(Int, length(C.base) + length(module_gens))

  for i = 1:length(module_gens)
    di = total_degree(module_gens[i])
    if di > d
      continue
    end
    if di == d
      push!(generators, module_gens[i])
      if remember_exponents
        for j = 1:length(C.base)
          coords[j] = 0
        end
        coords[length(C.base) + i] = 1
        exps[module_gens[i]] = copy(coords)
        coords[length(C.base) + i] = 0
      end
      continue
    end
    dd = d - di
    all_power_products_of_degree!(C, dd, remember_exponents)
    for j = 1:length(C.power_products[dd])
      f = C.power_products[dd][j]*module_gens[i]
      push!(generators, f)

      if remember_exponents
        coords[1:length(C.base)] = copy(C.exponent_vectors[C.power_products[dd][j]])
        coords[length(C.base) + i] = 1
        exps[f] = copy(coords)
        coords[length(C.base) + i] = 0
      end
    end
  end
  return generators, exps
end

# Adds f to the span of the polynomials in BasisOfPolynomials.
# Returns true iff the rank of the span increased, i.e. if f were not an element
# of the span.
function add_to_basis!(B::BasisOfPolynomials{T}, f::T) where {T <: MPolyElem}
  @assert B.R === parent(f)

  zero!.(B.v)
  c = ncols(B.M)
  for (a, m) in zip(coefficients(f), monomials(f))
    if !haskey(B.monomial_to_column, m)
      c += 1
      B.monomial_to_column[m] = c
      push!(B.v, deepcopy(a))
    else
      col = B.monomial_to_column[m]
      B.v[col] = deepcopy(a)
    end
  end
  new_cols = c - ncols(B.M)
  if !iszero(new_cols)
    append!(B.pivot_rows, zeros(Int, new_cols))
    B.M = hcat(B.M, zero_matrix(base_ring(B.M), nrows(B.M), new_cols))
  end
  if B.r == nrows(B.M)
    B.M = vcat(B.M, zero_matrix(base_ring(B.M), 1, ncols(B.M)))
  end
  b = Hecke._add_row_to_rref!(B.M, B.v, B.pivot_rows, B.r + 1)
  b && (B.r += 1)
  return b
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

  C = PowerProductCache(R, p_invars)
  for d = 1:sum( total_degree(f) - 1 for f in p_invars )
    Md = generators_for_given_degree!(C, s_invars, d, false)[1]

    B = BasisOfPolynomials(R, Md)

    for f in iterate_basis(RG, d)
      f = f.f
      if add_to_basis!(B, f)
        push!(s_invars, f)
      end
    end
  end
  return [ Rgraded(f) for f in s_invars ]
end

function irreducible_secondary_invariants_modular(RG::InvRing)
  Rgraded = polynomial_ring(RG)
  R = Rgraded.R
  K = coefficient_ring(R)

  s_invars = elem_type(R)[ f.f for f in secondary_invariants_modular(RG) ]
  maxd = maximum(total_degree(f) for f in s_invars)
  s_invars_sorted = Vector{Vector{elem_type(R)}}(undef, maxd)
  for i = 1:maxd
    s_invars_sorted[i] = elem_type(R)[]
  end
  for f in s_invars
    d = total_degree(f)
    if d == 0
      continue
    end
    push!(s_invars_sorted[d], f)
  end

  # secondary invariants of degree 1 are irreducible
  is_invars = s_invars_sorted[1]
  C = PowerProductCache(R, append!(elem_type(R)[ f.f for f in primary_invariants(RG) ], is_invars))
  for d = 2:maxd
    basis_d = all_power_products_of_degree!(C, d, false)
    B = BasisOfPolynomials(R, basis_d)
    for f in s_invars_sorted[d]
      if add_to_basis!(B, f)
        push!(is_invars, f)
        add_base_element!(C, f, false)
      end
    end
  end
  return is_invars
end

################################################################################
#
#  Non modular case
#
################################################################################

# DK15, Algorithm 3.7.2
function secondary_invariants_nonmodular(RG::InvRing)
  @assert !ismodular(RG)
  p_invars = primary_invariants(RG)
  I = ideal_of_primary_invariants(RG)
  LI = leading_ideal(I)

  fg = molien_series(RG)
  f = numerator(fg)
  g = denominator(fg)
  for p in p_invars
    d = total_degree(p.f)
    # multiply f by 1 - t^d
    f -= shift_left(f, d)
  end
  h = div(f, g)
  @assert g*h == f

  Rgraded = polynomial_ring(RG)
  R = Rgraded.R
  K = coefficient_ring(R)
  s_invars = elem_type(R)[ one(R) ] # secondary invariants
  is_invars = elem_type(R)[] # irreducible secondary invariants

  # sum(coefficients(h)) is the number of secondary invariants we need in total.
  B = BasisOfPolynomials(R, [ one(R) ], Int(numerator(sum(coefficients(h)))))

  # We try to construct as many invariants as possible as power products from
  # already computed ones using all_power_products_of_degree! .
  C = PowerProductCache(R, is_invars)
  for d = 1:degree(h)
    k = coeff(h, d) # number of invariants we need in degree d
    if iszero(k)
      continue
    end
    invars_found = 0

    invars = all_power_products_of_degree!(C, d, false)
    for f in invars
      nf = normal_form(f, I).f
      if add_to_basis!(B, nf)
        push!(s_invars, inv(leading_coefficient(f))*f)
        invars_found += 1
        invars_found == k && break
      end
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
      nf = normal_form(f, I).f
      if add_to_basis!(B, nf)
        f = inv(leading_coefficient(f))*f
        push!(s_invars, f)
        push!(is_invars, f)
        add_base_element!(C, f, false)
        invars_found += 1
        invars_found == k && break
      end
    end
  end
  return [ Rgraded(f) for f in s_invars], [ Rgraded(f) for f in is_invars ]
end

################################################################################
#
#  Singular
#
################################################################################

function secondary_invariants_via_singular(IR::InvRing)
  rey, mol = reynolds_molien_via_singular(IR)
  primary_invariants_via_successive_algo(IR)
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
  return s, is
end

@doc Markdown.doc"""
    secondary_invariants(IR::InvRing)

Return a system of secondary invariants for `IR` with respect to the currently
cached system of primary invariants for `IR` (if no system of primary invariants
for `IR` is cached, compute and cache such a system first).

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
    s, is = secondary_invariants_via_singular(IR)
    IR.secondary = s
    IR.irreducible_secondary = is
  end
  return copy(IR.secondary)
end

@doc Markdown.doc"""
    irreducible_secondary_invariants(IR::InvRing)

Return a system of irreducible secondary invariants for `IR` with respect to the
currently cached system of primary invariants for `IR` (if no system of primary
invariants for `IR` is cached, compute and cache such a system first).

If a system of irreducible secondary invariants is already cached, return the
cached system.
Otherwise, compute and cache such a system first.

NOTE: The irreducible secondary invariants are sorted by increasing degree.

# Examples
```jldoctest
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
    s, is = secondary_invariants_via_singular(IR)
    IR.secondary = s
    IR.irreducible_secondary = is
  end
  return copy(IR.irreducible_secondary)
end
