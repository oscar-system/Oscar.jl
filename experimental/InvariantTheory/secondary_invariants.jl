################################################################################
#
#  Helpers
#
################################################################################

# Computes all power products of elements of polys of total degree exactly d.
# This function is meant to be used in the case where one needs the result for
# different degrees d as results for lower degrees are reused.
# If one is only interested in all monomials of a certain degree, one should use
# all_monomials. Note that however also for a single degree d a naive evaluation
# of the result of all_monomials at the elements of polys is in general slower.
#
# database is a dictionary containing the already computed degree components.
# The value corresponding to a key d consists of:
# * a vector of polynomials: all power products of the polynomials of degree d.
# * a vector of vectors of integers: the representation of the products in the vector
#   of polynomials as an exponent vector.
# * a vector of integers: The first non-zero index of the exponent vector of the
#   corresponding element.
# If the key d in database exists, it is assumed that the value is correctly
# computed!
function all_power_products_of_degree!(polys::Vector{T}, d::Int, database::Dict{Int, Tuple{Vector{T}, Vector{Vector{Int}}, Vector{Int}}}) where {T <: MPolyElem}
  @assert all(!iszero, polys)
  # Degree d has already been computed?
  if haskey(database, d)
    return database[d][1]
  end

  coords = zeros(Int, length(polys))
  database[d] = (T[], Vector{Int}[], Int[])
  for i = 1:length(polys)
    di = total_degree(polys[i])

    # Larger degree, ignore
    if di > d
      continue
    end

    # Exactly matching degree
    if di == d
      push!(database[d][1], polys[i])
      coords[i] = 1;
      push!(database[d][2], copy(coords))
      coords[i] = 0;
      push!(database[d][3], i)
      continue
    end

    # Build all products with elements of degree d - total_degree(polys[i])
    dd = d - di
    all_power_products_of_degree!(polys, dd, database)
    for j = 1:length(database[dd][1])
      # We only need the products for which the second factor (of degree dd) does
      # not involve any factor of index < i. Otherwise we will get duplicates.
      if database[dd][3][j] < i
        continue
      end

      push!(database[d][1], database[dd][1][j]*polys[i])
      push!(database[d][2], copy(database[dd][2][j]))
      database[d][2][end][i] += 1
      push!(database[d][3], i)
    end
  end
  return database[d][1]
end

# Returns a generating set for the degree d component of a module
# \oplus_k K[f_1,\dots, f_n]*g_k.
# * algebra_gens is the list of the f_i
# * module_gens is the list of the g_k
# * database is a helper dictionary (possibly empty); see documentation for
#   all_power_products_of_degree!
# Returns
# * the generating set as a vector of polynomials
# * a vector of vectors of integers giving the representation of the polynomials
#   as products of the algebra_gens and module_gens
# The returned generating set should in general be a vector space basis, however,
# if e.g. the algebra_gens and module_gens are not disjoint, it might contain
# duplicates.
function generators_for_given_degree!(algebra_gens::Vector{T}, module_gens::Vector{T}, d::Int, database::Dict{Int, Tuple{Vector{T}, Vector{Vector{Int}}, Vector{Int}}}) where {T <: MPolyElem}
  @assert !isempty(algebra_gens) && !isempty(module_gens)
  @assert all(!iszero, algebra_gens) && all(!iszero, module_gens)

  R = parent(algebra_gens[1])
  generators = elem_type(R)[]
  exps = Vector{Int}[]
  coords = zeros(Int, length(algebra_gens) + length(module_gens))

  for i = 1:length(module_gens)
    di = total_degree(module_gens[i])
    if di > d
      continue
    end
    if di == d
      for j = 1:length(algebra_gens)
        coords[j] = 0
      end
      push!(generators, module_gens[i])
      coords[length(algebra_gens) + i] = 1
      push!(exps, copy(coords))
      coords[length(algebra_gens) + i] = 0
      continue
    end
    dd = d - di
    _ = all_power_products_of_degree!(algebra_gens, dd, database)
    for j = 1:length(database[dd][1])
      push!(generators, database[dd][1][j]*module_gens[i])

      coords[1:length(algebra_gens)] = database[dd][2][j]
      coords[length(algebra_gens) + i] = 1;
      push!(exps, copy(coords))
      coords[length(algebra_gens) + i] = 0;
    end
  end
  return generators, exps
end

# Adds f to the span of the polynomials in BasisOfPolynomials.
# Returns true iff the rank of the span increased, i.e. if f were not an element
# of the span.
function inspan!(B::BasisOfPolynomials{T}, f::T) where {T <: MPolyElem}
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

  database = Dict{Int, Tuple{Vector{elem_type(R)}, Vector{Vector{Int}}, Vector{Int}}}()
  for d = 1:sum( total_degree(f) - 1 for f in p_invars )
    Md = generators_for_given_degree!(p_invars, s_invars, d, database)[1]

    B = BasisOfPolynomials(R, Md)

    for f in iterate_basis(RG, d)
      f = f.f
      if inspan!(B, f)
        push!(s_invars, f)
      end
    end
  end
  return [ Rgraded(f) for f in s_invars ]
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

  R = polynomial_ring(RG)
  K = coefficient_ring(R)
  s_invars = elem_type(R)[ one(R) ] # secondary invariants

  # sum(coefficients(h)) is the number of secondary invariants we need in total.
  B = BasisOfPolynomials(R, [ one(R) ], Int(numerator(sum(coefficients(h)))))
  for d = 1:degree(h)
    k = coeff(h, d) # number of invariants we need in degree d
    if iszero(k)
      continue
    end
    invars_found = 0

    mons = all_monomials(polynomial_ring(RG), d)
    for m in mons
      f = reynolds_operator(RG, m)
      if iszero(f)
        continue
      end
      nf = normal_form(f, I)

      if inspan!(B, nf)
        push!(s_invars, inv(leading_coefficient(f))*f)
        invars_found += 1
        invars_found == k && break
      end
    end
  end
  return s_invars
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
