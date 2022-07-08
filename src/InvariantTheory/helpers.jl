################################################################################
#
#  PowerProductCache
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
  C.exponent_vectors_known[d] = remember_exponents

  if d == 0
    push!(C.power_products[d], one(C.ring))
    remember_exponents ? C.exponent_vectors[one(C.ring)] = zeros(Int, length(C.base)) : nothing
    C.last_factor[one(C.ring)] = 0
    return C.power_products[d]
  end

  for i = 1:length(C.base)
    extend!(C, d, i, remember_exponents)
  end

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

################################################################################
#
#  BasisOfPolynomials
#
################################################################################

# Adds f to the span of the polynomials in BasisOfPolynomials.
# Returns true iff the rank of the span increased, i.e. if f was not an element
# of the span.
function add_to_basis!(B::BasisOfPolynomials{T}, f::T) where {T <: MPolyElem}
  @assert B.R === parent(f)

  c = ncols(B.M)
  srow = sparse_row(coefficient_ring(B.R))
  for (a, m) in zip(coefficients(f), monomials(f))
    if !haskey(B.monomial_to_column, m)
      c += 1
      B.monomial_to_column[m] = c
      push!(srow.pos, c)
      push!(srow.values, deepcopy(a))
    else
      col = B.monomial_to_column[m]
      k = searchsortedfirst(srow.pos, col)
      insert!(srow.pos, k, col)
      insert!(srow.values, k, deepcopy(a))
    end
  end
  return Hecke._add_row_to_rref!(B.M, srow)
end

function enumerate_monomials(polys::Vector{PolyElemT}) where {PolyElemT <: MPolyElem}
  enum_mons = Dict{PolyElemT, Int}()
  c = 0
  for f in polys
    for m in monomials(f)
      if !haskey(enum_mons, m)
        c += 1
        enum_mons[m] = c
      end
    end
  end
  return enum_mons
end

function polys_to_smat(polys::Vector{PolyElemT}, monomial_to_column::Dict{PolyElemT, Int}; copy = true) where {PolyElemT <: MPolyElem}
  @assert !isempty(polys)
  K = coefficient_ring(parent(polys[1]))
  M = sparse_matrix(K)
  for i = 1:length(polys)
    srow = sparse_row(K)
    for (a, m) in zip(coefficients(polys[i]), monomials(polys[i]))
      col = get(monomial_to_column, m) do
        error("Monomial not found in the given dictionary")
      end
      k = searchsortedfirst(srow.pos, col)
      insert!(srow.pos, k, col)
      if copy
        insert!(srow.values, k, deepcopy(a))
      else
        insert!(srow.values, k, a)
      end
    end
    push!(M, srow)
  end
  return M
end
