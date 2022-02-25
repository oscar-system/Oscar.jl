################################################################################
#
#  Relations
#
################################################################################

# Relations between primary and irreducible secondary invariants, see
# [KS99, Section 17.5.5] or [DK15, Section 3.8.3].
function relations(RG::InvRing)
  Rgraded = polynomial_ring(RG)
  R = Rgraded.R
  K = coefficient_ring(R)

  p_invars = [ f.f for f in primary_invariants(RG) ]
  s_invars = [ f.f for f in secondary_invariants(RG) ]
  is_invars = [ f.f for f in irreducible_secondary_invariants(RG) ]
  s_invars_cache = RG.secondary

  S, t = grade(PolynomialRing(K, "t" => 1:(length(p_invars) + length(is_invars)))[1], append!([ total_degree(f) for f in p_invars ], [ total_degree(f) for f in is_invars ]))

  RtoS = Dict{elem_type(R), elem_type(S)}()
  for i = 1:length(p_invars)
    RtoS[p_invars[i]] = t[i]
  end
  for i = 1:length(s_invars)
    exps = append!(zeros(Int, length(p_invars)), s_invars_cache.sec_in_irred[i])
    g = set_exponent_vector!(one(S), 1, exps)
    RtoS[s_invars[i]] = g
  end

  # Assumes that s_invars and is_invars are sorted by degree
  maxd = total_degree(s_invars[end]) + total_degree(is_invars[end])
  products_sorted = Vector{Vector{Tuple{elem_type(R), elem_type(S)}}}(undef, maxd)
  for d = 1:maxd
    products_sorted[d] = Vector{Tuple{elem_type(R), elem_type(S)}}()
  end
  for i = 1:length(s_invars)
    for j = 1:length(s_invars)
      if !s_invars_cache.is_irreducible[j]
        continue
      end
      if s_invars_cache.is_irreducible[i] && i > j
        continue
      end
      m = RtoS[s_invars[i]]*RtoS[s_invars[j]]
      if m in values(RtoS)
        continue
      end
      f = s_invars[i]*s_invars[j]
      push!(products_sorted[total_degree(f)], (f, m))
    end
  end

  # TODO: In the modular case we need module syzygies!

  rels = elem_type(S)[]
  C = PowerProductCache(R, p_invars)
  p_and_s_invars = append!(copy(p_invars), s_invars)
  for d = 1:maxd
    if isempty(products_sorted[d])
      continue
    end
    gensd, expsd = generators_for_given_degree!(C, s_invars, d, true)

    monomial_to_column = Dict{elem_type(R), Int}()
    c = 0
    for f in gensd
      for m in monomials(f)
        if !haskey(monomial_to_column, m)
          c += 1
          monomial_to_column[m] = c
        end
      end
    end

    M = zero_matrix(K, length(gensd), c)
    for i = 1:length(gensd)
      f = gensd[i]
      for (a, m) in zip(coefficients(f), monomials(f))
        col = monomial_to_column[m]
        M[i, col] = deepcopy(a)
      end
    end

    N = zero_matrix(K, length(products_sorted[d]), c)
    for i = 1:length(products_sorted[d])
      f = products_sorted[d][i][1]
      for (a, m) in zip(coefficients(f), monomials(f))
        @assert haskey(monomial_to_column, m) "Monomial not found; this should not happen!"
        col = monomial_to_column[m]
        N[i, col] = deepcopy(a)
      end
    end

    fl, x = can_solve_with_solution(M, N, side = :left)
    @assert fl

    for i = 1:nrows(x)
      s = -products_sorted[d][i][2]
      for j = 1:ncols(x)
        m = S(x[i, j])
        for k = 1:length(p_and_s_invars)
          m *= RtoS[p_and_s_invars[k]]^expsd[gensd[j]][k]
        end
        s += m
      end
      push!(rels, s)
    end
  end
  return S, hom(S, Rgraded, append!(primary_invariants(RG), irreducible_secondary_invariants(RG))), rels
end
