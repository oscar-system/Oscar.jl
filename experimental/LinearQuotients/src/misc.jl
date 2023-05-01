# Try to 'update' the base_ring of G
function map_entries(K::Field, G::MatrixGroup)
  g = dense_matrix_type(K)[]
  for h in gens(G)
    push!(g, map_entries(K, h.elm))
  end
  return matrix_group(g)
end

function is_reflection(g::MatrixGroupElem)
  return rank(g.elm - one(parent(g)).elm) == 1
end

function subgroup_of_reflections(G::MatrixGroup)
  g = elem_type(G)[]
  for c in conjugacy_classes(G)
    if is_reflection(representative(c))
      append!(g, collect(c))
    end
  end
  return matrix_group(degree(G), base_ring(G), g)
end

# Check if G contains reflections
function is_small_group(G::MatrixGroup)
  return !any(c -> is_reflection(representative(c)), conjugacy_classes(G))
end

# Somehow one needs to look up powers of a given root of unity quite often...
# Assumes that zeta is a primitive l-th root of unity
function _powers_of_root_of_unity(zeta::FieldElem, l::Int)
  K = parent(zeta)
  powers_of_zeta = Dict{elem_type(K), Int}()
  t = one(K)
  for i in 0:l - 1
    powers_of_zeta[t] = i
    t *= zeta
  end
  @assert is_one(t) "Given element is not a $l-th root of unity"

  return powers_of_zeta
end

function is_subgroup_of_sl(G::MatrixGroup)
  return all(g -> is_one(det(g.elm)), gens(G))
end

# Mostly stolen from mpoly-graded.jl .
# If positive == false, this homogenizes 'negatively', that is, the resulting
# polynomial will have the minimal degree of the input.
# Assumes that the generators start_pos:start_pos + ngens(S.D) - 1 of S have
# degree [1, 0, ..., 0 ], [ 0, 1, 0, ..., 0 ], ..., [ 0, ..., 0, 1 ] if
# positive == true and with -1 if positive == false.
function homogenize(f::MPolyRingElem, S::MPolyDecRing, start_pos::Int, positive::Bool)
  @assert nvars(parent(f)) + ngens(S.D) == nvars(S)
  @assert is_free(S.D)

  exps = Vector{Int}[]
  for e in AbstractAlgebra.exponent_vectors(f)
    for j in 1:ngens(S.D)
      insert!(e, start_pos + j - 1, 0)
    end
    push!(exps, e)
  end

  degs_mons = Vector{Int}[]
  degs_vars = [ [ Int(S.d[i][j]) for j in 1:ngens(S.D) ] for i in 1:ngens(S) ]

  for e in exps
    push!(degs_mons, sum(e.*degs_vars))
  end

  if positive
    tdeg = maximum(hcat(degs_mons...), dims = 2)[:, 1]
  else
    ldeg = minimum(hcat(degs_mons...), dims = 2)[:, 1]
  end

  F = MPolyBuildCtx(S)
  for (i, (c, e)) in enumerate(zip(AbstractAlgebra.coefficients(f), exps))
    for j in 1:ngens(S.D)
      if positive
        e[start_pos + j - 1] = tdeg[j] - degs_mons[i][j]
      else
        e[start_pos + j - 1] = degs_mons[i][j] - ldeg[j]
      end
    end
    push_term!(F, c, e)
  end

  return finish(F)
end

# Homogenize the ideal with respect to the last variable of S.
# This might be either of degree 1 or -1.
# Uses a variant of Bayer's method, see [Sch23, Proposition 6.4.3].
function homogenize_at_last_variable(I::MPolyIdeal, S::MPolyDecRing)
  @assert nvars(S) == nvars(base_ring(I)) + 1
  @assert ngens(S.D) == 1

  # Permute the variables. The standard basis appears to be faster if the bigger
  # weights are in front.
  p = sortperm([ S.d[i].coeff[1] for i in 1:ngens(S) - 1 ], rev = true)

  push!(p, ngens(S))

  w_perm = ZZRingElem[ S.d[p[i]].coeff[1] for i in 1:ngens(S) ]
  @assert !iszero(w_perm[1]) "All given weights are zero"

  R = base_ring(I)
  Rp = polynomial_ring(coefficient_ring(R), "t" => 1:ngens(R))[1]
  RtoRp = hom(R, Rp, [ gen(Rp, findfirst(isequal(i), p)) for i in 1:ngens(Rp) ])

  Sp, _ = graded_polynomial_ring(coefficient_ring(S), ["t$i" for i in 1:ngens(S)], w_perm)
  SptoS = hom(Sp, S, [ gens(S)[p[i]] for i in 1:ngens(S) ])

  # Homogenize the generators
  positive = (Sp.d[ngens(Sp)].coeff[1] > 0)
  J = ideal(Sp, [ homogenize(RtoRp(f), Sp, ngens(Sp), positive) for f in gens(I) ])

  # Now saturate w.r.t. the last variable using a variant of Bayer's method

  # Build the matrix
  #      (       w_perm       )
  #      (  0  ... ...  0  -1 )
  # M := (  1               0 )
  #      (   . .     .      . )
  #      (  0       1   0   0 )
  # for the ordering
  M = zero_matrix(ZZ, nvars(S), nvars(S))
  for i in 1:nvars(Sp)
    M[1, i] = w_perm[i]
  end
  M[2, ncols(M)] = -1

  non_zero_weight = false
  for i in 2:nvars(Sp) - 1
    M[i + 1, i] = 1
  end

  o = matrix_ordering(gens(Sp), M)
  gb = standard_basis(J, ordering = o)

  t = gen(Sp, nvars(Sp))
  res = elem_type(Sp)[]
  for f in gb
    while iszero(mod(f, t))
      f = div(f, t)
    end
    push!(res, f)
  end
  return SptoS(ideal(Sp, res))
end

function ideal(S::MPolyDecRing, I::MPolyIdeal)
  @assert base_ring(I) === forget_grading(S)
  return ideal(S, [ S(f) for f in gens(I) ])
end
