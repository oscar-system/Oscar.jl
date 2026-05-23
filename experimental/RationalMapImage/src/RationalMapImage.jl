################################################################################
# RationalMapImage
#
# Translation of the Macaulay2 package TotalImage by Corey Harris,
# Mateusz Michalek, and Emre Sertöz, which accompanies the paper
# "Computing images of polynomial maps" (2017).
#
# The algorithm computes the image of a rational/polynomial map as a
# constructible set, represented as a tree of projective varieties.
################################################################################

module RationalMapImage

import ..Oscar:
  MPolyRing,
  MPolyRingElem,
  MPolyIdeal,
  base_ring,
  coefficient_ring,
  codomain,
  elem_type,
  ngens,
  gens,
  gen,
  symbols,
  polynomial_ring,
  hom,
  ideal,
  preimage,
  eliminate,
  saturation,
  minimal_primes,
  radical,
  dim,
  is_subset,
  is_zero,
  is_homogeneous,
  total_degree,
  coefficients,
  monomials

export partial_image,
       total_image,
       is_closed_image,
       print_image_tree

################################################################################
# Phase 1 – Core algebraic helpers
################################################################################

# Create the target polynomial ring k[b_0,...,b_m] from a list of polynomials.
function _generate_target_ring(L::Vector{<:MPolyRingElem})
  m = length(L) - 1
  R = parent(first(L))
  k = coefficient_ring(R)
  T, _ = polynomial_ring(k, [Symbol(:b_, i) for i in 0:m])
  return T
end

# Scheme-theoretic Zariski closure of the image of V(X) under the map defined
# by L.  Uses Oscar's built-in preimage via Singular.
function _scheme_theoretic_image(L::Vector{<:MPolyRingElem}, X::MPolyIdeal, T::MPolyRing)
  R = parent(first(L))
  @assert base_ring(X) === R "Domain ring mismatch"
  f = hom(T, R, L)
  return preimage(f, X)
end

# Rees-algebra blowup of the ideal I via elimination.
# Returns (Bl, W, x_vars, b_vars) where:
#   W  = K[domain_vars, target_vars] (combined ambient ring)
#   Bl = ideal of the blowup in W
#   x_vars = domain generators in W
#   b_vars = target/blowup generators in W
# The ordering is: x_1,...,x_n, b_0,...,b_{r-1} in W (no t variable).
function _blowup_of_ideal(I::MPolyIdeal, T::MPolyRing)
  R = base_ring(I)
  n = ngens(R)
  r = ngens(I)
  k = coefficient_ring(R)

  x_names = [Symbol(:_x, i) for i in 1:n]
  b_names = [Symbol(:_b, i) for i in 0:(r - 1)]

  # Extended ring with t first (for elimination ordering)
  Wt, Wt_vars = polynomial_ring(k, [:_t; x_names; b_names])
  t_var     = Wt_vars[1]
  x_vars_Wt = Wt_vars[2:(n + 1)]
  b_vars_Wt = Wt_vars[(n + 2):end]

  R_to_Wt = hom(R, Wt, x_vars_Wt)
  I_gens_Wt = [R_to_Wt(g) for g in gens(I)]

  J = ideal(Wt, [b_vars_Wt[j] - t_var * I_gens_Wt[j] for j in 1:r])
  Bl_in_Wt = eliminate(J, [t_var])

  # Drop t: create W = K[x_1,...,x_n, b_0,...,b_{r-1}]
  W, W_vars = polynomial_ring(k, [x_names; b_names])
  x_vars_W = W_vars[1:n]
  b_vars_W = W_vars[(n + 1):end]

  Wt_to_W = hom(Wt, W, [zero(W); W_vars])
  Bl = ideal(W, [Wt_to_W(g) for g in gens(Bl_in_Wt)])

  return (Bl, W, x_vars_W, b_vars_W)
end

# Random hyperplane section: returns a ring map R → R' that eliminates
# the first n_elim variables of R by substituting them with random linear
# combinations of the remaining variables.
function _hyperplane_section(X_ideal::MPolyIdeal, n_elim::Int)
  R = base_ring(X_ideal)
  N = ngens(R)
  @assert 0 < n_elim < N "n_elim must be between 1 and ngens(R)-1"

  lastvars_names = symbols(R)[(n_elim + 1):end]
  R_prime, R_prime_vars = polynomial_ring(coefficient_ring(R), collect(lastvars_names))

  # Build substitution: first n_elim variables → random linear combination
  # of the last (N - n_elim) variables
  sublist = [sum(rand(-200:200) * v for v in R_prime_vars) for _ in 1:n_elim]
  img = [sublist; R_prime_vars]

  return hom(R, R_prime, img)
end

# Monomial restriction: set variable i to 0, return ring map R → R'
function _monomial_restriction(R::MPolyRing, i::Int)
  var_syms = collect(symbols(R))
  keep_syms = [var_syms[j] for j in 1:length(var_syms) if j != i]
  R_prime, R_prime_vars = polynomial_ring(coefficient_ring(R), keep_syms)

  img = MPolyRingElem[]
  k = 1
  for j in 1:ngens(R)
    if j == i
      push!(img, zero(R_prime))
    else
      push!(img, R_prime_vars[k])
      k += 1
    end
  end
  return hom(R, R_prime, img)
end

# Binomial restriction: set variable i to a * variable j, return ring map R → R'
# The variable at position i is substituted; position j is kept.
# Indices are 1-based.
function _binomial_restriction(R::MPolyRing, i::Int, j::Int, a::Int)
  @assert i != j
  var_syms = collect(symbols(R))
  # New ring: keep all vars except position i
  keep_syms = [var_syms[k] for k in 1:length(var_syms) if k != i]
  R_prime, R_prime_vars = polynomial_ring(coefficient_ring(R), keep_syms)

  # Build map: x_k → R_prime var (with shift around i), except x_i → a * x_j_in_R_prime
  # j' = new position of j in R_prime (after removing i)
  j_prime = j < i ? j : j - 1

  img = MPolyRingElem[]
  k = 1
  for l in 1:ngens(R)
    if l == i
      push!(img, coefficient_ring(R)(a) * R_prime_vars[j_prime])
    else
      push!(img, R_prime_vars[k])
      k += 1
    end
  end
  return hom(R, R_prime, img)
end

# Apply a ring map (restriction) to an ideal: image of ideal generators under f.
function _apply_map_to_ideal(f, I::MPolyIdeal)
  S = codomain(f)
  return ideal(S, [f(g) for g in gens(I)])
end

# Pullback of exc_image (ideal in T) under the rational map defined by L.
# Scheme-theoretically: substitute b_i → L_i in generators of exc_image,
# add X (domain ideal), then saturate with the base locus ideal(L).
function _pullback_ideal(L::Vector{<:MPolyRingElem}, X::MPolyIdeal, exc_image::MPolyIdeal)
  R = parent(first(L))
  T = base_ring(exc_image)
  T_to_R = hom(T, R, L)
  pb = ideal(R, [T_to_R(g) for g in gens(exc_image)])
  return saturation(pb + X, ideal(R, L))
end

# Minimal prime decomposition of the pullback.
function _components_of_pullback(L::Vector{<:MPolyRingElem}, X::MPolyIdeal, exc_image::MPolyIdeal)
  return minimal_primes(_pullback_ideal(L, X, exc_image))
end

################################################################################
# Phase 2 – Main algorithm
################################################################################

# Apply a hyperplane restriction map r: R → R' to a polynomial map L.
function _restrict_map(r, L::Vector{<:MPolyRingElem})
  return [r(f) for f in L]
end

# Given a restriction ring map r: R → R', apply r to the domain ideal.
function _restrict_ideal(r, I::MPolyIdeal)
  R_prime = codomain(r)
  return ideal(R_prime, [r(g) for g in gens(I)])
end

# Attempt a single hyperplane restriction of type :monomial (set x_i = 0).
# Returns (new_L, new_X, r) or nothing if the restriction is not "good".
# "Good" means: dim(new_X) == dim(X) - 1 and image(new_X) ⊆ image_X.
function _try_monomial(i::Int, L::Vector{<:MPolyRingElem}, X::MPolyIdeal,
                       image_X::MPolyIdeal, T::MPolyRing)
  R = base_ring(X)
  ngens(R) < 2 && return nothing
  r = _monomial_restriction(R, i)
  new_X = _restrict_ideal(r, X)
  new_L = _restrict_map(r, L)
  # dim check
  dim(new_X) != dim(X) - 1 && return nothing
  # Image check: restricted image ⊆ image_X
  new_img = _scheme_theoretic_image(new_L, new_X, T)
  is_subset(new_img, image_X) || return nothing
  return (new_L, new_X, r)
end

# Attempt a single binomial restriction x_i → a * x_j.
function _try_binomial(i::Int, j::Int, a::Int, L::Vector{<:MPolyRingElem},
                       X::MPolyIdeal, image_X::MPolyIdeal, T::MPolyRing)
  R = base_ring(X)
  r = _binomial_restriction(R, i, j, a)
  new_X = _restrict_ideal(r, X)
  new_L = _restrict_map(r, L)
  dim(new_X) != dim(X) - 1 && return nothing
  new_img = _scheme_theoretic_image(new_L, new_X, T)
  is_subset(new_img, image_X) || return nothing
  return (new_L, new_X, r)
end

# Reduce the domain until dim(X) == dim(image_X) by finding hyperplane sections
# that preserve the image.  Returns (new_L, new_X, running_restriction_map).
function _shrink_X(L::Vector{<:MPolyRingElem}, X::MPolyIdeal,
                   image_X::MPolyIdeal, T::MPolyRing; tries::Int = 50)
  cur_L = L
  cur_X = X
  R_cur = base_ring(cur_X)

  dim_img = dim(image_X)

  while dim(cur_X) > dim_img
    R_cur = base_ring(cur_X)
    n = ngens(R_cur)
    found = false

    # Try monomial restrictions x_i = 0
    for i in 1:n
      res = _try_monomial(i, cur_L, cur_X, image_X, T)
      if res !== nothing
        (cur_L, cur_X, _) = res
        found = true
        break
      end
    end
    found && continue

    # Try binomial restrictions x_i = a * x_j
    for a in [1, -1, 2, -2, -3, 3, 5, 7], i in 1:(n - 1), j in (i + 1):n
      res = _try_binomial(i, j, a, cur_L, cur_X, image_X, T)
      if res !== nothing
        (cur_L, cur_X, _) = res
        found = true
        break
      end
      found && break
    end
    found && continue

    # Random hyperplane sections
    found_random = false
    for _ in 1:tries
      r = _hyperplane_section(cur_X, dim(cur_X) - dim_img)
      new_X = _restrict_ideal(r, cur_X)
      new_L = _restrict_map(r, cur_L)
      if dim(new_X) == dim_img
        new_img = _scheme_theoretic_image(new_L, new_X, T)
        if is_subset(new_img, image_X)
          cur_L = new_L
          cur_X = new_X
          found_random = true
          break
        end
      end
    end
    found_random && continue

    error("_shrink_X: could not find a good hyperplane section after $(tries) tries")
  end

  return (cur_L, cur_X)
end

@doc raw"""
    partial_image(L::Vector{<:MPolyRingElem}; verbose::Bool=false, verify::Bool=true,
                  pd::Bool=true, tries::Int=50)
    partial_image(L::Vector{<:MPolyRingElem}, X::MPolyIdeal; verbose::Bool=false,
                  verify::Bool=true, pd::Bool=true, tries::Int=50)

Compute the Zariski closure of the image of the rational map defined by `L`
on the projective variety `V(X)` (default: whole projective space), together
with a list of exceptional divisors where the image may fail to be closed.

Input:
- `L`: list of homogeneous polynomials `[f_0,...,f_m]` in a ring `R = k[x_0,...,x_n]`
  defining a rational map `f: Proj(R/X) --> P^m`
- `X`: ideal of the source variety (default: zero ideal, i.e., whole space)

Output: `(image_ideal, exceptional_loci)` where
- `image_ideal` is the ideal of the Zariski closure of the image in `k[b_0,...,b_m]`
- `exceptional_loci` is a vector of prime ideals in `k[b_0,...,b_m]` representing
  components of the complement where the image may not be dense

Options:
- `verbose`: print progress information
- `verify`: use deterministic hyperplane sections to verify correctness (slower)
- `pd`: decompose exceptional image into prime components (default true)
- `tries`: maximum number of random hyperplane section attempts
"""
function partial_image(L::Vector{<:MPolyRingElem}, X::MPolyIdeal;
                       verbose::Bool = false, verify::Bool = true,
                       pd::Bool = true, tries::Int = 50)
  R = parent(first(L))
  @assert base_ring(X) === R "Domain ideal must be over the same ring as L"

  T = _generate_target_ring(L)

  return _partial_image_with_target(L, X, T; verbose = verbose, verify = verify,
                                    pd = pd, tries = tries)
end

function partial_image(L::Vector{<:MPolyRingElem};
                       verbose::Bool = false, verify::Bool = true,
                       pd::Bool = true, tries::Int = 50)
  R = parent(first(L))
  return partial_image(L, ideal(R, elem_type(R)[]); verbose = verbose,
                       verify = verify, pd = pd, tries = tries)
end

# Internal version that accepts a pre-created target ring T (for reuse in tree_builder).
function _partial_image_with_target(L::Vector{<:MPolyRingElem}, X::MPolyIdeal, T::MPolyRing;
                                     verbose::Bool = false, verify::Bool = true,
                                     pd::Bool = true, tries::Int = 50)
  R = parent(first(L))

  verbose && println("partial_image: computing scheme-theoretic image ...")
  verbose && println("  (domain has Krull dim ", dim(X), ")")
  image_X = _scheme_theoretic_image(L, X, T)
  verbose && println("  (image has Krull dim ", dim(image_X), ")")

  base_locus = ideal(R, L)

  # If the base locus does not intersect V(X), the map is defined everywhere
  # on V(X) and the image is Zariski-closed; return early.
  if dim(base_locus + X) < 1
    verbose && println("  base locus empty on V(X); returning closed image")
    return (image_X, MPolyIdeal[])
  end

  dim_img_X = dim(image_X)
  fiber_dim = dim(X) - dim_img_X

  # Shrink domain to a linear subspace of relative dimension = fiber_dim
  if verify
    verbose && println("  shrinking domain ...")
    (cur_L, cur_X) = _shrink_X(L, X, image_X, T; tries = tries)
  else
    verbose && println("  using random hyperplane section (verify=false) ...")
    r = _hyperplane_section(X, fiber_dim)
    cur_X = _restrict_ideal(r, X)
    cur_L = _restrict_map(r, L)
  end

  R_cur = base_ring(cur_X)
  cur_base_locus = ideal(R_cur, cur_L)

  # Early exit if base locus does not intersect the restricted domain
  if dim(cur_base_locus + cur_X) < 1
    verbose && println("  base locus empty after restriction; returning")
    return (image_X, MPolyIdeal[])
  end

  verbose && println("  computing blowup ...")
  (Bl, W, x_vars_W, b_vars_W) = _blowup_of_ideal(cur_base_locus, T)

  R_to_W = hom(R_cur, W, x_vars_W)

  bl_locus_W = ideal(W, [R_to_W(g) for g in gens(cur_base_locus)])

  # Proper transform of cur_X in the blowup
  if is_zero(cur_X)
    proper_X = ideal(W, elem_type(W)[])
  else
    X_in_W = ideal(W, [R_to_W(g) for g in gens(cur_X)])
    proper_X = saturation(X_in_W + Bl, bl_locus_W)
    verbose && println("  computed proper transform")
  end

  # Exceptional fiber: the closure of the graph intersected with the fiber
  # over the base locus, inside the proper transform of X.
  # Combining all three ideals first and then saturating by the domain
  # irrelevant ideal removes the spurious affine-origin component that would
  # otherwise make the projected image equal to the whole target space.
  E_total = Bl + bl_locus_W + proper_X

  verbose && println("  saturating exceptional divisor (projective) ...")
  E_total = saturation(E_total, ideal(W, x_vars_W))

  verbose && println("  eliminating domain variables ...")
  image_E_in_W = eliminate(E_total, collect(x_vars_W))

  # Map back to T: send x variables to 0, b variables to T generators
  n_x = length(x_vars_W)
  W_to_T = hom(W, T, vcat([zero(T) for _ in 1:n_x], collect(gens(T))))
  image_E = ideal(T, [W_to_T(g) for g in gens(image_E_in_W)])

  if dim(image_E) < 1
    verbose && println("  exceptional image is empty")
    return (image_X, MPolyIdeal[])
  end

  verbose && println("  computing components of exceptional image ...")
  if pd
    exc_loci = minimal_primes(image_E)
  else
    exc_loci = [image_E]
  end
  verbose && println("  found ", length(exc_loci), " exceptional component(s)")

  return (image_X, exc_loci)
end

################################################################################
# Phase 3 – Tree utilities (pure Julia)
################################################################################

# Returns (children, parents) where children[i] = list of child indices of node i,
# and parents[i] = list of parent indices (should be at most one for a tree).
# Nodes are 1-based.
function _children_and_parents(N::Vector, E::Vector{<:Tuple{Int,Int}})
  n = length(N)
  children = [Int[] for _ in 1:n]
  parents  = [Int[] for _ in 1:n]
  for (p, c) in E
    push!(children[p], c)
    push!(parents[c], p)
  end
  return (children, parents)
end

# All descendant indices of node i (including i itself), sorted.
function _all_descendants(children::Vector{Vector{Int}}, i::Int)
  desc = [i]
  k = 1
  while k <= length(desc)
    append!(desc, children[desc[k]])
    k += 1
  end
  return sort(desc)
end

# Extract the subtree rooted at node i, reindexed starting from 1.
function _subtree(N::Vector, E::Vector{<:Tuple{Int,Int}}, i::Int)
  (children, _) = _children_and_parents(N, E)
  desc = _all_descendants(children, i)
  desc_set = Set(desc)
  rank = Dict(d => findfirst(==(d), desc) for d in desc)
  new_N = [N[d] for d in desc]
  new_E = [(rank[p], rank[c]) for (p, c) in E if p in desc_set && c in desc_set]
  return (new_N, new_E)
end

# Returns (Bool, child_index) if node i has a child equal to N[i], else (false, -1).
function _detect_duplicate(N::Vector, children::Vector{Vector{Int}}, i::Int)
  for c in children[i]
    N[i] == N[c] && return (true, c)
  end
  return (false, -1)
end

# Remove edge e = (p, c): merge c into p by re-parenting c's children to p.
# Node c becomes isolated and will be removed by _reindex_tree.
function _remove_edge(N::Vector, E::Vector{<:Tuple{Int,Int}}, e::Tuple{Int,Int})
  (p, c) = e
  (children, _) = _children_and_parents(N, E)
  # Remove edge (p, c), add (p, k) for each child k of c
  new_E = [(a, b) for (a, b) in E if (a, b) != (p, c)]
  for k in children[c]
    push!(new_E, (p, k))
  end
  return (N, new_E)
end

# Levels of the tree (BFS).
function _tree_levels(N::Vector, E::Vector{<:Tuple{Int,Int}})
  (children, _) = _children_and_parents(N, E)
  levels = [[1]]
  while true
    next = Int[]
    for nd in last(levels)
      append!(next, children[nd])
    end
    isempty(next) && break
    push!(levels, next)
  end
  return levels
end

# Remove nodes at odd depth that are equal to their parent.
function _remove_duplicates(N::Vector, E::Vector{<:Tuple{Int,Int}})
  lvls = _tree_levels(N, E)
  (children, _) = _children_and_parents(N, E)

  # Odd-depth levels (1-based level index: level 2, 4, ... are "odd depth" in 0-based)
  odd_depth_nodes = Int[]
  for li in 2:2:length(lvls)
    append!(odd_depth_nodes, lvls[li])
  end

  bad_edges = Tuple{Int,Int}[]
  for nd in odd_depth_nodes
    (dup, child_idx) = _detect_duplicate(N, children, nd)
    dup && push!(bad_edges, (nd, child_idx))
  end

  new_E = E
  for e in bad_edges
    (_, new_E) = _remove_edge(N, new_E, e)
  end

  return _subtree(N, new_E, 1)
end

# Get indices of leaf nodes (no children, has a parent).
function _get_leaves(N::Vector, E::Vector{<:Tuple{Int,Int}})
  (children, parents) = _children_and_parents(N, E)
  return [i for i in 1:length(N) if isempty(children[i]) && !isempty(parents[i])]
end

# Remove leaves that are dimension-less (dim < 1).
function _prune_leaves(N::AbstractVector{<:MPolyIdeal}, E::Vector{<:Tuple{Int,Int}})
  leaves = _get_leaves(N, E)
  empty_leaves = Set(l for l in leaves if dim(N[l]) < 1)
  new_E = [(p, c) for (p, c) in E if c ∉ empty_leaves]
  return (N, new_E)
end

# Clean tree: iteratively remove leaves whose parent has the same dimension,
# and also remove their parent if the parent becomes a leaf.
function _clean_tree(N::AbstractVector{<:MPolyIdeal}, E::Vector{<:Tuple{Int,Int}})
  (N2, E2) = _prune_leaves(N, E)
  (children2, parents2) = _children_and_parents(N2, E2)
  leaves2 = _get_leaves(N2, E2)

  bad_leaves = [l for l in leaves2 if !isempty(parents2[l]) &&
                dim(N2[parents2[l][1]]) == dim(N2[l])]
  bad_parents = [parents2[l][1] for l in bad_leaves]
  bad_edges = [(parents2[l][1], l) for l in bad_leaves]
  # Also remove edges from grandparent to bad parent
  for p in bad_parents
    if !isempty(parents2[p])
      push!(bad_edges, (parents2[p][1], p))
    end
  end

  new_E2 = [(a, b) for (a, b) in E2 if (a, b) ∉ Set(bad_edges)]
  new_E2 == E2 && return (N2, E2)
  return _clean_tree(N2, new_E2)
end

# Remove disconnected nodes (no parent except root) and reindex edges.
function _reindex_tree(N::Vector, E::Vector{<:Tuple{Int,Int}})
  (_, parents) = _children_and_parents(N, E)
  all_nodes = 1:length(N)
  dead = Set(i for i in 2:length(N) if isempty(parents[i]))
  live = [i for i in all_nodes if i ∉ dead]
  rank = Dict(v => k for (k, v) in enumerate(live))
  new_N = [N[i] for i in live]
  new_E = [(rank[p], rank[c]) for (p, c) in E if p ∉ dead && c ∉ dead]
  return (new_N, new_E)
end

# For affine mode: remove nodes contained in the hyperplane at infinity
# (first generator b_0 of the homogenized ring).
function _affine_prune(N::AbstractVector{<:MPolyIdeal}, E::Vector{<:Tuple{Int,Int}})
  isempty(N) && return (collect(N), E)
  R = base_ring(first(N))
  b0 = gen(R, 1)
  inf_ideal = ideal(R, [b0])

  bad = Set(i for i in 1:length(N) if is_subset(inf_ideal, N[i]))
  new_E = [(p, c) for (p, c) in E if c ∉ bad]

  # De-homogenize: send b_0 → 1 in remaining nodes
  affine_R, affine_vars = polynomial_ring(coefficient_ring(R), symbols(R)[2:end])
  to_affine = hom(R, affine_R, vcat([one(affine_R)], collect(affine_vars)))

  new_N = [ideal(affine_R, [to_affine(g) for g in gens(N[i])]) for i in 1:length(N)]
  (new_N, new_E) = _reindex_tree(new_N, new_E)
  return _subtree(new_N, new_E, 1)
end

################################################################################
# Phase 4 – Tree builder and public API
################################################################################

# BFS construction of the image tree.  T is the pre-created target ring.
# Returns (nodes::Vector{MPolyIdeal}, edges::Vector{Tuple{Int,Int}}) with 1-based indexing.
function _tree_builder(L::Vector{<:MPolyRingElem}, X::MPolyIdeal, T::MPolyRing;
                       verbose::Bool = false, verify::Bool = true,
                       tries::Int = 50, minimum_dim::Int = 0)
  # Each domain entry: (domain_ideal, parent_node_index)
  # parent_node_index == 0 means "root" (no parent edge)
  domains = [(X, 0)]

  nodes = MPolyIdeal[]
  edges = Tuple{Int,Int}[]

  j = 1
  while j <= length(domains)
    (D, parent_idx) = domains[j]
    j += 1

    verbose && println("tree_builder: processing domain $(j-1) of $(length(domains)) ...")

    (zariski_image, exc_loci) = _partial_image_with_target(
      L, D, T; verbose = verbose, verify = verify, pd = true, tries = tries
    )

    # Filter by minimum dimension
    exc_loci = [e for e in exc_loci if dim(e) > minimum_dim]

    # Compute pullback domains for each exceptional locus (pulled back to D, not X)
    exc_domains = [_components_of_pullback(L, D, e) for e in exc_loci]

    # Index new exceptional domains
    # The new nodes will be: zariski_image at position (length(nodes)+1),
    # then exceptional loci at positions (length(nodes)+2), ..., (length(nodes)+1+length(exc_loci))
    zariski_idx = length(nodes) + 1

    # Add edge from parent to this zariski image
    if parent_idx != 0
      push!(edges, (parent_idx, zariski_idx))
    end

    # Add edges from zariski image to each exceptional locus
    for k in 1:length(exc_loci)
      exc_idx = zariski_idx + k
      push!(edges, (zariski_idx, exc_idx))
    end

    push!(nodes, zariski_image)
    append!(nodes, exc_loci)

    # Add exceptional domains to the queue with their parent indices
    for (k, exc_domain_list) in enumerate(exc_domains)
      exc_locus_idx = zariski_idx + k
      for Z in exc_domain_list
        push!(domains, (Z, exc_locus_idx))
      end
    end
  end

  return (nodes, edges)
end

@doc raw"""
    total_image(L::Vector{<:MPolyRingElem}; verbose::Bool=false, clean::Bool=true,
                verify::Bool=true, affine::Bool=false, tries::Int=50,
                minimum_dim::Int=0)
    total_image(L::Vector{<:MPolyRingElem}, X::MPolyIdeal; verbose::Bool=false,
                clean::Bool=true, verify::Bool=true, affine::Bool=false,
                tries::Int=50, minimum_dim::Int=0)

Compute the constructible set which is the image of the rational map defined
by `L` on the projective variety `V(X)` (default: whole projective space).

The image is represented as a tree of varieties: the root is the Zariski closure
of the image, and each non-root node is an exceptional divisor where the image
is not closed.  The tree is returned as a pair `(nodes, edges)` where
- `nodes` is a vector of ideals in the target ring `k[b_0,...,b_m]`
- `edges` is a vector of pairs `(i, j)` indicating that node `j` is a child
  of node `i`

The image tree is also printed in human-readable form.

Options:
- `verbose`: print progress information
- `clean`: remove redundant nodes from the tree
- `verify`: use deterministic hyperplane sections (slower but more reliable)
- `affine`: treat the map as an affine map (also auto-activated if polynomials
  are not all homogeneous of the same degree)
- `tries`: maximum number of random hyperplane section attempts
- `minimum_dim`: ignore exceptional loci of dimension ≤ this value

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);
julia> L = [y*z, x*z, x*y];  # Cremona transformation
julia> (N, E) = total_image(L);
julia> sort(dim.(N))
[1, 1, 1, 1, 1, 1, 2, 2, 2, 3]
```
"""
function total_image(L::Vector{<:MPolyRingElem}, X::MPolyIdeal;
                     verbose::Bool = false, clean::Bool = true,
                     verify::Bool = true, affine::Bool = false,
                     tries::Int = 50, minimum_dim::Int = 0)
  R = parent(first(L))
  @assert base_ring(X) === R "Domain ideal must be over the same ring as L"

  do_affine = affine || !_all_same_degree(L)

  if do_affine
    verbose && println("total_image: homogenizing ...")
    (L, X, R) = _homogenize_input(L, X)
  end

  T = _generate_target_ring(L)

  tree = _tree_builder(L, X, T; verbose = verbose, verify = verify,
                        tries = tries, minimum_dim = minimum_dim)

  if clean
    verbose && println("total_image: cleaning tree ...")
    (nodes, edges) = tree
    (nodes, edges) = _remove_duplicates(nodes, edges)
    (nodes, edges) = _reindex_tree(nodes, edges)
    (nodes, edges) = _clean_tree(nodes, edges)
    if do_affine
      (nodes, edges) = _affine_prune(nodes, edges)
      # De-homogenization can reduce node dimensions; prune any new zero-dim leaves.
      (nodes, edges) = _prune_leaves(nodes, edges)
      (nodes, edges) = _reindex_tree(nodes, edges)
    end
    tree = (nodes, edges)
  end

  print_image_tree(tree[1], tree[2])
  return tree
end

function total_image(L::Vector{<:MPolyRingElem};
                     verbose::Bool = false, clean::Bool = true,
                     verify::Bool = true, affine::Bool = false,
                     tries::Int = 50, minimum_dim::Int = 0)
  R = parent(first(L))
  return total_image(L, ideal(R, elem_type(R)[]); verbose = verbose, clean = clean,
                     verify = verify, affine = affine, tries = tries,
                     minimum_dim = minimum_dim)
end

# Check if all polynomials in L are homogeneous and of the same degree.
function _all_same_degree(L::Vector{<:MPolyRingElem})
  isempty(L) && return true
  all(is_homogeneous, L) || return false
  d = total_degree(first(L))
  return all(f -> total_degree(f) == d, L)
end

# Manually homogenize polynomial f (already in S) to degree d using gen(S, zhom_idx).
function _homogenize_poly(f::MPolyRingElem, S::MPolyRing, zhom_idx::Int, d::Int)
  zhom = gen(S, zhom_idx)
  result = zero(S)
  for (c, m) in zip(coefficients(f), monomials(f))
    k = total_degree(m)
    result += c * m * zhom^(d - k)
  end
  return result
end

# Homogenize the input for affine mode: add variable zhom at position 1.
function _homogenize_input(L::Vector{<:MPolyRingElem}, X::MPolyIdeal)
  R = parent(first(L))
  k = coefficient_ring(R)
  new_names = vcat([:zhom], collect(symbols(R)))
  R_hom, R_hom_vars = polynomial_ring(k, new_names)
  zhom = R_hom_vars[1]
  R_to_Rhom = hom(R, R_hom, R_hom_vars[2:end])

  max_deg = maximum(total_degree(f) for f in L)

  # All components made degree max_deg; prepend zhom^max_deg as the new b_0
  new_L = vcat(
    [zhom^max_deg],
    [_homogenize_poly(R_to_Rhom(f), R_hom, 1, max_deg) for f in L]
  )
  new_X = ideal(R_hom, [R_to_Rhom(g) for g in gens(X)])

  return (new_L, new_X, R_hom)
end

@doc raw"""
    is_closed_image(L::Vector{<:MPolyRingElem}; verbose::Bool=false, tries::Int=50)
    is_closed_image(L::Vector{<:MPolyRingElem}, X::MPolyIdeal;
                    verbose::Bool=false, tries::Int=50)

Return `true` if the image of the rational map defined by `L` on `V(X)` is
(Zariski) closed, i.e., equal to its Zariski closure.

This is done by checking that the exceptional divisor produced by
`partial_image` has the same dimension as its pullback, and recursively
checking each prime component.
"""
function is_closed_image(L::Vector{<:MPolyRingElem}, X::MPolyIdeal;
                         verbose::Bool = false, tries::Int = 50)
  R = parent(first(L))
  T = _generate_target_ring(L)

  (image_X, exc_loci_vec) = _partial_image_with_target(
    L, X, T; verbose = verbose, pd = false, tries = tries
  )
  isempty(exc_loci_vec) && return true

  exc_image = only(exc_loci_vec)  # pd=false → at most one element
  dim(exc_image) < 1 && return true

  T_to_R = hom(T, R, L)
  pullback = ideal(R, [T_to_R(g) for g in gens(exc_image)]) + X
  pullback_image = _scheme_theoretic_image(L, saturation(pullback, ideal(R, L)), T)

  dim(exc_image) != dim(pullback_image) && return false

  verbose && println("is_closed_image: decomposing image of exceptional divisor ...")
  pd_exc = minimal_primes(exc_image)

  for p in pd_exc
    sat_pullback = saturation(ideal(R, [T_to_R(g) for g in gens(p)]) + X, ideal(R, L))
    is_closed_image(L, sat_pullback; verbose = verbose, tries = tries) || return false
  end

  return true
end

function is_closed_image(L::Vector{<:MPolyRingElem};
                         verbose::Bool = false, tries::Int = 50)
  R = parent(first(L))
  return is_closed_image(L, ideal(R, elem_type(R)[]); verbose = verbose, tries = tries)
end

################################################################################
# Phase 5 – Pretty printing
################################################################################

@doc raw"""
    print_image_tree(N::Vector{<:MPolyIdeal}, E::Vector{<:Tuple{Int,Int}})

Print the image tree returned by `total_image` in a human-readable format.
Each line is prefixed by `(d)` where `d` is the projective dimension of the
variety, followed by indentation showing the tree structure.
"""
function print_image_tree(N::AbstractVector{<:MPolyIdeal},
                          E::Vector{<:Tuple{Int,Int}})
  isempty(N) && return
  (children, _) = _children_and_parents(N, E)
  println()
  _print_node(N, children, 1, 0)
end

function _print_node(N::AbstractVector{<:MPolyIdeal},
                     children::Vector{Vector{Int}},
                     node::Int, level::Int)
  d = dim(N[node]) - 1  # projective dimension = Krull dim - 1
  prefix = "($d) "
  if level > 0
    indent = repeat("|    ", level - 1)
    branch = level % 2 == 0 ? " + " : " - "
    prefix = branch * prefix * indent * "|===="
  else
    prefix = "   " * prefix
  end
  println(prefix, N[node])
  for c in children[node]
    _print_node(N, children, c, level + 1)
  end
end

end # module RationalMapImage

using .RationalMapImage
