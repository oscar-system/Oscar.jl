# Translation of the Macaulay2 package TotalImage by Corey Harris,
# Mateusz Michalek, and Emre Sertöz, which accompanies the paper
# "Computing images of polynomial maps" (2017).
#
# The algorithm computes the image of a rational/polynomial map as a
# constructible set, represented as a tree of projective varieties.
################################################################################

module TotalImage

using ..Oscar

import AbstractAlgebra
import ..Oscar: edges, is_affine, vertices

export ConstructibleSetTree
export is_contained
export is_closed_image
export partial_image
export total_image

function __init__()
  add_verbosity_scope(:TotalImage)
end

@doc raw"""
    ConstructibleSetTree

A rooted tree representing the constructible image of a rational map.

The vertices are closed sets, represented by ideals in the target coordinate
ring. They are homogeneous for projective trees. The root is added to the
constructible set, its children are removed, their children are added, and so
on. Use `vertices` and `edges` to access copies of the vertex and edge lists.
Edge indices are one-based.
"""
struct ConstructibleSetTree{T<:MPolyIdeal}
  ideals::Vector{T}
  edge_list::Vector{Tuple{Int, Int}}
  affine::Bool
end

function ConstructibleSetTree(
    ideals::Vector{T}, edge_list::Vector{Tuple{Int, Int}}; affine::Bool=false
  ) where {T<:MPolyIdeal}
  @req !isempty(ideals) "An image tree must have a root"
  @req all(
    e -> 1 <= e[1] <= length(ideals) && 1 <= e[2] <= length(ideals), edge_list
  ) "Edge index out of range"
  return ConstructibleSetTree{T}(ideals, edge_list, affine)
end

@doc raw"""
    image_ideals(tree::ConstructibleSetTree)

Return a copy of the ideals at the vertices of `tree`, starting with its root.
"""
image_ideals(tree::ConstructibleSetTree) = copy(tree.ideals)

vertices(tree::ConstructibleSetTree) = image_ideals(tree)
edges(tree::ConstructibleSetTree) = copy(tree.edge_list)
is_affine(tree::ConstructibleSetTree) = tree.affine
Base.length(tree::ConstructibleSetTree) = length(tree.ideals)

function Base.show(io::IO, tree::ConstructibleSetTree)
  print(
    io,
    "Constructible set tree with ",
    length(tree),
    " vertices and ",
    length(tree.edge_list),
    " edges",
  )
end

function Base.show(io::IO, ::MIME"text/plain", tree::ConstructibleSetTree)
  output_tree(io, tree)
end

function _homogeneous_degree(f::MPolyRingElem)
  iszero(f) && return nothing
  degrees = unique(sum(e) for e in AbstractAlgebra.exponent_vectors(f))
  length(degrees) == 1 || return nothing
  return only(degrees)
end

function _is_homogeneous(f::MPolyRingElem)
  iszero(f) && return true
  return !isnothing(_homogeneous_degree(f))
end

function _same_degree(L::AbstractVector{<:MPolyRingElem})
  degree = nothing
  for f in L
    iszero(f) && continue
    d = _homogeneous_degree(f)
    isnothing(d) && return false
    if isnothing(degree)
      degree = d
    elseif degree != d
      return false
    end
  end
  return !isnothing(degree)
end

function _validate_input(L::AbstractVector{<:MPolyRingElem}, X::MPolyIdeal)
  @req !isempty(L) "The defining list must not be empty"
  R = parent(first(L))
  @req all(f -> parent(f) === R, L) "All defining polynomials must have the same parent"
  @req base_ring(X) === R "The source ideal and defining polynomials must have the same base ring"
  @req any(!iszero, L) "The defining polynomials must not all be zero"
  @req coefficient_ring(R) isa Field "The coefficient ring must be a field"
  return R
end

function _validate_projective_input(
    L::AbstractVector{<:MPolyRingElem}, X::MPolyIdeal, T::MPolyRing
  )
  R = _validate_input(L, X)
  @req _same_degree(L) "The defining polynomials must be homogeneous of the same degree"
  @req all(_is_homogeneous, gens(X)) "The source ideal must be homogeneous"
  @req ngens(T) == length(L) "The target ring must have one variable for each defining polynomial"
  same_coefficient_ring = coefficient_ring(T) === coefficient_ring(R)
  @req same_coefficient_ring "The source and target must have the same coefficient ring"
  return R
end

function _target(L::AbstractVector{<:MPolyRingElem})
  R = parent(first(L))
  return polynomial_ring(coefficient_ring(R), :b => 0:(length(L) - 1))[1]
end

function _map_ideal(phi, I::MPolyIdeal)
  return ideal(codomain(phi), phi.(gens(I)))
end

# Scheme-theoretic Zariski closure of the image of V(X) under the map defined
# by L. Uses Oscar's built-in preimage via Singular.
function _scheme_theoretic_image(
    L::Vector{<:MPolyRingElem}, X::MPolyIdeal, T::MPolyRing
  )
  f = hom(T, base_ring(X), L; check=false)
  return preimage(f, X)
end

_projectively_empty(I::MPolyIdeal) = dim(I) < 1

function _empty_exception(image::MPolyIdeal, decompose::Bool)
  if decompose
    return typeof(image)[]
  end
  T = base_ring(image)
  return ideal(T, [one(T)])
end

function _coordinate_section(R::MPolyRing, i::Int)
  variable_names = copy(symbols(R))
  deleteat!(variable_names, i)
  S, s = polynomial_ring(coefficient_ring(R), variable_names; cached=false)
  images = copy(s)
  insert!(images, i, zero(S))
  return hom(R, S, images; check=false)
end

function _binomial_section(R::MPolyRing, i::Int, j::Int, a::Int)
  @assert i < j
  variable_names = copy(symbols(R))
  deleteat!(variable_names, i)
  S, s = polynomial_ring(coefficient_ring(R), variable_names; cached=false)
  images = elem_type(S)[]
  for k in 1:ngens(R)
    if k == i
      push!(images, a*s[j - 1])
    elseif k < i
      push!(images, s[k])
    else
      push!(images, s[k - 1])
    end
  end
  return hom(R, S, images; check=false)
end

function _random_section(R::MPolyRing, number_of_sections::Int)
  @assert 0 <= number_of_sections < ngens(R)
  if iszero(number_of_sections)
    return hom(R, R, gens(R); check=false)
  end
  variable_names = symbols(R)[(number_of_sections + 1):end]
  S, s = polynomial_ring(coefficient_ring(R), variable_names; cached=false)
  images = elem_type(S)[]
  for _ in 1:number_of_sections
    linear_form = zero(S)
    for v in s
      linear_form += rand(-200:200)*v
    end
    push!(images, linear_form)
  end
  append!(images, s)
  return hom(R, S, images; check=false)
end

function _restrict(
    L::AbstractVector{<:MPolyRingElem}, X::MPolyIdeal, T::MPolyRing, phi
  )
  restricted_L = phi.(L)
  restricted_X = _map_ideal(phi, X)
  restricted_image = _scheme_theoretic_image(restricted_L, restricted_X, T)
  return restricted_L, restricted_X, restricted_image
end

function _good_section(
    L::AbstractVector{<:MPolyRingElem},
    X::MPolyIdeal,
    image::MPolyIdeal,
    T::MPolyRing,
    phi,
  )
  restricted_L, restricted_X, restricted_image = _restrict(L, X, T, phi)
  good_dimension = dim(restricted_X) == dim(X) - 1
  good_image = is_subset(restricted_image, image)
  return good_dimension && good_image, restricted_L, restricted_X
end

function _shrink(
    L::AbstractVector{<:MPolyRingElem},
    X::MPolyIdeal,
    image::MPolyIdeal,
    T::MPolyRing;
    tries::Int,
  )
  image_dimension = dim(image)
  while dim(X) > image_dimension
    found = false
    R = base_ring(X)
    for i in 1:ngens(R)
      phi = _coordinate_section(R, i)
      good, restricted_L, restricted_X = _good_section(L, X, image, T, phi)
      if good
        @vprintln :TotalImage "Found a coordinate hyperplane section (variable $i)"
        L, X = restricted_L, restricted_X
        found = true
        break
      end
    end
    found && continue

    R = base_ring(X)
    for a in (1, -1, 2, -2, -3, 3, 5, 7)
      for i in 1:(ngens(R) - 1)
        for j in (i + 1):ngens(R)
          phi = _binomial_section(R, i, j, a)
          good, restricted_L, restricted_X = _good_section(L, X, image, T, phi)
          if good
            @vprintln :TotalImage "Found a binomial hyperplane section (coefficient $a)"
            L, X = restricted_L, restricted_X
            found = true
            break
          end
        end
        found && break
      end
      found && break
    end
    found && continue

    number_of_sections = dim(X) - image_dimension
    for attempt in 1:tries
      phi = _random_section(base_ring(X), number_of_sections)
      restricted_L, restricted_X, restricted_image = _restrict(L, X, T, phi)
      if dim(restricted_X) == image_dimension && is_subset(restricted_image, image)
        @vprintln :TotalImage "Found a random linear section (attempt $attempt)"
        return restricted_L, restricted_X
      end
    end
    error("Could not find a linear section preserving the image after $tries tries")
  end
  return L, X
end

function _blowup_graph(L::AbstractVector{<:MPolyRingElem})
  R = parent(first(L))
  K = coefficient_ring(R)
  number_of_source_variables = ngens(R)
  number_of_target_variables = length(L)
  source_names = [Symbol("_graph_source_$i") for i in 1:number_of_source_variables]
  target_names = [Symbol("_graph_target_$i") for i in 1:number_of_target_variables]
  W, w = polynomial_ring(K, vcat([:_graph_parameter], source_names, target_names); cached=false)
  source_in_W = w[2:(number_of_source_variables + 1)]
  target_in_W = w[(number_of_source_variables + 2):end]
  R_to_W = hom(R, W, source_in_W; check=false)
  graph_equations = [
    target_in_W[i] - w[1]*R_to_W(L[i]) for i in 1:number_of_target_variables
  ]
  graph_with_parameter = ideal(W, graph_equations)
  eliminated_graph = eliminate(graph_with_parameter, [w[1]])

  U, u = polynomial_ring(K, vcat(source_names, target_names); cached=false)
  U_to_W = hom(U, W, w[2:end]; check=false)
  graph = preimage(U_to_W, eliminated_graph)
  return graph, U, u[1:number_of_source_variables], u[(number_of_source_variables + 1):end]
end

function _exceptional_image(
    L::AbstractVector{<:MPolyRingElem}, X::MPolyIdeal, T::MPolyRing
  )
  R = base_ring(X)
  graph, U, source_variables, target_variables = _blowup_graph(L)
  R_to_U = hom(R, U, source_variables; check=false)
  T_to_U = hom(T, U, target_variables; check=false)
  base_locus = ideal(R, collect(L))
  base_locus_in_graph = graph + _map_ideal(R_to_U, base_locus)

  proper_source = graph
  if !is_zero(X)
    source_in_graph = graph + _map_ideal(R_to_U, X)
    proper_source = saturation(source_in_graph, base_locus_in_graph)
  end
  @vprintln :TotalImage "Computed the proper transform of the domain"

  exceptional = proper_source + base_locus_in_graph
  exceptional = saturation(exceptional, ideal(U, target_variables))
  exceptional = saturation(exceptional, ideal(U, source_variables))
  projected_exceptional = eliminate(exceptional, source_variables)
  result = preimage(T_to_U, projected_exceptional)
  @vprintln :TotalImage "Computed the image of the exceptional divisor"
  return result
end

@doc raw"""
    partial_image(L[, X[, T]]; decompose=true, verify=true, tries=50)

Compute the Zariski closure of the image of the projective rational map defined
by the homogeneous polynomials in `L`, restricted to the projective domain
defined by `X`.

The second return value describes the possible complement of the actual image.
When `decompose` is `true`, it is the vector of minimal primes of the image of
the exceptional divisor. Otherwise it is that image as a single ideal. If `X`
is omitted, the whole projective space is used as the domain. If `T` is
omitted, a target ring with variables `b[0], ..., b[m]` is created.

When the generic fibers are positive dimensional, `verify=true` checks that a
linear section preserves the image before the blowup is computed. Up to
`tries` random sections are attempted after deterministic coordinate and
binomial sections.
"""
function partial_image(
    L::AbstractVector{<:MPolyRingElem};
    decompose::Bool=true,
    verify::Bool=true,
    tries::Int=50,
  )
  @req !isempty(L) "The defining list must not be empty"
  R = parent(first(L))
  X = ideal(R, [zero(R)])
  return partial_image(L, X; decompose, verify, tries)
end

function partial_image(
    L::AbstractVector{<:MPolyRingElem},
    X::MPolyIdeal;
    decompose::Bool=true,
    verify::Bool=true,
    tries::Int=50,
  )
  T = _target(L)
  return partial_image(L, X, T; decompose, verify, tries)
end

function partial_image(
    L::AbstractVector{<:MPolyRingElem},
    X::MPolyIdeal,
    T::MPolyRing;
    decompose::Bool=true,
    verify::Bool=true,
    tries::Int=50,
  )
  @req tries > 0 "The number of tries must be positive"
  _validate_projective_input(L, X, T)
  image = _scheme_theoretic_image(collect(L), X, T)
  @vprintln :TotalImage "Computed the Zariski closure of the image (dimension $(dim(image) - 1))"

  base_locus = ideal(base_ring(X), collect(L))
  if _projectively_empty(X + base_locus)
    return image, _empty_exception(image, decompose)
  end

  if dim(X) > dim(image)
    if verify
      L, X = _shrink(L, X, image, T; tries)
    else
      number_of_sections = dim(X) - dim(image)
      phi = _random_section(base_ring(X), number_of_sections)
      L, X, _ = _restrict(L, X, T, phi)
    end
    base_locus = ideal(base_ring(X), collect(L))
  end

  if _projectively_empty(X + base_locus)
    return image, _empty_exception(image, decompose)
  end

  exceptional_image = _exceptional_image(L, X, T)
  if _projectively_empty(exceptional_image)
    return image, _empty_exception(image, decompose)
  end
  if decompose
    return image, minimal_primes(exceptional_image)
  end
  return image, exceptional_image
end

function _components_of_pullback(
    L::AbstractVector{<:MPolyRingElem}, X::MPolyIdeal, exceptional_image::MPolyIdeal
  )
  R = base_ring(X)
  T = base_ring(exceptional_image)
  f = hom(T, R, collect(L); check=false)
  base_locus = ideal(R, collect(L))
  pullback = _map_ideal(f, exceptional_image) + X
  pullback = saturation(pullback, base_locus)
  _projectively_empty(pullback) && return typeof(X)[]
  return minimal_primes(pullback)
end

function _children_and_parents(number_of_nodes::Int, edge_list)
  children = [Int[] for _ in 1:number_of_nodes]
  parents = [Int[] for _ in 1:number_of_nodes]
  for (source, target) in edge_list
    push!(children[source], target)
    push!(parents[target], source)
  end
  return children, parents
end

function _is_contained_at(
    x::Vector, C::ConstructibleSetTree, children, node::Int
  )
  all(f -> iszero(evaluate(f, x)), gens(C.ideals[node])) || return false
  return !any(child -> _is_contained_at(x, C, children, child), children[node])
end

@doc raw"""
    is_contained(x::Vector, C::ConstructibleSetTree) -> Bool

Return `true` if the point `x` is contained in the constructible set `C`.

The point `x` must be given as a vector of coordinates matching the variables
of the target ring of `C`. For a projective map (the default), these are
homogeneous coordinates `[b_0 : b_1 : \cdots : b_m]`; for an affine map
(computed with `affine = true`) they are affine coordinates `[b_1, \ldots, b_m]`
(with `b_0 = 1` substituted).

The constructible set is the actual image (not just the Zariski closure). A
point lies in it if and only if it belongs to the Zariski closure of the image
(the root of the tree) and, whenever it falls in an exceptional subvariety,
it also lies in the sub-constructible set of that subvariety.
"""
function is_contained(x::Vector, C::ConstructibleSetTree)
  T = base_ring(first(C.ideals))
  @req length(x) == ngens(T) "The number of coordinates must match the target ring"
  @req C.affine || any(!iszero, x) "Projective coordinates must not all be zero"
  children, _ = _children_and_parents(length(C), C.edge_list)
  return _is_contained_at(x, C, children, 1)
end

function _reachable_tree(nodes, edge_list)
  children, _ = _children_and_parents(length(nodes), edge_list)
  reachable = falses(length(nodes))
  queue = [1]
  reachable[1] = true
  position = 1
  while position <= length(queue)
    current = queue[position]
    position += 1
    for child in children[current]
      if !reachable[child]
        reachable[child] = true
        push!(queue, child)
      end
    end
  end
  live = findall(reachable)
  new_index = zeros(Int, length(nodes))
  for (i, old_index) in enumerate(live)
    new_index[old_index] = i
  end
  new_nodes = nodes[live]
  new_edges = Tuple{Int, Int}[]
  for (source, target) in edge_list
    if reachable[source] && reachable[target]
      push!(new_edges, (new_index[source], new_index[target]))
    end
  end
  return new_nodes, new_edges
end

function _levels(number_of_nodes::Int, edge_list)
  children, _ = _children_and_parents(number_of_nodes, edge_list)
  levels = Vector{Vector{Int}}([[1]])
  seen = falses(number_of_nodes)
  seen[1] = true
  while true
    next_level = Int[]
    for node in last(levels)
      for child in children[node]
        if !seen[child]
          seen[child] = true
          push!(next_level, child)
        end
      end
    end
    isempty(next_level) && break
    push!(levels, next_level)
  end
  return levels
end

function _remove_duplicates(nodes, edge_list)
  while true
    children, parents = _children_and_parents(length(nodes), edge_list)
    levels = _levels(length(nodes), edge_list)
    duplicate_pair = nothing
    for level_index in 2:2:length(levels)
      for node in levels[level_index]
        duplicate = findfirst(child -> nodes[child] == nodes[node], children[node])
        if !isnothing(duplicate)
          duplicate_pair = (node, children[node][duplicate])
          break
        end
      end
      !isnothing(duplicate_pair) && break
    end
    isnothing(duplicate_pair) && return nodes, edge_list

    node, duplicate = duplicate_pair
    @assert length(parents[node]) == 1
    parent = only(parents[node])
    duplicate_children = copy(children[duplicate])
    filter!(e -> !(node in e || duplicate in e), edge_list)
    append!(edge_list, [(parent, child) for child in duplicate_children])
    nodes, edge_list = _reachable_tree(nodes, edge_list)
  end
end

function _prune_empty_leaves(nodes, edge_list)
  children, parents = _children_and_parents(length(nodes), edge_list)
  empty_leaves = [
    i for i in 2:length(nodes) if isempty(children[i]) && !isempty(parents[i]) &&
    _projectively_empty(nodes[i])
  ]
  filter!(e -> !(e[2] in empty_leaves), edge_list)
  return _reachable_tree(nodes, edge_list)
end

function _clean_tree(nodes, edge_list)
  nodes, edge_list = _remove_duplicates(nodes, edge_list)
  nodes, edge_list = _prune_empty_leaves(nodes, edge_list)
  while true
    children, parents = _children_and_parents(length(nodes), edge_list)
    bad_leaf = nothing
    for leaf in 2:length(nodes)
      isempty(children[leaf]) || continue
      length(parents[leaf]) == 1 || continue
      parent = only(parents[leaf])
      if parent != 1 && dim(nodes[parent]) == dim(nodes[leaf])
        bad_leaf = leaf
        break
      end
    end
    isnothing(bad_leaf) && return nodes, edge_list
    parent = only(parents[bad_leaf])
    grandparent = only(parents[parent])
    filter!(e -> e != (parent, bad_leaf) && e != (grandparent, parent), edge_list)
    nodes, edge_list = _reachable_tree(nodes, edge_list)
  end
end

function _tree_builder(
    L::AbstractVector{<:MPolyRingElem},
    X::MPolyIdeal;
    verify::Bool,
    tries::Int,
    minimum_dimension::Int,
  )
  T = _target(L)
  domains = Tuple{typeof(X), Int}[(X, 0)]
  nodes = MPolyIdeal[]
  edge_list = Tuple{Int, Int}[]
  domain_index = 1
  while domain_index <= length(domains)
    domain, parent_index = domains[domain_index]
    domain_index += 1
    image, exceptional_loci = partial_image(L, domain, T; decompose=true, verify, tries)
    filter!(E -> dim(E) > minimum_dimension, exceptional_loci)

    image_index = length(nodes) + 1
    exceptional_indices = (image_index + 1):(image_index + length(exceptional_loci))
    parent_index != 0 && push!(edge_list, (parent_index, image_index))
    append!(nodes, [image])
    append!(nodes, exceptional_loci)
    append!(edge_list, [(image_index, i) for i in exceptional_indices])

    for (exceptional_image, exceptional_index) in zip(exceptional_loci, exceptional_indices)
      components = _components_of_pullback(L, X, exceptional_image)
      append!(domains, [(component, exceptional_index) for component in components])
    end
  end
  typed_nodes = typeof(first(nodes))[node for node in nodes]
  return typed_nodes, edge_list
end

function _homogenize_affine(
    L::AbstractVector{<:MPolyRingElem}, X::MPolyIdeal
  )
  R = base_ring(X)
  H = homogenizer(R, :_homogenizing_coordinate; pos=1)
  homogeneous_X = H(X)
  decorated_ring = base_ring(homogeneous_X)
  h = gen(decorated_ring, 1)
  maximum_degree = maximum(total_degree(f) for f in L if !iszero(f))
  target_degree = maximum_degree + 1
  homogeneous_L = elem_type(forget_decoration(decorated_ring))[]
  push!(homogeneous_L, forget_decoration(h^target_degree))
  for f in L
    if iszero(f)
      push!(homogeneous_L, zero(forget_decoration(decorated_ring)))
    else
      push!(
        homogeneous_L,
        forget_decoration(h^(target_degree - total_degree(f))*H(f)),
      )
    end
  end
  return homogeneous_L, forget_decoration(homogeneous_X)
end

function _affine_prune(nodes, edge_list)
  T = base_ring(first(nodes))
  hyperplane_at_infinity = ideal(T, [gen(T, 1)])
  bad_nodes = [i for i in eachindex(nodes) if is_subset(hyperplane_at_infinity, nodes[i])]
  filter!(e -> !(e[2] in bad_nodes), edge_list)
  nodes, edge_list = _reachable_tree(nodes, edge_list)

  affine_target, affine_variables = polynomial_ring(
    coefficient_ring(T), symbols(T)[2:end]; cached=false
  )
  restriction = hom(T, affine_target, vcat([one(affine_target)], affine_variables); check=false)
  affine_nodes = [ideal(affine_target, restriction.(gens(node))) for node in nodes]
  return affine_nodes, edge_list
end

@doc raw"""
    total_image(L[, X]; clean=true, verify=true, affine=false, tries=50,
                minimum_dimension=0)

Compute the actual image of a rational map as a constructible set represented
by a [`ConstructibleSetTree`](@ref).

The map is defined by `L` on the domain cut out by `X`. If `X` is omitted, the
whole projective space is used as the domain. Projective input requires
homogeneous coordinates of one degree. Nonhomogeneous input is treated as
affine and homogenized automatically; set `affine=true` to force that
interpretation.

The root ideal is the Zariski closure of the image. At odd depths the ideals
are removed, and at even depths they are added back. Components of projective
dimension less than `minimum_dimension` are omitted. Set `clean=false` to keep
redundant branches.

Use `set_verbosity_level(:TotalImage, 1)` to print progress information.
"""
function total_image(
    L::AbstractVector{<:MPolyRingElem};
    clean::Bool=true,
    verify::Bool=true,
    affine::Bool=false,
    tries::Int=50,
    minimum_dimension::Int=0,
  )
  @req !isempty(L) "The defining list must not be empty"
  R = parent(first(L))
  X = ideal(R, [zero(R)])
  return total_image(L, X; clean, verify, affine, tries, minimum_dimension)
end

function total_image(
    L::AbstractVector{<:MPolyRingElem},
    X::MPolyIdeal;
    clean::Bool=true,
    verify::Bool=true,
    affine::Bool=false,
    tries::Int=50,
    minimum_dimension::Int=0,
  )
  _validate_input(L, X)
  @req tries > 0 "The number of tries must be positive"
  @req minimum_dimension >= 0 "The minimum dimension must be nonnegative"
  affine = affine || !_same_degree(L) || !all(_is_homogeneous, gens(X))
  if affine
    L, X = _homogenize_affine(L, X)
  end
  nodes, edge_list = _tree_builder(L, X; verify, tries, minimum_dimension)
  if clean
    nodes, edge_list = _clean_tree(nodes, edge_list)
  end
  if affine
    nodes, edge_list = _affine_prune(nodes, edge_list)
  end
  return ConstructibleSetTree(nodes, edge_list; affine)
end

@doc raw"""
    is_closed_image(L[, X]; verify=true, tries=50)

Return whether the actual image of the projective rational map defined by `L`
on the domain defined by `X` is Zariski closed.
"""
function is_closed_image(
    L::AbstractVector{<:MPolyRingElem};
    verify::Bool=true,
    tries::Int=50,
  )
  @req !isempty(L) "The defining list must not be empty"
  R = parent(first(L))
  return is_closed_image(L, ideal(R, [zero(R)]); verify, tries)
end

function is_closed_image(
    L::AbstractVector{<:MPolyRingElem},
    X::MPolyIdeal;
    verify::Bool=true,
    tries::Int=50,
  )
  image, exceptional_image = partial_image(L, X; decompose=false, verify, tries)
  _projectively_empty(exceptional_image) && return true
  R = base_ring(X)
  T = base_ring(image)
  f = hom(T, R, collect(L); check=false)
  pullback = _map_ideal(f, exceptional_image) + X
  image_of_pullback = preimage(f, pullback)
  dim(exceptional_image) == dim(image_of_pullback) || return false
  @vprintln :TotalImage "Decomposing the image of the exceptional divisor"
  base_locus = ideal(R, collect(L))
  for prime in minimal_primes(exceptional_image)
    restricted_source = saturation(_map_ideal(f, prime) + X, base_locus)
    _projectively_empty(restricted_source) && continue
    is_closed_image(L, restricted_source; verify, tries) || return false
  end
  return true
end

function _display_dimension(tree::ConstructibleSetTree, I::MPolyIdeal)
  d = dim(I)
  tree.affine && return d
  return d isa NegInf ? d : d - 1
end

function _print_node(
    io::IO, tree::ConstructibleSetTree, children, node::Int, depth::Int
  )
  dimension = _display_dimension(tree, tree.ideals[node])
  if iszero(depth)
    print(io, "   (", dimension, ") ")
  else
    sign = iseven(depth) ? '+' : '-'
    print(io, ' ', sign, " (", dimension, ") ", repeat("|    ", depth - 1), "|====")
  end
  show(io, tree.ideals[node])
  println(io)
  for child in children[node]
    _print_node(io, tree, children, child, depth + 1)
  end
end

@doc raw"""
    output_tree([io::IO,] tree::ConstructibleSetTree)

Print the alternating closed-set tree describing the constructible image.
Projective dimensions are shown for projective trees and Krull dimensions for
affine trees.
"""
function output_tree(io::IO, tree::ConstructibleSetTree)
  children, _ = _children_and_parents(length(tree), tree.edge_list)
  _print_node(io, tree, children, 1, 0)
  return nothing
end

output_tree(tree::ConstructibleSetTree) = output_tree(stdout, tree)

end
