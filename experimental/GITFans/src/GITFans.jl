"""
  GITFans is a collection of Julia functions for computations
  like those shown in the paper [BKR20](@cite).
"""
module GITFans

# the necessary Julia packages
using Oscar

export git_fan

#############################################################################

# utility functions

## a <= b for two arrays `a`, `b` of bitsets,
## defined lexicographically
function less_or_equal_array_bitlist(a::Vector{BitSet}, b::Vector{BitSet})
    for i in 1:length(a)
        if a[i].bits < b[i].bits
            return true
        elseif a[i].bits > b[i].bits
            return false
        end
    end
    return true
end

##  Return a new bitset, consisting of the images of `bitset` under the
##  permutation `perm`.
bitlist_oper(bitset::BitSet, perm::Vector{Int}) = BitSet(perm[i] for i in bitset)

function bitlist_oper_tuple(bitset_tuple::Vector{BitSet}, perm_tuple::Vector{Vector{Int}})
    return [bitlist_oper(bitset_tuple[i], perm_tuple[i]) for i in 1:length(bitset_tuple)]
end

function find_smallest_orbit_element(elem, ggens, action::Function, comparator::Function, leq::Function)
    current_orbit = orbit(elem, ggens, action, comparator)
    sorted_orbit = sort(current_orbit; lt=leq)
    return sorted_orbit[1]
end

function rewrite_action_to_orbits(homs::Vector{GAPGroupHomomorphism})
    G = domain(homs[1])
    n = ngens(G)
    generators_new_perm = [Vector{Int}[] for x in 1:n]

    for hom in homs
      for j in 1:n
        img = Vector{Int}(image(hom, gen(G, j)))
        if length(img) == 0
          img = Int[1]
        end
        push!(generators_new_perm[j], img)
      end
    end

    return generators_new_perm
end


#############################################################################

# user functions

"""
    is_monomial_free(I::MPolyIdeal, vars_to_zero::Vector{Int} = Int[])

See Prop. 3.1 in [BKR20](@cite).
"""
function is_monomial_free(I::MPolyIdeal, vars_to_zero::Vector{Int} = Int[])
    Oscar.singular_assure(I)
    SingI = I.gens.S
    R = base_ring(SingI)
    Rgens = gens(R)
    nr_variables = length(Rgens)
    poly_list = [evaluate(f, vars_to_zero, fill(R(0), length(vars_to_zero))) for f in gens(SingI)]

    # perm is an n-cycle
    perm = collect(2:nr_variables)
    push!(perm, 1)

    Rgens_permuted = Rgens[perm]
    si = Singular.Ideal(R, poly_list)

    for i in 1:nr_variables
        if !(nr_variables in vars_to_zero)
            si = Singular.satstd( si, Singular.MaximalIdeal(R, 1))
            if Singular.ngens(si) == 1 && si[1] == R(1)
                return false
            elseif reduce(R(1), si) == 0
                return false
            end
        end
        vars_to_zero = perm[vars_to_zero]
        permuted_gens = [evaluate(j, Rgens_permuted) for j in gens(si)]
        si = Singular.Ideal(R, permuted_gens)
    end

    return true
end

"""
    orbit_cones(I::MPolyIdeal, Q::Matrix{Int}, G::PermGroup = symmetric_group(1))

Return orbit representatives of `G` on the set of those cones
whose defining rays are given by subsets `S` of the rows of `Q`,
such that the matrix `S` has full rank and such that `I` is
monomial-free (see [`is_monomial_free`](@ref)) w.r.t. the variables `x_i`
for which the `i`-th row of `Q` is not contained in `S`.
"""
function orbit_cones(I::MPolyIdeal, Q::Matrix{Int}, G::PermGroup = symmetric_group(1))
    nr_variables, projected_dimension = size(Q)

    collector_cones = Cone{QQFieldElem}[]

    # We need not consider sets of smaller size because of the rank condition.
    for k in projected_dimension:nr_variables
        for orb in orbits(gset(G, on_sets, subsets(nr_variables, k)))
            i = representative(orb)
            current_mat = Q[i,:]
            if rank(current_mat) == projected_dimension &&
               is_monomial_free(I, setdiff(1:nr_variables, i))
                cone = positive_hull(current_mat)
                # When comparing two cones it is fastest to solve a bunch of
                # LP's to figure out whether the rays of one are contained in
                # the other and vice versa. It would be even easier to just do
                # some matrix multiplication if one already has the facets and
                # rays. However to get the facets and rays is very hard
                # compared to the LP's. But in this case we are comparing cones
                # a lot, so the cost of the LP's outweigh the heavy lifting of
                # computing rays and facets.
                facets(cone)
                rays(cone)
                if ! any(j -> j == cone,
                      collector_cones)
                    push!(collector_cones, cone)
                end
            end
        end
    end

    return collector_cones
end
#T what if some projections lie in the same orbit?
#T later we expand the orbits, do we want to check this here?

@doc raw"""
    action_on_target(Q::Matrix{Int}, G::PermGroup)

Let `Q` be an $n \times m$ Q-matrix, and `G` be a permutation group
on $n$ points that describes an action on the rows of `Q`.
The function returns the group homomorphism $\rho$ from `G`
to its induced matrix action on the $m$-dimensional space over the Rationals.

# Examples

```jldoctest
julia> Q = [
        1  1   0   0   0
        1  0   1   1   0
        1  0   1   0   1
        1  0   0   1   1
        0  1   0   0  -1
        0  1   0  -1   0
        0  1  -1   0   0
        0  0   1   0   0
        0  0   0   1   0
        0  0   0   0   1 ];

julia> n = nrows(Q)
10

julia> perms_list = [[1,3,2,4,6,5,7,8,10,9], [5,7,1,6,9,2,8,4,10,3]];

julia> sym10 = symmetric_group(n);

julia> permgens = [sym10(x) for x in perms_list];

julia> G, _ = sub(sym10, permgens);

julia> Oscar.GITFans.action_on_target(Q, G)
Group homomorphism
  from permutation group of degree 10
  to matrix group of degree 5 over QQ

```

Note that $\rho(\pi)$, for $\pi \in$ `G`, satisfies
`Q`$^\pi$ = `Q` * $\rho(\pi^{-1})$, where `Q`$^\pi$ is the matrix obtained
from `Q` by permuting its rows by $\pi$.

Let $Q\{J\}$ denote the matrix whose rows are those rows of $Q$ indexed by
the subset $J$ of $\{ 1, 2, \ldots, n \}$.
The cone defined by $Q{J}$ is mapped to
$Q{J} \cdot \rho(\pi^{-1}) = (Q \cdot \rho(\pi^{-1}))\{J\} = Q^\pi\{J\}$,
which is equal to $Q\{J^{\pi^{-1}}\}$.
"""
function action_on_target(Q::Matrix{Int}, G::PermGroup)

    # For each permutation generator describing the action on the rows of Q,
    # compute the induced action on the column space,
    # by solving a linear equation system.
    m, n = size(Q)
    mat = matrix(QQ, Q)
    permgens = gens(G)
    matgens = typeof(mat)[]
    for ppi in permgens
      matimg = mat[Vector{Int}(ppi), 1:n]  # permute the rows with `ppi`
      push!(matgens, solve(mat, matimg; side = :right))
    end

    # Create the matrix group.
    matgroup = matrix_group(QQ, n, matgens)

    # Create the group homomorphism.
    return hom(G, matgroup, permgens, gens(matgroup))
end


"""
    orbit(point, generators, action::Function, compare_func::Function)

Return the orbit of the point `point` under the action of the group that is
generated by the elements in the array `generators`.
The function `action` defines the action (pt,g) -> pt^g on some set.
The function `compare_func` compares two elements of the set
and returns `true` if the two objects are considered as equal,
and `false` otherwise.
The elements of the returned array are in general not sorted.
"""
function orbit(point, generators, action::Function, compare_func::Function)
    orb = [point]
    # Note that in Julia (like in GAP),
    # the 'for' loop runs also over entries that get added
    # inside the loop.
    for b in orb
        for g in generators
            c = action(b, g)
            if ! any(i->compare_func(i, c), orb)
                push!(orb, c)
            end
        end
    end
    return orb
end

"""
    as_permutation(element, set, action::Function, compare::Function)

Return the permutation (as a group element in the symmetric group of degree
`length(set)`)
induced by the action of the group element `element` on the array `set`
via the function `action`: (pt,element) -> pt^element.
The equality of points is decided via the binary function `compare`.
"""
function as_permutation(element, set, action::Function, compare::Function)
    n = length(set)
    permlist = Vector{Int64}(undef, n)
    for i in 1:length(set)
        img = action(set[i], element)
        pos = findfirst(j -> compare(j, img), set)
        if pos === nothing
          error("the set is not invariant under the given action")
        end
        permlist[i] = pos
    end
    return perm(permlist)
end

"""
    matrix_action_on_cones(C::Cone{QQFieldElem}, M::QQMatrix)

Return the cone defined by transforming the rays and lineality conditions
of `C` by multiplication with `M`.

```jldoctest
julia> R = [0 1];

julia> L = [1 0];

julia> M = matrix(QQ, [0 1; -1 0])
[ 0   1]
[-1   0]

julia> C1 = positive_hull([0 1]);

julia> img = Oscar.GITFans.matrix_action_on_cones(C1, M);

julia> data = rays_modulo_lineality(img);

julia> data.rays_modulo_lineality
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [-1, 0]

julia> data.lineality_basis
0-element SubObjectIterator{RayVector{QQFieldElem}}

julia> C2 = positive_hull([0 1], [1 0]);

julia> img = Oscar.GITFans.matrix_action_on_cones(C2, M);

julia> data = rays_modulo_lineality(img);

julia> data.rays_modulo_lineality
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [-1, 0]

julia> data.lineality_basis
1-element SubObjectIterator{RayVector{QQFieldElem}}:
 [0, 1]
```
"""
function matrix_action_on_cones(C::Cone{QQFieldElem}, M::QQMatrix)
    data = rays_modulo_lineality(C)
    rayscone = matrix(QQ, data.rays_modulo_lineality)
    linescone = matrix(QQ, data.lineality_basis)
    return positive_hull(rayscone * M, linescone * M)
end

function orbit_cone_orbits(cones::Vector{Cone{T}}, ghom::GAPGroupHomomorphism) where T
    matgens = [image(ghom, g).elm for g in gens(domain(ghom))]
    act = matrix_action_on_cones
   
    result = Vector{Vector{Cone{T}}}([])
    for c in cones
        # When comparing two cones it is fastest to solve a bunch of LP's to
        # figure out whether the rays of one are contained in the other and
        # vice versa. It would be even easier to just do some matrix
        # multiplication if one already has the facets and rays. However to get
        # the facets and rays is very hard compared to the LP's. But in this
        # case we are comparing cones a lot, so the cost of the LP's outweigh
        # the heavy lifting of computing rays and facets.
        rays(c)
        facets(c)
        if all(o -> all(x -> c != x, o), result)
            push!(result, orbit(c, matgens, act, Base.:(==)))
        end
    end

    return result
end

function action_on_orbit_cone_orbits(orbs::Vector{Vector{Cone{T}}}, ghom::GAPGroupHomomorphism) where T
#T Prescribing `GAPGroupHomomorphism` is quite technical;
#T do we have a better way to request a "map that is a group homomorphism"?
    G = domain(ghom)

    if ngens(G) == 0
        m = id_hom(G)
        return GAPGroupHomomorphism[m for x in orbs]
    end

    Ggens = gens(G)
    matgens = [image(ghom, g).elm for g in Ggens]
    act = matrix_action_on_cones
    
    res = GAPGroupHomomorphism[]
    for orb in orbs
        list = [as_permutation(gen, orb, act, Base.:(==)) for gen in matgens]
        push!(res, hom(G, sub(list...)[1], Ggens, list))
    end

    return res
end



"""
    compute_bit_list(orbs::Vector{Vector{Cone{T}}}, point:AbstractVector{T}) where T

Return an array of bitsets,
containing at position i the bitset for the i-th entry of `orbs`,
that is, `true` at the positions of those cones in `orbs[i]`
that contain the point `point`.
"""
function compute_bit_list(orbs::Vector{Vector{Cone{T}}}, point::AbstractVector{T}) where T
    bitset_list = [BitSet() for orb in orbs]
    for (i, orb) in enumerate(orbs)
        for (j, cone) in enumerate(orb)
            if point in cone
                push!(bitset_list[i], j)
            end
        end
    end
    return bitset_list
end

"""
    cones_from_bitlist(cone_list::Vector{Vector{Cone{T}}}, bit_list_tuple::Vector{BitSet})

Return the array of all cones at `true` positions in the bitsets.
"""
function cones_from_bitlist(cone_list::Vector{Vector{Cone{T}}}, bit_list_tuple::Vector{BitSet}) where T
  return_list = Cone{T}[]
    for i in 1:length(cone_list)
        for j in bit_list_tuple[i]
            push!(return_list, cone_list[i][j])
        end
    end
    return return_list
end


hash_to_cone(orbit_list::Vector{Vector{Cone{T}}}, hash::Vector{BitSet}) where T = intersect(cones_from_bitlist(orbit_list, hash)...)


"""
    get_neighbor_hash(orbits::Vector{Vector{Cone{T}}}, facet_point::Vector{QQFieldElem}, inner_normal_vector::Vector{QQFieldElem})

Return the list of bitsets describing the cones adjacent to the point
`facet_point` in the direction of `inner_normal_vector`.
"""
function get_neighbor_hash(orbits::Vector{Vector{Cone{T}}}, facet_point::Vector{QQFieldElem}, inner_normal_vector::Vector{QQFieldElem}) where T
    lambda = 1024
    facet_point_bl = compute_bit_list(orbits, facet_point)
    while true
        current_point = lambda * facet_point - inner_normal_vector
## FIXME: compute only necessary part of BL
        current_bl = compute_bit_list(orbits, current_point)
        if all(i->issubset(current_bl[i], facet_point_bl[i]), 1:length(facet_point_bl))
            return current_bl
        end
        lambda *= 2
    end
end


#T is not needed anymore, thanks to `relative_interior_point`
function get_interior_point_weighted(cone::Cone)
    weights = rand(-100:100, n_rays(cone))
    return sum(rays(cone) .* weights)
end


"""
    fan_traversal(orbit_list::Vector{Vector{Cone{T}}}, q_cone::Cone{T}, perm_actions::Vector{GAPGroupHomomorphism}) where T

Return the pair `(hash_list, edges)` where `hash_list` is an array that
encodes orbit representatives of the maximal cones of the GIT fan described
by `orbit_list`, `q_cone`, and `perm_actions`,
and `edges` encodes the `Set` of edges of the incidence graph of the orbits.
"""
function fan_traversal(orbit_list::Vector{Vector{Cone{T}}}, q_cone::Cone{T}, perm_actions::Vector{GAPGroupHomomorphism}) where T
    # the induced actions on each of the orbits
    generators_new_perm = rewrite_action_to_orbits(perm_actions)

    q_cone_int_point = relative_interior_point(q_cone)
 
    start_hash = compute_bit_list(orbit_list, q_cone_int_point)
    orbit_start_hash_smallest = find_smallest_orbit_element(
        start_hash, generators_new_perm, bitlist_oper_tuple, ==,
        less_or_equal_array_bitlist)
    hash_list = [orbit_start_hash_smallest]
    edges = Set(Vector{Int}[])

    for (current_pos, current_hash) in enumerate(hash_list)

        # note that we run also over elements added inside the loop
        current_cone_list = cones_from_bitlist(orbit_list, current_hash)
        intersected_cone = intersect(current_cone_list)
        ic_facets = -linear_inequality_matrix(facets(intersected_cone))
        facet_points = [T[relative_interior_point(f)...] for f in facets(Cone, intersected_cone)]

        neighbor_hashes = []
        for i in 1:length(facet_points)
          !any(f->facet_points[i] in f, facets(Cone, q_cone)) || continue
          push!(neighbor_hashes, get_neighbor_hash(orbit_list, facet_points[i], Vector{T}([ic_facets[i:i, :]...])))
        end

        neighbor_hashes = [find_smallest_orbit_element(i, generators_new_perm, bitlist_oper_tuple, ==, less_or_equal_array_bitlist) for i in neighbor_hashes]
        for i in neighbor_hashes
            if i in hash_list
                # perhaps we have found a new incidence
                push!(edges, sort!([findfirst(x->x == i, hash_list), current_pos]))
            else
                # new representative found
                push!(hash_list, i)
                # new incidence
                push!(edges, [current_pos, length(hash_list)])
            end
        end
    end

    return (hash_list, edges)
end


orbits_of_maximal_GIT_cones(orbit_list::Vector{Vector{Cone{T}}}, hash_list::Vector{Vector{BitSet}}, matrix_action::GAPGroupHomomorphism) where T =
    orbit_cone_orbits([hash_to_cone(orbit_list, x) for x in hash_list], matrix_action)


function hashes_to_polyhedral_fan(orbit_list::Vector{Vector{Cone{T}}}, hash_list::Vector{Vector{BitSet}}, ghom::GAPGroupHomomorphism) where T
    # translate the descriptions of the orbit repres. of maximal cones
    # to cone objects
    result_cones = [hash_to_cone(orbit_list, x) for x in hash_list]

    # expand their orbits
    expanded = orbit_cone_orbits(result_cones, ghom)
    maxcones = vcat(expanded...)

    # the defining rays for all maximal cones
    rays_maxcones = [[Vector{QQFieldElem}(x) for x in rays(cone)]
                      for cone in maxcones]

    # the defining rays for the orbit representatives of maximal cones
    rays_result_cones = [[Vector{QQFieldElem}(x) for x in rays(cone)]
                      for cone in result_cones]

    # the set of rays
    allrays = sort!(unique(vcat(rays_maxcones...)))

    # the indices of rays that belong to each maximal cone (0-based)
    index_maxcones = [sort([findfirst(x -> x == v, allrays)-1
                            for v in rays])
                      for rays in rays_maxcones]

    # the indices of rays that belong to each repres. cone
    index_result_cones = [sort([findfirst(x -> x == v, allrays)
                            for v in rays])
                      for rays in rays_result_cones]

    # The fan object without information about the symmetry group
    # could now be constructed as follows.
    # # return Polymake.fan.PolyhedralFan(INPUT_RAYS = hcat(allrays...)',
    # #                                   INPUT_CONES = index_maxcones)

    # In order to create a fan object with symmetry information, prepare
    # - the matrix of rays,
    # - the permutation generators of the symmetry group G acting on the rays,
    #   computed from the matrix action,
    # - representatives of the G-orbits on the maximal cones,
    #   each given by its defining rays, via the row indices (zero based)
    #   in the matrix of rays.
    matgens = [image(ghom, x) for x in gens(domain(ghom))]
    mats_transp = [transpose(x.elm) for x in matgens]

    # Note that the matrices have been transposed because we use
    # matrix times vector multiplication.
    actfun = function(ray, mat)
      ray = mat * ray
      # normalize the image vector
      for i in 1:length(ray)
        if ray[i] != 0
          return ray // ray[i]
        end
      end
      error("a zero vector should not appear")
    end

    rewr = [as_permutation(x, allrays, actfun, Base.:(==)) for x in mats_transp]

    pf = polyhedral_fan_from_rays_action(allrays, IncidenceMatrix(index_result_cones), rewr)

    return pf
end


"""
    git_fan(a::MPolyIdeal, Q::Matrix{Int}, G::PermGroup)

Return the object that represents the polyhedral fan given by
the ideal `a`, the grading matrix `Q`, and the symmetry group `G`.
"""
function git_fan(a::MPolyIdeal, Q::Matrix{Int}, G::PermGroup)
    collector_cones = orbit_cones(a, Q, G)
    matrix_action = action_on_target(Q, G)
    orbit_list = orbit_cone_orbits(collector_cones, matrix_action)
    perm_actions = action_on_orbit_cone_orbits(orbit_list, matrix_action)
    q_cone = positive_hull(Q)

    (hash_list, edges) = fan_traversal(orbit_list, q_cone, perm_actions)
    return hashes_to_polyhedral_fan(orbit_list, hash_list, matrix_action)
end


"""
    edges_intersection_graph(maxcones::Vector{Cone{T}}, inter_dim::Int) where T

Return the array of those pairs `[i, j]` such that `i < j` and such that
the intersection of the cones `maxcones[i]` and `maxcones[j]` has dimension
`inter_dim`.

If `maxcones` is an array of cones of dimension `inter_dim + 1`
then the returned array describes the edges of the intersection graph.
"""
function edges_intersection_graph(maxcones::Vector{Cone{T}}, inter_dim::Int) where T
    edges = Vector{Int}[]
    for j in 1:length(maxcones)
        for i in 1:(j-1)
            if dim(intersect(
                   maxcones[i], maxcones[j])) == inter_dim
                push!(edges, [i, j])
            end
        end
    end

    return edges
end

end
