function PermGroup_to_polymake_array(G::PermGroup)
   generators = gens(G)
   d = degree(G)
   result = Polymake.Array{Polymake.Array{Polymake.to_cxx_type(Int)}}(length(generators))
   i = 1
   for g in generators
      array = Polymake.Array{Polymake.to_cxx_type(Int)}(d)
      for j in 1:d
         array[j] = g(j)-1
      end
      result[i] = array
      i = i+1
   end
   return result
end

function _pm_group_to_oscar_group(G)
    gens = _pm_arr_arr_to_group_generators(G.GENERATORS)
    return _gens_to_group(gens)
end

function _gens_to_group(gens::Vector{PermGroupElem})
    S = parent(gens[1])
    return sub(S, gens)[1]
end

function _pm_arr_arr_to_group_generators(M)
    n=length(M[1])
    S=symmetric_group(n)
    perm_bucket = Vector{Oscar.BasicGAPGroupElem{PermGroup}}()
    for g in M
        push!(perm_bucket,perm(S, Polymake.to_one_based_indexing(g)))
    end
    return perm_bucket
end

#TODO: rename this. It is the automorphism group of the vertex facet incidence
#also applies to tests and exports
function vf_group(P::Polyhedron)
    if pm_object(P).BOUNDED
        _pm_group_to_oscar_group(Polymake.group.automorphism_group(pm_object(P).VERTICES_IN_FACETS).PERMUTATION_ACTION)
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end
end



@doc Markdown.doc"""
    combinatorial_symmetries(P::Polyhedron)

Compute the combinatorial symmetries (i.e., automorphisms of the face lattice)
of a given polytope $P$.  The result is given as permutations of the vertices
(or rather vertex indices) of $P$. This group contains the `linear_symmetries`
as a subgroup.

# Example
The quadrangle one obtains from moving one vertex of the square out along the
diagonal has eight combinatorial symmetries, but only two linear symmetries.
```jldoctest
julia> quad = convex_hull([0 0; 1 0; 2 2; 0 1])
A polyhedron in ambient dimension 2

julia> G = combinatorial_symmetries(quad)
Group([ (1,2)(3,4), (1,3) ])

julia> length(elements(G))
8

julia> G = linear_symmetries(quad)
Group([ (2,4) ])

julia> length(elements(G))
2
```
"""
function combinatorial_symmetries(P::Polyhedron)
    result = automorphism_group(P; type=:combinatorial)
    return result[:vertex_action]
end

@doc Markdown.doc"""
    linear_symmetries(P::Polyhedron)

Get the group of linear symmetries on the vertices of a polyhedron. These are
morphisms of the form $x\mapsto Ax+b$,with $A$ a matrix and $b$ a vector, that
preserve the polyhedron $P$. The result is given as permutations of the
vertices (or rather vertex indices) of $P$.

# Examples
The 3-dimensional cube has 48 linear symmetries.
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> G = linear_symmetries(c)
Group([ (3,5)(4,6), (2,3)(6,7), (1,2)(3,4)(5,6)(7,8) ])

julia> length(elements(G))
48
```

The quadrangle one obtains from moving one vertex of the square out along the
diagonal has two linear symmetries.
```jldoctest
julia> quad = convex_hull([0 0; 1 0; 2 2; 0 1])
A polyhedron in ambient dimension 2

julia> G = linear_symmetries(quad)
Group([ (2,4) ])

julia> length(elements(G))
2
```
"""
function linear_symmetries(P::Polyhedron)
    result = automorphism_group(P; type=:linear)
    return result[:vertex_action]
end


@doc Markdown.doc"""
    automorphism_group_generators(P::Polyhedron; type = :combinatorial)

Compute generators of the group of automorphisms of a polyhedron.

The optional parameter `type` takes two values:
- `:combinatorial` (default) -- Return the combinatorial symmetries, the
    automorphisms of the face lattice.
- `:linear` -- Return the linear automorphisms.

The return value is a `Dict{Symbol, Vector{PermGroupElem}}` with two entries,
one for the key `:vertex_action` containing the generators for the action
permuting the vertices, and `:facet_action` for the facets.


# Examples
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> automorphism_group_generators(c)
Dict{Symbol, Vector{PermGroupElem}} with 2 entries:
  :vertex_action => [(3,5)(4,6), (2,3)(6,7), (1,2)(3,4)(5,6)(7,8)]
  :facet_action  => [(3,5)(4,6), (1,3)(2,4), (1,2)]
```
"""
function automorphism_group_generators(P::Polyhedron; type = :combinatorial)
    if type == :combinatorial
        gens = Polymake.graph.automorphisms(vertex_indices(facets(P)))
        facet_action = _pm_arr_arr_to_group_generators([first(g) for g in gens])
        vertex_action = _pm_arr_arr_to_group_generators([last(g) for g in gens])
    elseif type == :linear
        gp = Polymake.polytope.linear_symmetries(pm_polytope(P))
        facet_action = _pm_arr_arr_to_group_generators(gp.FACETS_ACTION.GENERATORS)
        vertex_action = _pm_arr_arr_to_group_generators(gp.VERTICES_ACTION.GENERATORS)
    else
        throw(ArgumentError("Action type $(type) not supported."))
    end
    return Dict{Symbol, Vector{PermGroupElem}}(:vertex_action => vertex_action,
            :facet_action => facet_action)
end



@doc Markdown.doc"""
    automorphism_group(P::Polyhedron; type = :combinatorial)

Compute the group of automorphisms of a polyhedron.

The optional parameter `type` takes two values:
- `:combinatorial` (default) -- Return the combinatorial symmetries, the
    automorphisms of the face lattice.
- `:linear` -- Return the linear automorphisms.

The return value is a `Dict{Symbol, PermGroup}` with two entries, one for the
key `:vertex_action` containing the group permuting the vertices, and
`:facet_action` for the facets.

# Examples
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> automorphism_group(c)
Dict{Symbol, PermGroup} with 2 entries:
  :vertex_action => Group([ (3,5)(4,6), (2,3)(6,7), (1,2)(3,4)(5,6)(7,8) ])
  :facet_action  => Group([ (3,5)(4,6), (1,3)(2,4), (1,2) ])
```
"""
function automorphism_group(P::Polyhedron; type = :combinatorial)
    result = automorphism_group_generators(P; type = type)
    va = _gens_to_group(result[:vertex_action])
    fa = _gens_to_group(result[:facet_action])
    return Dict{Symbol, PermGroup}(:vertex_action => va,
            :facet_action => fa)
end
