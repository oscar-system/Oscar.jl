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
function _gens_to_group(gens::Dict{Symbol, Vector{PermGroupElem}})
    return Dict{Symbol, PermGroup}([k => _gens_to_group(v) for (k,v) in gens])
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


@deprecate vf_group(P::Polyhedron) automorphism_group(P; action = :on_facets)


@doc Markdown.doc"""
    combinatorial_symmetries(P::Polyhedron)

Compute the combinatorial symmetries (i.e., automorphisms of the face lattice)
of a given polytope $P$.  The result is given as permutations of the vertices
(or rather vertex indices) of $P$. This group contains the `linear_symmetries`
as a subgroup.

# Examples
The quadrangle one obtains from moving one vertex of the square out along the
diagonal has eight combinatorial symmetries, but only two linear symmetries.
```jldoctest
julia> quad = convex_hull([0 0; 1 0; 2 2; 0 1])
A polyhedron in ambient dimension 2

julia> G = combinatorial_symmetries(quad)
Group([ (2,4), (1,2)(3,4) ])

julia> order(G)
8

julia> G = linear_symmetries(quad)
Group([ (2,4) ])

julia> order(G)
2
```
"""
function combinatorial_symmetries(P::Polyhedron)
    return automorphism_group(P; type=:combinatorial, action=:on_vertices)
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

julia> order(G)
48
```

The quadrangle one obtains from moving one vertex of the square out along the
diagonal has two linear symmetries.
```jldoctest
julia> quad = convex_hull([0 0; 1 0; 2 2; 0 1])
A polyhedron in ambient dimension 2

julia> G = linear_symmetries(quad)
Group([ (2,4) ])

julia> order(G)
2
```
"""
function linear_symmetries(P::Polyhedron)
    return automorphism_group(P; type=:linear, action=:on_vertices)
end


@doc Markdown.doc"""
    automorphism_group_generators(P::Polyhedron; type = :combinatorial, action = :all)

Compute generators of the group of automorphisms of a polyhedron.

The optional parameter `type` takes two values:
- `:combinatorial` (default) -- Return the combinatorial symmetries, the
    automorphisms of the face lattice.
- `:linear` -- Return the linear automorphisms.

The optional parameter `action` takes three values:
- `:all` (default) -- Return the generators of the permutation action on both
    vertices and facets as a Dict{Symbol, Vector{PermGroupElem}}.
- `:on_vertices` -- Only return generators of the permutation action on the
    vertices.
- `:on_facets` -- Only return generators of the permutation action on the
    facets.

The return value is a `Dict{Symbol, Vector{PermGroupElem}}` with two entries,
one for the key `:on_vertices` containing the generators for the action
permuting the vertices, and `:on_facets` for the facets.


# Examples
Compute the automorphisms of the 3dim cube:
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> automorphism_group_generators(c)
Dict{Symbol, Vector{PermGroupElem}} with 2 entries:
  :on_vertices => [(3,5)(4,6), (2,3)(6,7), (1,2)(3,4)(5,6)(7,8)]
  :on_facets   => [(3,5)(4,6), (1,3)(2,4), (1,2)]

julia> automorphism_group_generators(c; action = :on_vertices)
3-element Vector{PermGroupElem}:
 (3,5)(4,6)
 (2,3)(6,7)
 (1,2)(3,4)(5,6)(7,8)

julia> automorphism_group_generators(c; action = :on_facets)
3-element Vector{PermGroupElem}:
 (3,5)(4,6)
 (1,3)(2,4)
 (1,2)
```

Compute the automorphisms of a non-quadratic quadrangle. Since it has less
symmetry than the square, it has less linear symmetries.
```jldoctest
julia> quad = convex_hull([0 0; 1 0; 2 2; 0 1])
A polyhedron in ambient dimension 2

julia> automorphism_group_generators(quad)
Dict{Symbol, Vector{PermGroupElem}} with 2 entries:
  :on_vertices => [(2,4), (1,2)(3,4)]
  :on_facets   => [(1,2)(3,4), (1,3)]

julia> automorphism_group_generators(quad; type = :combinatorial)
Dict{Symbol, Vector{PermGroupElem}} with 2 entries:
  :on_vertices => [(2,4), (1,2)(3,4)]
  :on_facets   => [(1,2)(3,4), (1,3)]

julia> automorphism_group_generators(quad; type = :linear)
Dict{Symbol, Vector{PermGroupElem}} with 2 entries:
  :on_vertices => [(2,4)]
  :on_facets   => [(1,2)(3,4)]
```
"""
function automorphism_group_generators(P::Polyhedron; type = :combinatorial, action = :all)
    is_bounded(P) || throw(ArgumentError("Automorphism groups not supported for unbounded polyhedra."))
    if type == :combinatorial
        IM = vertex_indices(facets(P))
        if action == :all
            result = automorphism_group_generators(IM; action = action)
            return Dict{Symbol, Vector{PermGroupElem}}(
                    :on_vertices => result[:on_cols],
                    :on_facets => result[:on_rows])
        elseif action == :on_vertices
            return automorphism_group_generators(IM; action = :on_cols)
        elseif action == :on_facets
            return automorphism_group_generators(IM; action = :on_rows)
        else
            throw(ArgumentError("Action $(action) not supported."))
        end
    elseif type == :linear
        return _linear_symmetries_generators(P; action = action)
    else
        throw(ArgumentError("Action type $(type) not supported."))
    end
    return Dict{Symbol, Vector{PermGroupElem}}(:on_vertices => vertex_action,
            :on_facets => facet_action)
end


@doc Markdown.doc"""
    automorphism_group_generators(IM::IncidenceMatrix; action = :all)

Compute the generators of the group of automorphisms of an IncidenceMatrix. 

The optional parameter `action` takes three values:
- `:all` (default) -- Return the generators of the permutation action on both
    columns and rows as a Dict{Symbol, Vector{PermGroupElem}}.
- `:on_cols` -- Only return generators of the permutation action on the
    columns.
- `:on_rows` -- Only return generators of the permutation action on the
    rows.

# Examples
Compute the automorphisms of the incidence matrix of the 3dim cube:
```jldoctest
julia> c = cube(3)
A polyhedron in ambient dimension 3

julia> IM = vertex_indices(facets(c))
6×8 IncidenceMatrix
[1, 3, 5, 7]
[2, 4, 6, 8]
[1, 2, 5, 6]
[3, 4, 7, 8]
[1, 2, 3, 4]
[5, 6, 7, 8]


julia> automorphism_group_generators(IM)
Dict{Symbol, Vector{PermGroupElem}} with 2 entries:
  :on_cols => [(3,5)(4,6), (2,3)(6,7), (1,2)(3,4)(5,6)(7,8)]
  :on_rows => [(3,5)(4,6), (1,3)(2,4), (1,2)]

julia> automorphism_group_generators(IM; action = :on_rows)
3-element Vector{PermGroupElem}:
 (3,5)(4,6)
 (1,3)(2,4)
 (1,2)

julia> automorphism_group_generators(IM; action = :on_cols)
3-element Vector{PermGroupElem}:
 (3,5)(4,6)
 (2,3)(6,7)
 (1,2)(3,4)(5,6)(7,8)
```
"""
function automorphism_group_generators(IM::IncidenceMatrix; action = :all)
    gens = Polymake.graph.automorphisms(IM)
    rows_action = _pm_arr_arr_to_group_generators([first(g) for g in gens])
    if action == :on_rows
        return rows_action
    end
    cols_action = _pm_arr_arr_to_group_generators([last(g) for g in gens])
    if action == :on_cols
        return cols_action
    elseif action == :all
        return Dict{Symbol, Vector{PermGroupElem}}(:on_rows => rows_action,
                :on_cols => cols_action)
    else
        throw(ArgumentError("Action $(action) not supported."))
    end
end
    


function _linear_symmetries_generators(P::Polyhedron; action = :all)
    if is_bounded(P)
        gp = Polymake.polytope.linear_symmetries(vertices(P))
        vgens = gp.PERMUTATION_ACTION.GENERATORS
        if action == :on_vertices
            return _pm_arr_arr_to_group_generators(vgens)
        end
        hgens = Polymake.group.induced_permutations(vgens, vertex_indices(facets(P)))
        if action == :on_facets
            return _pm_arr_arr_to_group_generators(hgens)
        end
        return Dict{Symbol, Vector{PermGroupElem}}(
                :on_vertices => _pm_arr_arr_to_group_generators(vgens),
                :on_facets => _pm_arr_arr_to_group_generators(hgens)
            )
    else
        throw(ArgumentError("Linear symmetries supported for bounded polyhedra only"))
    end
end





@doc Markdown.doc"""
    automorphism_group(P::Polyhedron; type = :combinatorial, action = :all)

Compute the group of automorphisms of a polyhedron. The parameters and return
values are the same as for [`automorphism_group_generators(P::Polyhedron; type
= :combinatorial, action = :all)`](@ref) except that groups are returned
instead of generators of groups.
"""
function automorphism_group(P::Polyhedron; type = :combinatorial, action = :all)
    result = automorphism_group_generators(P; type = type, action = action)
    return _gens_to_group(result)
end


@doc Markdown.doc"""
    automorphism_group(IM::IncidenceMatrix; action = :all)

Compute the group of automorphisms of an IncidenceMatrix. The parameters and
return values are the same as for
[`automorphism_group_generators(IM::IncidenceMatrix; action = :all)`](@ref)
except that groups are returned instead of generators of groups.
"""
function automorphism_group(IM::IncidenceMatrix; action = :all)
    result = automorphism_group_generators(P; action = action)
    return _gens_to_group(result)
end
