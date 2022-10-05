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

function pm_group_to_oscar_group(G)
    pm_arr_arr_to_group(G.GENERATORS)
end

function pm_arr_arr_to_group(M)
    n=length(M[1])
    S=symmetric_group(n)
    perm_bucket = Vector{Oscar.BasicGAPGroupElem{PermGroup}}()
    for g in M
        push!(perm_bucket,perm(S, Polymake.to_one_based_indexing(g)))
    end
    H=sub(S,perm_bucket)[1]
    return H
end

#TODO: rename this. It is the automorphism group of the vertex facet incidence
#also applies to tests and exports
function vf_group(P::Polyhedron)
    if pm_object(P).BOUNDED
        pm_group_to_oscar_group(Polymake.group.automorphism_group(pm_object(P).VERTICES_IN_FACETS).PERMUTATION_ACTION)
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
    if pm_object(P).BOUNDED
        pm_group_to_oscar_group(Polymake.polytope.combinatorial_symmetries(pm_object(P)))
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end
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
    if pm_object(P).BOUNDED
        pm_group_to_oscar_group(Polymake.polytope.linear_symmetries(pm_object(P).VERTICES).PERMUTATION_ACTION)
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end
end
