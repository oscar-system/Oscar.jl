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
        push!(perm_bucket,perm(S,[a+1 for a in g]))
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
function combinatorial_symmetries(P::Polyhedron)
    if pm_object(P).BOUNDED
        pm_group_to_oscar_group(Polymake.polytope.combinatorial_symmetries(pm_object(P)))
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end
end

"""
    linear_symmetries(P::Polyhedron)

Get the group of linear symmetries on the vertices of a polyhedron.
"""
function linear_symmetries(P::Polyhedron)
    if pm_object(P).BOUNDED
        pm_group_to_oscar_group(Polymake.polytope.linear_symmetries(pm_object(P).VERTICES).PERMUTATION_ACTION)
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end
end
