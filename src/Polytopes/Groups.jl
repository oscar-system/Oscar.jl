
function pm_group_to_oscar_group(G)
    pm_arr_arr_to_group(G.GENERATORS)
end

function pm_arr_arr_to_group(M)
    n=length(M[1])
    S=symmetric_group(n)
    perm_bucket = Array{Oscar.BasicGAPGroupElem{PermGroup},1}()
    for g in M
        push!(perm_bucket,perm(S,[a+1 for a in g]))
    end
    H=sub(S,perm_bucket)[1]
    return(H)
end

#TODO: rename this. It is the automorphism group of the vertex facet incidence
function vf_group(P::Polyhedron)
    if P.pm_polytope.BOUNDED
        pm_group_to_oscar_group(Polymake.group.automorphism_group(P.pm_polytope.VERTICES_IN_FACETS).PERMUTATION_ACTION)
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end
end
function combinatorial_symmetries(P::Polyhedron)
    if P.pm_polytope.BOUNDED
        pm_group_to_oscar_group(Polymake.polytope.combinatorial_symmetries(P.pm_polytope))
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end
end

function linear_symmetries(P::Polyhedron)
    if P.pm_polytope.BOUNDED
        pm_group_to_oscar_group(Polymake.polytope.linear_symmetries(P.pm_polytope.VERTICES).PERMUTATION_ACTION)
    else
        throw(ArgumentError("Symmetry groups currently supported for bounded polyhedra only"))
    end  
end
