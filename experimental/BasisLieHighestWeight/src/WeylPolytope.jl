function orbit_weylgroup(type::String, rank::Int, lie_algebra::GAP.Obj, weight_vector::Vector{Int})
    """
    operates weyl-group of type type and rank rank on vector weight_vector and returns list of vectors in orbit
    input and output weights in terms of w_i
    """
    # initialization
    weyl_group = GAP.Globals.WeylGroup(GAP.Globals.RootSystem(lie_algebra))
    orbit_iterator = GAP.Globals.WeylOrbitIterator(weyl_group, GAP.Obj(weight_vector))
    vertices = []
    
    # operate with the weylgroup on weight_vector
    GAP.Globals.IsDoneIterator(orbit_iterator)
    while !(GAP.Globals.IsDoneIterator(orbit_iterator))
        w = GAP.Globals.NextIterator(orbit_iterator) 
        push!(vertices, Vector{Int}(w))
    end

    # return
    vertices = convert(Vector{Vector{Int}}, vertices)
    return vertices
end

function get_points_polytope(polytope)
    """
    returns all points (interior and vertices) of a polytope in regular (i.e. not homogenoues coordinates).
    """
    interior_points = convert(Matrix{Int64}, polytope.INTERIOR_LATTICE_POINTS)
    vertices_points = convert(Matrix{Int64}, polytope.VERTICES)
    points = [interior_points; vertices_points][:, 2:end]
    return points
end


function convert_lattice_points_to_monomials(ZZx, lattice_points_weightspace)
    return [finish(push_term!(MPolyBuildCtx(ZZx), ZZ(1), convert(Vector{Int}, convert(Vector{Int64}, lattice_point)))) 
              for lattice_point in lattice_points_weightspace]
end

function get_lattice_points_of_weightspace(weights, weight, type)
    """
    calculates all lattice points in a given weightspace for a lie algebra of type type
    input:
    weights: the operator weights in eps_i
    weight: lambda - mu

    output: all lattice points with weight weight
    """
    if type in ["A", "G"]
        return get_lattice_points_of_weightspace_A_G_n(weights, weight)
    else
        return get_lattice_points_of_weightspace_Xn(weights, weight)
    end
end

function get_lattice_points_of_weightspace_A_G_n(weights, weight)
    """
    calculates all monomials in a given weightspace for lie algebras that have type A or G
    input:
    weights: the operator weights in eps_i
    weight: lambda - mu

    output: all monomials with weight weight
    
    works by calculating all integer solutions to the following linear program:
    [ 1     |              |    ]       [   x   ]      
    [ 1 weights[1]... weights[k]]   *   [   |   ]   =   weight 
    [...    |              |    ]       [  res  ] 
    [ 1     |              |    ]       [   |   ]
    where res[i] >= 0 for all i

    example:
    weights = [[1, 0, 2], [-1, 1, 1], [0, -1, 0]] (i.e. a_1 = eps_1 - eps_2, a_2 = eps_2 - eps_3, a_12 = eps_1 - eps_3)
    weight = [2, 1, 0]
    -> poly = polytope.Polytope(INEQUALITIES=[0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1], 
                                  EQUATIONS=[-2 1 1 0 2; -1 1 -1 1 1; 0 1 0 -1 0])
    => returns [[1 0 0], [1 1 0]]
    """
    # build linear (in-)equalities
    weights = [reshape(w, 1, :) for w in weights]
    n = length(weights)
    ineq = zeros(Int64, n, n+2)
    for i in 1:n
        ineq[i, 2+i] = 1
    end
    equ = cat([-i for i in vec(weight)], [1 for i=1:length(weight)], dims = (2,2))
    equ = cat(equ, [transpose(w) for w in weights] ..., dims = (2,2))

    # find integer solutions of linear (in-)equation as lattice points of polytope
    poly = polytope.Polytope(INEQUALITIES=ineq, EQUATIONS=equ)

    # convert lattice-points to Oscar monomials
    lattice_points_weightspace = lattice_points(Polyhedron(poly))
    lattice_points_weightspace = [lattice_point[2:end] for lattice_point in lattice_points_weightspace]
    return lattice_points_weightspace

end

function get_lattice_points_of_weightspace_Xn(weights, weight)
    """
    calculates all lattice points in a given weightspace for lie algebras that don't have type A
    input:
    weights: the operator weights in eps_i
    weight: lambda - mu

    output: all lattice points with weight weight
    
    works by calculating all integer solutions to the following linear program:
    [   |              |    ]       [   x   ]      
    [weights[1]...weights[k]]   *   [   |   ]   =   weight 
    [   |              |    ]       [  res  ] 
    [   |              |    ]       [   |   ]
    where res[i] >= 0 for all i

    example:
    weights = [[1, 0, 2], [-1, 1, 1], [0, -1, 0]] (i.e. a_1 = eps_1 - eps_2, a_2 = eps_2 - eps_3, a_12 = eps_1 - eps_3)
    weight = [2, 1, 0]
    -> poly = polytope.Polytope(INEQUALITIES=[0 1 0 0; 0 0 1 0; 0 0 0 1], EQUATIONS=[-2 1 0 2; -1 -1 1 1; 0 0 -1 0])
    => returns 
    """
    # build linear (in-)equalities
    weights = [reshape(w, 1, :) for w in weights]
    n = length(weights)
    ineq = zeros(Int64, n, n+1)
    for i in 1:n
        ineq[i, 1+i] = 1
    end
    equ = [-i for i in vec(weight)]
    equ = cat(equ, [transpose(w) for w in weights] ..., dims = (2,2))

    # find integer solutions of linear (in-)equation as lattice points of polytope
    poly = polytope.Polytope(INEQUALITIES=ineq, EQUATIONS=equ)
    
    # convert lattice-points to Oscar monomials
    lattice_points_weightspace = lattice_points(Polyhedron(poly))
    return lattice_points_weightspace
end
