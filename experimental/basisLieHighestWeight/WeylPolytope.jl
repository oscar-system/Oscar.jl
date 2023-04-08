#using Gapjm
using Oscar

include("./RootConversion.jl")

fromGap = Oscar.GAP.gap_to_julia

#########################
#     weyl-polytope     #
#########################

function weylpolytope(t::String, n::Int, hw::Vector{Int})
    """
    returns weyl-polytope in homogeneous coordinates, i.e. convex hull of orbit of weyl-group of type t,n on highest weight vector hw 
    """
    vertices = orbit_weylgroup(t, n, hw)
    #println("vertices: ", vertices)
    println(vertices)
    vertices_hom = [ones(Int64, size(vertices)[1]) vertices]# homogeneous coordinates
    println(vertices_hom)
    #println("vertices_hom: ", vertices_hom)
    #weylpoly = convex_hull(vertices) # normalized to first coordinate 1? See documentation in PolyMake
    weylpoly = polytope.Polytope(POINTS=vertices_hom)
    #println("weylpoly: ", weylpoly)
    #println("contains?:", contains(weylpoly, [0, 0]))
    return weylpoly
end

function orbit_weylgroup(t::String, n::Int, hw)
    """
    operates weyl-group of type t,n on highest weight vector hw and returns list of vector
    input in terms of w_i, output in eps_i (fix!!)
    """
    # also possible with polymake orbit_polytope(Vector input_point, Group g), root_system(String type) to save the equations, constant summand missing
    # initialization
    L, CH = lieAlgebra(t, n)    
    W = GAP.Globals.WeylGroup(GAP.Globals.RootSystem(L))
    orb = GAP.Globals.WeylOrbitIterator(W, GAP.Obj(hw))
    vertices = []
    
    # operate with the weylgroup on hw
    GAP.Globals.IsDoneIterator(orb)
    while !(GAP.Globals.IsDoneIterator(orb))
        w = GAP.Globals.NextIterator(orb) 
        push!(vertices, fromGap(w))
    end

    # return result
    vertices = transpose(hcat(vertices ...))
    vertices = [w_to_eps(t, n, w) for w in eachrow(vertices)]
    vertices = transpose(hcat(vertices ...))
    #println("wts_eps ", [w_to_eps(t, n, w) for w in eachrow(vertices)])
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


#############################################
#     compute monomials for weightspace     #
#############################################

function convert_lattice_points_to_monomials(ZZx, lattice_points_weightspace)
    return [finish(push_term!(MPolyBuildCtx(ZZx), ZZ(1), convert(Vector{Int}, convert(Vector{Int64}, lattice_point)))) for lattice_point in lattice_points_weightspace]
end

function get_lattice_points_of_weightspace(wts, weight, t)
    """
    calculates all lattice points in a given weightspace for lie algebras of type t
    input:
    wts: the operator weights in eps_i
    weight: lambda - mu

    output: all lattice points with weight weight
    """
    if t in ["A", "G"]
        return get_lattice_points_of_weightspace_A_G_n(wts, weight)
    else
        return get_lattice_points_of_weightspace_Xn(wts, weight)
    end
end

function get_lattice_points_of_weightspace_A_G_n(wts, weight)
    """
    calculates all monomials in a given weightspace for lie algebras that have type A or G
    input:
    wts: the operator weights in eps_i
    weight: lambda - mu

    output: all monomials with weight weight
    
    works by calculating all integer solutions to the following linear program:
    [ 1     |              |    ]       [   x   ]      
    [ 1   wts[1]   ...   wts[k] ]   *   [   |   ]   =   weight 
    [...    |              |    ]       [  res  ] 
    [ 1     |              |    ]       [   |   ]
    where res[i] >= 0 for all i

    example:
    wts = [[1, 0, 2], [-1, 1, 1], [0, -1, 0]] (i.e. a_1 = eps_1 - eps_2, a_2 = eps_2 - eps_3, a_12 = eps_1 - eps_3)
    weight = [2, 1, 0]
    -> poly = polytope.Polytope(INEQUALITIES=[0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1], EQUATIONS=[-2 1 1 0 2; -1 1 -1 1 1; 0 1 0 -1 0])
    => returns [[1 0 0], [1 1 0]]
    """
    # build linear (in-)equations
    wts = [reshape(w, 1, :) for w in wts]
    n = length(wts)
    ineq = zeros(Int64, n, n+2)
    for i in 1:n
        ineq[i, 2+i] = 1
    end
    equ = cat([-i for i in vec(weight)], [1 for i=1:length(weight)], dims = (2,2))
    equ = cat(equ, [transpose(w) for w in wts] ..., dims = (2,2))

    # find integer solutions of linear (in-)equation as lattice points of polytope
    poly = polytope.Polytope(INEQUALITIES=ineq, EQUATIONS=equ)

    # convert lattice-points to Oscar monomials
    lattice_points_weightspace = lattice_points(Polyhedron(poly))
    lattice_points_weightspace = [lattice_point[2:end] for lattice_point in lattice_points_weightspace]
    return lattice_points_weightspace

end

function get_lattice_points_of_weightspace_Xn(wts, weight)
    """
    calculates all lattice points in a given weightspace for lie algebras that don't have type A
    input:
    wts: the operator weights in eps_i
    weight: lambda - mu

    output: all lattice points with weight weight
    
    works by calculating all integer solutions to the following linear program:
    [  |              |    ]       [   x   ]      
    [wts[1]   ...   wts[k] ]   *   [   |   ]   =   weight 
    [  |              |    ]       [  res  ] 
    [  |              |    ]       [   |   ]
    where res[i] >= 0 for all i

    example:
    wts = [[1, 0, 2], [-1, 1, 1], [0, -1, 0]] (i.e. a_1 = eps_1 - eps_2, a_2 = eps_2 - eps_3, a_12 = eps_1 - eps_3)
    weight = [2, 1, 0]
    -> poly = polytope.Polytope(INEQUALITIES=[0 1 0 0; 0 0 1 0; 0 0 0 1], EQUATIONS=[-2 1 0 2; -1 -1 1 1; 0 0 -1 0])
    => returns 
    """
    wts = [reshape(w, 1, :) for w in wts]
    n = length(wts)
    ineq = zeros(Int64, n, n+1)
    for i in 1:n
        ineq[i, 1+i] = 1
    end
    #equ = cat([-i for i in vec(weight)], [1 for i=1:length(weight)], dims = (2,2))
    equ = [-i for i in vec(weight)]
    equ = cat(equ, [transpose(w) for w in wts] ..., dims = (2,2))
    poly = polytope.Polytope(INEQUALITIES=ineq, EQUATIONS=equ)
    lattice_points_weightspace = lattice_points(Polyhedron(poly))
    return lattice_points_weightspace
end
