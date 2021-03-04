
import Markdown
import Base: ==

const AnyVecOrMat = Union{MatElem, AbstractVecOrMat}

export Cone,
    PolyhedralFan,
    Polyhedron,
    IncidenceMatrix,
    ambient_dim,
    codim,
    convex_hull,
    cross,
    cube,
    dim,
    faces,
    facets,
    facets_as_halfspace_matrix_pair,
    facets_as_point_matrix,
    face_fan,
    feasible_region,
    f_vector,
    hilbert_basis,
    intersect,
    isbounded,
    iscomplete,
    isfeasible,
    isfulldimensional,
    isnormal,
    ispointed,
    isregular,
    issmooth,
    lattice_points,
    lineality_space,
    linear_symmetries,
    LinearProgram,
    maximal_cones,
    maximal_cones_as_incidence_matrix,
    maximal_value,
    maximal_vertex,
    minimal_value,
    minimal_vertex,
    minkowski_sum,
    newton_polytope,
    normalized_volume,
    normal_fan,
    nfacets,
    nmaximal_cones,
    nrays,
    nvertices,
    objective_function,
    orbit_polytope,
    recession_cone,
    simplex,
    solve_lp,
    support_function,
    pm_cone,
    pm_polytope,
    positive_hull,
    rays,
    rays_as_point_matrix,
    vertices,
    vertices_as_point_matrix,
    visual,
    volume

include("helpers.jl")
include("Polyhedron.jl")
include("Cone.jl")
include("LinearProgram.jl")
include("PolyhedralFan.jl")
include("Groups.jl")
