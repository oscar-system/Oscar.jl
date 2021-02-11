
import Markdown
import Base: ==

export Cone,
    Polyhedron,
    ambient_dim,
    codim,
    convex_hull,
    cube,
    dim,
    facets,
    feasible_region,
    lineality_space,
    LinearProgram,
    maximal_value,
    maximal_vertex,
    minimal_value,
    minimal_vertex,
    newton_polytope,
    normalized_volume,
    objective_function,
    recession_cone,
    solve_lp,
    support_function,
    pm_polytope,
    positive_hull,
    rays,
    vertices,
    visual,
    volume

include("helpers.jl")
include("Polyhedron.jl")
include("Cone.jl")
include("LinearProgram.jl")
