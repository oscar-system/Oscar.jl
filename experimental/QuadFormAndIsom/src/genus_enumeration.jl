###############################################################################
#
#  General interface for enumeration of genera for definite ZZLat
#
###############################################################################

#=====TODO
- Set up a general interface for flexbility on user usage;
- Setup a default "invariants function";
- Re-adapt neighbours methods for ZZLat and use randomness;
- Implement a serialisation process to save and upload partial results;
- Include the possibility to use isometry enumeration;
- Implement a type for neighbour multigraph
======#

#======User interface
Input:
- known -> finite list of known isometry classes (always non-empty by starting from a single lattice)
- alg_type -> how to enumerate neighbours: all of them (:exhaustive), orbits of them (:orbit), a random number (:random)
Optional:
- rand_neigh -> for random enumeration, how many randome neighbours are computed
- distinct -> if the lattices in "known" are known to be pairwise non-isometric
- invariant_func -> functions to compute isometry invariants for comparing lattices
- save_partial -> whether one wants to save iteratively new isometry classes (for instance for large genera)
- save_path -> in the case "save_partial" is true, where to save the lattices
- use_mass -> whether to use the mass formula as termination condition
- missing_mass -> if "use_mass" and "distinct" are true, and the partial mass of "known" is known, mention what is the part of the mass missing
- neigh_graph -> whether one wants to compute the neighbour multigraph (so remembering the edges)
- known_graph -> if "neigh_graph" is true and a part of the multigraph is known, one can mention it here
- prime_graph -> if "neigh_graph" is true, one can mention at which prime number construct the graphs.
======#

function enumerate_genus(
    known::Vector{ZZLat},
    alg_type::Symbol = :exhaustive;
    rand_neigh::{Int, Nothing} = nothing
    distinct::Bool = true,
    invariant_func::Function = default_func,
    save_partial::Bool = false,
    save_path::Union{IO, String, Nothing} = nothing,
    use_mass::Bool = true,
    missing_mass::Union{QQFieldElem, Nothing} = nothing,
    neigh_graph::Bool = false,
    known_graph::Union{KneserGraph, Nothing} = nothing,
    prime_graph::Union{Vector{ZZRingElem}, Nothing} = nothing
  )
end

function enumerate_genus(
    G::ZZGenus,
    alg_type::Symbol = :exhaustive;
    rand_neigh::{Int, Nothing} = nothing
    invariant_func::Function = default_func,
    save_partial::Bool = false,
    save_path::Union{IO, String, Nothing} = nothing,
    use_mass::Bool = true,
    neigh_graph::Bool = false,
    prime_graph::Union{Vector{ZZRingElem}, Nothing} = nothing
  )
end

function enumerate_genus(
    G::ZZLat,
    alg_type::Symbol = :exhaustive;
    rand_neigh::{Int, Nothing} = nothing
    invariant_func::Function = default_func,
    save_partial::Bool = false,
    save_path::Union{IO, String, Nothing} = nothing,
    use_mass::Bool = true,
    missing_mass::Union{QQFieldElem, Nothing} = nothing,
    neigh_graph::Bool = false,
    prime_graph::Union{Vector{ZZRingElem}, Nothing} = nothing
  )
end
