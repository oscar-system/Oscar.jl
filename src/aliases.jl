include(joinpath(pathof(AbstractAlgebra), "..", "Aliases.jl"))

# HACK/FIXME: remove these aliases once we have them in AA/Nemo/Hecke
@alias characteristic_polynomial charpoly  # FIXME
@alias minimal_polynomial minpoly  # FIXME

import Nemo: is_cyclo_type
import Nemo: is_maxreal_type
import Nemo: ZZModRing  # FIXME: remove if/once Nemo exports this
import Nemo: zzModRing  # FIXME: remove if/once Nemo exports this
import Nemo: FpField  # FIXME: remove if/once Nemo exports this
import Nemo: fpField  # FIXME: remove if/once Nemo exports this
include(joinpath(pathof(Nemo), "..", "Aliases.jl"))

#import Hecke: quadratic_genera
#import Hecke: has_algebra
#import Hecke: has_embedding
#import Hecke: has_matrix_action
#import Hecke: has_root
#import Hecke: ...
#include(joinpath(pathof(Hecke), "..", "Aliases.jl"))

# make some Julia names compatible with our naming conventions
@alias is_subset issubset
@alias is_valid isvalid

# add some shorthands for our functions
function number_of_cones end
function number_of_edges end
function number_of_facets end
function number_of_maximal_cells end
function number_of_maximal_cones end
function number_of_maximal_polyhedra end
function number_of_partitions end
function number_of_patches end
function number_of_points end
function number_of_polyhedra end
function number_of_positive_roots end               # from experimental/LieAlgebras
function number_of_rays end
function number_of_roots end                        # from experimental/LieAlgebras
function number_of_simple_roots end                 # from experimental/LieAlgebras
function number_of_vertices end

# these are kept for compatibility with Graphs.jl / GraphsBase.jl
@alias ne number_of_edges
@alias nv number_of_vertices

# aliases should (in general) also follow the usual snake case style
@alias n_cones number_of_cones
@alias n_edges number_of_edges
@alias n_facets number_of_facets
@alias n_maximal_cells number_of_maximal_cells
@alias n_maximal_cones number_of_maximal_cones
@alias n_maximal_polyhedra number_of_maximal_polyhedra
@alias n_partitions number_of_partitions
@alias n_patches number_of_patches
@alias n_points number_of_points
@alias n_polyhedra number_of_polyhedra
@alias n_positive_roots number_of_positive_roots
@alias n_rays number_of_rays
@alias n_roots number_of_roots
@alias n_simple_roots number_of_simple_roots
@alias n_vertices number_of_vertices
