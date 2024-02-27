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

# predeclare some functions to allow defining aliases
function n_cones end
function n_connected_components end
function n_edges end
function n_facets end
function n_maximal_cells end
function n_maximal_cones end
function n_maximal_polyhedra end
function n_points end
function n_polyhedra end
function n_rays end
function n_vertices end

function number_of_partitions end
function number_of_patches end
function number_of_positive_roots end               # from experimental/LieAlgebras
function number_of_roots end                        # from experimental/LieAlgebras
function number_of_simple_roots end                 # from experimental/LieAlgebras

# these are kept for compatibility with Graphs.jl / GraphsBase.jl
@alias ne n_edges
@alias nv n_vertices

# aliases should (in general) also follow the usual snake case style
@alias number_of_cones n_cones
@alias number_of_connected_components n_connected_components
@alias number_of_edges n_edges
@alias number_of_facets n_facets
@alias number_of_maximal_cells n_maximal_cells
@alias number_of_maximal_cones n_maximal_cones
@alias number_of_maximal_polyhedra n_maximal_polyhedra
@alias number_of_points n_points
@alias number_of_polyhedra n_polyhedra
@alias number_of_rays n_rays
@alias number_of_vertices n_vertices

@alias n_partitions number_of_partitions
@alias n_patches number_of_patches
@alias n_positive_roots number_of_positive_roots
@alias n_roots number_of_roots
@alias n_simple_roots number_of_simple_roots

# aliases for consistency with oscar style
@alias n_columns ncols
@alias n_digits ndigits
@alias n_generators ngens
@alias n_rows nrows
@alias n_variables nvars
