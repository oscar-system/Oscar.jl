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
@alias ncones number_of_cones
@alias nedges number_of_edges
@alias nfacets number_of_facets
@alias nmaxcells number_of_maximal_cells            # decided to have but not used anywhere
@alias nmaxcones number_of_maximal_cones
@alias nmaxpolyhedra number_of_maximal_polyhedra    # decided to have but not used anywhere
@alias npartitions number_of_partitions             # decided to have but not used anywhere
@alias npatches number_of_patches
@alias npoints number_of_points                     # decided to have but not used anywhere
@alias npolyhedra number_of_polyhedra               # decided to have but not used anywhere
@alias nposroots number_of_positive_roots
@alias nrays number_of_rays
@alias nroots number_of_roots
@alias nsimpleroots number_of_simple_roots
@alias nvertices number_of_vertices

# these are kept for compatibility with Graphs.jl / GraphsBase.jl
@alias ne number_of_edges
@alias nv number_of_vertices
