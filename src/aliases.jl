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
function number_of_maximal_cells end
function number_of_maximal_cones end
function number_of_maximal_polyhedra end
@alias ncones number_of_cones
@alias nedges number_of_edges
@alias nmaxcells number_of_maximal_cells            # decided to have but not used anywhere
@alias nmaxcones number_of_maximal_cones
@alias nmaxpolyhedra number_of_maximal_polyhedra    # decided to have but not used anywhere
