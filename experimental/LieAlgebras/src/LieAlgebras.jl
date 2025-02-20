module LieAlgebras

using ..Oscar

using Oscar:
  _root_system_type_string,
  _vec,
  GAPWrap,
  IntegerUnion,
  MapHeader,
  set_root_system_type!

import Random

# not importet in Oscar
using AbstractAlgebra:
  ProductIterator, _number_of_direct_product_factors, ordinal_number_string

using AbstractAlgebra.PrettyPrinting

# functions with new methods
import ..Oscar:
  FPGroup,
  PermGroup,
  _is_exterior_power,
  _is_tensor_product,
  _iso_oscar_gap,
  action,
  basis_matrix,
  basis,
  canonical_injection,
  canonical_injections,
  canonical_projection,
  canonical_projections,
  cartan_matrix,
  center,
  centralizer,
  change_base_ring,
  character,
  characteristic,
  check_parent,
  coeff,
  coefficient_ring,
  coefficients,
  compose,
  demazure_character,
  derived_series,
  dim,
  direct_sum,
  dual,
  elem_type,
  expressify,
  exterior_power,
  gen,
  gens,
  hom,
  hom_direct_sum,
  hom_tensor,
  ideal,
  identity_map,
  image,
  induced_map_on_exterior_power,
  inner_direct_product,
  inv,
  is_abelian,
  is_isomorphism,
  is_nilpotent,
  is_perfect,
  is_semisimple,
  is_simple,
  is_solvable,
  is_welldefined,
  isomorphism,
  kernel,
  lower_central_series,
  matrix,
  normalizer,
  number_of_generators,
  ngens,
  parent_type,
  rank,
  root_system,
  structure_constant_table,
  sub,
  symbols,
  symmetric_power,
  tensor_product,
  zero_map,
  ⊕,
  ⊗

Oscar.@import_all_serialization_functions

import Base: getindex, deepcopy_internal, hash, issubset, iszero, parent, zero

include("exports.jl")

# The following export statements export things *into* module Oscar;
# they are however for internal use and thus are not exported *from* Oscar.
export _is_direct_sum
export _is_dual
export _is_exterior_power
export _is_standard_module
export _is_symmetric_power
export _is_tensor_power
export _is_tensor_product

include("Types.jl")
include("Combinatorics.jl")
include("Util.jl")

include("CoxeterGroup.jl")
include("DynkinDiagram.jl")
include("WeylGroup.jl")

include("LieAlgebra.jl")
include("AbstractLieAlgebra.jl")
include("LinearLieAlgebra.jl")
include("DirectSumLieAlgebra.jl")

include("LieSubalgebra.jl")
include("LieAlgebraIdeal.jl")
include("LieAlgebraHom.jl")

include("LieAlgebraModule.jl")
include("LieAlgebraModuleHom.jl")

include("iso_oscar_gap.jl")
include("iso_gap_oscar.jl")
include("GapWrapper.jl")

include("SSLieAlgebraModule.jl")

include("serialization.jl")

end # module LieAlgebras

using .LieAlgebras

include("exports.jl")

export n_positive_roots  # alias lives in a submodule
export n_roots           # alias lives in a submodule
export n_simple_roots    # alias lives in a submodule
