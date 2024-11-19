isdefined(Oscar, :word) || function word end

module LieAlgebras

using ..Oscar

import Oscar: GAPWrap, IntegerUnion, MapHeader

import Random

# not importet in Oscar
using AbstractAlgebra:
  ProductIterator, _number_of_direct_product_factors, ordinal_number_string

using AbstractAlgebra.PrettyPrinting

# functions with new methods
import ..Oscar:
  _is_exterior_power,
  _is_tensor_product,
  _iso_oscar_gap,
  _vec,
  action,
  add!,
  addmul!,
  basis_matrix,
  basis,
  canonical_injection,
  canonical_injections,
  canonical_projection,
  canonical_projections,
  center,
  centralizer,
  character,
  characteristic,
  coeff,
  coefficient_ring,
  coefficients,
  compose,
  demazure_character,
  derived_series,
  dim,
  direct_sum,
  dot,
  dual,
  elem_type,
  expressify,
  exterior_power,
  fp_group,
  gen,
  gens,
  height,
  hom,
  hom_direct_sum,
  hom_tensor,
  ideal,
  identity_map,
  image,
  induced_map_on_exterior_power,
  inv,
  is_abelian,
  is_finite,
  is_gen,
  is_isomorphism,
  is_nilpotent,
  is_perfect,
  is_simple,
  is_solvable,
  is_welldefined,
  isomorphism,
  kernel,
  lower_central_series,
  matrix,
  mul!,
  neg!,
  normalizer,
  number_of_generators,
  ngens,
  order,
  parent_type,
  permutation_group,
  rank,
  root,
  roots,
  sub,
  sub!,
  symbols,
  symmetric_power,
  tensor_product,
  weyl_vector,
  word,
  zero!,
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

# Aliases
function number_of_positive_roots end
function number_of_roots end
function number_of_simple_roots end

@alias n_positive_roots number_of_positive_roots
@alias n_roots number_of_roots
@alias n_simple_roots number_of_simple_roots

include("Types.jl")
include("Combinatorics.jl")
include("Util.jl")

include("CartanMatrix.jl")
include("CoxeterGroup.jl")
include("RootSystem.jl")
include("WeightLattice.jl")
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

include("serialization.jl")

end # module LieAlgebras

using .LieAlgebras

include("exports.jl")

export n_positive_roots  # alias lives in a submodule
export n_roots           # alias lives in a submodule
export n_simple_roots    # alias lives in a submodule
