# Add your new types, functions, and methods here.

module FiniteRingUnits

import ..Oscar:
  @attr,
  @req,
  @vprint,
  @vprintln,
  Dedent,
  DirectProductGroup,
  FPGroup,
  FinField,
  FinGenAbGroup,
  FiniteRing,
  FiniteRingElem,
  GAPGroup,
  GL,
  Hecke,
  Hecke.AbstractAssociativeAlgebra,
  Hecke.FiniteRings.OnePlusIdeal,
  Hecke.FiniteRings.OnePlusIdealModuloOnePlusIdeal,
  Hecke.FiniteRings.OnePlusIdealModuloOnePlusIdealElem,
  Hecke.FiniteRings._ideal_zgens,
  Hecke.FiniteRings._zgens,
  Hecke.FiniteRings.data,
  Hecke.FiniteRings.decompose_into_p_rings,
  Hecke.FiniteRings.decompose_semisimple_p_ring,
  Hecke.FiniteRings.underlying_abelian_group,
  Indent,
  Lowercase,
  Map,
  Oscar,
  StructureConstantAlgebra,
  ZZ,
  abelian_group,
  absolute_degree,
  base_ring,
  characteristic,
  codomain,
  coefficients,
  decompose,
  decompose_into_indecomposable_rings,
  degree,
  dim,
  direct_product,
  domain,
  elem_type,
  free_group,
  gen,
  gens,
  get_attribute!,
  has_preimage_with_preimage,
  hom,
  id_hom,
  ideal,
  identity_matrix,
  image,
  is_one,
  is_prime,
  is_unit,
  is_zero,
  isomorphism,
  kernel,
  lift,
  map_word,
  matrix,
  ngens,
  nrows,
  order,
  preimage,
  preimage,
  pretty,
  quo,
  radical,
  relators,
  set_order,
  simplified_fp_group,
  snf,
  unit_group,
  zero_matrix,
  finite_ring


include("Types.jl")
include("EffectivePresentation.jl")
include("Units.jl")
include("GL.jl")
include("Map.jl")
include("Ktheory.jl")

export k1
export abelianization_of_unit_group

end

import .FiniteRingUnits:
  k1,
  abelianization_of_unit_group
