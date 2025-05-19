module QuantumGroups

using ..Oscar

import Base: hash, deepcopy_internal

import AbstractAlgebra.Generic:
  FracFieldElem,
  LaurentPolyWrap

import Nemo:
  set!

import ..Oscar:
  add!,
  addmul!,
  coeff,
  coefficient_ring,
  div,
  div!,
  elem_type,
  exponent_vector,
  expressify,
  gen,
  gens,
  inv,
  is_domain_type,
  isone,
  iszero,
  length,
  mul!,
  ngens,
  nvars,
  one,
  one!,
  parent,
  parent_type,
  q_binomial,
  q_factorial,
  root_system,
  setcoeff!,
  sub!,
  symbols,
  zero,
  zero!

include("exports.jl")

include("Types.jl")

include("CanonicalBasis.jl")
include("MPolyRing.jl")
include("PBWAlgebra.jl")
include("PBWAlgebraHom.jl")
include("QuantumField.jl")
include("QuantumGroup.jl")
include("QuantumGroupHom.jl")

end

using .QuantumGroups

include("exports.jl")
