module QuantumGroups

using ..Oscar

import Base: hash, deepcopy_internal, rand

import Random

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
  characteristic,
  div,
  div!,
  divexact,
  divexact!,
  divexact_right,
  elem_type,
  exponent_vector,
  expressify,
  gen,
  gens,
  image,
  inv,
  inv!,
  is_exact_type,
  is_domain_type,
  is_unit,
  isone,
  iszero,
  length,
  monomial,
  mul!,
  neg!,
  ngens,
  nvars,
  one,
  one!,
  parent,
  parent_type,
  q_integer,
  q_binomial,
  q_factorial,
  root_system,
  setcoeff!,
  sub!,
  submul!,
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
