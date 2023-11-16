# standard packages
using Pkg
using Random
using RandomExtensions
using UUIDs

# our packages
import AbstractAlgebra
import AlgebraicSolving
# we currently need to load Polymake before GAP to avoid the crashe mentioned in
# https://github.com/oscar-system/Oscar.jl/pull/1902
# Once there is a GAP_pkg_browse that links to the correct ncurses we might
# switch this back.
import Polymake
import GAP
import Hecke
import Nemo
import Singular

# import stuff from Base for which we want to provide extra methods
import Base:
  +,
  -,
  *,
  ^,
  ==,
  conj,
  convert,
  deepcopy_internal,
  eltype,
  exponent,
  getindex,
  hash,
  intersect,
  inv,
  isfinite,
  issubset,
  iterate,
  length,
  mod,
  one,
  parent,
  print,
  reduce,
  show,
  sum,
  union,
  values,
  Vector,
  zero

import AbstractAlgebra:
  @alias,
  @attr,
  @attributes,
  @show_name,
  @show_special,
  addeq!,
  base_ring,
  canonical_unit,
  codomain,
  data,
  Dedent,
  degree,
  dim,
  domain,
  elem_type,
  evaluate,
  expressify,
  Field,
  FieldElem,
  force_coerce,
  force_op,
  gen,
  Generic,
  Generic.finish,
  Generic.MPolyBuildCtx,
  Generic.MPolyCoeffs,
  Generic.MPolyExponentVectors,
  Generic.push_term!,
  gens,
  get_attribute,
  get_attribute!,
  Ideal,
  Indent,
  Lowercase,
  LowercaseOff,
  map,
  Map,
  MatElem,
  matrix,
  MatSpace,
  MPolyRing,
  MPolyRingElem,
  NCRing,
  NCRingElem,
  ngens,
  nvars,
  ordering,
  parent_type,
  polynomial_ring,
  PolyRing,
  PolyRingElem,
  pretty,
  Ring,
  RingElem,
  RingElement,
  set_attribute!,
  SetMap,
  symbols,
  total_degree

import AbstractAlgebra.GroupsCore
import AbstractAlgebra.GroupsCore:
  hasgens,
  isfiniteorder,
  istrivial

import GAP:
  @gapattribute,
  GapInt,
  GapObj

import Nemo:
  bell,
  binomial,
  denominator,
  divexact,
  divides,
  divisor_sigma,
  euler_phi,
  factorial,
  fibonacci,
  fits,
  fqPolyRepFieldElem,
  fraction_field,
  height,
  is_prime,
  is_probable_prime,
  is_square,
  is_unit,
  isqrtrem,
  jacobi_symbol,
  matrix_space,
  moebius_mu,
  number_of_partitions,
  numerator,
  primorial,
  QQ,
  QQField,
  QQFieldElem,
  QQMatrix,
  rising_factorial,
  root,
  unit,
  ZZ,
  ZZMatrix,
  ZZRing,
  ZZRingElem

let exclude_hecke = [
    :change_uniformizer,
    :coefficients,
    :exponent_vectors,
    :leading_coefficient,
    :leading_monomial,
    :leading_term,
    :monomials,
    :narrow_class_group,
    :Partition,
    :perm,
    :QQBar,
    :SymmetricGroup,
    :tail,
    :terms,
  ]
  for i in names(Hecke)
    (i in exclude_hecke || !isdefined(Hecke, i)) && continue
    @eval import Hecke: $i
    @eval export $i
  end
end

import Hecke:
  conjugate,
  expand,
  field_extension,
  hensel_qf,
  IntegerUnion,
  MapHeader,
  multiplicative_jordan_decomposition,
  primitive_element,
  QQBar

import cohomCalg_jll
