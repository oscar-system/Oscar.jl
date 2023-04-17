# standard packages
using Pkg
using Random
using RandomExtensions
using Test

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
    *,
    ^,
    ==,
    conj,
    convert,
    eltype,
    exponent,
    getindex,
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
    Map,
    map,
    MatElem,
    matrix,
    MatSpace,
    MPolyRingElem,
    MPolyRing,
    NCRing,
    NCRingElem,
    ngens,
    nvars,
    ordering,
    parent_type,
    PolyRingElem,
    polynomial_ring,
    PolyRing,
    Ring,
    RingElem,
    RingElement,
    set_attribute!,
    SetMap,
    symbols,
    total_degree

# FIXME/TODO: clean up the following once AbstractAlgebra provides the new name
if isdefined(AbstractAlgebra, :MPolyRingElem)
  import AbstractAlgebra: MPolyRingElem
else
  @alias MPolyRingElem MPolyRingElem
end

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
    ZZRing,
    QQField,
    QQFieldElem,
    QQMatrix,
    ZZRingElem,
    ZZMatrix,
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
    rising_factorial,
    root,
    unit,
    ZZ

exclude = [:Nemo, :AbstractAlgebra, :Rational, :change_uniformizer,
    :genus_symbol, :data, :narrow_class_group, :perm, :SymmetricGroup,
    :coefficients, :leading_coefficient, :coefficients_and_exponents,
    :exponents, :exponent_vectors, :monomials, :leading_monomial, :terms,
    :leading_term, :tail, :Partition]

for i in names(Hecke)
  (i in exclude || !isdefined(Hecke, i)) && continue
  @eval import Hecke: $i
  @eval export $i
end

import Hecke:
    _two_adic_normal_forms,
    _block_indices_vals,
    _jordan_2_adic,
    _jordan_odd_adic,
    _min_val,
    _normalize,
    _rational_canonical_form_setup,
    _solve_X_ker,
    _val,
    @req,
    abelian_group,
    automorphism_group,
    center,
    cokernel,
    compose,
    defining_polynomial,
    derived_series,
    det,
    direct_product,
    elements,
    field_extension,
    FinField,
    FinFieldElem,
    fqPolyRepField,
    free_abelian_group,
    gens,
    gram_matrix,
    gram_matrix_quadratic,
    haspreimage,
    hensel_qf,
    hom,
    id_hom,
    image,
    index,
    IntegerUnion,
    inv!,
    is_abelian,
    is_bijective,
    is_characteristic,
    is_conjugate,
    is_cyclic,
    is_injective,
    is_invertible,
    is_isomorphic,
    is_normal,
    is_primitive,
    is_regular,
    is_simple,
    is_subgroup,
    is_surjective,
    kernel,
    Map,
    MapHeader,
    mul,
    mul!,
    multiplicative_jordan_decomposition,
    normal_closure,
    nrows,
    one!,
    order,
    preimage,
    primitive_element,
    quo,
    radical,
    refine_for_jordan,
    representative,
    small_group,
    sub,
    subsets,
    subgroups,
    TorQuadModule,
    TorQuadModuleElem,
    TorQuadModuleMor,
    tr,
    trace

import cohomCalg_jll
