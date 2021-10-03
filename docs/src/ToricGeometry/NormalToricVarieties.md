```@meta
CurrentModule = JToric
```

```@contents
Pages = ["NormalToricVarieties.md"]
```

# Normal Toric Varieties

## Introduction

We introduce two main types of normal toric varieties, distinguishing between
the affine and non-affine case:
- `AffineNormalToricVariety` is the toric variety associated to a cone $\sigma$, denoted by $U_{\sigma}$ in [CLS11](@cite)
- `NormalToricVariety` is the toric variety associated to a polyhedral fan $\Sigma$, denoted by $X_{\Sigma}$ in [CLS11](@cite)

!!! warning
    The lattice is always assumed to be the standard lattice $\mathbb{Z}^n$.
    Transformations for non-standard lattices will have to be done by the user.


## Constructors

```@docs
AffineNormalToricVariety(C::Cone)
NormalToricVariety(C::Cone)
NormalToricVariety(PF::PolyhedralFan)
NormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )
projective_space
hirzebruch_surface
del_pezzo
```

## Conversion among GAP and Polymake toric varieties

```@docs
ntv_gap2polymake
ntv_polymake2gap
```

## Properties of toric varieties

```@docs
isnormal
isaffine
isprojective
issmooth
iscomplete
has_torusfactor
is_orbifold
issimplicial
is_isomorphic_to_projective_space
is_direct_product_of_projective_spaces
is_gorenstein
is_q_gorenstein
is_fano
```

## Attributes of toric varieties

```@docs
affine_open_covering
cox_ring
list_of_variables_of_cox_ring
class_group
torus_invariant_divisor_group
map_from_character_to_principal_divisor
map_from_weil_divisors_to_class_group
dim
dim_of_torusfactor
coordinate_ring_of_torus
list_of_variables_of_coordinate_ring_of_torus
is_product_of
factors
character_lattice
torus_invariant_prime_divisors
irrelevant_ideal
stanley_reisner_ideal
morphism_from_cox_variety
cox_variety
fan_of_variety
fan
cartier_torus_invariant_divisor_group
picard_group
name_of_variety
zariski_cotangent_sheaf
cotangent_sheaf
euler_characteristic
weil_divisors_of_variety
zariski_cotangent_sheaf_via_euler_sequence
zariski_cotangent_sheaf_via_poincare_residue_map
nef_cone
mori_cone
toric_ideal_binomial_generators(antv::AffineNormalToricVariety)
```


## Methods of toric varieties

```@docs
set_name_of_variety
coordinate_ring_of_torus( v::NormalToricVariety, names::Vector{String} )
cox_ring( v::NormalToricVariety, name::String )
Base.:*( v::NormalToricVariety, w::NormalToricVariety )
character_to_rational_function( l::Vector{Int}, v::NormalToricVariety )
blowup_on_ith_minimal_torus_orbit( v::NormalToricVariety, i::Int )
ith_betti_number( v::NormalToricVariety, i::Int )
nr_of_q_rational_points( v::NormalToricVariety, i::Int )
```

