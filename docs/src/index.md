# JToric -- toric geometry in Julia

```@contents
```

```@meta
CurrentModule = JToric
```

## Goal

The goal of the package `JToric` is to make computations on toric geometry accessible in the [OSCAR computer algebra system](https://oscar.computeralgebra.de/). As such, we follow [the book](https://www.ams.org/publications/authors/books/postpub/gsm-124) `Toric Varieties` by `David A. Cox`, `John B. Little` and `Henry K. Schenck`]. Currently, this project is work-in-progress. Our goal is to ensure that one can perform all computations of `Appendix B`.


### Generalities

```@docs
version
```

## Toric Varieties

### Constructors

```@docs
NormalToricVariety(PF::PolyhedralFan)
NormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )
projective_space
hirzebruch_surface
del_pezzo
```

### Conversion among GAP and Polymake toric varieties

```@docs
ntv_gap2polymake
ntv_polymake2gap
```

### Properties of toric varieties

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

### Attributes of toric varieties

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
```


### Methods of toric varieties

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


## Toric Divisors

### Constructors

```@docs
create_divisor
divisor_of_character
divisor_of_class
```

### Properties of toric divisors

```@docs
is_cartier
is_principal
is_primedivisor
is_basepoint_free
is_ample
is_very_ample
is_nef
```
