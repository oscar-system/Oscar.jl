# Origamis

## Basic construction of origamis

```@docs
origami(h::PermGroupElem, v::PermGroupElem; check=true)
horizontal_perm(o::Origami)
vertical_perm(o::Origami)
degree(o::Origami)
perm_group(o::Origami)
are_equivalent(o1::Origami, o2::Origami)
```

## Computing attributes of origamis

```@docs
stratum(o::Origami)
genus(o::Origami)
index_monodromy_group(o::Origami)
sum_of_lyapunov_exponents(o::Origami)
translations(o::Origami)
point_reflections(o::Origami)
is_hyperelliptic(o::Origami)
cylinder_structure(o::Origami)
```

## Normal forms of origamis

```@docs
normal_form(o::Origami)
```

## ${\rm SL}_2(\mathbb{Z})$ actions on origamis

```@docs
action_s(o::Origami)
action_t(o::Origami)
action_s_inv(o::Origami)
action_t_inv(o::Origami)
action_sl2(A::ZZMatrix,o::Origami)
```

## Cyclic torus covers

```@docs
generalized_cyclic_torus_cover(n::Int, d::Int, vslits::Vector{Int64}, hslits::Vector{Int64})
comb_origami(n::Int, x::Int, y::Int)
cyclic_torus_cover_origamiS(n::Int, d::Int, v::Vector{Int64})
cyclic_torus_cover_origamiL(n::Int, d::Int, v::Vector{Int64})
base_change_l_to_s(n::Int)
```
