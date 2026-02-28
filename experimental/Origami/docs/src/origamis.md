```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Origamis

## Basic construction of origamis

```@docs
origami(h::PermGroupElem, v::PermGroupElem)
origami_disconnected(h::PermGroupElem, v::PermGroupElem, d::Integer)
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
