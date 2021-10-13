```@meta
CurrentModule = Oscar
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
NormalToricVariety(P::Polyhedron)
NormalToricVariety( r::Matrix{Int}, c::Vector{Vector{Int}} )
projective_space
hirzebruch_surface
del_pezzo
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
is_gorenstein
is_q_gorenstein
is_fano
```


## Attributes of toric varieties

```@docs
dim
dim_of_torusfactor
picard_group
nef_cone
mori_cone
toric_ideal_binomial_generators(antv::AffineNormalToricVariety)
```
