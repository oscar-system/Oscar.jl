```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["tg_ntv.md"]
```

# Normal Toric Varities


## Introduction

We introduce two main types of normal toric varieties, distinguishing between
the affine and non-affine case:
- `AffineNormalToricVariety` is the toric variety associated to a cone $\sigma$, denoted by $U_{\sigma}$ in [CLS11](@cite)
- `NormalToricVariety` is the toric variety associated to a polyhedral fan $\Sigma$, denoted by $X_{\Sigma}$ in [CLS11](@cite)

!!! warning
    The lattice is always assumed to be the standard lattice $\mathbb{Z}^n$.
    Transformations for non-standard lattices will have to be done by the user.

## Construction

```@docs
AffineNormalToricVariety(C::Cone)
NormalToricVariety(C::Cone)
NormalToricVariety(PF::PolyhedralFan)
NormalToricVariety(P::Polyhedron)
```

## Auxiliary functions
```@docs
isaffine(ntv::NormalToricVarieties)
iscomplete(ntv::NormalToricVarieties)
isnormal_variety(ntv::NormalToricVarieties)
isprojective(ntv::NormalToricVarieties)
issimplicial(ntv::NormalToricVarieties)
issmooth(ntv::NormalToricVarieties)
toric_ideal_binomial_generators(ntv::AffineNormalToricVariety)
```
