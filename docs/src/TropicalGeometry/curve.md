```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["curve.md"]
```

# Curves


## Introduction

An *abstract tropical curve* is a finite, loopless, mulitgraph.
It is defined by an incidence matrix with vertices as columns and edges as rows.

A *divisor* on an abstract tropical curve is a formal linear combination of the vertices with integer coefficients. The degree of a divisor is the sum of its coefficients. A divisor is effective if all its coefficients are nonnegative.

For more definitions on the theory of divisors and linear sisyems on abstract tropical curve, we refer to [BN07](@cite).

## Construction

```@docs
TropicalCurve{M}(PC::PolyhedralComplex) where {M}
DivisorOnTropicalCurve(tc::TropicalCurve{M, EMB}, coeffs::Vector{Int}) where {M, EMB}
StructureTropicalJacobian(TC::TropicalCurve)
```

## Auxiliary functions
```@docs
graph(tc::TropicalCurve{M, EMB}) where {M, EMB}
n_nodes(tc::TropicalCurve{M, EMB}) where {M, EMB}
coefficients(dtc::DivisorOnTropicalCurve{M, EMB}) where {M, EMB}
degree(dtc::DivisorOnTropicalCurve{M, EMB}) where {M, EMB}
is_effective(dtc::DivisorOnTropicalCurve{M, EMB}) where {M, EMB}
chip_firing_move(dtc::DivisorOnTropicalCurve{M, EMB}, position::Int) where {M, EMB}
v_reduced(dtc::DivisorOnTropicalCurve{M, EMB}, vertex::Int) where {M, EMB}
is_linearly_equivalent(dtc1::DivisorOnTropicalCurve{M, EMB}, dtc2::DivisorOnTropicalCurve{M, EMB}) where {M, EMB}
```
