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

The *tropical Jacobian* of an abstract tropical curve is the group of divisors of degree zero modulo the subgroup of principal divisors. Here a principal divisor is the divisor associated to a piecewise-linear function on the vertices by the Laplacian operator. The tropical Jacobian is a finite abelian group, with order equal to the number of maximal spanning trees in the graph. It is isomorphic to $\prod{\mathbb{Z}/n_{i}\mathbb{Z}}$, where the $n_{i}$ are the nonzero elementary divisors of the Laplacian matrix. For more details, see [BN07](@cite).

## Construction

```@docs
TropicalCurve{M}(PC::PolyhedralComplex) where {M}
DivisorOnTropicalCurve(tc::TropicalCurve{M, EMB}, coeffs::Vector{Int}) where {M, EMB}
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
structure_tropical_jacobian(tc::TropicalCurve{M, EMB}) where {M, EMB}
```
