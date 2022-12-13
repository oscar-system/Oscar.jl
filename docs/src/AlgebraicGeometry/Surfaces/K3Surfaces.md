```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["K3Surfaces.md"]
```
# Algebraic Surfaces

## K3 surfaces

A complex K3 surface is a compact complex surface $S$
with vanishing irregularity $h^1(S, \mathcal{O}_S)=0$
and trivial canonical bundle $\mathcal{O}_S\cong \omega_S$.

Much of the theory of (complex) K3 surfaces is governed by
its Hodge structure and the $\mathbb{Z}$-lattices
$NS(X) \subseteq H^2(X, \mathbb{Z})$.
See [Huy16](@cite) for the theory of K3 surfaces.

### Automorphisms

```@docs
K3_surface_automorphism_group(S::ZLat)
borcherds_method
K3Chamber
chamber(data::BorcherdsCtx, weyl_vector::fmpz_mat, parent_wall::fmpz_mat=zero_matrix(ZZ, 0, 0))
weyl_vector(D::K3Chamber)
walls(::K3Chamber)
inner_point(::K3Chamber)
rays(::K3Chamber)
aut(::K3Chamber)
hom(::K3Chamber,::K3Chamber)
adjacent_chamber(D::K3Chamber, v)
separating_hyperplanes(S::ZLat, v::fmpq_mat, h::fmpq_mat, d)
has_zero_entropy
```
