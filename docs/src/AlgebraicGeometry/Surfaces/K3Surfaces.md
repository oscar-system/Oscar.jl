```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Automorphism Groups of  K3 surfaces

A complex K3 surface is a compact complex surface $X$
with vanishing irregularity $h^1(X, \mathcal{O}_X)=0$
and trivial canonical bundle $\mathcal{O}_X\cong \omega_X$.

Much of the theory of (complex) K3 surfaces is governed by
its Hodge structure and the $\mathbb{Z}$-lattices
$NS(X) \subseteq H^2(X, \mathbb{Z})$.

See [Huy16](@cite) for the theory of K3 surfaces.

## Automorphisms

```@docs
K3_surface_automorphism_group(S::ZZLat)
borcherds_method
K3Chamber
chamber(data::BorcherdsCtx, weyl_vector::ZZMatrix, parent_wall::ZZMatrix=zero_matrix(ZZ, 0, 0))
weyl_vector(D::K3Chamber)
walls(::K3Chamber)
inner_point(::K3Chamber)
rays(::K3Chamber)
aut(::K3Chamber)
hom(::K3Chamber,::K3Chamber)
adjacent_chamber(D::K3Chamber, v::ZZMatrix)
separating_hyperplanes(S::ZZLat, v::QQMatrix, h::QQMatrix, d)
has_zero_entropy
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en).
* [Matthias Zach](https://math.rptu.de/en/wgs/agag/people/members),

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
