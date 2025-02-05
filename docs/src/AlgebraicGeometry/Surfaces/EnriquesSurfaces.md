```@meta
CurrentModule = Oscar
```

# Automorphism Groups of Enriques surfaces

An Enriques surface is a smooth, proper surface $Y$ over a field $k$ 
such that $Y$ has second étale Betti-number $b_2(Y)=10$ and numerically 
trivial canonical bundle. 

If the characteristic of $k$ is not $2$, then canonical bundle induces a ``2:1``-étale cover 
$\pi \colon X \to Y$ where $X$ is a K3 surface. 

Let now $k$ be algebraically closed. We denote by $S_Y = \mathrm{Num}(Y)$
and $S_X = \mathrm{Num}(X)$ the numerical lattices of divisors on $Y$ and $X$. 

Under certain genericity assumptions on $Y$, Oscar can compute the image of the map
```math
\mathrm{Aut}(Y) \to O(S_Y)
```
and further invariants. 

[BS, BSR, BS](@cite)
## Automorphisms

```@docs
K3_surface_automorphism_group(S::ZZLat)
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
