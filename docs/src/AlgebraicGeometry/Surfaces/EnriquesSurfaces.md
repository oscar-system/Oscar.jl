```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```

# Borcherds' method for Enriques surfaces

An Enriques surface is a smooth, proper surface $Y$ over a field $k$ 
such that $Y$ has second étale Betti-number $b_2(Y)=10$ and numerically 
trivial canonical bundle. 

Let now $k$ be algebraically closed of characteristic not $2$. 
Under certain genericity assumptions on $Y$ a version of Borcherds' method allows to compute 
the image of the map
```math
\mathrm{Aut}(Y) \to O(\mathrm{Num}(Y)),
```
a fundamental domain of the action of ``\mathrm{Aut}(Y)`` on the nef cone of ``Y``,
and a complete set of representatives of the ``\mathrm{Aut}(Y)``-orbits of 
smooth rational curves, elliptic fibrations and polarizations.

See [BS22](@cite), [BS22*1](@cite), [BRS23](@cite) for the underlying algorithms and theory.

## Automorphisms
The two main entry points are the following:
```@docs
generic_enriques_surface(n::Int)
enriques_surface_automorphism_group(SY2::ZZLat, SX::ZZLat)
```
## The surface
```@docs
EnriquesBorcherdsCtx
borcherds_method(Y::EnriquesBorcherdsCtx; max_nchambers=-1)
splitting_roots_mod2(Y::EnriquesBorcherdsCtx)
root_invariant(Y::EnriquesBorcherdsCtx)
mass(ECtx::EnriquesBorcherdsCtx)
numerical_lattice(Y::EnriquesBorcherdsCtx)
numerical_lattice_of_K3_cover(Y::EnriquesBorcherdsCtx)
invariant_lattice_of_K3_cover(Y::EnriquesBorcherdsCtx)
```
## Chambers
```@docs
rays(D::EnriquesChamber)
isotropic_rays(D::EnriquesChamber)
walls(D::EnriquesChamber)
hom(D1::EnriquesChamber, D2::EnriquesChamber)
adjacent_chamber(D::EnriquesChamber, v::ZZMatrix)
chamber_invariants(Y::EnriquesBorcherdsCtx)
```
## Orbits of Nef divisors
```@docs
isomorphism_classes_polarizations(Y::EnriquesBorcherdsCtx, h::ZZMatrix)
isomorphism_classes_elliptic_fibrations(Y::EnriquesBorcherdsCtx)
reducible_fibers(Y::EnriquesBorcherdsCtx, fbar::TorQuadModuleElem)
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Simon Brandhorst](https://www.math.uni-sb.de/ag/brandhorst/index.php?lang=en).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
