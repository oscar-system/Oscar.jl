```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Sheaves on Projective Space

We present two algorithms for computing sheaf cohomology over projective $n$-space.
The algorithms are based on Tate resolutions via the **B**ernstein-**G**elfand-**G**elfand-correspondence
as introduced in [EFS03](@cite) and on local cohomology (see [Eis98](@cite)), respectively. While the first
algorithm makes use of syzygy computations over the exterior algebra, the second algorithm is based on
syzygy computations over the symmetric algebra (see [DE02](@cite) for a tutorial). Thus, in most examples,
the first algorithm is much faster.

```@docs
sheaf_cohomology(M::ModuleFP{T}, l::Int, h::Int; algorithm::Symbol = :bgg) where {T <: MPolyDecRingElem}
```



