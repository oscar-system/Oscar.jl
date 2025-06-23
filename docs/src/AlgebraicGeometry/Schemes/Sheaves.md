# Sheaves on covered schemes

Oscar supports modeling sheaves by means of a covering by affine charts.

## Presheaves
```@docs
AbsPreSheaf
PreSheafOnScheme
```

## Structure sheaves
```@docs
StructureSheafOfRings
```

## Ideal sheaves 
```@docs
AbsIdealSheaf
IdealSheaf
PrimeIdealSheafFromChart
SumIdealSheaf
ProductIdealSheaf
SimplifiedIdealSheaf
PullbackIdealSheaf
RadicalOfIdealSheaf
ToricIdealSheafFromCoxRingIdeal
SingularLocusIdealSheaf
```

## Coherent sheaves of modules 

These are some types for coherent sheaves.
```@docs
SheafOfModules
HomSheaf
PushforwardSheaf
PullbackSheaf
DirectSumSheaf
```

We provide some common constructors.
```@docs
twisting_sheaf(IP::AbsProjectiveScheme{<:Field}, d::Int)
tautological_bundle(IP::AbsProjectiveScheme{<:Field})
cotangent_sheaf(X::AbsCoveredScheme)
free_module(R::StructureSheafOfRings, n::Int)
dual(M::SheafOfModules)
tangent_sheaf(X::AbsCoveredScheme)
```

```@docs
projectivization(E::AbsCoherentSheaf; var_names::Vector{String}=Vector{String}(), check::Bool=true)
```
