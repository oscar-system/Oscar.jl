```@meta
CurrentModule = Oscar
```

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
```

## Coherent sheaves of modules 

```@docs
SheafOfModules
twisting_sheaf(IP::AbsProjectiveScheme{<:Field}, d::Int)
tautological_bundle(IP::AbsProjectiveScheme{<:Field})
cotangent_sheaf(X::AbsCoveredScheme)
free_module(R::StructureSheafOfRings, n::Int)
projectivization(E::AbsCoherentSheaf; var_names::Vector{String}=Vector{String}(), check::Bool=true)
```
