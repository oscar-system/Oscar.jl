```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Tropical curves

## Introduction
A tropical curve is a graph with multiplicities on its edges.  If embedded, it is a polyhedral complex of dimension (at most) one.

#### Note:
- The type `TropicalCurve` can be thought of as subtype of `TropicalVariety` in the sense that it should have all properties and features of the latter.

## Construction
In addition to converting from `TropicalVariety`, objects of type `TropicalCurve` can be constructed from:
```@docs
```

## Properties
In addition to the properties inherited from `TropicalVariety`, objects of type `TropicalCurve` have the following exclusive properties:
```@docs
graph(tc::TropicalCurve)
```
