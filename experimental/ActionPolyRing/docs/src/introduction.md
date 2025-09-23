```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

This project aims to provide functionality for what we call action polynomial rings. This phrase refers to a generalized
framework that allows for an algorithmic treatment of both difference polynomial rings and differential polynomial rings.
In the future, further similar algebraic structure might be covered as well.

## Content

This project is structured as follows:

- The section [Action polynomial rings](@ref actionpolyring) contains generic methods that work for all subtypes of the action polynomial ring interface.
- The sections [Difference polynomial rings](@ref differencepolyring) and [Differential polynomial rings](@ref differentialpolyring) contain functionality that is unique to difference polynomial rings and differential polynomial rings respectively.
- The section [Rankings](@ref actionpolyranking) contains the functionality for rankings of jet variables of action polynomial rings.

## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what this means.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Tobias Braun](https://www.math.rwth-aachen.de/homes/Tobias.Braun/)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

