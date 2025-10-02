```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# OscarDB

This module provides a general database framework which, conceptually, works for all OSCAR types.

There are two main ingredients:
- the database backend is [MongoDB](https://www.mongodb.com);
- communication with the database relies on our serialization.

The overall design is inspired by polymake's [PolyDB](https://polydb.org/).

## Collections

The database is organized into collections.
Within each collection the datasets are uniform.

### Vertex-transitive combinatorial manifolds [TransitiveSimplicialComplexes]
The OSCAR DB provides access to [Frank Lutz' collection of vertex transitive triangulations](https://www3.math.tu-berlin.de/IfM/Nachrufe/Frank_Lutz/stellar/vertex-transitive-triangulations.html).
It comprises all vertex-transitive combinatorial manifolds with up to 15 vertices in dimensions $d=2,3,9,10,11,12$.
In the remaining dimensions below 12, the enumeration is complete up to 13 vertices.

See also [Lut08](@cite).

### Leech pairs [LeechPairs]
The OSCAR DB provides access to the Leech pairs computed by G. Höhn and G. Mason [HM16](@ref).

## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what this means.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Antony Della Vecchia](https://antonydellavecchia.github.io/)
* Benjamin Lorenz

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

```@meta
CurrentModule = Oscar.OscarDB
```

```@docs
get_db
find_one
find
length
getindex(db::Database, name::AbstractString)
get_collection_names
```

