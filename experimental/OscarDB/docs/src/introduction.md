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


## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what this means.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Antony Della Vecchia](https://antonydellavecchia.github.io/)
* Benjamin Lorenz

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

```@docs
get_db
find
find_one
```

