```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```


# Introduction

Algebraic shifting is a widely applicable method for converting a uniform hypergraph $S$ into another hypergraph, $\Delta(S)$, which is somehow simpler but still retains key properties of $S$.
Often this is applied to simplicial complexes, where each layer of $(k-1)$-faces forms a $k$-uniform hypergraph.

Here we focus on full and partial algebraic shifting in the exterior algebra.
For general background on the subject we refer to [Kal02](@cite); for details concerning the implementation see [D-VJL24](@cite).


## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what this means.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Antony Della Vecchia](https://antonydellavecchia.github.io)
* [Michael Joswig](https://page.math.tu-berlin.de/~joswig/)


You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
