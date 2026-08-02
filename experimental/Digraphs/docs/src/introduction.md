```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Digraphs

This part of Oscar provides directed graphs (digraphs), backed by the GAP
package [Digraphs](https://digraphs.github.io/Digraphs/). A `Digraph` is a
set of vertices, identified with `1, ..., n`, together with a list of
out-neighbours for every vertex, so that multiple edges and loops are
supported. The functionality mirrors the GAP Digraphs manual: constructing
digraphs, operators and attributes, properties, and isomorphisms and
homomorphisms. Digraphs can also be converted to and from Oscar's directed,
undirected and mixed `Graph` types; see [Conversions with Oscar graphs](@ref).

```jldoctest
julia> D = cycle_digraph(5)
Digraph with 5 vertices, 5 edges

julia> is_isomorphic(D, digraph_reverse(D))
true

julia> order(automorphism_group(D))
5
```

## Status

This part of OSCAR is in an experimental state; please see [Adding new
projects to experimental](@ref) for what this means.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Zhi Ming Long](mailto:longgoforit@163.com)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
