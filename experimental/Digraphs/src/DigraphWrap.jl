module DigraphWrap

using GAP

# Use the @wrap macro to set up optimized versions of frequently used Digraphs
# GAP functions. Instead of writing e.g. GAP.Globals.DigraphNrVertices(x)
# you would write DigraphWrap.DigraphNrVertices(x). The former always performs
# a variable lookup in the GAP kernel and is not type stable. The latter
# accesses the underlying GAP function object directly and also has information
# about the return type.
#
# Note that the macro GAP.@wrap has a similar purpose as @gapattribute has,
# but works on a much lower level on purpose. We may actually phase out
# use of @gapattribute in the future.
#
# This list is sorted according to LC_COLLATE=C sort -f
GAP.@wrap AdjacencyMatrix(x::GapObj)::GapObj
GAP.@wrap AutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap BlissAutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap BlissCanonicalDigraph(x::GapObj)::GapObj
GAP.@wrap BlissCanonicalLabelling(x::GapObj)::GapObj
GAP.@wrap ChromaticNumber(x::GapObj)::GapInt
GAP.@wrap CliqueNumber(x::GapObj)::GapInt
GAP.@wrap CompleteBipartiteDigraph(x::Int, y::Int)::GapObj
GAP.@wrap CompleteDigraph(x::Int)::GapObj
GAP.@wrap CycleDigraph(x::Int)::GapObj
GAP.@wrap DigraphBicomponents(x::GapObj)::GapObj
GAP.@wrap DigraphByAdjacencyMatrix(x::GapObj)::GapObj
GAP.@wrap DigraphByEdges(x::GapObj)::GapObj
GAP.@wrap DigraphCartesianProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphConnectedComponents(x::GapObj)::GapObj
GAP.@wrap DigraphCons(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphDiameter(x::GapObj)::GapInt
GAP.@wrap DigraphDirectProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphDisjointUnion(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphDual(x::GapObj)::GapObj
GAP.@wrap DigraphEdgeUnion(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphGroup(x::GapObj)::GapObj
GAP.@wrap DigraphJoin(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphLexicographicProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphMaximalCliques(x::GapObj)::GapObj
GAP.@wrap DigraphNrEdges(x::GapObj)::GapInt
GAP.@wrap DigraphNrVertices(x::GapObj)::GapInt
GAP.@wrap DigraphReverse(x::GapObj)::GapObj
GAP.@wrap DigraphShortestDistances(x::GapObj)::GapObj
GAP.@wrap DigraphSinks(x::GapObj)::GapObj
GAP.@wrap DigraphSources(x::GapObj)::GapObj
GAP.@wrap DigraphStronglyConnectedComponents(x::GapObj)::GapObj
GAP.@wrap DigraphTopologicalSort(x::GapObj)::GapObj
GAP.@wrap Graph6String(x::GapObj)::GapObj
GAP.@wrap InDegrees(x::GapObj)::GapObj
GAP.@wrap InducedSubdigraph(x::GapObj, y::GapObj)::GapObj
GAP.@wrap InNeighbours(x::GapObj)::GapObj
GAP.@wrap IsAcyclicDigraph(x::GapObj)::Bool
GAP.@wrap IsAntisymmetricDigraph(x::GapObj)::Bool
GAP.@wrap IsBiconnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsBipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsBridgelessDigraph(x::GapObj)::Bool
GAP.@wrap IsCayleyDigraph(x::GapObj)::Bool
GAP.@wrap IsCompleteBipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsCompleteDigraph(x::GapObj)::Bool
GAP.@wrap IsConnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsDigraph(x::GapObj)::Bool
GAP.@wrap IsDirectedTree(x::GapObj)::Bool
GAP.@wrap IsEdgeTransitive(x::GapObj)::Bool
GAP.@wrap IsEmptyDigraph(x::GapObj)::Bool
GAP.@wrap IsEulerianDigraph(x::GapObj)::Bool
GAP.@wrap IsFunctionalDigraph(x::GapObj)::Bool
GAP.@wrap IsHamiltonianDigraph(x::GapObj)::Bool
GAP.@wrap IsImmutableDigraph(x::GapObj)::Bool
GAP.@wrap IsIsomorphicDigraph(x::GapObj, y::GapObj)::Bool
GAP.@wrap IsMultiDigraph(x::GapObj)::Bool
GAP.@wrap IsMutableDigraph(x::GapObj)::Bool
GAP.@wrap IsNullDigraph(x::GapObj)::Bool
GAP.@wrap IsomorphismDigraphs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap IsPlanarDigraph(x::GapObj)::Bool
GAP.@wrap IsReflexiveDigraph(x::GapObj)::Bool
GAP.@wrap IsRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsStronglyConnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsSymmetricDigraph(x::GapObj)::Bool
GAP.@wrap IsTournament(x::GapObj)::Bool
GAP.@wrap IsTransitiveDigraph(x::GapObj)::Bool
GAP.@wrap IsVertexTransitive(x::GapObj)::Bool
GAP.@wrap LaplacianMatrix(x::GapObj)::GapObj
GAP.@wrap MaximalSymmetricSubdigraph(x::GapObj)::GapObj
GAP.@wrap ModularProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap NullDigraph(x::Int)::GapObj
GAP.@wrap OutDegrees(x::GapObj)::GapObj
GAP.@wrap OutNeighbours(x::GapObj)::GapObj
GAP.@wrap QuotientDigraph(x::GapObj, y::GapObj)::GapObj
GAP.@wrap ReadDigraphs(x::GapObj)::GapObj
GAP.@wrap ReducedDigraph(x::GapObj)::GapObj
GAP.@wrap Sparse6String(x::GapObj)::GapObj
GAP.@wrap StrongProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK23(x::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK33(x::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK4(x::GapObj)::GapObj
GAP.@wrap WriteDigraphs(x::GapObj, y::GapObj)::GapObj

function __init__()
    GAP.Packages.load("Digraphs") || error("cannot load GAP package Digraphs")
end

end
