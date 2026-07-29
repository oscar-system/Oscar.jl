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
GAP.@wrap ArticulationPoints(x::GapObj)::GapObj
GAP.@wrap AutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap BlissAutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap BlissCanonicalDigraph(x::GapObj)::GapObj
GAP.@wrap BlissCanonicalLabelling(x::GapObj)::GapObj
GAP.@wrap BooleanAdjacencyMatrix(x::GapObj)::GapObj
GAP.@wrap Bridges(x::GapObj)::GapObj
GAP.@wrap CayleyDigraph(filt::GapObj, G::GapObj, gens::GapObj)::GapObj
GAP.@wrap CayleyDigraph(filt::GapObj, G::GapObj)::GapObj
GAP.@wrap CharacteristicPolynomial(x::GapObj)::GapObj
GAP.@wrap ChromaticNumber(x::GapObj)::GapInt
GAP.@wrap CliqueNumber(x::GapObj)::GapInt
GAP.@wrap CompleteBipartiteDigraph(filt::GapObj, x::Int, y::Int)::GapObj
GAP.@wrap CompleteDigraph(filt::GapObj, x::Int)::GapObj
GAP.@wrap CycleDigraph(filt::GapObj, x::Int)::GapObj
GAP.@wrap CycleDigraph(x::Int)::GapObj
GAP.@wrap ChainDigraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap DegreeMatrix(x::GapObj)::GapObj
GAP.@wrap DigraphAllSimpleCircuits(x::GapObj)::GapObj
GAP.@wrap DigraphAbsorptionProbabilities(x::GapObj)::GapObj
GAP.@wrap DigraphAbsorptionExpectedSteps(x::GapObj)::GapObj
GAP.@wrap DigraphConnectedComponent(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphCycleBasis(x::GapObj)::GapObj
GAP.@wrap DigraphDegeneracy(x::GapObj)::GapInt
GAP.@wrap DigraphDegeneracyOrdering(x::GapObj)::GapObj
GAP.@wrap DigraphDijkstra(x::GapObj, s::Int)::GapObj
GAP.@wrap DigraphEdges(x::GapObj)::GapObj
GAP.@wrap DigraphFloydWarshall(x::GapObj)::GapObj
GAP.@wrap DigraphGirth(x::GapObj)::GapInt
GAP.@wrap DigraphOddGirth(x::GapObj)::GapInt
GAP.@wrap DigraphUndirectedGirth(x::GapObj)::GapInt
GAP.@wrap DigraphHash(x::GapObj)::GapInt
GAP.@wrap DigraphLayers(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphLoops(x::GapObj)::GapObj
GAP.@wrap DigraphLongestSimpleCircuit(x::GapObj)::GapObj
GAP.@wrap DigraphNrAdjacencies(x::GapObj)::GapInt
GAP.@wrap DigraphNrAdjacenciesWithoutLoops(x::GapObj)::GapInt
GAP.@wrap DigraphNrLoops(x::GapObj)::GapInt
GAP.@wrap DigraphPath(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphPeriod(x::GapObj)::GapInt
GAP.@wrap DigraphRandomWalk(x::GapObj, len::Int, start::Int)::GapObj
GAP.@wrap DigraphShortestDistance(x::GapObj, s::Int, t::Int)::GapInt
GAP.@wrap DigraphShortestPath(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphStronglyConnectedComponent(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphVertexConnectivity(x::GapObj)::GapInt
GAP.@wrap DigraphVertices(x::GapObj)::GapObj
GAP.@wrap EmptyDigraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap EdgeWeights(x::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraph(D::GapObj, weights::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraphMinimumSpanningTree(x::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraphShortestPaths(x::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraphShortestPaths(x::GapObj, s::Int)::GapObj
GAP.@wrap EdgeWeightedDigraphShortestPath(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap EdgeWeightedDigraphTotalWeight(x::GapObj)::GapObj
GAP.@wrap GeneratorsOfCayleyDigraph(x::GapObj)::GapObj
GAP.@wrap GroupOfCayleyDigraph(x::GapObj)::GapObj
GAP.@wrap HamiltonianPath(x::GapObj)::GapObj
GAP.@wrap IsReachable(x::GapObj, s::Int, t::Int)::Bool
GAP.@wrap JohnsonDigraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap Digraph(filt::GapObj, adj::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, G::GapObj, list::GapObj, act::GapObj, rel::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, labels::GapObj, source::GapObj, range::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, list::GapObj, func::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, n::Int64, source::GapObj, range::GapObj)::GapObj
GAP.@wrap Digraph(obj::GapObj)::GapObj
GAP.@wrap DigraphBicomponents(x::GapObj)::GapObj
GAP.@wrap DigraphByAdjacencyMatrix(filt::GapObj, x::GapObj)::GapObj
GAP.@wrap DigraphByAdjacencyMatrix(x::GapObj)::GapObj
GAP.@wrap DigraphByEdges(filt::GapObj, x::GapObj, n::Int64)::GapObj
GAP.@wrap DigraphByEdges(filt::GapObj, x::GapObj)::GapObj
GAP.@wrap DigraphByInNeighbours(filt::GapObj, x::GapObj)::GapObj
GAP.@wrap DigraphByInNeighbours(x::GapObj)::GapObj
GAP.@wrap DigraphCartesianProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphConnectedComponents(x::GapObj)::GapObj
GAP.@wrap DigraphCons(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphCore(x::GapObj)::GapObj
GAP.@wrap DigraphDiameter(x::GapObj)::GapInt
GAP.@wrap DigraphDirectProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphDisjointUnion(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphDual(x::GapObj)::GapObj
GAP.@wrap DigraphEdgeUnion(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphGroup(x::GapObj)::GapObj
GAP.@wrap DigraphHasLoops(x::GapObj)::Bool
GAP.@wrap DigraphHasNoVertices(x::GapObj)::Bool
GAP.@wrap DigraphJoin(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphJoinTable(x::GapObj)::GapObj
GAP.@wrap DigraphLexicographicProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphMaximalCliques(x::GapObj)::GapObj
GAP.@wrap DigraphMaximalMatching(x::GapObj)::GapObj
GAP.@wrap DigraphMaximumFlow(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphMaximumMatching(x::GapObj)::GapObj
GAP.@wrap DigraphMinimumCut(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphMinimumCutSet(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphMeetTable(x::GapObj)::GapObj
GAP.@wrap DigraphNrEdges(x::GapObj)::GapInt
GAP.@wrap DigraphNrConnectedComponents(x::GapObj)::GapInt
GAP.@wrap DigraphNrVertices(x::GapObj)::GapInt
GAP.@wrap DigraphNrStronglyConnectedComponents(x::GapObj)::GapInt
GAP.@wrap DigraphReverse(x::GapObj)::GapObj
GAP.@wrap DigraphShortestDistances(x::GapObj)::GapObj
GAP.@wrap DigraphSinks(x::GapObj)::GapObj
GAP.@wrap DigraphSources(x::GapObj)::GapObj
GAP.@wrap DigraphStronglyConnectedComponents(x::GapObj)::GapObj
GAP.@wrap DigraphTopologicalSort(x::GapObj)::GapObj
GAP.@wrap EdgeOrbitsDigraph(G::GapObj, edges::GapObj, n::Int)::GapObj
GAP.@wrap EdgeOrbitsDigraph(G::GapObj, edges::GapObj)::GapObj
GAP.@wrap Graph6String(x::GapObj)::GapObj
GAP.@wrap InDegrees(x::GapObj)::GapObj
GAP.@wrap InducedSubdigraph(x::GapObj, y::GapObj)::GapObj
GAP.@wrap InNeighbours(x::GapObj)::GapObj
GAP.@wrap Is2EdgeTransitive(x::GapObj)::Bool
GAP.@wrap IsAcyclicDigraph(x::GapObj)::Bool
GAP.@wrap IsAntisymmetricDigraph(x::GapObj)::Bool
GAP.@wrap IsAperiodicDigraph(x::GapObj)::Bool
GAP.@wrap IsBiconnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsBipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsBridgelessDigraph(x::GapObj)::Bool
GAP.@wrap IsCayleyDigraph(x::GapObj)::Bool
GAP.@wrap IsChainDigraph(x::GapObj)::Bool
GAP.@wrap IsCograph(x::GapObj)::Bool
GAP.@wrap IsCompleteBipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsCompleteDigraph(x::GapObj)::Bool
GAP.@wrap IsCompleteMultipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsConnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsCycleDigraph(x::GapObj)::Bool
GAP.@wrap IsDigraph(x::GapObj)::Bool
GAP.@wrap IsDigraphCore(x::GapObj)::Bool
GAP.@wrap IsDirectedForest(x::GapObj)::Bool
GAP.@wrap IsDistanceRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsDistributiveLatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsDirectedTree(x::GapObj)::Bool
GAP.@wrap IsEdgeTransitive(x::GapObj)::Bool
GAP.@wrap IsEmptyDigraph(x::GapObj)::Bool
GAP.@wrap IsEquivalenceDigraph(x::GapObj)::Bool
GAP.@wrap IsEulerianDigraph(x::GapObj)::Bool
GAP.@wrap IsFunctionalDigraph(x::GapObj)::Bool
GAP.@wrap IsHamiltonianDigraph(x::GapObj)::Bool
GAP.@wrap IsInRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsJoinSemilatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsImmutableDigraph(x::GapObj)::Bool
GAP.@wrap IsIsomorphicDigraph(x::GapObj, y::GapObj)::Bool
GAP.@wrap IsMeetSemilatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsLowerSemimodularDigraph(x::GapObj)::Bool
GAP.@wrap IsModularLatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsMultiDigraph(x::GapObj)::Bool
GAP.@wrap IsNegativeEdgeWeightedDigraph(x::GapObj)::Bool
GAP.@wrap IsNonemptyDigraph(x::GapObj)::Bool
GAP.@wrap IsMutableDigraph(x::GapObj)::Bool
GAP.@wrap IsNullDigraph(x::GapObj)::Bool
GAP.@wrap IsOuterPlanarDigraph(x::GapObj)::Bool
GAP.@wrap IsOutRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsPartialOrderDigraph(x::GapObj)::Bool
GAP.@wrap IsPermutationDigraph(x::GapObj)::Bool
GAP.@wrap IsPreorderDigraph(x::GapObj)::Bool
GAP.@wrap IsomorphismDigraphs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap IsPlanarDigraph(x::GapObj)::Bool
GAP.@wrap IsReflexiveDigraph(x::GapObj)::Bool
GAP.@wrap IsRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsStronglyConnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsSymmetricDigraph(x::GapObj)::Bool
GAP.@wrap IsTournament(x::GapObj)::Bool
GAP.@wrap IsTransitiveDigraph(x::GapObj)::Bool
GAP.@wrap IsUndirectedTree(x::GapObj)::Bool
GAP.@wrap IsUndirectedForest(x::GapObj)::Bool
GAP.@wrap IsUpperSemimodularDigraph(x::GapObj)::Bool
GAP.@wrap IsVertexTransitive(x::GapObj)::Bool
GAP.@wrap LaplacianMatrix(x::GapObj)::GapObj
GAP.@wrap MaximalSymmetricSubdigraph(x::GapObj)::GapObj
GAP.@wrap MinimalCyclicEdgeCut(x::GapObj)::GapObj
GAP.@wrap ModularProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap NrSpanningTrees(x::GapObj)::GapInt
GAP.@wrap NonLowerSemimodularPair(x::GapObj)::GapObj
GAP.@wrap NonUpperSemimodularPair(x::GapObj)::GapObj
GAP.@wrap NullDigraph(filt::GapObj, x::Int)::GapObj
GAP.@wrap OutDegrees(x::GapObj)::GapObj
GAP.@wrap OutNeighbours(x::GapObj)::GapObj
GAP.@wrap QuotientDigraph(x::GapObj, y::GapObj)::GapObj
GAP.@wrap RandomDigraph(filt::GapObj, n::Int, p::GapObj)::GapObj
GAP.@wrap RandomDigraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap RandomLattice(n::Int)::GapObj
GAP.@wrap RandomMultiDigraph(n::Int, m::Int)::GapObj
GAP.@wrap RandomMultiDigraph(n::Int)::GapObj
GAP.@wrap RandomTournament(filt::GapObj, n::Int)::GapObj
GAP.@wrap ReadDigraphs(x::GapObj)::GapObj
GAP.@wrap ReducedDigraph(x::GapObj)::GapObj
GAP.@wrap Sparse6String(x::GapObj)::GapObj
GAP.@wrap StrongProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap StrongOrientation(x::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK23(x::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK33(x::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK4(x::GapObj)::GapObj
GAP.@wrap WriteDigraphs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap AndrasfaiGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap BananaTree(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap BinaryTree(filt::GapObj, n::Int)::GapObj
GAP.@wrap BinomialTreeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap BishopsGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap BondyGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap BookGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap BurntPancakeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap CirculantGraph(filt::GapObj, n::Int, par::GapObj)::GapObj
GAP.@wrap CompleteMultipartiteDigraph(filt::GapObj, orders::GapObj)::GapObj
GAP.@wrap CycleGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap GearGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap GeneralisedPetersenGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap HaarGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HalvedCubeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HanoiGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HelmGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HypercubeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap KellerGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap KingsGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap KneserGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap KnightsGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap LindgrenSousselierGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap LollipopGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap MobiusLadderGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap MycielskiGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap OddGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap PancakeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap PathGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap PermutationStarGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap PrismGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap QueensGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap RooksGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap SquareGridGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap StackedBookGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap StackedPrismGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap StarGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap TadpoleGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap TriangularGridGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap WalshHadamardGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap WebGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap WheelGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap WindmillGraph(filt::GapObj, n::Int, m::Int)::GapObj
function __init__()
    GAP.Packages.load("Digraphs") || error("cannot load GAP package Digraphs")
end

end



