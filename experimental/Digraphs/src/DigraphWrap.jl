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
GAP.@wrap AdjacencyMatrixMutableCopy(x::GapObj)::GapObj
GAP.@wrap AndrasfaiGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap ArticulationPoints(x::GapObj)::GapObj
GAP.@wrap AsBinaryRelation(d::GapObj)::GapObj
GAP.@wrap AsDigraph(filt::GapObj, d::GapObj)::GapObj
GAP.@wrap AsDigraph(filt::GapObj, f::GapObj, n::Int64)::GapObj
GAP.@wrap AsGraph(d::GapObj)::GapObj
GAP.@wrap AsMonoid(filt::GapObj, digraph::GapObj)::GapObj
GAP.@wrap AsSemigroup(filt::GapObj, digraph::GapObj)::GapObj
GAP.@wrap AsSemigroup(filt::GapObj, Y::GapObj, gps::GapObj, homs::GapObj)::GapObj
GAP.@wrap AsTransformation(d::GapObj)::GapObj
GAP.@wrap AutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap AutomorphismGroup(x::GapObj, colours::GapObj)::GapObj
GAP.@wrap AutomorphismGroup(x::GapObj, vert_colours::GapObj, edge_colours::GapObj)::GapObj
GAP.@wrap BananaTree(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap BinaryTree(filt::GapObj, n::Int)::GapObj
GAP.@wrap BinomialTreeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap BipartiteDoubleDigraph(x::GapObj)::GapObj
GAP.@wrap BishopsGraph(filt::GapObj, color::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap BishopsGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap BlissAutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap BlissAutomorphismGroup(x::GapObj, colours::GapObj)::GapObj
GAP.@wrap BlissAutomorphismGroup(x::GapObj, vert_colours::GapObj, edge_colours::GapObj)::GapObj
GAP.@wrap BlissCanonicalDigraph(x::GapObj)::GapObj
GAP.@wrap BlissCanonicalLabelling(x::GapObj)::GapObj
GAP.@wrap BlissCanonicalLabelling(x::GapObj, colours::GapObj)::GapObj
GAP.@wrap BondyGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap BookGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap BooleanAdjacencyMatrix(x::GapObj)::GapObj
GAP.@wrap BooleanAdjacencyMatrixMutableCopy(x::GapObj)::GapObj
GAP.@wrap Bridges(x::GapObj)::GapObj
GAP.@wrap BurntPancakeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap CayleyDigraph(filt::GapObj, G::GapObj)::GapObj
GAP.@wrap CayleyDigraph(filt::GapObj, G::GapObj, gens::GapObj)::GapObj
GAP.@wrap ChainDigraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap CharacteristicPolynomial(x::GapObj)::GapObj
GAP.@wrap ChromaticNumber(x::GapObj)::GapInt
GAP.@wrap CirculantGraph(filt::GapObj, n::Int, par::GapObj)::GapObj
GAP.@wrap CliqueNumber(x::GapObj)::GapInt
GAP.@wrap CompleteBipartiteDigraph(filt::GapObj, x::Int, y::Int)::GapObj
GAP.@wrap CompleteDigraph(filt::GapObj, x::Int)::GapObj
GAP.@wrap CompleteMultipartiteDigraph(filt::GapObj, orders::GapObj)::GapObj
GAP.@wrap ConormalProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap CycleDigraph(filt::GapObj, x::Int)::GapObj
GAP.@wrap CycleDigraph(x::Int)::GapObj
GAP.@wrap CycleGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap DegreeMatrix(x::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, adj::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, G::GapObj, list::GapObj, act::GapObj, rel::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, labels::GapObj, source::GapObj, range::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, list::GapObj, func::GapObj)::GapObj
GAP.@wrap Digraph(filt::GapObj, n::Int64, source::GapObj, range::GapObj)::GapObj
GAP.@wrap Digraph(obj::GapObj)::GapObj
GAP.@wrap DigraphAbsorptionExpectedSteps(x::GapObj)::GapObj
GAP.@wrap DigraphAbsorptionProbabilities(x::GapObj)::GapObj
GAP.@wrap DigraphAddAllLoops(x::GapObj)::GapObj
GAP.@wrap DigraphAddAllLoopsAttr(x::GapObj)::GapObj
GAP.@wrap DigraphAddEdge(x::GapObj, e::GapObj)::GapObj
GAP.@wrap DigraphAddEdge(x::GapObj, src::Int, ran::Int)::GapObj
GAP.@wrap DigraphAddEdgeOrbit(x::GapObj, e::GapObj)::GapObj
GAP.@wrap DigraphAddEdges(x::GapObj, edges::GapObj)::GapObj
GAP.@wrap DigraphAddVertex(x::GapObj)::GapObj
GAP.@wrap DigraphAddVertex(x::GapObj, label::GapObj)::GapObj
GAP.@wrap DigraphAddVertices(x::GapObj, labels::GapObj)::GapObj
GAP.@wrap DigraphAddVertices(x::GapObj, m::Int)::GapObj
GAP.@wrap DigraphAdjacencyFunction(x::GapObj)::GapObj
GAP.@wrap DigraphAllChordlessCycles(x::GapObj)::GapObj
GAP.@wrap DigraphAllChordlessCyclesOfMaximalLength(x::GapObj, maxLength::Int)::GapObj
GAP.@wrap DigraphAllSimpleCircuits(x::GapObj)::GapObj
GAP.@wrap DigraphAllUndirectedSimpleCircuits(x::GapObj)::GapObj
GAP.@wrap DigraphBicomponents(x::GapObj)::GapObj
GAP.@wrap DigraphByAdjacencyMatrix(filt::GapObj, x::GapObj)::GapObj
GAP.@wrap DigraphByAdjacencyMatrix(x::GapObj)::GapObj
GAP.@wrap DigraphByEdges(filt::GapObj, x::GapObj)::GapObj
GAP.@wrap DigraphByEdges(filt::GapObj, x::GapObj, n::Int64)::GapObj
GAP.@wrap DigraphByInNeighbours(filt::GapObj, x::GapObj)::GapObj
GAP.@wrap DigraphByInNeighbours(x::GapObj)::GapObj
GAP.@wrap DigraphCartesianProduct(list::GapObj)::GapObj
GAP.@wrap DigraphCartesianProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphCartesianProductProjections(x::GapObj)::GapObj
GAP.@wrap DigraphClosure(x::GapObj, k::Int)::GapObj
GAP.@wrap DigraphColouring(x::GapObj, n::Int)::GapObj
GAP.@wrap DigraphConnectedComponent(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphConnectedComponents(x::GapObj)::GapObj
GAP.@wrap DigraphCons(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphContractEdge(x::GapObj, e::GapObj)::GapObj
GAP.@wrap DigraphContractEdge(x::GapObj, src::Int, ran::Int)::GapObj
GAP.@wrap DigraphCopy(x::GapObj)::GapObj
GAP.@wrap DigraphCopySameMutability(d::GapObj)::GapObj
GAP.@wrap DigraphCore(x::GapObj)::GapObj
GAP.@wrap DigraphCycleBasis(x::GapObj)::GapObj
GAP.@wrap DigraphDegeneracy(x::GapObj)::GapInt
GAP.@wrap DigraphDegeneracyOrdering(x::GapObj)::GapObj
GAP.@wrap DigraphDiameter(x::GapObj)::GapInt
GAP.@wrap DigraphDijkstra(x::GapObj, s::Int)::GapObj
GAP.@wrap DigraphDijkstra(x::GapObj, source::Int, target::Int)::GapObj
GAP.@wrap DigraphDirectProduct(list::GapObj)::GapObj
GAP.@wrap DigraphDirectProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphDirectProductProjections(x::GapObj)::GapObj
GAP.@wrap DigraphDisjointUnion(list::GapObj)::GapObj
GAP.@wrap DigraphDisjointUnion(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphDistanceSet(x::GapObj, v::Int, d::Int)::GapObj
GAP.@wrap DigraphDistanceSet(x::GapObj, v::Int, dists::GapObj)::GapObj
GAP.@wrap DigraphDual(x::GapObj)::GapObj
GAP.@wrap DigraphDualAttr(x::GapObj)::GapObj
GAP.@wrap DigraphEdgeLabel(x::GapObj, i::Int, j::Int)::GapInt
GAP.@wrap DigraphEdgeLabels(x::GapObj)::Any
GAP.@wrap DigraphEdges(x::GapObj)::GapObj
GAP.@wrap DigraphEdgeUnion(list::GapObj)::GapObj
GAP.@wrap DigraphEdgeUnion(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphEmbedding(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphEpimorphism(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphFloydWarshall(x::GapObj, func::Any, nopath::Any, edge::Any)::GapObj
GAP.@wrap DigraphGirth(x::GapObj)::GapInt
GAP.@wrap DigraphGreedyColouring(x::GapObj)::GapObj
GAP.@wrap DigraphGreedyColouring(x::GapObj, order::GapObj)::GapObj
GAP.@wrap DigraphGroup(x::GapObj)::GapObj
GAP.@wrap DigraphHasAVertex(x::GapObj)::Bool
GAP.@wrap DigraphHash(x::GapObj)::GapInt
GAP.@wrap DigraphHasLoops(x::GapObj)::Bool
GAP.@wrap DigraphHasNoVertices(x::GapObj)::Bool
GAP.@wrap DigraphHomomorphism(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphImmutableCopy(d::GapObj)::GapObj
GAP.@wrap DigraphImmutableCopyIfImmutable(d::GapObj)::GapObj
GAP.@wrap DigraphImmutableCopyIfMutable(x::GapObj)::GapObj
GAP.@wrap DigraphInEdges(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphIsKing(x::GapObj, v::Int, k::Int)::Bool
GAP.@wrap DigraphJoin(list::GapObj)::GapObj
GAP.@wrap DigraphJoin(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphJoinTable(x::GapObj)::GapObj
GAP.@wrap DigraphKings(x::GapObj, n::Int)::GapObj
GAP.@wrap DigraphLayers(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphLexicographicProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphLongestDistanceFromVertex(x::GapObj, v::Int)::GapInt
GAP.@wrap DigraphLongestSimpleCircuit(x::GapObj)::GapObj
GAP.@wrap DigraphLoops(x::GapObj)::GapObj
GAP.@wrap DigraphMaximalCliques(x::GapObj)::GapObj
GAP.@wrap DigraphMaximalMatching(x::GapObj)::GapObj
GAP.@wrap DigraphMaximumFlow(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphMaximumMatching(x::GapObj)::GapObj
GAP.@wrap DigraphMeetTable(x::GapObj)::GapObj
GAP.@wrap DigraphMinimumCut(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphMinimumCutSet(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphMonomorphism(x::GapObj, y::GapObj)::GapObj
GAP.@wrap DigraphMutableCopy(d::GapObj)::GapObj
GAP.@wrap DigraphMutableCopyIfImmutable(x::GapObj)::GapObj
GAP.@wrap DigraphMutableCopyIfMutable(x::GapObj)::GapObj
GAP.@wrap DigraphMycielskian(x::GapObj)::GapObj
GAP.@wrap DigraphMycielskianAttr(x::GapObj)::GapObj
GAP.@wrap DigraphNrAdjacencies(x::GapObj)::GapInt
GAP.@wrap DigraphNrAdjacenciesWithoutLoops(x::GapObj)::GapInt
GAP.@wrap DigraphNrConnectedComponents(x::GapObj)::GapInt
GAP.@wrap DigraphNrEdges(x::GapObj)::GapInt
GAP.@wrap DigraphNrLoops(x::GapObj)::GapInt
GAP.@wrap DigraphNrStronglyConnectedComponents(x::GapObj)::GapInt
GAP.@wrap DigraphNrVertices(x::GapObj)::GapInt
GAP.@wrap DigraphOddGirth(x::GapObj)::GapInt
GAP.@wrap DigraphOrbitReps(x::GapObj)::GapObj
GAP.@wrap DigraphOrbits(x::GapObj)::GapObj
GAP.@wrap DigraphOutEdges(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphPath(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphPeriod(x::GapObj)::GapInt
GAP.@wrap DigraphRandomWalk(x::GapObj, v::Int, t::Int)::GapObj
GAP.@wrap DigraphReflexiveTransitiveClosure(x::GapObj)::GapObj
GAP.@wrap DigraphReflexiveTransitiveClosureAttr(x::GapObj)::GapObj
GAP.@wrap DigraphReflexiveTransitiveReduction(x::GapObj)::GapObj
GAP.@wrap DigraphReflexiveTransitiveReductionAttr(x::GapObj)::GapObj
GAP.@wrap DigraphRemoveAllMultipleEdges(x::GapObj)::GapObj
GAP.@wrap DigraphRemoveAllMultipleEdgesAttr(x::GapObj)::GapObj
GAP.@wrap DigraphRemoveEdge(x::GapObj, e::GapObj)::GapObj
GAP.@wrap DigraphRemoveEdge(x::GapObj, src::Int, ran::Int)::GapObj
GAP.@wrap DigraphRemoveEdgeOrbit(x::GapObj, e::GapObj)::GapObj
GAP.@wrap DigraphRemoveEdges(x::GapObj, edges::GapObj)::GapObj
GAP.@wrap DigraphRemoveLoops(x::GapObj)::GapObj
GAP.@wrap DigraphRemoveLoopsAttr(x::GapObj)::GapObj
GAP.@wrap DigraphRemoveVertex(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphRemoveVertices(x::GapObj, verts::GapObj)::GapObj
GAP.@wrap DigraphReverse(x::GapObj)::GapObj
GAP.@wrap DigraphReverseAttr(x::GapObj)::GapObj
GAP.@wrap DigraphReverseEdge(x::GapObj, e::GapObj)::GapObj
GAP.@wrap DigraphReverseEdge(x::GapObj, src::Int, ran::Int)::GapObj
GAP.@wrap DigraphReverseEdges(x::GapObj, edges::GapObj)::GapObj
GAP.@wrap DigraphSchreierVector(x::GapObj)::GapObj
GAP.@wrap DigraphShortestDistance(x::GapObj, list1::GapObj, list2::GapObj)::GapInt
GAP.@wrap DigraphShortestDistance(x::GapObj, list::GapObj)::GapInt
GAP.@wrap DigraphShortestDistance(x::GapObj, s::Int, t::Int)::GapInt
GAP.@wrap DigraphShortestDistances(x::GapObj)::GapObj
GAP.@wrap DigraphShortestPath(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap DigraphShortestPathSpanningTree(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphSinks(x::GapObj)::GapObj
GAP.@wrap DigraphSources(x::GapObj)::GapObj
GAP.@wrap DigraphStabilizer(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphStronglyConnectedComponent(x::GapObj, v::Int)::GapObj
GAP.@wrap DigraphStronglyConnectedComponents(x::GapObj)::GapObj
GAP.@wrap DigraphSymmetricClosure(x::GapObj)::GapObj
GAP.@wrap DigraphSymmetricClosureAttr(x::GapObj)::GapObj
GAP.@wrap DigraphsRespectsColouring(src::GapObj, ran::GapObj, x::GapObj, col1::GapObj, col2::GapObj)::Bool
GAP.@wrap DigraphsUseBliss()::Nothing
GAP.@wrap DigraphsUseNauty()::Nothing
GAP.@wrap DigraphTopologicalSort(x::GapObj)::GapObj
GAP.@wrap DigraphTransitiveClosure(x::GapObj)::GapObj
GAP.@wrap DigraphTransitiveClosureAttr(x::GapObj)::GapObj
GAP.@wrap DigraphTransitiveReduction(x::GapObj)::GapObj
GAP.@wrap DigraphTransitiveReductionAttr(x::GapObj)::GapObj
GAP.@wrap DigraphUndirectedGirth(x::GapObj)::GapInt
GAP.@wrap DigraphVertexConnectivity(x::GapObj)::GapInt
GAP.@wrap DigraphVertexLabel(x::GapObj, i::Int)::GapInt
GAP.@wrap DigraphVertexLabels(x::GapObj)::Any
GAP.@wrap DigraphVertices(x::GapObj)::GapObj
GAP.@wrap DigraphWelshPowellOrder(x::GapObj)::GapObj
GAP.@wrap DistanceDigraph(x::GapObj, i::Int)::GapObj
GAP.@wrap DistanceDigraph(x::GapObj, list::GapObj)::GapObj
GAP.@wrap Dominators(x::GapObj, root::Int)::GapObj
GAP.@wrap DominatorTree(x::GapObj, root::Int)::GapObj
GAP.@wrap DoubleDigraph(x::GapObj)::GapObj
GAP.@wrap DualPlanarGraph(x::GapObj)::GapObj
GAP.@wrap EdgeDigraph(x::GapObj)::GapObj
GAP.@wrap EdgeOrbitsDigraph(G::GapObj, edges::GapObj)::GapObj
GAP.@wrap EdgeOrbitsDigraph(G::GapObj, edges::GapObj, n::Int)::GapObj
GAP.@wrap EdgeUndirectedDigraph(x::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraph(D::GapObj, weights::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraphMinimumSpanningTree(x::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraphShortestPath(x::GapObj, s::Int, t::Int)::GapObj
GAP.@wrap EdgeWeightedDigraphShortestPaths(x::GapObj)::GapObj
GAP.@wrap EdgeWeightedDigraphShortestPaths(x::GapObj, s::Int)::GapObj
GAP.@wrap EdgeWeightedDigraphTotalWeight(x::GapObj)::Any
GAP.@wrap EdgeWeights(x::GapObj)::GapObj
GAP.@wrap EdgeWeightsMutableCopy(x::GapObj)::GapObj
GAP.@wrap EmbeddingsDigraphs(D1::GapObj, D2::GapObj)::GapObj
GAP.@wrap EmbeddingsDigraphsRepresentatives(D1::GapObj, D2::GapObj)::GapObj
GAP.@wrap EmptyDigraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap EpimorphismsDigraphs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap EpimorphismsDigraphsRepresentatives(x::GapObj, y::GapObj)::GapObj
GAP.@wrap FacialWalks(x::GapObj, list::GapObj)::GapObj
GAP.@wrap GearGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap GeneralisedPetersenGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap GeneratorsOfCayleyDigraph(x::GapObj)::GapObj
GAP.@wrap GeneratorsOfEndomorphismMonoid(x::GapObj)::GapObj
GAP.@wrap GeneratorsOfEndomorphismMonoid(x::GapObj, colours::GapObj)::GapObj
GAP.@wrap GeneratorsOfEndomorphismMonoid(x::GapObj, colours::GapObj, limit::Int)::GapObj
GAP.@wrap GeneratorsOfEndomorphismMonoidAttr(x::GapObj)::GapObj
GAP.@wrap Graph(d::GapObj)::GapObj
GAP.@wrap Graph6String(x::GapObj)::GapObj
GAP.@wrap GroupOfCayleyDigraph(x::GapObj)::GapObj
GAP.@wrap HaarGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HalvedCubeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HamiltonianPath(x::GapObj)::GapObj
GAP.@wrap HanoiGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HelmGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap HomomorphicProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap HomomorphismDigraphsFinder(D1::GapObj, D2::GapObj, hook::GapObj, user_param::GapObj, max_results, hint, injective::Int, image::GapObj, partial_map::GapObj, colors1::GapObj, colors2::GapObj)::GapObj
GAP.@wrap HomomorphismDigraphsFinder(D1::GapObj, D2::GapObj, hook::GapObj, user_param::GapObj, max_results, hint, injective::Int, image::GapObj, partial_map::GapObj, colors1::GapObj, colors2::GapObj, order::GapObj)::GapObj
GAP.@wrap HomomorphismDigraphsFinder(D1::GapObj, D2::GapObj, hook::GapObj, user_param::GapObj, max_results, hint, injective::Int, image::GapObj, partial_map::GapObj, colors1::GapObj, colors2::GapObj, order::GapObj, aut_grp::GapObj)::GapObj
GAP.@wrap HomomorphismsDigraphs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap HomomorphismsDigraphsRepresentatives(x::GapObj, y::GapObj)::GapObj
GAP.@wrap HypercubeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap InDegreeOfVertex(x::GapObj, v::Int)::GapInt
GAP.@wrap InDegrees(x::GapObj)::GapObj
GAP.@wrap InDegreeSequence(x::GapObj)::GapObj
GAP.@wrap InDegreeSet(x::GapObj)::GapObj
GAP.@wrap InducedSubdigraph(x::GapObj, y::GapObj)::GapObj
GAP.@wrap InNeighbors(x::GapObj)::GapObj
GAP.@wrap InNeighborsMutableCopy(x::GapObj)::GapObj
GAP.@wrap InNeighborsOfVertex(x::GapObj, v::Int)::GapObj
GAP.@wrap InNeighbours(x::GapObj)::GapObj
GAP.@wrap InNeighboursMutableCopy(x::GapObj)::GapObj
GAP.@wrap InNeighboursOfVertex(x::GapObj, v::Int)::GapObj
GAP.@wrap Is2EdgeTransitive(x::GapObj)::Bool
GAP.@wrap IsAcyclicDigraph(x::GapObj)::Bool
GAP.@wrap IsAntisymmetricDigraph(x::GapObj)::Bool
GAP.@wrap IsAperiodicDigraph(x::GapObj)::Bool
GAP.@wrap IsBiconnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsBipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsBridgelessDigraph(x::GapObj)::Bool
GAP.@wrap IsCayleyDigraph(x::GapObj)::Bool
GAP.@wrap IsChainDigraph(x::GapObj)::Bool
GAP.@wrap IsCompleteBipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsCompleteDigraph(x::GapObj)::Bool
GAP.@wrap IsCompleteMultipartiteDigraph(x::GapObj)::Bool
GAP.@wrap IsConnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsCycleDigraph(x::GapObj)::Bool
GAP.@wrap IsDigraph(x::GapObj)::Bool
GAP.@wrap IsDigraphAutomorphism(d::GapObj, x::GapObj)::Bool
GAP.@wrap IsDigraphAutomorphism(d::GapObj, x::GapObj, col::GapObj)::Bool
GAP.@wrap IsDigraphColouring(d::GapObj, list::GapObj)::Bool
GAP.@wrap IsDigraphCore(x::GapObj)::Bool
GAP.@wrap IsDigraphEdge(x::GapObj, e::GapObj)::Bool
GAP.@wrap IsDigraphEdge(x::GapObj, u::Int, v::Int)::Bool
GAP.@wrap IsDigraphEmbedding(src::GapObj, ran::GapObj, x::GapObj)::Bool
GAP.@wrap IsDigraphEmbedding(src::GapObj, ran::GapObj, x::GapObj, col1::GapObj, col2::GapObj)::Bool
GAP.@wrap IsDigraphEndomorphism(d::GapObj, x::GapObj)::Bool
GAP.@wrap IsDigraphEndomorphism(d::GapObj, x::GapObj, col::GapObj)::Bool
GAP.@wrap IsDigraphEpimorphism(src::GapObj, ran::GapObj, x::GapObj)::Bool
GAP.@wrap IsDigraphEpimorphism(src::GapObj, ran::GapObj, x::GapObj, col1::GapObj, col2::GapObj)::Bool
GAP.@wrap IsDigraphHomomorphism(src::GapObj, ran::GapObj, x::GapObj)::Bool
GAP.@wrap IsDigraphHomomorphism(src::GapObj, ran::GapObj, x::GapObj, col1::GapObj, col2::GapObj)::Bool
GAP.@wrap IsDigraphIsomorphism(src::GapObj, ran::GapObj, x::GapObj)::Bool
GAP.@wrap IsDigraphIsomorphism(src::GapObj, ran::GapObj, x::GapObj, col1::GapObj, col2::GapObj)::Bool
GAP.@wrap IsDigraphMonomorphism(src::GapObj, ran::GapObj, x::GapObj)::Bool
GAP.@wrap IsDigraphMonomorphism(src::GapObj, ran::GapObj, x::GapObj, col1::GapObj, col2::GapObj)::Bool
GAP.@wrap IsDigraphPath(x::GapObj, list::GapObj)::Bool
GAP.@wrap IsDigraphPath(x::GapObj, v::GapObj, a::GapObj)::Bool
GAP.@wrap IsDigraphWithAdjacencyFunction(x::GapObj)::Bool
GAP.@wrap IsDirectedForest(x::GapObj)::Bool
GAP.@wrap IsDirectedTree(x::GapObj)::Bool
GAP.@wrap IsDistanceRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsDistributiveLatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsEdgeTransitive(x::GapObj)::Bool
GAP.@wrap IsEmptyDigraph(x::GapObj)::Bool
GAP.@wrap IsEquivalenceDigraph(x::GapObj)::Bool
GAP.@wrap IsEulerianDigraph(x::GapObj)::Bool
GAP.@wrap IsFunctionalDigraph(x::GapObj)::Bool
GAP.@wrap IsHamiltonianDigraph(x::GapObj)::Bool
GAP.@wrap IsImmutableDigraph(x::GapObj)::Bool
GAP.@wrap IsInRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsIsomorphicDigraph(x::GapObj, y::GapObj)::Bool
GAP.@wrap IsIsomorphicDigraph(x::GapObj, y::GapObj, col1::GapObj, col2::GapObj)::Bool
GAP.@wrap IsJoinSemilatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsLatticeEmbedding(L1::GapObj, L2::GapObj, map::GapObj)::Bool
GAP.@wrap IsLatticeEndomorphism(L::GapObj, map::GapObj)::Bool
GAP.@wrap IsLatticeEpimorphism(L1::GapObj, L2::GapObj, map::GapObj)::Bool
GAP.@wrap IsLatticeHomomorphism(L1::GapObj, L2::GapObj, map::GapObj)::Bool
GAP.@wrap IsLatticeMonomorphism(L1::GapObj, L2::GapObj, map::GapObj)::Bool
GAP.@wrap IsLowerSemimodularDigraph(x::GapObj)::Bool
GAP.@wrap IsMatching(x::GapObj, list::GapObj)::Bool
GAP.@wrap IsMaximalMatching(x::GapObj, list::GapObj)::Bool
GAP.@wrap IsMaximumMatching(x::GapObj, list::GapObj)::Bool
GAP.@wrap IsMeetSemilatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsModularLatticeDigraph(x::GapObj)::Bool
GAP.@wrap IsMultiDigraph(x::GapObj)::Bool
GAP.@wrap IsMutableDigraph(x::GapObj)::Bool
GAP.@wrap IsNegativeEdgeWeightedDigraph(x::GapObj)::Bool
GAP.@wrap IsNonemptyDigraph(x::GapObj)::Bool
GAP.@wrap IsNullDigraph(x::GapObj)::Bool
GAP.@wrap IsomorphismDigraphs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap IsomorphismDigraphs(x::GapObj, y::GapObj, col1::GapObj, col2::GapObj)::GapObj
GAP.@wrap IsOrderFilter(x::GapObj, subset::GapObj)::Bool
GAP.@wrap IsOrderIdeal(x::GapObj, subset::GapObj)::Bool
GAP.@wrap IsOuterPlanarDigraph(x::GapObj)::Bool
GAP.@wrap IsOutRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsPartialOrderDigraph(x::GapObj)::Bool
GAP.@wrap IsPerfectMatching(x::GapObj, list::GapObj)::Bool
GAP.@wrap IsPermutationDigraph(x::GapObj)::Bool
GAP.@wrap IsPlanarDigraph(x::GapObj)::Bool
GAP.@wrap IsPreorderDigraph(x::GapObj)::Bool
GAP.@wrap IsQuasiorderDigraph(x::GapObj)::Bool
GAP.@wrap IsReachable(x::GapObj, s::Int, t::Int)::Bool
GAP.@wrap IsReflexiveDigraph(x::GapObj)::Bool
GAP.@wrap IsRegularDigraph(x::GapObj)::Bool
GAP.@wrap IsStronglyConnectedDigraph(x::GapObj)::Bool
GAP.@wrap IsSubdigraph(x::GapObj, y::GapObj)::Bool
GAP.@wrap IsSymmetricDigraph(x::GapObj)::Bool
GAP.@wrap IsTournament(x::GapObj)::Bool
GAP.@wrap IsTransitiveDigraph(x::GapObj)::Bool
GAP.@wrap IsUndirectedForest(x::GapObj)::Bool
GAP.@wrap IsUndirectedSpanningForest(x::GapObj, y::GapObj)::Bool
GAP.@wrap IsUndirectedSpanningTree(x::GapObj, y::GapObj)::Bool
GAP.@wrap IsUndirectedTree(x::GapObj)::Bool
GAP.@wrap IsUpperSemimodularDigraph(x::GapObj)::Bool
GAP.@wrap IsVertexTransitive(x::GapObj)::Bool
GAP.@wrap IteratorOfPaths(x::GapObj, u::Int, v::Int)::GapObj
GAP.@wrap JohnsonDigraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap KellerGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap KingsGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap KneserGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap KnightsGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap KuratowskiOuterPlanarSubdigraph(x::GapObj)::GapObj
GAP.@wrap KuratowskiPlanarSubdigraph(x::GapObj)::GapObj
GAP.@wrap LaplacianMatrix(x::GapObj)::GapObj
GAP.@wrap LatticeDigraphEmbedding(L1::GapObj, L2::GapObj)::GapObj
GAP.@wrap LexicographicProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap LindgrenSousselierGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap LineDigraph(x::GapObj)::GapObj
GAP.@wrap LineUndirectedDigraph(x::GapObj)::GapObj
GAP.@wrap ListNamedDigraphs(s::GapObj)::GapObj
GAP.@wrap ListNamedDigraphs(s::GapObj, level::Int)::GapObj
GAP.@wrap LollipopGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap MaximalAntiSymmetricSubdigraph(d::GapObj)::GapObj
GAP.@wrap MaximalAntiSymmetricSubdigraphAttr(x::GapObj)::GapObj
GAP.@wrap MaximalCommonSubdigraph(D1::GapObj, D2::GapObj)::GapObj
GAP.@wrap MaximalSymmetricSubdigraph(x::GapObj)::GapObj
GAP.@wrap MaximalSymmetricSubdigraphAttr(d::GapObj)::GapObj
GAP.@wrap MaximalSymmetricSubdigraphWithoutLoops(d::GapObj)::GapObj
GAP.@wrap MaximalSymmetricSubdigraphWithoutLoopsAttr(d::GapObj)::GapObj
GAP.@wrap MinimalCommonSuperdigraph(D1::GapObj, D2::GapObj)::GapObj
GAP.@wrap MinimalCyclicEdgeCut(x::GapObj)::GapObj
GAP.@wrap MobiusLadderGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap ModularProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap MonomorphismsDigraphs(x::GapObj, y::GapObj)::GapObj
GAP.@wrap MonomorphismsDigraphsRepresentatives(x::GapObj, y::GapObj)::GapObj
GAP.@wrap MycielskiGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap NautyAutomorphismGroup(x::GapObj)::GapObj
GAP.@wrap NautyAutomorphismGroup(x::GapObj, colours::GapObj)::GapObj
GAP.@wrap NautyCanonicalDigraph(x::GapObj)::GapObj
GAP.@wrap NautyCanonicalLabelling(x::GapObj)::GapObj
GAP.@wrap NautyCanonicalLabelling(x::GapObj, colours::GapObj)::GapObj
GAP.@wrap NonLowerSemimodularPair(x::GapObj)::GapObj
GAP.@wrap NonUpperSemimodularPair(x::GapObj)::GapObj
GAP.@wrap NrSpanningTrees(x::GapObj)::GapInt
GAP.@wrap NullDigraph(filt::GapObj, x::Int)::GapObj
GAP.@wrap OddGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap OnDigraphs(x::GapObj, p::GapObj)::GapObj
GAP.@wrap OnMultiDigraphs(x::GapObj, pair::GapObj)::GapObj
GAP.@wrap OnMultiDigraphs(x::GapObj, p1::GapObj, p2::GapObj)::GapObj
GAP.@wrap OnSetsDigraphs(list::GapObj, p::GapObj)::GapObj
GAP.@wrap OnTuplesDigraphs(list::GapObj, p::GapObj)::GapObj
GAP.@wrap OutDegreeOfVertex(x::GapObj, v::Int)::GapInt
GAP.@wrap OutDegrees(x::GapObj)::GapObj
GAP.@wrap OutDegreeSequence(x::GapObj)::GapObj
GAP.@wrap OutDegreeSet(x::GapObj)::GapObj
GAP.@wrap OuterPlanarEmbedding(x::GapObj)::GapObj
GAP.@wrap OutNeighbors(x::GapObj)::GapObj
GAP.@wrap OutNeighborsMutableCopy(x::GapObj)::GapObj
GAP.@wrap OutNeighborsOfVertex(x::GapObj, v::Int)::GapObj
GAP.@wrap OutNeighbours(x::GapObj)::GapObj
GAP.@wrap OutNeighboursMutableCopy(x::GapObj)::GapObj
GAP.@wrap OutNeighboursOfVertex(x::GapObj, v::Int)::GapObj
GAP.@wrap PancakeGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap PartialOrderDigraphJoinOfVertices(x::GapObj, u::Int, v::Int)::GapInt
GAP.@wrap PartialOrderDigraphMeetOfVertices(x::GapObj, u::Int, v::Int)::GapInt
GAP.@wrap PathGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap PermutationStarGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap PetersenGraph(filt::GapObj)::GapObj
GAP.@wrap PlanarEmbedding(x::GapObj)::GapObj
GAP.@wrap PrismGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap QueensGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap QuotientDigraph(x::GapObj, y::GapObj)::GapObj
GAP.@wrap RandomDigraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap RandomDigraph(filt::GapObj, n::Int, p::GapObj)::GapObj
GAP.@wrap RandomDigraph(n::Int)::GapObj
GAP.@wrap RandomDigraph(n::Int, p::GapObj)::GapObj
GAP.@wrap RandomLattice(filt::GapObj, n::Int)::GapObj
GAP.@wrap RandomLattice(n::Int)::GapObj
GAP.@wrap RandomMultiDigraph(n::Int)::GapObj
GAP.@wrap RandomMultiDigraph(n::Int, m::Int)::GapObj
GAP.@wrap RandomTournament(filt::GapObj, n::Int)::GapObj
GAP.@wrap RandomTournament(n::Int)::GapObj
GAP.@wrap RandomUniqueEdgeWeightedDigraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap RandomUniqueEdgeWeightedDigraph(filt::GapObj, n::Int, p::GapObj)::GapObj
GAP.@wrap ReadDigraphs(x::GapObj)::GapObj
GAP.@wrap ReducedDigraph(x::GapObj)::GapObj
GAP.@wrap ReducedDigraphAttr(x::GapObj)::GapObj
GAP.@wrap RepresentativeOutNeighbours(x::GapObj)::GapObj
GAP.@wrap RooksGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap SemigroupOfCayleyDigraph(x::GapObj)::GapObj
GAP.@wrap SetDigraphEdgeLabel(x::GapObj, i::Int, j::Int, obj::GapObj)::Nothing
GAP.@wrap SetDigraphEdgeLabels(x::GapObj, labels::GapObj)::Nothing
GAP.@wrap SetDigraphVertexLabel(x::GapObj, i::Int, obj::GapObj)::Nothing
GAP.@wrap SetDigraphVertexLabels(x::GapObj, list::GapObj)::Nothing
GAP.@wrap Sparse6String(x::GapObj)::GapObj
GAP.@wrap SquareGridGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap StackedBookGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap StackedPrismGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap StarGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap StrongOrientation(x::GapObj)::GapObj
GAP.@wrap StrongOrientationAttr(x::GapObj)::GapObj
GAP.@wrap StrongProduct(x::GapObj, y::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK23(x::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK33(x::GapObj)::GapObj
GAP.@wrap SubdigraphHomeomorphicToK4(x::GapObj)::GapObj
GAP.@wrap SubdigraphsMonomorphisms(x::GapObj, y::GapObj)::GapObj
GAP.@wrap SubdigraphsMonomorphismsRepresentatives(x::GapObj, y::GapObj)::GapObj
GAP.@wrap TadpoleGraph(filt::GapObj, m::Int, n::Int)::GapObj
GAP.@wrap TriangularGridGraph(filt::GapObj, n::Int, k::Int)::GapObj
GAP.@wrap UndirectedSpanningForest(x::GapObj)::GapObj
GAP.@wrap UndirectedSpanningForestAttr(x::GapObj)::GapObj
GAP.@wrap UndirectedSpanningTree(d::GapObj)::GapObj
GAP.@wrap UndirectedSpanningTreeAttr(x::GapObj)::GapObj
GAP.@wrap VerticesReachableFrom(x::GapObj, list::GapObj)::GapObj
GAP.@wrap VerticesReachableFrom(x::GapObj, root::Int)::GapObj
GAP.@wrap WalshHadamardGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap WebGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap WheelGraph(filt::GapObj, n::Int)::GapObj
GAP.@wrap WindmillGraph(filt::GapObj, n::Int, m::Int)::GapObj
GAP.@wrap WriteDigraphs(x::GapObj, y::GapObj)::GapObj

function __init__()
    GAP.Packages.load("Digraphs") || error("cannot load GAP package Digraphs")
end

end
