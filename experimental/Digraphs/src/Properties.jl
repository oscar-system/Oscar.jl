is_digraph(::Digraph) = true

is_connected(D::Digraph) = DigraphWrap.IsConnectedDigraph(GapObj(D))::Bool

is_strongly_connected(D::Digraph) = DigraphWrap.IsStronglyConnectedDigraph(GapObj(D))::Bool

is_acyclic(D::Digraph) = DigraphWrap.IsAcyclicDigraph(GapObj(D))::Bool

is_bipartite(D::Digraph) = DigraphWrap.IsBipartiteDigraph(GapObj(D))::Bool

is_complete(D::Digraph) = DigraphWrap.IsCompleteDigraph(GapObj(D))::Bool

is_complete_bipartite(D::Digraph) = DigraphWrap.IsCompleteBipartiteDigraph(GapObj(D))::Bool

is_tournament(D::Digraph) = DigraphWrap.IsTournament(GapObj(D))::Bool

is_symmetric(D::Digraph) = DigraphWrap.IsSymmetricDigraph(GapObj(D))::Bool

is_antisymmetric(D::Digraph) = DigraphWrap.IsAntisymmetricDigraph(GapObj(D))::Bool

is_transitive(D::Digraph) = DigraphWrap.IsTransitiveDigraph(GapObj(D))::Bool

is_reflexive(D::Digraph) = DigraphWrap.IsReflexiveDigraph(GapObj(D))::Bool

is_empty(D::Digraph) = DigraphWrap.IsEmptyDigraph(GapObj(D))::Bool

is_eulerian(D::Digraph) = DigraphWrap.IsEulerianDigraph(GapObj(D))::Bool

is_hamiltonian(D::Digraph) = DigraphWrap.IsHamiltonianDigraph(GapObj(D))::Bool

is_regular(D::Digraph) = DigraphWrap.IsRegularDigraph(GapObj(D))::Bool

is_multi(D::Digraph) = DigraphWrap.IsMultiDigraph(GapObj(D))::Bool

is_functional(D::Digraph) = DigraphWrap.IsFunctionalDigraph(GapObj(D))::Bool

is_biconnected(D::Digraph) = DigraphWrap.IsBiconnectedDigraph(GapObj(D))::Bool

is_bridgeless(D::Digraph) = DigraphWrap.IsBridgelessDigraph(GapObj(D))::Bool

is_cayley(D::Digraph) = DigraphWrap.IsCayleyDigraph(GapObj(D))::Bool

is_directed_tree(D::Digraph) = DigraphWrap.IsDirectedTree(GapObj(D))::Bool

is_edge_transitive(D::Digraph) = DigraphWrap.IsEdgeTransitive(GapObj(D))::Bool

is_vertex_transitive(D::Digraph) = DigraphWrap.IsVertexTransitive(GapObj(D))::Bool

is_planar(D::Digraph) = DigraphWrap.IsPlanarDigraph(GapObj(D))::Bool
