@testset "Quiver Grassmannians" begin
    G = graph_from_edges(Directed, [[1,2]])
    A = matrix(QQ,[1 0 0 0; 0 1 0 0])
    Q = quiver_representation(G, [2,4], [A], QQ)
    Qsr = quiver_grassmannian(Q, [1,2])
    @test length(gens(Qsr.defining_ideal)) == 5 
end
