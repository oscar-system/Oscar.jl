@testset "LieAlgebras.SimpleLieAlgebra" begin
    L = lie_algebra(QQ, "A2")
    a = L()
    @test parent(a) == L
    @test parent_type(a) == SimpleLieAlgebra{QQFieldElem}
    @test elem_type(L) == SimpleLieAlgebraElem{QQFieldElem}
    @test dim(L) == 52
    @test base_ring(L) == QQ
    @test base_ring(a) == QQ
end