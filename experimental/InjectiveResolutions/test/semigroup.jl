@testset "construct affine semigroups" begin
    #Q = ZZ^2
    Q = affine_semigroup([1 0;0 1])
    @test is_pointed(Q)
    C = cone(Q)
    @test ambient_dim(C) == ambient_dimension(Q)
    H = hyperplanes(Q)
    @test length(H) == 2
    @test dim(zonotope(Q)[1]) == rank(Q) 

    #non-normal example
    Q = affine_semigroup([2 3 0 0;0 0 2 3])
    @test is_pointed(Q)

    kQ = monoid_algebra(Q,QQ)
    @test !is_normal(kQ)


    #non-pointed example
    Q = affine_semigroup([[-1],[2]])
    @test !is_pointed(Q) 
end
