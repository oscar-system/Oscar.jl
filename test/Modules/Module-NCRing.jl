@testset "modules over nc rings" begin
    #Free Module tests
    E,x = exterior_algebra(QQ, 3)
    M = free_module(E, 3)
    @test ngens(M) == 3
    @test  gens(M)[1] in M
    @test parent(M[1]) == M
    v = [x[1], x[1] + x[2], x[1]*x[2]]
    @test v == Vector(M(v))

    #SubModuleOfFreeModule tests
    N = Oscar.SubModuleOfFreeModule(M, gens(M)[1:2])
    @test gens(N) == gens(M)[1:2]
    @test N[1] in M
    #fails! #NEEDS A 'DEFAULT ORDERING' on the PBWAlgQuo ie on E
    #@test M(v) in N
    
    #FreeModuleHom tests
    W = free_module(E,2)
    G = hom(W,M, gens(M)[1:2])
    @test ngens(image(G)[1]) == 2

    #SubQuoModule tests
    Q1 = SubquoModule(M, [x[1]*M[1]])
    Q2 = SubquoModule(M, [x[2]*M[1]])
    #fail! #NEEDS A 'DEFAULT ORDERING' on the PBWAlgQuo ie on E
    #@test !is_canonically_isomorphic(Q2,Q1)
    #simplify(Q[1])
end