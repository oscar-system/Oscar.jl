@testset "modules over PBWAlgQuo" begin
    #Free Module tests
    E,x = exterior_algebra(QQ, 3)
    M = FreeMod(E, 3)
    @test ngens(M) == 3
    @test parent(M[1]) === M
    v = [x[1], x[1] + x[2], 5*x[1]*x[2]]
    @test M(2*v) == 2*M(v)

    #SubModuleOfFreeModule tests
    N = Oscar.SubModuleOfFreeModule(M, gens(M)[1:2])
    @test gens(N) == gens(M)[1:2]
    @test N[1] in M
    #fails! #NEEDS A 'DEFAULT ORDERING' on the PBWAlgQuo ie on E
    #@test M(v) in N

    #FreeModuleHom tests
    W = FreeMod(E,2)
    G = hom(W,M, gens(M)[1:2])
    @test ngens(image(G)[1]) == 2

    #SubQuoModule tests
    Q1 = SubquoModule(M, [x[1]*M[1]])
    Q2 = SubquoModule(M, [x[2]*M[1]])
    #fail! #NEEDS A 'DEFAULT ORDERING' on the PBWAlgQuo ie on E
    #@test !is_canonically_isomorphic(Q2,Q1)
    #simplify(Q[1])
end

@testset "modules over PBWAlgRing" begin
    R, (x, y, z) = QQ[:x, :y, :z]
    L = [x*y, x*z, y*z + 1]
    REL = strictly_upper_triangular_matrix(L)
    A, (x, y, z) = pbw_algebra(R, REL, deglex(gens(R)))
    M = FreeMod(A, 3)
    @test ngens(M) == 3
    @test  gens(M)[1] in M
    @test parent(M[1]) == M
    v = [x, x+y, x*y]
    @test v == Vector(M(v))

    #SubModuleOfFreeModule tests
    N = Oscar.SubModuleOfFreeModule(M, gens(M)[1:2])
    @test gens(N) == gens(M)[1:2]
    @test N[1] in M
    #fails! #NEEDS A 'DEFAULT ORDERING' on the PBWAlgQuo ie on E
    #@test M(v) in N

    #FreeModuleHom tests
    W = FreeMod(A,2)
    G = hom(W,M, gens(M)[1:2])
    @test ngens(image(G)[1]) == 2

    #SubQuoModule tests
    Q1 = SubquoModule(M, [x*M[1]])
    Q2 = SubquoModule(M, [y*M[1]])
    #fail! #NEEDS A 'DEFAULT ORDERING' on the PBWAlgQuo ie on E
    #@test !is_canonically_isomorphic(Q2,Q1)
    #simplify(Q[1])
end
