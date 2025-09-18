@testset "Invariant Theory of SL_m" begin
    S, z = polynomial_ring(QQ, :z=> (1:2, 1:2))
    G = linearly_reductive_group(:SL,2,S)
    @test group_type(G) == :SL
    @test group_dim(G) == 2
    @test group_ideal(G) == ideal([z[1, 1]*z[2, 2] - z[2, 1]*z[1, 2] - 1])
    @test ncols(canonical_representation(G)) == 2
    @test natural_representation(G) == canonical_representation(G)

    rep1 = representation_on_forms(G,2)
    @test ncols(representation_matrix(rep1)) == 3
    @test group(rep1) == G
    @test vector_space_dim(rep1) == 3

    R_rep1 = invariant_ring(rep1)
    FI_rep1 = fundamental_invariants(R_rep1)
    @test length(FI_rep1) == 1

    #same rep mat as in rep1, only without multinomial coefficients.
    M = matrix(S,3,3,[z[1,1]^2 z[1,1]*z[2,1] z[2,1]^2 ; 2*z[1,1]*z[1,2] z[1,1]*z[2,2] + z[2,1]*z[1,2] 2*z[2,1]*z[2,2]; z[1,2]^2 z[1,2]*z[2,2] z[2,2]^2])
    rep2 = representation_reductive_group(G, M)
    R_rep2 = invariant_ring(R_rep1.poly_ring, rep2)
    FI_rep2 = fundamental_invariants(R_rep2) 
    x = gens(R_rep1.poly_ring)
    @test FI_rep2 == [-4*x[1]*x[3] + x[2]^2]

    #direct sum of representation_linearly_reductive_group
    D = direct_sum(rep1, rep2)
    @test ncols(representation_matrix(D)) == 6
    T1 = tensor([rep1, rep2])
    T2 = tensor(rep1, rep2)
    @test representation_matrix(T1) == representation_matrix(T2)
    @test ncols(representation_matrix(T1)) == 9

    #ternary cubics
    T, X = graded_polynomial_ring(QQ, :X=>1:10)
    g = linearly_reductive_group(:SL, 3, QQ)
    rep3 = representation_on_forms(g, 3)
    R_rep3 = invariant_ring(T, rep3)
    f = X[1]*X[4]*X[8]*X[10] - X[1]*X[4]*X[9]^2 - X[1]*X[5]*X[7]*X[10] + X[1]*X[5]*X[8]*X[9] + X[1]*X[6]*X[7]*X[9] - X[1]*X[6]*X[8]^2 - X[2]^2*X[8]*X[10] + X[2]^2*X[9]^2 + X[2]*X[3]*X[7]*X[10] - X[2]*X[3]*X[8]*X[9] + X[2]*X[4]*X[5]*X[10] - X[2]*X[4]*X[6]*X[9] - 2*X[2]*X[5]^2*X[9] + 3*X[2]*X[5]*X[6]*X[8] - X[2]*X[6]^2*X[7] - X[3]^2*X[7]*X[9] + X[3]^2*X[8]^2 - X[3]*X[4]^2*X[10] + 3*X[3]*X[4]*X[5]*X[9] - X[3]*X[4]*X[6]*X[8] - 2*X[3]*X[5]^2*X[8] + X[3]*X[5]*X[6]*X[7] + X[4]^2*X[6]^2 - 2*X[4]*X[5]^2*X[6] + X[5]^4
    @test reynolds_operator(R_rep3, f) == f

    #tori
    #example in the derksen book
    T = torus_group(QQ,1)
    r = representation_from_weights(T, [-3, -1, 1, 2])
    I = invariant_ring(r)
    f = fundamental_invariants(I)
    @test length(f) == 6

    #example from Macaulay2
    T = torus_group(QQ,2)
    r = representation_from_weights(T, [1 0; 0 1; -1 -1; -1 1])
    I = invariant_ring(r)
    R = polynomial_ring(I)
    X = gens(R)
    f = fundamental_invariants(I)
    @test f == [X[1]*X[2]*X[3], X[1]^2*X[3]*X[4]]

    #another example, with affine algebra computation
    T = torus_group(QQ,2)
    r = representation_from_weights(T, [-1 1; -1 1; 2 -2; 0 -1])
    RT = invariant_ring(r)
    A, _ = affine_algebra(RT)
    @test ngens(modulus(A)) == 1

    #no reynolds operator required/available:

    #SL(2) over symmetric forms of degree 2
    g = linearly_reductive_group(:SL,2,QQ)
    r = representation_on_forms(g,2)
    M = representation_matrix(r)
    ringg = parent(M[1,1])
    z = gens(ringg)
    f = z[1]*z[4] - z[2]*z[3] - 1
    G = linearly_reductive_group(ideal([f]))
    R = representation_reductive_group(G,M)
    RG = invariant_ring(R)
    F = fundamental_invariants(RG)
    @test length(F) == 1
    X = polynomial_ring(RG)
    g = gens(X)
    @test F[1] == g[1]*g[3] - g[2]^2

    #SL(2) over symmetric forms of degree 4:
    #Here we see that without the Reynolds operator we get many more generators.
    g = linearly_reductive_group(:SL,2,QQ)
    r = representation_on_forms(g,4)
    #with reynolds operator we get 2 fundamental invariants:
    rg = invariant_ring(r)
    FF = fundamental_invariants(rg)
    @test length(FF) == 2
    M = representation_matrix(r)
    ringg = parent(M[1,1])
    z = gens(ringg)
    f = z[1]*z[4] - z[2]*z[3] - 1
    G = linearly_reductive_group(ideal([f]))
    R = representation_reductive_group(G,M)
    RG = invariant_ring(R)
    F = fundamental_invariants(RG)
    @test length(F) == 27
end
