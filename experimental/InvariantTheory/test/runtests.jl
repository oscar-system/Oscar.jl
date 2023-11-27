@testset "Invariant Theory of SL_m" begin
    S, z = polynomial_ring(QQ, "z"=> (1:2, 1:2))
    G = reductive_group(:SL,2,S)
    rep1 = representation_on_forms(G,2)
    @test ncols(representation_matrix(rep1)) == 3

    R_rep1 = invariant_ring(rep1)
    FI_rep1 = fundamental_invariants(R_rep1)
    @test length(FI_rep1) == 1

    #same rep mat as in rep1, only without multinomial coefficients.
    M = matrix(S,3,3,[z[1,1]^2 z[1,1]*z[2,1] z_[2,1^2] ; 2*z[1,1]*z[1,2] z[1,1]*z[2,2] + z[2,1]*z[1,2] 2*z[2,1]*z[2,2]; z[1,2]^2 z[1,2]*z[2,2] z[2,2]^2])
    rep2 = representation_reductive_group(G, M)
    poly_ring_rep2,x = polynomial_ring(QQ, "x"=>1:3)
    R_rep2 = invariant_ring(rep2,poly_ring_rep2)
    FI_rep2 = fundamental_invariants(R_rep2) 
    @test FI_rep2 == -4*x[1]*x[3] + x[2]^2
    

end