@testset "hartshorne example" begin
    # this is Hartshorne's example from Section 20.5 in 24HLC (24 hours of local cohohomology)
    kQ = monoid_algebra([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
    a, b, c, d = gens(kQ)

    # M = k[Q] (as a k[Q]-module)
    I_M = ideal(kQ, [])
    M = Oscar.quotient_ring_as_module(I_M)
    I = ideal(kQ, [a, b])

    # cohomological degree 0
    H0 = Oscar.zeroth_local_cohomology(M, I)
    @test is_zero(H0)

    # cohomological degrees 1 and 2 from one injective resolution
    H = local_cohomology_all(I_M, I, 2)
    H1_sectors = [h for h in sectors(H[1]) if !is_zero(h)] #sectors with non-zero local cohomology
    @test !is_zero(H[1])
    @test length(H1_sectors) == 1

    H2_sectors = [h for h in sectors(H[2]) if !is_zero(h)] #sectors with non-zero local cohomology
    @test !is_zero(H[2])
    @test length(H2_sectors) == 1

    # k[Q] is Cohen-Macaulay of dimension 3, so its Bass numbers at the maximal
    # ideal vanish below degree 3. The shift for a resolution up to degree 1
    # must nevertheless use them up to degree 1 + dim Q, otherwise J^1 is lost.
    for shift in (:bound, :helm_miller, :milp)
      res = injective_resolution(M, 1; shift)
      @test res.upto == 1
      @test length(indecomposable_injectives(injective_modules(res)[2])) == 4
    end
    @test is_minimal(injective_resolution(M, 1; shift=:helm_miller))
end
