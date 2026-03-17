@testset "simplicial_cohomology" begin
  function ConformanceTests.generate_element(R::Oscar.DGAlgCohRing{ZZRingElem})
    n_degrees = rand(0:5)
    x = zero(R)
    for i=1:n_degrees
      x = x+Oscar.generate_homogeneous_element(R)
    end
    return x
  end

  K = torus()
  CZ = Oscar.SimplicialCochainComplex(ZZ, K)

  @test typeof(CZ) <: Oscar.SimplicialCochainComplex

  AKZ = Oscar.DGAlgCohRing(CZ)

  @test AKZ isa Oscar.DGAlgCohRing

  @test base_ring(AKZ) == ZZ

  CF_23 = Oscar.SimplicialCochainComplex(GF(23), K)

  AK_23 = Oscar.DGAlgCohRing(CF_23)

  @test AK_23 isa Oscar.DGAlgCohRing

  @test base_ring(AK_23) === GF(23)

  K = simplicial_complex([[1,2,3],[2,3,4]])
  Cmplx_p = Oscar.SimplicialCochainComplex(padic_field(7, precision = 30),K)
  R_p = Oscar.DGAlgCohRing(Cmplx_p)
  @test is_exact_type(typeof(one(R_p))) == false
  @test is_domain_type(typeof(one(R_p))) == false

  Cmplx = Oscar.SimplicialCochainComplex(ZZ,K)
  R = Oscar.DGAlgCohRing(Cmplx)
  @test is_exact_type(typeof(one(R))) == true
  @test is_domain_type(typeof(one(R))) == false

  ConformanceTests.test_NCRing_interface(R);
  #ConformanceTests.test_NCRing_interface(R_p);
end

