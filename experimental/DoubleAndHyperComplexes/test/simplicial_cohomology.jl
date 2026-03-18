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
end

@testset "integration" begin
  T = torus()
  C = Oscar.SimplicialCochainComplex(ZZ, T)
  A = Oscar.DGAlgCohRing(C)
  top = Oscar.graded_part(A, 2)
  v = first(filter(!is_zero, gens(top)))
  v = Oscar.DGAlgCohRingElem(A, 2, v)
  Oscar.set_volume_form!(A, v)
  # double assignment is forbidden, as it would create inconsistencies. 
  @test_throws ErrorException Oscar.set_volume_form!(A, v)
  
  # automated version which seeks out a single generator for the top cohomology
  A = Oscar.DGAlgCohRing(C)
  Oscar.set_volume_form!(A)
  v = Oscar.volume_form(A)
  @test is_one(integral(v))
  @test length(small_generating_set(A, 1)) == 2
  g = small_generating_set(A, 1)
  # compose the gram matrix for the middle cohomology degree
  G = matrix_space(ZZ, 2, 2)([integral(a*b) for a in g for b in g])
  @assert is_one(det(G))
end

