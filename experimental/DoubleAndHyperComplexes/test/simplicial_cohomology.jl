using Test

@testset "simplicial_cohomology" begin
    K = torus()
    CZ = Oscar.SimplicialCoComplex(ZZ, K)

    @test typeof(CZ) <: Oscar.SimplicialCoComplex

    AKZ = Oscar.SimplicialCohomologyRing(CZ)

    @test AKZ <: Oscar.SimplicialCohomologyRing

    @test base_ring(AKZ) == ZZ

    CF_23 = Oscar.SimplicialCoComplex(GF(23), K)

    AK_23 = Oscar.SimplicialCohomologyRing(CZ)

    @test AK_23 <: Oscar.SimplicialCohomologyRing

    @test base_ring(AK_23) == GF(23)

    K = simplicial_complex([[1,2,3],[2,3,4]])
    Cmplx_p = Oscar.SimplicialCoComplex(padic_field(7, precision = 30),K)
    R_p = Oscar.SimplicialCohomologyRing(Cmplx_p)
    @test is_exact_type(typeof(one(R_p))) == false
    @test is_domain_type(typeof(one(R_p))) == false
    
    Cmplx = Oscar.SimplicialCoComplex(ZZ,K)
    R = Oscar.SimplicialCohomologyRing(Cmplx)
    @test is_exact_type(typeof(one(R))) == true
    @test is_domain_type(typeof(one(R))) == false

end

function generate_homogeneous_element(R::Oscar.SimplicialCohomologyRing{ZZRingElem})
    degree = rand(1:dim(Oscar.simplicial_complex(R))+1) # indexing starts at 1 (for degree 0)
    n_gens = rand(0:2*length(gens(Oscar.graded_parts(R)[degree])))
    x = zero(R)
    for i=1:n_gens
        n = rand(-100:100)
        x = x+n*R[degree-1,rand(1:length(gens(Oscar.graded_parts(R)[degree])))]
    end
    return x
end

function ConformanceTests.generate_element(R::Oscar.SimplicialCohomologyRing{ZZRingElem})
    n_degrees = rand(0:5)
    x = zero(R)
    for i=1:n_degrees
        x = x+Oscar.generate_homogeneous_element(R)
    end
    return x
end






