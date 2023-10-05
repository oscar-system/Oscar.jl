@testset "hypersurface" begin
  for addition in (max, min)
    T = TropicalSemiring(addition)
    R,(x,y,z) = polynomial_ring(T,3)
    b = [T(10),T(3//2095), T(-5//2)]
    m = [[1,0,1],[0,3,4],[1,0,7]]
    f = R(b,m)
    hyp = TropicalHypersurface(f)
    coeffs, exps = Oscar.homogenize_and_convert_to_pm(f)
    @test hyp isa TropicalHypersurface{addition}
    @test coeffs == get_attribute(hyp, :polymake_bigobject).COEFFICIENTS
    @test exps == matrix(ZZ, get_attribute(hyp, :polymake_bigobject).MONOMIALS)
  end
end
