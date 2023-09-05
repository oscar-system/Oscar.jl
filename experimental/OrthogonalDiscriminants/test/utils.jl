@testset "order_omega_mod_N" begin
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 168) == (false, false)
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(6, 2, 168) == (true, false)
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 10) == (false, true)
  @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(4, 2, 12) == (true, true)
  for (d, q) in [(4, 2), (6, 3), (12, 5), (14, 7)]
    @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(d, q, 60)[2]
    @test Oscar.OrthogonalDiscriminants.order_omega_mod_N(d, q, ZZ(60))[2]
  end
end

@testset "reduce_mod_squares" begin
  F, z = cyclotomic_field(7)
  for val in [F(0), F(1), -F(5), 1 + 2*z + z^5]
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(val) == val
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(12 * val) == 3 * val
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(27 * val) == 3 * val
    @test Oscar.OrthogonalDiscriminants.reduce_mod_squares(val // 12) == 3 * val
  end
end
