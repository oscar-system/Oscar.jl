
@testset "Test RootConversion" begin
  w_to_alpha = BasisLieHighestWeight.w_to_alpha
  alpha_to_w = BasisLieHighestWeight.alpha_to_w

  function test_inverse_alpha_w(L, weight)
    @test w_to_alpha(L, alpha_to_w(L, weight)) == weight # alpha -> w -> alpha
    @test alpha_to_w(L, w_to_alpha(L, weight)) == weight # w -> alpha -> w
  end

  @testset "Dynkin type $dynkin" for dynkin in (:A, :B, :C, :D, :E, :F, :G)
    @testset "n = $n" for n in 1:10
      if (
        !(dynkin == :B && n < 2) &&
        !(dynkin == :C && n < 2) &&
        !(dynkin == :D && n < 4) &&
        !(dynkin == :E && !(n == 6 || n == 7 || n == 8)) &&
        !(dynkin == :F && n != 4) &&
        !(dynkin == :G && (n != 2))
      )
        weight = [rand(QQ, -10:10) for _ in 1:n]
        L = BasisLieHighestWeight.lie_algebra(dynkin, n)
        test_inverse_alpha_w(L, weight)
      end
    end
  end
end
