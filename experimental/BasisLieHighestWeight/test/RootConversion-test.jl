include("RootConversion-test-data.jl")

@testset "Test RootConversion" begin
  w_to_alpha = BasisLieHighestWeight.w_to_alpha
  alpha_to_w = BasisLieHighestWeight.alpha_to_w
  w_to_eps = BasisLieHighestWeight.w_to_eps
  eps_to_w = BasisLieHighestWeight.eps_to_w
  alpha_to_eps = BasisLieHighestWeight.alpha_to_eps
  eps_to_alpha = BasisLieHighestWeight.eps_to_alpha

  function test_inverse_alpha_w(dynkin, n, weight)
    lie_type = string(dynkin)
    @test isequal(w_to_alpha(lie_type, n, alpha_to_w(lie_type, n, weight)), weight) # alpha -> w -> alpha
    @test isequal(alpha_to_w(lie_type, n, w_to_alpha(lie_type, n, weight)), weight) # w -> alpha -> w
  end

  function test_inverse_eps_w(dynkin, n, weight)
    lie_type = string(dynkin)
    if lie_type in ["A", "G"]
      weight_representative = w_to_eps(lie_type, n, eps_to_w(lie_type, n, weight))
      weight_representative .-= weight_representative[end]
      pop!(weight_representative)
      @test isequal(weight_representative, weight) # eps -> w -> eps
    else
      @test isequal(w_to_eps(lie_type, n, eps_to_w(lie_type, n, weight)), weight) # eps -> w -> eps
    end
    @test isequal(eps_to_w(lie_type, n, w_to_eps(lie_type, n, weight)), weight) # w -> eps -> w
  end

  function test_inverse_eps_alpha(dynkin, n, weight)
    lie_type = string(dynkin)
    if lie_type in ["A", "G"]
      weight_representative = alpha_to_eps(lie_type, n, eps_to_alpha(lie_type, n, weight))
      weight_representative .-= weight_representative[end]
      pop!(weight_representative)
      @test isequal(weight_representative, weight) # eps -> alpha -> eps
    else
      @test isequal(alpha_to_eps(lie_type, n, eps_to_alpha(lie_type, n, weight)), weight) # eps -> alpha -> eps
    end
    @test isequal(eps_to_alpha(lie_type, n, alpha_to_eps(lie_type, n, weight)), weight) # alpha -> eps -> alpha
  end

  @testset "Dynkin type $dynkin" for dynkin in ('A', 'B', 'C', 'D', 'E', 'F', 'G')
    @testset "n = $n" for n in 1:10
      if (
        !(dynkin == 'B' && n < 2) &&
        !(dynkin == 'C' && n < 2) &&
        !(dynkin == 'D' && n < 4) &&
        !(dynkin == 'E' && !(n == 6 || n == 7 || n == 8)) &&
        !(dynkin == 'F' && n != 4) &&
        !(dynkin == 'G' && (n != 2))
      )
        if dynkin == 'E' && n in [6, 7]
          @test_broken false # TODO: fix these test cases
          continue
        end
        for weight in test_sample_weights[n]
          print(".")
          test_inverse_alpha_w(dynkin, n, weight)
          test_inverse_eps_w(dynkin, n, weight)
          test_inverse_eps_alpha(dynkin, n, weight)
        end
      end
    end
  end
end
