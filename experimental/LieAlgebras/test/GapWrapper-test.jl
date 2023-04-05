num_random_tests = 10

@testset "LieAlgebras.GapWrapper" begin
  @testset "AbstractLieAlgebra" begin
    function sl2_struct_consts(R::Ring)
      sc = zeros(R, 3, 3, 3)
      sc[1, 2, 3] = R(1)
      sc[2, 1, 3] = R(-1)
      sc[3, 1, 1] = R(2)
      sc[1, 3, 1] = R(-2)
      sc[3, 2, 2] = R(-2)
      sc[2, 3, 2] = R(2)
      return sc
    end

    baserings = [QQ, cyclotomic_field(4)[1]]
    @testset for R in baserings
      lie_algebras = [
        lie_algebra(R, sl2_struct_consts(R), ["e", "f", "h"]),
        lie_algebra(R, ('A', 3)),
        lie_algebra(R, ('B', 2)),
      ]

      @testset for L in lie_algebras
        iso = Oscar.LieAlgebras._iso_oscar_gap(L)
        @test domain(iso) == L
        @test L == lie_algebra(AbstractLieAlgebra, codomain(iso), symbols(L))

        for _ in 1:num_random_tests
          x = L(rand(-10:10, dim(L)))
          y = L(rand(-10:10, dim(L)))
          x_ = iso(x)
          y_ = iso(y)

          # self-inverse
          @test preimage(iso, x_) == x
          @test preimage(iso, y_) == y

          # f homomorphic
          @test iso(2 * x) == 2 * iso(x)
          @test iso(x + y) == iso(x) + iso(y)
          @test iso(x * y) == iso(x) * iso(y)

          # f^-1 homomorphic
          @test preimage(iso, 2 * x_) == 2 * preimage(iso, x_)
          @test preimage(iso, x_ + y_) == preimage(iso, x_) + preimage(iso, y_)
          @test preimage(iso, x_ * y_) == preimage(iso, x_) * preimage(iso, y_)
        end
      end
    end
  end

  @testset "LinearLieAlgebra" begin
    baserings = [QQ, cyclotomic_field(4)[1]]
    @testset for R in baserings
      lie_algebras = [
        general_linear_lie_algebra(R, 2),
        special_linear_lie_algebra(R, 3),
        special_orthogonal_lie_algebra(R, 4),
      ]

      @testset for L in lie_algebras
        iso = Oscar.LieAlgebras._iso_oscar_gap(L)
        @test domain(iso) == L
        @test L == lie_algebra(LinearLieAlgebra, codomain(iso), symbols(L))

        for _ in 1:num_random_tests
          x = L(rand(-10:10, dim(L)))
          y = L(rand(-10:10, dim(L)))
          x_ = iso(x)
          y_ = iso(y)

          # self-inverse
          @test preimage(iso, x_) == x
          @test preimage(iso, y_) == y

          # f homomorphic
          @test iso(2 * x) == 2 * iso(x)
          @test iso(x + y) == iso(x) + iso(y)
          @test iso(x * y) == iso(x) * iso(y)

          # f^-1 homomorphic
          @test preimage(iso, 2 * x_) == 2 * preimage(iso, x_)
          @test preimage(iso, x_ + y_) == preimage(iso, x_) + preimage(iso, y_)
          @test preimage(iso, x_ * y_) == preimage(iso, x_) * preimage(iso, y_)
        end
      end
    end
  end
end
