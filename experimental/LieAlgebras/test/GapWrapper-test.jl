num_random_tests = 10

@testset "LieAlgebras.GapWrapper" begin
  @testset "Oscar -> GAP" begin
    baserings = [QQ, cyclotomic_field(4)[1]]

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

      @testset for R in baserings
        lie_algebras = [
          lie_algebra(R, sl2_struct_consts(R), ["e", "f", "h"]),
          lie_algebra(R, ('A', 3)),
          lie_algebra(R, ('B', 2)),
        ]

        @testset for L in lie_algebras
          iso = Oscar.LieAlgebras._iso_oscar_gap(L)
          @test domain(iso) == L
          @test GAP.Globals.IsLieAlgebra(codomain(iso))

          for _ in 1:num_random_tests
            x = L(rand(-10:10, dim(L)))
            y = L(rand(-10:10, dim(L)))
            x_ = iso(x)
            y_ = iso(y)

            # mapping non-zero to non-zero
            @test iszero(x) == iszero(x_)
            @test iszero(y) == iszero(y_)

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
      @testset for R in baserings
        lie_algebras = [
          general_linear_lie_algebra(R, 2),
          special_linear_lie_algebra(R, 3),
          special_orthogonal_lie_algebra(R, 4),
        ]

        @testset for L in lie_algebras
          iso = Oscar.LieAlgebras._iso_oscar_gap(L)
          @test domain(iso) == L
          @test GAP.Globals.IsLieAlgebra(codomain(iso))

          for _ in 1:num_random_tests
            x = L(rand(-10:10, dim(L)))
            y = L(rand(-10:10, dim(L)))
            x_ = iso(x)
            y_ = iso(y)

            # mapping non-zero to non-zero
            @test iszero(x) == iszero(x_)
            @test iszero(y) == iszero(y_)

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

  @testset "GAP -> Oscar" begin
    baserings = [GAP.Globals.Rationals, GAP.Globals.CyclotomicField(4)]
    @testset for RG in baserings
      @testset "AbstractLieAlgebra" begin
        lie_algebras = [
          GAP.Globals.SimpleLieAlgebra(GAP.Obj("A"), 3, RG),
          GAP.Globals.SimpleLieAlgebra(GAP.Obj("B"), 2, RG),
          GAP.Globals.LieAlgebraByStructureConstants(
            RG,
            GAP.evalstr(
              "[ [ [ [  ], [  ] ], [ [ 5 ], [ -1 ] ], [ [ 6 ], [ -1 ] ], [ [ 7 ], [ -1 ] ], [ [ 2 ], [ 1 ] ], [ [ 3 ], [ 1 ] ], [ [ 4 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ] ], [ [ [ 5 ], [ 1 ] ], [ [  ], [  ] ], [ [ 8 ], [ -1 ] ], [ [ 9 ], [ -1 ] ], [ [ 1 ], [ -1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 3 ], [ 1 ] ], [ [ 4 ], [ 1 ] ], [ [  ], [  ] ] ], [ [ [ 6 ], [ 1 ] ], [ [ 8 ], [ 1 ] ], [ [  ], [  ] ], [ [ 10 ], [ -1 ] ], [ [  ], [  ] ], [ [ 1 ], [ -1 ] ], [ [  ], [  ] ], [ [ 2 ], [ -1 ] ], [ [  ], [  ] ], [ [ 4 ], [ 1 ] ] ], [ [ [ 7 ], [ 1 ] ], [ [ 9 ], [ 1 ] ], [ [ 10 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 1 ], [ -1 ] ], [ [  ], [  ] ], [ [ 2 ], [ -1 ] ], [ [ 3 ], [ -1 ] ] ], [ [ [ 2 ], [ -1 ] ], [ [ 1 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 8 ], [ -1 ] ], [ [ 9 ], [ -1 ] ], [ [ 6 ], [ 1 ] ], [ [ 7 ], [ 1 ] ], [ [  ], [  ] ] ], [ [ [ 3 ], [ -1 ] ], [ [  ], [  ] ], [ [ 1 ], [ 1 ] ], [ [  ], [  ] ], [ [ 8 ], [ 1 ] ], [ [  ], [  ] ], [ [ 10 ], [ -1 ] ], [ [ 5 ], [ -1 ] ], [ [  ], [  ] ], [ [ 7 ], [ 1 ] ] ], [ [ [ 4 ], [ -1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 1 ], [ 1 ] ], [ [ 9 ], [ 1 ] ], [ [ 10 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 5 ], [ -1 ] ], [ [ 6 ], [ -1 ] ] ], [ [ [  ], [  ] ], [ [ 3 ], [ -1 ] ], [ [ 2 ], [ 1 ] ], [ [  ], [  ] ], [ [ 6 ], [ -1 ] ], [ [ 5 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 10 ], [ -1 ] ], [ [ 9 ], [ 1 ] ] ], [ [ [  ], [  ] ], [ [ 4 ], [ -1 ] ], [ [  ], [  ] ], [ [ 2 ], [ 1 ] ], [ [ 7 ], [ -1 ] ], [ [  ], [  ] ], [ [ 5 ], [ 1 ] ], [ [ 10 ], [ 1 ] ], [ [  ], [  ] ], [ [ 8 ], [ -1 ] ] ], [ [ [  ], [  ] ], [ [  ], [  ] ], [ [ 4 ], [ -1 ] ], [ [ 3 ], [ 1 ] ], [ [  ], [  ] ], [ [ 7 ], [ -1 ] ], [ [ 6 ], [ 1 ] ], [ [ 9 ], [ -1 ] ], [ [ 8 ], [ 1 ] ], [ [  ], [  ] ] ], -1, 0 ]",
            ),
          ),
        ]

        @testset for LG in lie_algebras
          iso = Oscar.LieAlgebras._iso_gap_oscar(LG)
          @test domain(iso) == LG
          LO = codomain(iso)
          @test LO isa AbstractLieAlgebra

          for _ in 1:num_random_tests
            x = LO(rand(-10:10, dim(LO)))
            y = LO(rand(-10:10, dim(LO)))
            x_ = preimage(iso, x)
            y_ = preimage(iso, y)

            # mapping non-zero to non-zero
            @test iszero(x) == iszero(x_)
            @test iszero(y) == iszero(y_)

            # self-inverse
            @test iso(x_) == x
            @test iso(y_) == y

            # f homomorphic
            @test iso(2 * x_) == 2 * iso(x_)
            @test iso(x_ + y_) == iso(x_) + iso(y_)
            @test iso(x_ * y_) == iso(x_) * iso(y_)

            # f^-1 homomorphic
            @test preimage(iso, 2 * x) == 2 * preimage(iso, x)
            @test preimage(iso, x + y) == preimage(iso, x) + preimage(iso, y)
            @test preimage(iso, x * y) == preimage(iso, x) * preimage(iso, y)
          end
        end
      end

      @testset "LinearLieAlgebra" begin
        lie_algebras = [
          GAP.Globals.LieAlgebra(
            RG,
            GAP.evalstr(
              "[ [ [ 1, 0 ], [ 0, -1 ] ], [ [ 0, 1 ], [ 0, 0 ] ], [ [ 0, 0 ], [ 1, 0] ] ]"
            ),
          ),
          GAP.Globals.LieAlgebra(
            RG,
            GAP.evalstr(
              "[ [ [ 1, 0 ], [ 0, -1 ] ], [ [ 0, 1 ], [ 0, 0 ] ], [ [ 0, 0 ], [ 1, 0] ] ]"
            ),
            GAP.Obj("basis"),
          ),
        ]

        @testset for LG in lie_algebras
          iso = Oscar.LieAlgebras._iso_gap_oscar(LG)
          @test domain(iso) == LG
          LO = codomain(iso)
          @test LO isa LinearLieAlgebra

          for _ in 1:num_random_tests
            x = LO(rand(-10:10, dim(LO)))
            y = LO(rand(-10:10, dim(LO)))
            x_ = preimage(iso, x)
            y_ = preimage(iso, y)

            # mapping non-zero to non-zero
            @test iszero(x) == iszero(x_)
            @test iszero(y) == iszero(y_)

            # self-inverse
            @test iso(x_) == x
            @test iso(y_) == y

            # f homomorphic
            @test iso(2 * x_) == 2 * iso(x_)
            @test iso(x_ + y_) == iso(x_) + iso(y_)
            @test iso(x_ * y_) == iso(x_) * iso(y_)

            # f^-1 homomorphic
            @test preimage(iso, 2 * x) == 2 * preimage(iso, x)
            @test preimage(iso, x + y) == preimage(iso, x) + preimage(iso, y)
            @test preimage(iso, x * y) == preimage(iso, x) * preimage(iso, y)
          end
        end
      end
    end
  end
end
