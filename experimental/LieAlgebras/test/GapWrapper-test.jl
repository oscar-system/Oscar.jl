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
        @test L == lie_algebra(Oscar.LieAlgebras._gap_object(L), symbols(L))
      end
    end
  end

  @testset "LinearLieAlgebra" begin
    baserings = [QQ, cyclotomic_field(4)[1]]
    @testset for R in baserings
      lie_algebras = [
        general_linear_lie_algebra(R, 3),
        special_linear_lie_algebra(R, 4),
        special_orthogonal_lie_algebra(R, 4),
      ]
      @testset for L in lie_algebras
        @test L ==
          lie_algebra(LinearLieAlgebra, Oscar.LieAlgebras._gap_object(L), symbols(L))
      end
    end
  end
end
