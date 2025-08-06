import Oscar: GAPWrap

function test_iso_oscar_gap(LO::LieAlgebra; num_random_tests::Int=25)
  iso = Oscar.iso_oscar_gap(LO)
  @test domain(iso) === LO
  LG = codomain(iso)
  @test GAPWrap.IsLieAlgebra(LG)

  @test iso === Oscar.iso_oscar_gap(LO) # test caching

  for _ in 1:num_random_tests
    x = LO(rand(-10:10, dim(LO)))
    y = LO(rand(-10:10, dim(LO)))
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

@testset "LieAlgebras.iso_oscar_gap" begin
  baserings = [QQ, cyclotomic_field(4)[1]]

  @testset for RO in baserings
    @testset "AbstractLieAlgebra" begin
      function sl2_struct_consts(R::Field)
        sc = zeros(R, 3, 3, 3)
        sc[1, 2, 3] = R(1)
        sc[2, 1, 3] = R(-1)
        sc[3, 1, 1] = R(2)
        sc[1, 3, 1] = R(-2)
        sc[3, 2, 2] = R(-2)
        sc[2, 3, 2] = R(2)
        return sc
      end

      lie_algebras = [
        lie_algebra(RO, sl2_struct_consts(RO), ["e", "f", "h"]),
        lie_algebra(RO, :A, 3),
        lie_algebra(RO, :B, 2),
        (L0 = lie_algebra(RO, :B, 5); root_system(L0); L0),
      ]

      @testset for LO in lie_algebras
        test_iso_oscar_gap(LO)
      end
    end

    @testset "LinearLieAlgebra" begin
      lie_algebras = [
        general_linear_lie_algebra(RO, 2),
        special_linear_lie_algebra(RO, 3),
        special_orthogonal_lie_algebra(RO, 4),
        symplectic_lie_algebra(RO, 6),
        (LO = special_linear_lie_algebra(RO, 4); root_system(LO); LO),
      ]

      @testset for LO in lie_algebras
        test_iso_oscar_gap(LO)
      end
    end
  end
end
