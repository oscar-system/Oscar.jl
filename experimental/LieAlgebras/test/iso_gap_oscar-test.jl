import Oscar: GAPWrap

function test_iso_gap_oscar(LG, oscarT; num_random_tests::Int=10)
  iso = Oscar.iso_gap_oscar(LG)
  @test domain(iso) === LG
  LO = codomain(iso)
  @test LO isa oscarT

  @test iso === Oscar.iso_gap_oscar(LG) # test caching

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

@testset "LieAlgebras.iso_gap_oscar" begin
  baserings = [GAP.Globals.Rationals, GAPWrap.CF(4)]

  @testset for RG in baserings
    @testset "AbstractLieAlgebra" begin
      lie_algebras = [
        GAP.Globals.SimpleLieAlgebra(GAP.Obj("A"), 3, RG),
        GAP.Globals.SimpleLieAlgebra(GAP.Obj("B"), 2, RG),
        GAP.Globals.LieAlgebraByStructureConstants(
          RG,
          GAP.evalstr(
            "[ [ [ [  ], [  ] ], [ [ 5 ], [ -1 ] ], [ [ 6 ], [ -1 ] ], [ [ 7 ], [ -1 ] ], [ [ 2 ], [ 1 ] ], [ [ 3 ], [ 1 ] ], [ [ 4 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ] ], [ [ [ 5 ], [ 1 ] ], [ [  ], [  ] ], [ [ 8 ], [ -1 ] ], [ [ 9 ], [ -1 ] ], [ [ 1 ], [ -1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 3 ], [ 1 ] ], [ [ 4 ], [ 1 ] ], [ [  ], [  ] ] ], [ [ [ 6 ], [ 1 ] ], [ [ 8 ], [ 1 ] ], [ [  ], [  ] ], [ [ 10 ], [ -1 ] ], [ [  ], [  ] ], [ [ 1 ], [ -1 ] ], [ [  ], [  ] ], [ [ 2 ], [ -1 ] ], [ [  ], [  ] ], [ [ 4 ], [ 1 ] ] ], [ [ [ 7 ], [ 1 ] ], [ [ 9 ], [ 1 ] ], [ [ 10 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 1 ], [ -1 ] ], [ [  ], [  ] ], [ [ 2 ], [ -1 ] ], [ [ 3 ], [ -1 ] ] ], [ [ [ 2 ], [ -1 ] ], [ [ 1 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 8 ], [ -1 ] ], [ [ 9 ], [ -1 ] ], [ [ 6 ], [ 1 ] ], [ [ 7 ], [ 1 ] ], [ [  ], [  ] ] ], [ [ [ 3 ], [ -1 ] ], [ [  ], [  ] ], [ [ 1 ], [ 1 ] ], [ [  ], [  ] ], [ [ 8 ], [ 1 ] ], [ [  ], [  ] ], [ [ 10 ], [ -1 ] ], [ [ 5 ], [ -1 ] ], [ [  ], [  ] ], [ [ 7 ], [ 1 ] ] ], [ [ [ 4 ], [ -1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 1 ], [ 1 ] ], [ [ 9 ], [ 1 ] ], [ [ 10 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 5 ], [ -1 ] ], [ [ 6 ], [ -1 ] ] ], [ [ [  ], [  ] ], [ [ 3 ], [ -1 ] ], [ [ 2 ], [ 1 ] ], [ [  ], [  ] ], [ [ 6 ], [ -1 ] ], [ [ 5 ], [ 1 ] ], [ [  ], [  ] ], [ [  ], [  ] ], [ [ 10 ], [ -1 ] ], [ [ 9 ], [ 1 ] ] ], [ [ [  ], [  ] ], [ [ 4 ], [ -1 ] ], [ [  ], [  ] ], [ [ 2 ], [ 1 ] ], [ [ 7 ], [ -1 ] ], [ [  ], [  ] ], [ [ 5 ], [ 1 ] ], [ [ 10 ], [ 1 ] ], [ [  ], [  ] ], [ [ 8 ], [ -1 ] ] ], [ [ [  ], [  ] ], [ [  ], [  ] ], [ [ 4 ], [ -1 ] ], [ [ 3 ], [ 1 ] ], [ [  ], [  ] ], [ [ 7 ], [ -1 ] ], [ [ 6 ], [ 1 ] ], [ [ 9 ], [ -1 ] ], [ [ 8 ], [ 1 ] ], [ [  ], [  ] ] ], -1, 0 ]"
          ),
        ),
      ]

      @testset for LG in lie_algebras
        test_iso_gap_oscar(LG, AbstractLieAlgebra)
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
        test_iso_gap_oscar(LG, LinearLieAlgebra)
      end
    end

    @testset "Infinite dimensional Lie algebra" begin
      lie_algebras = [GAP.Globals.FreeLieAlgebra(RG, 1), GAP.Globals.FreeLieAlgebra(RG, 2)]

      @testset for LG in lie_algebras
        @test_throws ErrorException Oscar.iso_gap_oscar(LG)
      end
    end
  end
end
