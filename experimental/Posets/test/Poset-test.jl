@testset "Posets" begin
  # covering relations, A2 adjoint rep
  a2_adj_cov = [
    0 1 1 0 0 0
    0 0 0 2 1 0
    0 0 0 1 2 0
    0 0 0 0 0 1
    0 0 0 0 0 1
    0 0 0 0 0 0
  ]

  # general relations, A2 adjoint rep
  a2_adj_rel = BitMatrix(
    [
      0 1 1 1 1 1
      0 0 0 1 1 1
      0 0 0 1 1 1
      0 0 0 0 0 1
      0 0 0 0 0 1
      0 0 0 0 0 0
    ]
  )

  # covering relations, B2 adjoint rep
  b2_adj_cov = [
    0 1 1 0 0 0 0 0
    0 0 0 3 1 0 0 0
    0 0 0 1 2 0 0 0
    0 0 0 0 0 2 1 0
    0 0 0 0 0 1 3 0
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 1
    0 0 0 0 0 0 0 0
  ]

  # general relations, B2 adjoint rep
  b2_adj_rel = BitMatrix(
    [
      0 1 1 1 1 1 1 1
      0 0 0 1 1 1 1 1
      0 0 0 1 1 1 1 1
      0 0 0 0 0 1 1 1
      0 0 0 0 0 1 1 1
      0 0 0 0 0 0 0 1
      0 0 0 0 0 0 0 1
      0 0 0 0 0 0 0 0
    ],
  )

  @testset "<(x::PosetElem, y::PosetElem)" begin
    # rel: expected relations
    function test_poset(cov::Matrix{Int}, rel::BitMatrix)
      sz = ncols(cov)
      for _ in 1:10
        ps = poset(cov)
        for _ in 1:5
          x = rand(1:sz)
          y = rand(1:sz)

          @test rel[x, y] == (ps(x) < ps(y))
          @test all(
            ps.set[i, j] == false || ps.rel[i, j] == rel[i, j] for i in 1:sz for
            j in (i + 1):sz
          )
        end
      end
    end

    test_poset(a2_adj_cov, a2_adj_rel)
    test_poset(b2_adj_cov, b2_adj_rel)
  end

  @testset "iterate(::MaximalChainsIterator, ::Vector{Int})" begin
    ps = poset(a2_adj_cov)
    @test collect(maximal_chains(ps)) ==
      [[1, 2, 4, 6], [1, 2, 5, 6], [1, 3, 4, 6], [1, 3, 5, 6]]

    ps = poset(b2_adj_cov)
    @test collect(maximal_chains(ps)) == [
      [1, 2, 4, 6, 8],
      [1, 2, 4, 7, 8],
      [1, 2, 5, 6, 8],
      [1, 2, 5, 7, 8],
      [1, 3, 4, 6, 8],
      [1, 3, 4, 7, 8],
      [1, 3, 5, 6, 8],
      [1, 3, 5, 7, 8],
    ]
  end
end
