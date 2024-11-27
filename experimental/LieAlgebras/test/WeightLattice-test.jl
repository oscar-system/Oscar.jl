@testset "LieAlgebras.WeightLattice" begin
  @testset "WeightLatticeElem" begin
    R = root_system(:A, 2)
    w = WeightLatticeElem(R, [2, 2])

    @test root_system(w) === R
  end

  @testset "conjugate_dominant_weight_with_*_elem(w::WeightLatticeElem)" begin
    for (R, vec) in [
      (root_system(:A, 5), [1, -1, 2, 0, 2]),
      (root_system(:B, 3), [1, 1, 1]),
      (root_system(:C, 4), [2, 1, 0, 1]),
      (root_system(:D, 5), [-1, 2, 2, -1, -1]),
      (root_system(:E, 6), [1, 2, 0, 0, 2, 1]),
      (root_system(:F, 4), [1, 2, 3, 4]),
      (root_system(:G, 2), [-1, -1]),
    ]
      wt = WeightLatticeElem(R, vec)
      d, x = conjugate_dominant_weight_with_elem(wt)
      @test is_dominant(d)
      @test x * wt == d
      @test wt * inv(x) == d
    end
  end
end
