@testset "Total space of direct sum of line bundles on toric space" begin
  @testset "Test that some vector bundles on P1 are Calabi-Yau" begin
    P1 = projective_space(NormalToricVariety, 1)
    for a in 0:5, b in 0:5
      la = toric_line_bundle(P1, [a])
      lb = toric_line_bundle(P1, [-b])
      X = total_space(la, lb)
      @test is_smooth(X) == true
      @test !is_fano(X)
      @test !is_complete(X)
      @test torsion_free_rank(picard_group_with_map(X)[1]) == 1
      @test dim(X) == 3
      @test (degree(canonical_bundle(X)) == 0) == (a - b == -2)
    end
  end

  @testset "Test that omega_S is Calabi-Yau for any smooth surface S" begin
    @testset "S is Hirzebruch" begin
      for a in 0:10
        S = hirzebruch_surface(NormalToricVariety, a)
        X = total_space(canonical_bundle(S))
        @test is_smooth(X) == true
        @test !is_fano(X)
        @test !is_complete(X)
        @test torsion_free_rank(picard_group_with_map(X)[1]) ==
          torsion_free_rank(picard_group_with_map(S)[1])
        @test dim(X) == 3
        @test degree(canonical_bundle(X)) == 0
      end
    end

    @testset "S is del Pezzo" begin
      for a in 0:3
        S = del_pezzo_surface(NormalToricVariety, a)
        X = total_space(canonical_divisor(S))
        @test is_smooth(X) == true
        @test !is_fano(X)
        @test torsion_free_rank(picard_group_with_map(X)[1]) ==
          torsion_free_rank(picard_group_with_map(S)[1])
        @test dim(X) == 3
        @test degree(canonical_bundle(X)) == 0
      end
    end
  end
end
