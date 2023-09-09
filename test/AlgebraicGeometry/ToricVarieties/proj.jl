@testset "Projectivization of direct sum of line bundles on toric space (set_attributes = $set_attributes)" for set_attributes in [true, false]
  P1 = projective_space(NormalToricVariety, 1; set_attributes)
  @testset "Test of the Hirzebruch surfaces from 0 to 10" begin
    for a in 0:10
      Da = toric_divisor(P1, [a,0])
      D0 = toric_divisor(P1, [0,0])
      X = proj(D0, Da)
      @test is_smooth(X) == true
      @test (a < 2) == is_fano(X)
      @test rank(picard_group(X)) == 2
      @test integrate(cohomology_class(anticanonical_divisor(X))^dim(X)) == integrate(cohomology_class(anticanonical_divisor(hirzebruch_surface(NormalToricVariety, a)))^2)
    end
  end

  # We will test the rank of the Picard group and the constant (-K_X)^3
  # Our reference is https://www.fanography.info/toric
  
  # Let us start from the projective bundle over P2
  P2 = projective_space(NormalToricVariety, 2; set_attributes)
  D0 = toric_divisor(P2, [0,0,0])
  X1_P2 = proj(D0, D0)
  D1 = toric_divisor(P2, [0,0,1])
  X2_P2 = proj(D0, D1)
  D2 = toric_divisor(P2, [0,2,0])
  X3_P2 = proj(D0, D2)
  @testset "Test of some Fano projective bundles of dimension 3 over P2" begin
    @test rank(picard_group(X1_P2)) == 2
    @test is_fano(X1_P2) == true
    @test integrate(cohomology_class(anticanonical_divisor(X1_P2))^dim(X1_P2)) == 54
    @test rank(picard_group(X2_P2)) == 2
    @test is_fano(X2_P2) == true
    @test integrate(cohomology_class(anticanonical_divisor(X2_P2))^dim(X2_P2)) == 56
    @test rank(picard_group(X3_P2)) == 2
    @test is_fano(X3_P2) == true
    @test integrate(cohomology_class(anticanonical_divisor(X3_P2))^dim(X3_P2)) == 62
  end
  
  # Projective bundle over F1 = hirzebruch_surface(1)
  D0 = toric_divisor(P1, [0,0])
  D1 = toric_divisor(P1, [-1,0])
  F1 = proj(D0, D1)
  X1_F1 = proj(toric_divisor(F1, [1,1,0,0]), toric_divisor(F1, [1,1,0,0]))
  l1_F1 = toric_line_bundle(F1, [0,1])
  X2_F1 = proj(toric_line_bundle(F1, [0,0]), l1_F1)
  @testset "Test of some Fano projective bundles of dimension 3 over F1" begin
    @test rank(picard_group(X1_F1)) == 3
    @test is_fano(X1_F1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X1_F1))^dim(X1_F1)) == 48  
    @test integrate(cohomology_class(l1_F1)^2) == 1
    @test rank(picard_group(X2_F1)) == 3
    @test is_fano(X2_F1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X2_F1))^dim(X2_F1)) == 50
  end
  
  P1xP1 = proj(D0, D0)
  O00 = toric_line_bundle(P1xP1, [0,0])
  O10 = toric_line_bundle(P1xP1, [1,0])
  O01 = toric_line_bundle(P1xP1, [0,1])
  O11 = toric_line_bundle(P1xP1, [1,1])
  X1_P1xP1 = proj(O10, O01)
  X2_P1xP1 = proj(O10, O10)
  X3_P1xP1 = proj(O00, O11)
  @testset "Test of some Fano projective bundles of dimension 3 over P1 * P1" begin
    @test rank(picard_group(X1_P1xP1)) == 3
    @test is_fano(X1_P1xP1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X1_P1xP1))^dim(X1_P1xP1)) == 44
    @test rank(picard_group(X2_P1xP1)) == 3
    @test is_fano(X2_P1xP1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X2_P1xP1))^dim(X2_P1xP1)) == 48
    @test rank(picard_group(X3_P1xP1)) == 3
    @test is_fano(X3_P1xP1) == true
    @test integrate(cohomology_class(anticanonical_divisor(X3_P1xP1))^dim(X3_P1xP1)) == 52
  end
  
end

@testset "Total space of direct sum of line bundles on toric space (set_attributes = $set_attributes)" for set_attributes in [true, false]

	@testset "Test that some vb on P1 are Calabi-Yau" begin
		P1 = projective_space(NormalToricVariety, 1; set_attributes)
		for a in 0:5
			for b in 0:5
				la = toric_line_bundle(P1, [a])
				lb = toric_line_bundle(P1, [-b])
				X = total_space(la, lb)
				@test is_smooth(X) == true
				@test !is_fano(X)
				@test !is_complete(X)
				@test rank(picard_group(X)) == 1
				@test dim(X) == 3
				@test (degree(canonical_bundle(X)) == 0) == (a - b == -2)
			end
		end
	end

	@testset "Test that omega_S is Calabi-Yau for any S smooth surface" begin
		@testset "S is Hirzebruch" begin
			for a in 0:10
				S = hirzebruch_surface(NormalToricVariety, a; set_attributes)
				X = total_space(canonical_bundle(S))
				@test is_smooth(X) == true
				@test !is_fano(X)
        @test !is_complete(X)
				@test rank(picard_group(X)) == rank(picard_group(S))
				@test dim(X) == 3
				@test degree(canonical_bundle(X)) == 0
			end
		end

		@testset "S is del Pezzo" begin
			for a in 0:3
				S = del_pezzo_surface(NormalToricVariety, a; set_attributes)
				X = total_space(canonical_divisor(S))
				@test is_smooth(X) == true
				@test !is_fano(X)
				@test rank(picard_group(X)) == rank(picard_group(S))
				@test dim(X) == 3
				@test degree(canonical_bundle(X)) == 0
			end
		end
	end

end
