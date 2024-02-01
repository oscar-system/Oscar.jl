@testset "lifting of maps" begin
  # We create two Koszul complexes and lift the natural map induced 
  # on the associated quotients

  n = 4
  d1 = 3
  R, x = polynomial_ring(QQ, n)

  R1 = FreeMod(R, n)

  v = sum(x^d1*e for (x, e) in zip(gens(R), gens(R1)))

  w = (x[1] + x[2])*R1[1] + (x[2]^2 - x[3])*R1[2] - x[1]^2*R1[3] + (x[4])*R1[4]
  
  K = Oscar.SimpleComplexWrapper(koszul_complex(w)) # The first Koszul complex

  K2 = Oscar.SimpleComplexWrapper(koszul_complex(v)) # The second Koszul complex

  phi = hom(K2[0], K[0], [K[0][1]]) # The identity on degree zero

  Phi = Oscar.lift_map(K2, K, phi) # This creates a lazy morphism of complexes 

  prev = Phi[0]
  for i in 1:n
    next = Phi[i] # We can call for the lifting maps by getindex. It is computed iteratively on the spot.
    @test compose(map(K2, i), prev) == compose(next, map(K, i)) # Squares should commute
    prev = next
  end

  @test !Oscar.can_compute_index(Phi, 10) # An error is thrown if we call an index out of bounds
end
