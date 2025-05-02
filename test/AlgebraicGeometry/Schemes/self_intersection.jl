@testset "self intersection" begin
  IP3 = projective_space(QQ, [:x, :y, :z, :w])
  S = homogeneous_coordinate_ring(IP3)
  (x, y, z, w) = gens(S)
  # The following line produces an example where the current heuristic 
  # for moving around divisors to general position still runs into an 
  # infinite loop. We would like to solve this problem in the long run.
  #I = ideal(S, x^3 + y*z*w)
  # This example is fine, though.
  I = ideal(S, x^2 + y*z)
  IPX, inc_IPX = sub(IP3, I)
  X = covered_scheme(IPX)
  J = ideal(homogeneous_coordinate_ring(IPX), [x, y, z])
  JJ = ideal_sheaf(IPX, J)
  bl = blow_up(JJ)
  Y = domain(bl)
  @test is_smooth(Y)
  E = exceptional_divisor(bl)
  E_weil = irreducible_decomposition(weil_divisor(E))
  @test integral(Oscar._intersect(E, E_weil)) == -2
  mov = Oscar.move_divisor(E_weil)
  mov = Oscar.move_divisor(mov)
  # The following line will also produce an infinite loop. TODO: Fix this!
  #mov = Oscar.move_divisor(Oscar.move_divisor(Oscar.move_divisor(E_weil)))
  @test integral(Oscar._intersect(E, mov)) == -2
end

