@testset "du Val Tester" begin
  R,(x,y,z,w) = QQ[:x, :y, :z, :w]
  I = ideal(R,[w,(x^2+1)^2*x^3+(y^2+2)^3*y^2+z^4])
  I2 = ideal(R,[w,x^2+y^2*z + z^7])
  X = spec(quo(R,I)[1])
  X2 = spec(quo(R,I2)[1])
  J1 = ideal(R,[x,y,z,w])
  J2 = ideal(R,[x^2+1,y,z,w])
  J3 = ideal(R,[x,y^2+2,z,w])
  J4 = ideal(R,[x^2+1,y^2+2,z,w])
  @test is_du_val_singularity(X,J1) == true
  @test is_du_val_singularity(X,J3) == false
  @test decide_du_val_singularity(X,J2)[1][3] == (:A, 3)
  @test decide_du_val_singularity(X,J1)[1][3] == (:E, 6)
  @test decide_du_val_singularity(X,J3)[1][1] == false
  @test has_du_val_singularities(X) == false
  @test decide_du_val_singularity(X2,J1)[1][3] == (:D, 8)

  Rg,(x,y,z,w) = graded_polynomial_ring(QQ,[:x,:y,:z,:w])
  X = proj(Rg)
  I = ideal(Rg,[x^2*w+y^2*z+z^3])
  Y = subscheme(X,I)
  @test has_du_val_singularities(Y) == true
end

