@testset "du Val Tester" begin
  R,(x,y,z,w) = QQ["x","y","z","w"]
  I = ideal(R,[w,(x^2+1)^2*x^3+(y^2+2)^3*y^2+z^4])
  I2 = ideal(R,[w,x^2+y^2*z + z^7])
  X = Spec(quo(R,I)[1])
  X2 = Spec(quo(R,I2)[1])
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

  Rg,(x,y,z,w) = graded_polynomial_ring(QQ,["x","y","z","w"])
  X = projective_scheme(Rg)
  I = ideal(Rg,[x^2*w+y^2*z+z^3])
  Y = subscheme(X,I)
  @test has_du_val_singularities(Y) == true
end

@testset "vector_space_dimension and vector_space_basis" begin
  R,(x,y,z,w) = QQ["x","y","z","w"]
  U=complement_of_point_ideal(R,[0,0,0,0])
  RL, loc_map = localization(R,U)
  Floc = free_module(RL,2)
  v = gens(Floc)
  Mloc, _=quo(Floc,[x*v[1],y*v[1],z*v[1],w^2*(w+1)^3*v[1],v[2]])
  @test vector_space_dimension(Mloc) == 2
  @test vector_space_dimension(Mloc,1) == 1
  @test length(vector_space_basis(Mloc)) == 2
  @test length(vector_space_basis(Mloc,0)) == 1
end
