@testset "flattenings of modules I" begin
  kk = GF(40009);
  B, (t,) = polynomial_ring(kk, [:t]);
  IP5 = projective_space(B, [:u, :v, :w, :x, :y, :z]);
  S = homogeneous_coordinate_ring(IP5);
  (u, v, w, x, y, z) = gens(S);
  f = x*y - t^2*v^2 - w^2;
  g = z*w - t*u*v;
  I = ideal(S, [f, g]);
  X, inc_X = sub(IP5, I);
  Om1X = Oscar.relative_cotangent_module(X);
  SX = homogeneous_coordinate_ring(X);
  F1X = graded_free_module(SX,[0]);
  TX,_ = hom(Om1X,F1X);
  TXamb,_ = pushforward(inc_X,TX);
  W1X,_ = hom(TX,F1X)
  M, b = Oscar.simplify(W1X)
  a = get_attribute(b, :inverse)
  @test is_isomorphism(a)
  @test is_isomorphism(b)
  @test compose(a, b) == id_hom(domain(a))
  @test compose(b, a) == id_hom(domain(b))
end

@testset "flattenings of modules II" begin
  # The same as above, but over a quotient ring as coefficient ring
  kk = GF(40009);
  B, (t,) = polynomial_ring(kk, [:t]);
  B, _ = quo(B, ideal(B, [t^3]))
  IP5 = projective_space(B, [:u, :v, :w, :x, :y, :z]);
  S = homogeneous_coordinate_ring(IP5);
  (u, v, w, x, y, z) = gens(S);
  f = x*y - t^2*v^2 - w^2;
  g = z*w - t*u*v;
  I = ideal(S, [f, g]);
  X, inc_X = sub(IP5, I);
  Om1X = Oscar.relative_cotangent_module(X);
  M, b = Oscar.simplify(Om1X)
  a = get_attribute(b, :inverse)
  @test is_isomorphism(a)
  @test is_isomorphism(b)
  @test compose(a, b) == id_hom(domain(a))
  @test compose(b, a) == id_hom(domain(b))
end
