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

@testset "flattening of a graded free module" begin
  R, (x,) = polynomial_ring(QQ, [:x])
  S, (y,) = polynomial_ring(R, [:y])
  Sg, _ = grade(S)
  F = graded_free_module(Sg, 2)
  flat_map = Oscar.RingFlattening(Sg)
  # manual invocation of `change_base_ring` should not clutter 
  # the cache of `flat_map`
  FS, phi = change_base_ring(flat_map, F)
  @test !isdefined(flat_map, :flat_counterparts) # cache not even created
  @test ngens(FS) == 2
  @test is_graded(FS)
  @test phi(F[1]) == FS[1]
  # If we want modules and elements to be memoized, we use `flat_map` directly.
  FF, iso, iso_inv = flat_map(F)
  @test F in keys(Oscar.flat_counterparts(flat_map)) && Oscar.flat_counterparts(flat_map)[F] === (FF, iso, iso_inv)
  @test FF !== FS
  I = sub_object(F, [x*F[1]])
  IS, phi = change_base_ring(flat_map, I)
  @test !(I in keys(Oscar.flat_counterparts(flat_map))) # no caching
  II, iso, iso_inv = flat_map(I)
  @test I in keys(Oscar.flat_counterparts(flat_map)) # caching when using `flat_map` directly
end

