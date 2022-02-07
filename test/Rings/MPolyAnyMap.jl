@testset begin "MPolyAnyMap"
  A, (x,y) = ZZ["x", "y"]
  f = hom(A, A, [2*x, 5*y])
  R, (u, v) = A["u", "v"]
  h = hom(R, R, f, [u+v, u*v])
  @test (@inferred h(x*u)) == 2*x*u + 2*x*v
end
