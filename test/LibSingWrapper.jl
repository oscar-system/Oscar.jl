@testset "LibSing wrapper" begin
  R, (x,y) = QQ["x","y"]
  f = x^5 + y^5 + x^3*y^3
  shift = hom(R, R, [x-1, y-2])
  (gb, mu, g) = milnor(f)
  shift_f = shift(f)
  (gbs, mus, gs) = milnor(shift_f, point=[1, 2])
  @test mu == mus
  @test mu == 16

  (gb, tau, g) = tjurina(f)
  (gbs, taus, gs) = tjurina(shift_f, point=[1, 2])
  @test tau == 15
  @test tau == taus

  h = (y-x^2)*(y+x^2-2)
  @test global_milnor_number(h) == 3
  @test global_tjurina_number(h) == 2

  g = y^3 - (x-1)^2*(x+1)^2
  @test global_milnor_number(g) == 6
  @test global_tjurina_number(g) == 4

  f = f*(y-2)
  @test global_milnor_number(f) == 30
  @test global_tjurina_number(f) == 20

  R, (x,y,z) = QQ["x", "y", "z"]
  f = [x^2 + y^2 + z^2, x*y]
  @test milnor(f) == 5
end
