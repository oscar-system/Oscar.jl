dir = joinpath(Oscar.oscardir, "experimental", "DrawingCurves", "test")
T = matrix(QQ, [2052//2055 -111//2055; 111//2055 2052//2055])

@testset "Robinson6" begin
  R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  r6 =
    19 * x^6 − 20 * x^4 * y^2 − 20 * x^2 * y^4 + 19 * y^6 − 20 * x^4 * z^2 +
    60 * x^2 * y^2 * z^2 − 20 * y^4 * z^2 − 20 * x^2 * z^4 − 20 * y^2 * z^4 + 19 * z^6
  zchart_ring, (zx, zy) = polynomial_ring(QQ, [:x, :y])
  zchart_hom = hom(R, zchart_ring, [zx, zy, 1])
  IG = Oscar._compute_isotopy_graph(zchart_hom(r6), T, 1)[1]
  desiredG = load(joinpath(dir, "robinson6.graph"))
  computedG = IG.G
  @test Polymake.graph.isomorphic(desiredG.pm_graph, computedG.pm_graph)
  @test IG.singularNodes == load(joinpath(dir, "robinson6.sn"))
  @test IG.ytangentNodes == load(joinpath(dir, "robinson6.yn"))
  @test IG.node2pair == load(joinpath(dir, "robinson6.n2p"))
  fn = tempname()
  io = open(fn, "w")
  Oscar.draw_graph_tikz(IG, io)
  close(io)
  desiredfn = joinpath(dir, "robinson6_graph.tikz")
  @test success(`cmp --quiet $fn $desiredfn`)
end

@testset "Harnack6" begin
  h6 = Polymake.tropical.harnack_curve(6)
  f = Oscar._patchworks2polynomial(h6; t=QQ(1, 16))[1]
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  zchart = hom(parent(f), R, [x, y, 1])
  fn = tempname()
  desiredfn = joinpath(dir, "harnack6_graph.tikz")
  draw_curve_tikz(fn, zchart(f); graph=true)
  IG = Oscar._compute_isotopy_graph(zchart(f), T, 1)[1]
  # Number of affine components should always be the same.
  # This is a better check than the file comparison later, nevertheless the file
  # check is also necessary. But the file might change a lot, so this check
  # gives some confidence whether we can update it.
  @test length(connected_components(IG.G)) == 16
  @test success(`cmp --quiet $fn $desiredfn`)
end
