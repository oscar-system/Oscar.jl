let dir, T
  dir = @__DIR__
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
    fna = tempname()
    ioa = open(fna, "w")
    Oscar.draw_graph_tikz_aligned(IG, ioa)
    close(ioa)
    desiredfn = joinpath(dir, "robinson6_graph.tikz")
    desiredfna = joinpath(dir, "robinson6_graph_aligned.tikz")
    @test success(`cmp --quiet $fn $desiredfn`)
    @test success(`cmp --quiet $fna $desiredfna`)
  end

  @testset "Harnack6" begin
    h6 = Polymake.tropical.harnack_curve(6)
    f = Oscar._patchworks2polynomial(h6; t=QQ(1, 16))[1]
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    zchart = hom(parent(f), R, [x, y, 1])
    fn = tempname()
    desiredfn = joinpath(dir, "harnack6_graph.tikz")
    desiredfna = joinpath(dir, "harnack6_graph_aligned.tikz")
    draw_curve_tikz(fn, zchart(f); graph=true, overwrite=true, aligned=false)
    fna = tempname()
    draw_curve_tikz(fna, zchart(f); graph=true, overwrite=true, aligned=true)
    IG = Oscar._compute_isotopy_graph(zchart(f), T, 1)[1]
    # Number of affine components should always be the same.
    # This is a better check than the file comparison later, nevertheless the file
    # check is also necessary. But the file might change a lot, so this check
    # gives some confidence whether we can update it.
    @test length(connected_components(IG.G)) == 16
    @test success(`cmp --quiet $fn $desiredfn`)
    @test success(`cmp --quiet $fna $desiredfna`)
  end

  @testset "Shustin's Q(C_4)" begin
    fin = load(joinpath(dir, "QC4.mrdi"))
    Rh = parent(fin)
    B = transpose(matrix(QQ, ([1 0 1; 0 1 1; 0 0 1])))
    f = fin((inv(B) * gens(Rh))...)
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    zchart = hom(Rh, R, [x, y, 1])
    fn = tempname()
    desiredfn = joinpath(dir, "QC4_graph.tikz")
    desiredfna = joinpath(dir, "QC4_graph_aligned.tikz")
    draw_curve_tikz(fn, zchart(f); graph=true, overwrite=true, aligned=false)
    fna = tempname()
    draw_curve_tikz(fna, zchart(f); graph=true, overwrite=true, aligned=true)
    IG = Oscar._compute_isotopy_graph(zchart(f), T, 1)[1]
    @test length(connected_components(IG.G)) == 4
    @test !Oscar._compute_isotopy_graph(zchart(f), identity_matrix(QQ, 2), 1)[3]
    @test success(`cmp --quiet $fn $desiredfn`)
    @test success(`cmp --quiet $fna $desiredfna`)
  end

  @testset "Not implemented yet" begin
    # The following testsets contain situations we are currently unable to draw.
    # Some of them will be easy to deal with and we will get back to those.
    # @testset "Parabola" begin
    #   R, (x, y) = polynomial_ring(QQ, [:x, :y])
    #   f = y - x^2
    #   fn = tempname()
    #   @test_throws NotImplementedError draw_curve_tikz(fn, f; transform=identity_matrix(QQ,2))
    # end

    @testset "Empty" begin
      R, (x, y) = polynomial_ring(QQ, [:x, :y])
      f = x^2 + y^2 + 1
      fn = tempname()
      @test_throws NotImplementedError draw_curve_tikz(fn, f)
    end

    @testset "Lines" begin
      R, (x, y) = polynomial_ring(QQ, [:x, :y])
      f1 = y - 2 * x
      f2 = f1 * (y - 2 * x - 1)
      @test_throws NotImplementedError draw_curve_tikz(tempname(), f1)
      @test_throws NotImplementedError draw_curve_tikz(tempname(), f2)
    end
  end
end
