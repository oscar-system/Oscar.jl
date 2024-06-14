@testset "Triangulations" begin
  @testset "FTheory" begin
    pts = matrix(ZZ, [0 0 0 0; -3 -4 -8 2; 0 0 1 -1; 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
    triang = Oscar._find_full_star_triangulation(pts; seed=994)
    sop = subdivision_of_points(pts, triang)
    @test Oscar._is_full_triangulation(sop)
    @test Oscar._is_star_triangulation(sop)
  end

  @testset "Square" begin
    pts = matrix(ZZ, [0 0; 1000 0; 0 1000; 501 501])
    pts = 1000000000 * pts
    triang = Oscar._find_full_star_triangulation(pts; seed=47235)
    sop = subdivision_of_points(pts, triang)
    @test Oscar._is_full_triangulation(sop)
    @test Oscar._is_star_triangulation(sop)
  end
end
