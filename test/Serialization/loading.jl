@testset "loading" begin
  @testset "loading file format paper example" begin
    F = GF(7, 2)
    o = gen(F)
    Fyz, (y, z) = F[:x, :y]
    load(joinpath(@__DIR__,"polynomial-example.mrdi");)
    loaded = load(joinpath(@__DIR__,"polynomial-example.mrdi"); params=Fyz)
    @test loaded == 2*y^3*z^4 + 5*o*y + (o + 3)*z^2 + 1
  end

  @testset "loading Vector{LinearProgram}" begin
    c = cube(3)
    LP0 = linear_program(c, [2,2,-3])
    LP1 = linear_program(c, [2,2,4])
    v = [LP0, LP1]
    loaded = load(joinpath(@__DIR__,"vlp.mrdi"))
    @test length(v) == length(loaded)
    @test feasible_region(loaded[1]) == feasible_region(loaded[2])
    @test feasible_region(loaded[1]) == feasible_region(LP0)
    @test objective_function(loaded[1]) == objective_function(v[1])
    @test objective_function(loaded[2]) == objective_function(v[2])
    @test optimal_value(loaded[1]) == optimal_value(v[1])
    @test optimal_value(loaded[2]) == optimal_value(v[2])
  end
end
