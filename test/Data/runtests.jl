using Oscar
using Test

@testset "Data Loading" begin
  data_dir = joinpath(dirname(dirname(@__DIR__)), "data")

  @testset "Surfaces" begin
    surfaces_dir = joinpath(data_dir, "Surfaces")
    for filename in readdir(surfaces_dir)
      surface = load(joinpath(surfaces_dir, filename))
      @test surface isa MPolyIdeal
    end
  end
end
