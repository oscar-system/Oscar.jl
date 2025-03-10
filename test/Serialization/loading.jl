using LazyArtifacts
artifact_toml = LazyArtifacts.find_artifacts_toml(Oscar.oscardir)
LazyArtifacts.ensure_artifact_installed("version-1-3-0-files", artifact_toml)
_hash = artifact_hash("version-1-3-0-files", artifact_toml)

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

  @testset "Loading all files from v1.3.0" begin
    dir = artifact_path(_hash)
    type_folders = joinpath(dir, "version_1_3_0_files")
    for dir_name in readdir(type_folders)
      type_str = split(dir_name, "-")[1]
      T = Oscar.reverse_type_map[type_str]
      for filename in readdir(joinpath(type_folders, dir_name))
        @testset "$T $filename" begin
          @test load(joinpath(type_folders, dir_name, filename)) isa T
        end
      end
    end
  end
end
