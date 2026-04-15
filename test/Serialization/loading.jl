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

@testset "saving and loading Gzip'ed filed" begin
  @testset "loading file format paper example (Gzip'ed)" begin
    F = GF(7, 2)
    o = gen(F)
    Fyz, (y, z) = F[:x, :y]
    load(joinpath(@__DIR__,"polynomial-example.mrdi.gz"))
    loaded = load(joinpath(@__DIR__,"polynomial-example.mrdi.gz"); params=Fyz)
    @test loaded == 2*y^3*z^4 + 5*o*y + (o + 3)*z^2 + 1
  end

  @testset "saving and loading Gzip'ed file" begin
    R, x = polynomial_ring(QQ, :x)
    f = x^3 + 2x + 1
    g = x^2 + 3
    mktempdir() do path
      filename = joinpath(path, "poly.mrdi.gz")
      save(filename, [f,g]; compression=:gzip)
      loaded = load(filename)
      @test loaded == [f,g]
      Oscar.reset_global_serializer_state()
      loaded = load(filename; params=R)
      @test loaded == [f,g]
    end
  end
end

@testset "pretty printing" begin
  mktempdir() do path
    filename = joinpath(path, "pretty.mrdi")
    save(filename, [[1, 2], [3, 4], [5, 6]]; pretty_print=true)
    str = read(filename, String)
    version_info = Oscar.Serialization.get_oscar_serialization_version()[:Oscar][2]
    cmp_str = "{\n  \"_ns\":{\n    \"Oscar\":[\n      \"https://github.com/oscar-system/Oscar.jl\",\n      \"" * version_info * "\"\n    ]\n  },\n  \"_type\":{\n    \"name\":\"Vector\",\n    \"params\":{\n      \"name\":\"Vector\",\n      \"params\":\"Base.Int\"\n    }\n  },\n  \"data\":[\n    [\n      \"1\",\n      \"2\"\n    ],\n    [\n      \"3\",\n      \"4\"\n    ],\n    [\n      \"5\",\n      \"6\"\n    ]\n  ]\n}"

    @test str == cmp_str
  end
end
