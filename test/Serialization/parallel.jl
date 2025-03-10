using Distributed

procs = addprocs(1)
@everywhere using Oscar

@everywhere begin
  struct SampleDataStruct{T} <: Oscar.ParallelTask where T
    elems::T
  end

  # The following line communicates that this struct is available for serialization.
  Oscar.@register_serialization_type SampleDataStruct

  function Oscar._compute(ds::SampleDataStruct)
    return true, gcd(ds.elems...)
  end
end

@testset "easy parallelization" begin
  R, (x, y) = QQ[:x, :y]
  l = [x^2, x*y, y^2]
  
  a = [SampleDataStruct([a for (i, a) in enumerate(l) if i != k]) for k in 1:length(l)]
  success, res1 = Oscar.parallel_all(a)
  @test success
  @test all(p in [one(R), x, y] for p in res1)
  success, _, res2 = Oscar.parallel_any(a)
  @test success
  @test res2 in [one(R), x, y]
end

@testset "sending dictionaries" begin
  R, (x, y) = QQ[:x, :y]
  l = [x^2, x*y, y^2]
  a = [SampleDataStruct([Dict(:a => R[b;]) for (i, b) in enumerate(l) if i != k]) for k in 1:length(l)]
  @everywhere begin
    function Oscar._compute(ds::SampleDataStruct)
      return true, prod([c[:a] for c in ds.elems])
    end
  end
  success, res1 = Oscar.parallel_all(a)
  @test success
  @test res1 == [x*y^3, x^2 * y^2, x^3 * y]
end

map(rmprocs, procs)

