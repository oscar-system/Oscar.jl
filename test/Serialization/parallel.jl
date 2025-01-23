using Distributed

procs = addprocs(1)
@everywhere using Oscar

@everywhere begin
  struct SampleDataStruct{T} <: Oscar.ParallelTask where T <: Vector{<:RingElem}
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
  res1 = Oscar.wait_all_parallel(a)
  @test all(p in [(true, one(R)), (true, y), (true, x)] for p in res1)

  res2 = Oscar.wait_first_parallel(a)
  @test res2 in [one(R), x, y]
  map(rmprocs, procs)
end
