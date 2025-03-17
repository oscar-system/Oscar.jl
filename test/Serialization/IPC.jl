using Distributed

procs = addprocs(1)
@everywhere procs using Oscar

@everywhere [myid(), procs...] begin
  struct SampleDataStruct{T} <: Oscar.ParallelTask where T
    elems::T
  end

  # The following line communicates that this struct is available for serialization.
  Oscar.@register_serialization_type SampleDataStruct
end

@testset "easy parallelization" begin
  R, (x, y) = QQ[:x, :y]
  l = [x^2, x*y, y^2]
  a = [SampleDataStruct([a for (i, a) in enumerate(l) if i != k]) for k in 1:length(l)]

  @everywhere procs begin
    function Oscar._compute(ds::SampleDataStruct{Vector{MPolyRingElem}})
      return true, gcd(ds.elems...)
    end
  end

  success, res1 = Oscar.parallel_all(a; workers=procs)
  @test success
  @test all(p in [one(R), x, y] for p in res1)

  success, _, res2 = Oscar.parallel_any(a; workers=procs)
  @test success
  @test res2 in [one(R), x, y]
end

@testset "sending dictionaries" begin
  R, (x, y) = QQ[:x, :y]
  l = [x^2, x*y, y^2]
  a = [SampleDataStruct([Dict(:a => R[b;]) for (i, b) in enumerate(l) if i != k]) for k in 1:length(l)]
  @everywhere procs begin
    function Oscar._compute(ds::SampleDataStruct{Vector{Dict{Symbol, MatElem}}})
      return true, prod([c[:a] for c in ds.elems])
    end
  end
  success, res1 = Oscar.parallel_all(a; workers=procs)
  @test success
  @test res1 == [x*y^3, x^2 * y^2, x^3 * y]
end

@testset "sending matrices" begin
  Qx, x = QQ[:x]
  F, a = number_field(x^2 + x + 1)
  MR = matrix_space(F, 2, 2)
  c = [SampleDataStruct(MR([a^i F(1); a a + 1])) for i in 1:5]
  
  @everywhere procs begin
    Oscar._compute(ds::SampleDataStruct{MatElem}) = true, det(ds.elems)
  end
  
  success, dets = Oscar.parallel_all(c; workers=procs)
  @test success
  total = reduce(*, dets)
  @test total == F(4)
end

map(rmprocs, procs)
