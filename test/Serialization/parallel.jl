@testset "easy parallelization" begin
  using Distributed
  procs = addprocs(1)
  @everywhere using Oscar

  R, (x, y) = QQ[:x, :y]

  l = [x^2, x*y, y^2]
  a = [Oscar.SampleDataStruct([a for (i, a) in enumerate(l) if i != k]) for k in 1:length(l)]
  res1 = Oscar.wait_all_parallel(a)
  @test all(p in [(true, one(R)), (true, y), (true, x)] for p in res1)

  res2 = Oscar.wait_first_parallel(a)
  @test res2 in [one(R), x, y]
  map(rmprocs, procs)
end


