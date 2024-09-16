using Distributed

process_ids = addprocs(1)

@everywhere using Oscar

@testset "Interprocess Serialization" begin
  Qx, x = QQ["x"]
  F, a = number_field(x^2 + x + 1)
  MR = matrix_space(F, 2, 2)
  
  c = [MR([a^i F(1); a a + 1]) for i in 1:5]
  dets = pmap(det, c)
  total = reduce(*, dets)

  @test total == F(4)
end

map(rmprocs, process_ids)


