using Distributed

process_ids = addprocs(1)

@everywhere using Oscar

@testset "Interprocess Serialization" begin
  (_, _, params_channels) = set_channels(MatElem, FieldElem, Union{Ring, MatSpace})

  Qx, x = QQ["x"]
  F, a = number_field(x^2 + x + 1)
  MR = matrix_space(F, 2, 2)

  put_params(params_channels, Qx)
  put_params(params_channels, F)
  put_params(params_channels, MR)
  
  c = [MR([a^i F(1); a a + 1]) for i in 1:5]
  @everywhere function print_det(m)
    w = myid()
    println("working on worker $w")
    return det(m)
  end
  dets = pmap(print_det, c)
  total = reduce(*, dets)

  @test total == F(4)
end

map(rmprocs, process_ids)

