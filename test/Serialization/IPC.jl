using Distributed

process_ids = addprocs(5)

@everywhere using Oscar

@testset "Interprocess Serialization" begin

  chnnls = set_channels(MatElem, FieldElem, Tuple{Ring, Field, MatSpace})

  Qx, x = QQ["x"]
  F, a = number_field(x^2 + x + 1)
  MR = matrix_space(F, 2, 2)

  put_params(chnnls[3], (Qx, F, MR))
  w_pool = WorkerPool(workers())
  c = [[a^i F(1); a a + 1] for i in 1:10]

  function t1()
    for m in c
      put!(chnnls[1], MR(m))
    end
  end
      
  errormonitor(@async t1)
  
  for m in c
    w = take!(w_pool)
    remote_do(det, w, chnnls[1], chnnls[2])
  end
  errormonitor(task2)

  total = F(1)
  task3 = @async for m in c
    determinant = take!(chnnls[2])
    total *= determinant
  end

  errormonitor(task3)
  @test total == 

end

map(rmprocs, process_ids)

