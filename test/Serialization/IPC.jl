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

  function put_collection()
    for m in c
      put!(chnnls[1], MR(m))
    end
  end

  put_collection()
  
  function wrap_take(f)
    function g(input_channel, output_channel)
      println("g")
      args = take!(input_channel)
      println(args)
      result = f(args)
      put!(output_channel, result)
    end
    return g
  end
  
  function remote_det()
    for w in workers()
      remote_do(w, chnnls[1], chnnls[2]) do args
        println(args)
      end
    end
  end


  total = F(1)
  function t3()
    for w in workers()
      println("$w")
      determinant = take!(chnnls[2])
      total *= determinant
      yield()
    end
  end
  
  
  @sync begin
    
    @async t2;
    @async t3;
  end



  @test total == 



end

map(rmprocs, process_ids)

