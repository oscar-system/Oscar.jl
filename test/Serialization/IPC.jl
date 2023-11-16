using Distributed

process_ids = addprocs(1)

@everywhere begin
    using Oscar
end

@testset "Interprocess Serialization" begin

  @everywhere function do_work(rings, jobs, results) # define work function everywhere
    Qx, F, MR = take!(rings)
    
    for i in 1:10
      m = take!(jobs)
      put!(results, det(m))
    end
  end

  rings = Distributed.RemoteChannel(() -> Channel{Tuple{Ring, Field, MatSpace}}(10));
  jobs = Distributed.RemoteChannel(() -> Channel{MatElem}(10));
  results = Distributed.RemoteChannel(() -> Channel{FieldElem}(32));

  Qx, x = QQ["x"]
  F, a = number_field(x^2 + x + 1)
  MR = matrix_space(F, 2, 2)
  
  for p in workers()
    put!(rings, (Qx, F, MR))
  end

  n = 3
  for i in 1:n
    put!(jobs, MR([a^i F(1); F(0) F(1)]))
  end

  for p in workers() # start tasks on the workers to process requests in parallel
    remote_do(do_work, p, rings, jobs, results)
  end

  total = F(1)
  
  while n > 0
    determinant = take!(results)
    n = n - 1
    total *= determinant
  end

  @test total == a^6

end

map(rmprocs, process_ids)

