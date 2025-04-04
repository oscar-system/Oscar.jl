@testset "Parallel" begin
  using Distributed
  # Start a worker pool with Oscar running on the workers.
  oscar_worker_pool(1) do wp
    @testset "pmap" begin
      R, (x, y, z) = GF(101)[:x, :y, :z]
      @test pmap(f->f^2, wp, gens(R)) == [x^2, y^2, z^2]
    end

    @testset "compute distributed" begin
      # Define an auxiliary context object for the distributed computation
      mutable struct CtxIO
        a::Vector
        b::Dict{Int, Any}
        function CtxIO(a::Vector)
          return new(a, Dict{Int, Any}())
        end
      end

      # Implement methods for `pop_task!` and `process_result!` for this 
      # type as explained in their docstrings. 
      Oscar.pop_task!(ctx::CtxIO) = is_empty(ctx.a) ? nothing : (length(ctx.a), x->x^2, pop!(ctx.a))
      Oscar.process_result!(ctx::CtxIO, id::Int, result) = ctx.b[id] = result
      
      # Create a context for the computation.
      R, (x, y, z) = QQ[:x, :y, :z]
      ctx = CtxIO(gens(R))
      
      # Trigger distributed execution.
      Oscar.compute_distributed!(ctx, wp)
      
      # Verify that the result is correct. 
      @test all(x in keys(ctx.b) for x in [1, 2, 3])
      @test all(a^2 in values(ctx.b) for a in gens(R))
    end
  end
end
  
