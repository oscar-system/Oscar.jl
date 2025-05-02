@testset "group recognition, basic functions" begin
  @testset for G in [cyclic_group(PermGroup, 1),
                     symmetric_group(4),
                     symmetric_group(100),
                     general_linear_group(3, 2),
                     general_linear_group(4, 9),
                     special_linear_group(4, 7),
                     ]
    res = recognize(G)
    @test is_ready(res)
    x = rand(G)
    @test x in res
    if G isa MatrixGroup
      @test matrix(x) in res
    end
    slp = straight_line_program(res, x)
    ng = nice_gens(res)
    @test evaluate(slp, ng) == x
  end

  @test is_leaf(recognize(symmetric_group(4)))
  @test !is_leaf(recognize(general_linear_group(3, 2)))

  @test_throws MethodError recognize(small_group(12, 1))
end
