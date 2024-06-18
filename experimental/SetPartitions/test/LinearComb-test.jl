@testset "LinearCombinations of SetPartitions" begin
    @testset "LinearSetPartition Constructor" begin
        S, d = polynomial_ring(QQ, "d")
        @test get_term(linear_partition(Dict(set_partition([1, 2], [1, 1]) => S(4), set_partition([1, 1], [1, 1]) => 8*d))) == Dict(set_partition([1, 2], [1, 1]) => S(4), set_partition([1, 1], [1, 1]) => 8*d)
        @test get_term(linear_partition(Dict(set_partition([1, 2], [1, 1]) => S(0), set_partition([1, 1], [1, 1]) => 8*d))) == Dict(set_partition([1, 1], [1, 1]) => 8*d)
        @test get_term(linear_partition([(set_partition([1, 1], [1, 1]), S(5)), (set_partition([1, 1], [1, 1]), 4*d)])) == Dict(set_partition([1, 1], [1, 1]) => 5 + 4*d)
    end

    @testset "LinearPartition Operations" begin
        S, d = polynomial_ring(QQ, "d")
        a = linear_partition([(set_partition([1, 2], [1, 1]), S(5)), (set_partition([1, 1], [1, 1]), 5*d)])
        @test add(a, a) == a + a == linear_partition(Dict(set_partition([1, 2], [1, 1]) => S(10), set_partition([1, 1], [1, 1]) => 10*d))
        @test add(a, linear_partition([(set_partition([1, 1], [1, 1]), S(1)), (set_partition([1, 1], [1, 1]), 2*d)])) == linear_partition([(set_partition([1, 2], [1, 1]), S(5)), (set_partition([1, 1], [1, 1]), 7*d + 1)])

        @test scale(S(2), linear_partition(Dict(set_partition([1, 2], [1, 1]) => S(10), set_partition([1, 1], [1, 1]) => 8*d))) == S(2) * linear_partition(Dict(set_partition([1, 2], [1, 1]) => S(10), set_partition([1, 1], [1, 1]) => 8*d)) == linear_partition(Dict(set_partition([1, 2], [1, 1]) => S(20), set_partition([1, 1], [1, 1]) => 16*d))

        a = linear_partition([(set_partition([1, 2], [1, 1]), S(4)), (set_partition([1, 1], [1, 1]), 4*d)])
        @test linear_composition(a, a, d) == linear_partition(Dict(set_partition([1, 2], [1, 1]) => 16*d + 16, set_partition([1, 1], [1, 1]) => 16*d^2 + 16*d))
    
        @test linear_tensor_product(a, a) == linear_partition(Dict(set_partition([1, 1, 2, 2], [1, 1, 2, 2]) => 16*d^2, set_partition([1, 2, 3, 3], [1, 1, 3, 3]) => 16*d, set_partition([1, 2, 3, 4], [1, 1, 3, 3]) => S(16), set_partition([1, 1, 2, 3], [1, 1, 2, 2]) => 16*d))
        
        @test subtract(a, a) == a - a == linear_partition(Dict(set_partition([1, 1], [1, 1]) => S(0)))
    end
end
