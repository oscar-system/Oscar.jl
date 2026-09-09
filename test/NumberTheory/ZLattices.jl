@testset "root_overlattices" begin
 @test length(Oscar.root_overlattices(4))==6
end

@testset "invariant functions for integer quadratic lattices" begin
  A = [[2 -1 0 0 0 0; -1 2 -1 0 0 0; 0 -1 2 -1 0 0; 0 0 -1 2 -1 0; 0 0 0 -1 2 0; 0 0 0 0 0 20], [2 0 0 0 -1 -1; 0 2 0 -1 0 -1; 0 0 2 -1 1 0; 0 -1 -1 4 1 2; -1 0 1 1 4 1; -1 -1 0 2 1 4], [2 -1 1 0 0 0; -1 2 -1 0 0 0; 1 -1 2 0 0 0; 0 0 0 2 0 0; 0 0 0 0 2 1; 0 0 0 0 1 8], [2 1 -1 -1 0 0; 1 2 -1 -1 0 0; -1 -1 2 1 0 0; -1 -1 1 2 0 0; 0 0 0 0 2 0; 0 0 0 0 0 12], [2 -1 0 0 0 -1; -1 2 0 0 0 0; 0 0 2 0 1 0; 0 0 0 2 1 0; 0 0 1 1 4 0; -1 0 0 0 0 4], [2 -1 1 0 -1 -1; -1 2 -1 0 1 1; 1 -1 2 0 0 0; 0 0 0 2 0 0; -1 1 0 0 4 1; -1 1 0 0 1 6], [2 -1 1 1 -1 0; -1 2 -1 -1 0 0; 1 -1 2 0 0 0; 1 -1 0 2 -1 0; -1 0 0 -1 2 0; 0 0 0 0 0 30]]
  L = [integer_lattice(gram=matrix(QQ,6,6,i),cached=false) for i in A]
  @test 7 == length(unique!(Oscar.invariant_function_graph_hash.(L)))
  @test 7 == length(unique!(Oscar.oscar_invariant_function.(L)))
end 
