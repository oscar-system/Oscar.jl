@testset "Examples.Permutations" begin
  p = Oscar.Permutations.@perm (1,2)(3,4)(5,6)
  @test p  == cperm([1,2],[3,4],[5,6])
  p = Oscar.Permutations.@perm (1,2)(3,4,5)(7,8,9)
  @test p == cperm([1,2],[3,4,5],[7,8,9])

  a=1;b=2;f=x->3*x + 2;
  p = Oscar.Permutations.@perm (a,f(a),b)(a+1,b*2)
  @test p == cperm([1,5,4,2])
end


