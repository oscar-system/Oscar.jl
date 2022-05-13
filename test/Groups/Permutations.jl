@testset "Examples.Permutations" begin
  p = @perm (1,2)(3,4)(5,6)
  @test p == cperm([1,2],[3,4],[5,6])
  p = @perm (1,2)(3,4,5)(7,8,9)
  @test p == cperm([1,2],[3,4,5],[7,8,9])

  a=1;b=2;f=x->3*x + 2;
  p = @perm (a,f(a),b)(a+1,b*2)
  @test p == cperm([1,5,4,2])
  
  @test_throws ErrorException @perm (-1, 1)
  @test_throws LoadError @eval @perm "bla"
  @test_throws LoadError @eval @perm 1 + 1
  
  gens = @perm 14 [
         (1,10)
        (2,11)
        (3,12)
        (4,13)
        (5,14)
        (6,8)
        (7,9)
        (1,2,3,4,5,6,7)(8,9,10,11,12,13,14)
        (1,2)(10,11)
       ]
  p = Vector{PermGroupElem}(undef,9)
  p[1] = cperm(symmetric_group(14),[1,10])
  p[2] = cperm(symmetric_group(14),[2,11])
  p[3] = cperm(symmetric_group(14),[3,12])
  p[4] = cperm(symmetric_group(14),[4,13])
  p[5] = cperm(symmetric_group(14),[5,14])
  p[6] = cperm(symmetric_group(14),[6,8])
  p[7] = cperm(symmetric_group(14),[7,9])
  p[8] = cperm(symmetric_group(14),[1,2,3,4,5,6,7],[8,9,10,11,12,13,14])
  p[9] = cperm(symmetric_group(14),[1,2],[10,11])
  @test gens == p
  
  @test_throws ArgumentError @perm 10 [(1,11)]
  
  G = sub(symmetric_group(14),gens)[1]
  @test order(G) == 645120
end


