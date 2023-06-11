@testset "Examples.Permutations" begin
  p = @perm (1,2)(3,4)(5,6)
  @test p == cperm([1,2],[3,4],[5,6])
  p = @perm (1,2)(3,4,5)(7,8,9)
  @test p == cperm([1,2],[3,4,5],[7,8,9])

  a=1;b=2;f=x->3*x + 2;
  p = @perm (a,f(a),b)(a+1,b*2)
  @test p == cperm([1,5,4,2])
  p = @perm ()
  @test p == cperm()
  
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

@testset "permutation_group" begin
  g = permutation_group(2, PermGroupElem[])
  @test order(g) == 1
  @test degree(g) == 2

  g = permutation_group(2, [cperm()])
  @test order(g) == 1
  @test degree(g) == 2

  @test order(permutation_group(3, [cperm(1:3)])) == 3
  @test order(permutation_group(5, [cperm(1:3, 4:5), cperm(1:4)])) == 120
end

@testset "@permutation_group" begin
  g = @permutation_group(2)
  @test order(g) == 1
  @test degree(g) == 2

  g = @permutation_group(2, ())
  @test order(g) == 1
  @test degree(g) == 2

  @test order(@permutation_group(3, (1,2,3))) == 3
  @test order(@permutation_group(5, (1,2,3)(4,5), (1,2,3,4))) == 120

  @test_throws ArgumentError @permutation_group(1, (1,2))
  @test_throws ErrorException @permutation_group(1, (1,0))
end

@testset "constructors for projective groups" begin
   @testset for n in 4:5, q in [2, 3, 4]
      G = projective_general_linear_group(n, q)
      S = projective_special_linear_group(n, q)
      qq = ZZRingElem(q)
      @test order(G) == prod([qq^n-qq^i for i in 0:(n-1)])÷(qq-1)
      @test index(G, S) == gcd(qq-1, n)
      @test_throws ArgumentError projective_general_linear_group(2, 6)
      @test_throws ArgumentError projective_special_linear_group(2, 6)
   end

   @testset for n in [4,6], q in [2, 3, 4]
      G = projective_symplectic_group(n, q)
      qq = ZZRingElem(q)
      m = div(n, 2)
      @test order(G) == qq^(m^2) * prod([qq^i-1 for i in 2:2:n])÷gcd(qq-1, 2)
      @test_throws ArgumentError projective_symplectic_group(3, 4)
      @test_throws ArgumentError projective_symplectic_group(2, 6)
   end

   @testset for n in 1:3, q in [2, 3, 4]
      G = projective_unitary_group(n, q)
      S = projective_special_unitary_group(n, q)
      qq = ZZRingElem(q)
      @test order(G) == qq^((n*(n-1))÷2) * prod([qq^i-(-1)^i for i in 2:n])
      @test index(G, S) == gcd(qq+1, n)
      @test_throws ArgumentError projective_unitary_group(2, 6)
      @test_throws ArgumentError projective_special_unitary_group(2, 6)
   end

   @testset for q in [2, 3, 4]
      qq = ZZRingElem(q)
      @testset for n in [2, 4, 6], e in [+1, -1]
         m = div(n, 2)
         G = projective_orthogonal_group(e, n, q)
         S = projective_special_orthogonal_group(e, n, q)
         O = projective_omega_group(e, n, q)
         N = qq^(m*(m-1)) * (q^m-e) * prod([qq^i-1 for i in 2:2:(2*m-2)], init = 1)
         @test order(G) == div(2 * N, gcd(2, qq^m-e))
         @test order(G) == gcd(2, qq^m-e) * order(S)
         @test order(G) == div(2 * gcd(4, qq^m-e), gcd(2, qq^m-e)) * order(O)
      end
      @testset for n in [1, 3, 5]
         e = 0
         m = div(n - 1, 2)
         G = projective_orthogonal_group(e, n, q)
         S = projective_special_orthogonal_group(e, n, q)
         O = projective_omega_group(e, n, q)
         @test G == projective_orthogonal_group(n, q)
         @test S == projective_special_orthogonal_group(n, q)
         @test O == projective_omega_group(n, q)
         N = qq^(m^2) * prod([qq^i-1 for i in 2:2:(2*m)], init = 1)
         @test order(G) == N
         @test order(S) == N
         @test N == 1 || order(S) == gcd(2, q-1) * order(O)
      end
   end

   @test_throws ArgumentError projective_orthogonal_group(1, 3, 5)
   @test_throws ArgumentError projective_orthogonal_group(1, 3, 6)
   @test_throws ArgumentError projective_orthogonal_group(0, 2, 5)
   @test_throws ArgumentError projective_special_orthogonal_group(1, 3, 5)
   @test_throws ArgumentError projective_special_orthogonal_group(1, 3, 6)
   @test_throws ArgumentError projective_special_orthogonal_group(0, 2, 5)
   @test_throws ArgumentError projective_omega_group(1, 3, 5)
   @test_throws ArgumentError projective_omega_group(1, 3, 6)
   @test_throws ArgumentError projective_omega_group(0, 2, 5)
end

@testset "CycleType" begin
  @testset "degree and order" begin
    g = cperm(1:3, 4:5, 6:7, 8:10, 11:15)
    c = cycle_structure(g)
    @test degree(g) == degree(c)
    @test order(g) == order(c)
    @test order(Int, g) == order(Int, c)

    for g in symmetric_group(6)
      c = cycle_structure(g)
      @test degree(g) == degree(c)
      @test order(g) == order(c)
      @test order(Int, g) == order(Int, c)
    end
  end
end
