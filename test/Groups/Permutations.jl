@testset "Examples.Permutations" begin
  p = @perm (1,2)(3,4)(5,6)
  @test p == cperm([1,2],[3,4],[5,6])
  p = @perm (1,2)(3,4,5)(7,8,9)
  @test p == cperm([1,2],[3,4,5],[7,8,9])

  a=1;b=2;f=x->3*x + 2;
  p = @perm (a,f(a),b)(a+2,b*2)
  @test p == cperm([1,5,2],[3,4])
  p = @perm ()
  @test p == cperm()
    
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
  G = symmetric_group(14)
  p = Vector{PermGroupElem}(undef,9)
  p[1] = cperm(G,[1,10])
  p[2] = cperm(G,[2,11])
  p[3] = cperm(G,[3,12])
  p[4] = cperm(G,[4,13])
  p[5] = cperm(G,[5,14])
  p[6] = cperm(G,[6,8])
  p[7] = cperm(G,[7,9])
  p[8] = cperm(G,[1,2,3,4,5,6,7],[8,9,10,11,12,13,14])
  p[9] = cperm(G,[1,2],[10,11])
  @test gens == p
    
  G = sub(G,gens)[1]
  @test order(G) == 645120
end

@testset "@perm" begin
  @testset "@perm <...> cycles" begin
    @testset "@perm cycles" begin
      p = @perm (1,2,6)(4,5)
      @test p isa PermGroupElem
      G = parent(p)
      @test degree(G) == 6
      @test is_natural_symmetric_group(G)
      @test order(p) == 6
  
      p = @perm ()
      @test p isa PermGroupElem
      G = parent(p)
      @test degree(G) == 1
      @test is_natural_symmetric_group(G)
      @test order(p) == 1
  
      @test_throws ArgumentError @perm (-1,1)
    end

    @testset "@perm n cycles" begin
      n = 7

      p = @perm n (1,2,6)(4,5)
      @test p isa PermGroupElem
      G = parent(p)
      @test degree(G) == n
      @test is_natural_symmetric_group(G)
      @test order(p) == 6

      p = @perm n ()
      @test p isa PermGroupElem
      G = parent(p)
      @test degree(G) == n
      @test is_natural_symmetric_group(G)
      @test order(p) == 1

      @test_throws ArgumentError @perm 10 (1,11)
      @test_throws ArgumentError @perm 5 (-1,1)
    end

    @testset "@perm G cycles" begin
      S = symmetric_group(7)
      p = @perm S (1,2,6)(4,5)
      @test p isa PermGroupElem
      @test parent(p) == S
      @test order(p) == 6
      p = @perm S ()
      @test p isa PermGroupElem
      @test parent(p) == S
      @test order(p) == 1

      @test_throws ArgumentError @perm S (-1,1)
      @test_throws ArgumentError @perm S (2,5,8)

      A = alternating_group(7)
      p = @perm A (1,6)(4,5)
      @test p isa PermGroupElem
      @test parent(p) == A
      @test order(p) == 2
      p = @perm A ()
      @test p isa PermGroupElem
      @test parent(p) == A
      @test order(p) == 1

      # perms with negative sign
      @test_throws ArgumentError @perm A (3,5) 
      @test_throws ArgumentError @perm A (1,2,3,4) 
      @test_throws ArgumentError @perm A (2,5)(4,1,3) 
    end
  end

  @testset "@perm <...> vector" begin
    let n = 7
      S = symmetric_group(n)

      @test_throws ArgumentError @macroexpand @perm []
      @test_throws ArgumentError @macroexpand @perm n []
      @test_throws ArgumentError @macroexpand @perm S []
    end

    @testset "@perm vector" begin
      ps = @perm [(1,4)(2,3)(6,5)]
      @test ps isa Vector{PermGroupElem}
      G = parent(ps[1])
      @test degree(G) == 6
      @test is_natural_symmetric_group(G)
      @test ps == [@perm(G, (1,4)(2,3)(6,5))]
      
      ps = @perm [(1,2,3)(4,5), (6,2,3,1), (), (2,)(1,4)]
      @test ps isa Vector{PermGroupElem}
      @test allequal(parent, ps)
      G = parent(ps[1])
      @test degree(G) == 6
      @test is_natural_symmetric_group(G)
      @test ps == [@perm(G, (1,2,3)(4,5)), @perm(G, (6,2,3,1)), @perm(G, ()), @perm(G, (2,)(1,4))]
  
      @test_throws ArgumentError @perm [(1,-2)]
    end

    @testset "@perm n vector" begin
      n = 7
      ps = @perm n [(1,4)(2,3)(6,5)]
      @test ps isa Vector{PermGroupElem}
      G = parent(ps[1])
      @test degree(G) == n
      @test is_natural_symmetric_group(G)
      @test ps == [@perm(n, (1,4)(2,3)(6,5))]
      
      ps = @perm n [(1,2,3)(4,5), (6,2,3,1), (), (2,)(1,4)]
      @test ps isa Vector{PermGroupElem}
      @test allequal(parent, ps)
      G = parent(ps[1])
      @test degree(G) == n
      @test is_natural_symmetric_group(G)
      @test ps == [@perm(n, (1,2,3)(4,5)), @perm(n, (6,2,3,1)), @perm(n, ()), @perm(n, (2,)(1,4))]

      @test_throws ArgumentError @perm 10 [(1,11)]
      @test_throws ArgumentError @perm 5 [(1,-2)]
    end

    @testset "@perm G vector" begin
      S = symmetric_group(7)
      A = alternating_group(7)
      
      ps = @perm S [(1,4)(2,3)(6,5)]
      @test ps isa Vector{PermGroupElem}
      @test parent(ps[1]) == S
      @test ps == [@perm(S, (1,4)(2,3)(6,5))]
      
      ps = @perm S [(1,2,3)(4,5), (6,2,3,1), (), (2,)(1,4)]
      @test ps isa Vector{PermGroupElem}
      @test all(p -> parent(p) == S, ps)
      @test ps == [@perm(S, (1,2,3)(4,5)), @perm(S, (6,2,3,1)), @perm(S, ()), @perm(S, (2,)(1,4))]

      @test_throws ArgumentError @perm S [(1,11)]
      @test_throws ArgumentError @perm S [(1,-2)]

      ps = @perm A [(1,4)(6,5)]
      @test ps isa Vector{PermGroupElem}
      @test parent(ps[1]) == A
      @test ps == [@perm(A, (1,4)(6,5))]
      
      ps = @perm A [(1,2,3)(4,5,7), (6,2,3,1,5), (), (2,)(1,4,5)]
      @test ps isa Vector{PermGroupElem}
      @test all(p -> parent(p) == A, ps)
      @test ps == [@perm(A, (1,2,3)(4,5,7)), @perm(A, (6,2,3,1,5)), @perm(A, ()), @perm(A, (2,)(1,4,5))]

      @test_throws ArgumentError @perm A [(1,4)(2,3)(6,5)]
    end
  end

  @testset "@perm <...> tuple" begin
    # the empty tuple does not throw an error as it denotes the identity permutation
    @testset "@perm tuple" begin
      ps = @perm ((1,4)(2,3)(6,5),)
      @test ps isa NTuple{1, PermGroupElem}
      G = parent(ps[1])
      @test degree(G) == 6
      @test is_natural_symmetric_group(G)
      @test ps == (@perm(G, (1,4)(2,3)(6,5)),)
      
      ps = @perm ((1,2,3)(4,5), (6,2,3,1), (), (2,)(1,4))
      @test ps isa NTuple{4, PermGroupElem}
      @test allequal(parent, ps)
      G = parent(ps[1])
      @test degree(G) == 6
      @test is_natural_symmetric_group(G)
      @test ps == (@perm(G, (1,2,3)(4,5)), @perm(G, (6,2,3,1)), @perm(G, ()), @perm(G, (2,)(1,4)))

      @test_throws ArgumentError @perm ((1,-2),)
    end

    @testset "@perm n tuple" begin
      n = 7

      ps = @perm n ((1,4)(2,3)(6,5),)
      @test ps isa NTuple{1, PermGroupElem}
      G = parent(ps[1])
      @test degree(G) == n
      @test is_natural_symmetric_group(G)
      @test ps == (@perm(n, (1,4)(2,3)(6,5)),)
      
      ps = @perm n ((1,2,3)(4,5), (6,2,3,1), (), (2,)(1,4))
      @test ps isa NTuple{4, PermGroupElem}
      @test allequal(parent, ps)
      G = parent(ps[1])
      @test degree(G) == n
      @test is_natural_symmetric_group(G)
      @test ps == (@perm(n, (1,2,3)(4,5)), @perm(n, (6,2,3,1)), @perm(n, ()), @perm(n, (2,)(1,4)))

      @test_throws ArgumentError @perm 10 ((1,11),)
      @test_throws ArgumentError @perm 5 ((1,-2),)
    end

    @testset "@perm G tuple" begin
      S = symmetric_group(7)
      A = alternating_group(7)

      ps = @perm S ((1,4)(2,3)(6,5),)
      @test ps isa NTuple{1, PermGroupElem}
      @test parent(ps[1]) == S
      @test ps == (@perm(S, (1,4)(2,3)(6,5)),)
      
      ps = @perm S ((1,2,3)(4,5), (6,2,3,1), (), (2,)(1,4))
      @test ps isa NTuple{4, PermGroupElem}
      @test all(p -> parent(p) == S, ps)
      @test ps == (@perm(S, (1,2,3)(4,5)), @perm(S, (6,2,3,1)), @perm(S, ()), @perm(S, (2,)(1,4)))

      @test_throws ArgumentError @perm S ((1,11),)
      @test_throws ArgumentError @perm S ((1,-2),)

      ps = @perm A ((1,4)(6,5),)
      @test ps isa NTuple{1, PermGroupElem}
      @test parent(ps[1]) == A
      @test ps == (@perm(A, (1,4)(6,5)),)
      
      ps = @perm A ((1,2,3)(4,5,7), (6,2,3,5,1), (), (2,)(1,4,5))
      @test ps isa NTuple{4, PermGroupElem}
      @test all(p -> parent(p) == A, ps)
      @test ps == (@perm(A, (1,2,3)(4,5,7)), @perm(A, (6,2,3,5,1)), @perm(A, ()), @perm(A, (2,)(1,4,5)))

      @test_throws ArgumentError @perm A ((1,4)(2,3)(6,5),)
    end
  end


  @test_throws MethodError @macroexpand @perm "bla"
  @test_throws ErrorException @macroexpand @perm 1 + 1
  @test_throws ErrorException @macroexpand [@perm (1,2), @perm (3,4), @perm (5,6)]
  @test_throws ErrorException @macroexpand [@perm (1,2)(3,4), @perm (5,6)(7,8)]
  @test_throws ErrorException @macroexpand p1, p2, p3 = @perm (1,2), @perm (3,4), @perm (5,6)
  @test_throws ErrorException @macroexpand p1, p2 = @perm (1,2)(3,4), @perm (5,6)(7,8)
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
  @test_throws ArgumentError @permutation_group(1, (1,0))
end

@testset "parent coercion for permutation groups" begin
  g = symmetric_group(5)
  s = stabilizer(g, 5)[1]
  g4 = symmetric_group(4)
  @test s != g4   # different degree
  s4 = g4(s)
  @test s4 == g4  # same degree

  @test_throws ArgumentError g4(g)
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
  @testset "constructors" begin
    @test Oscar.CycleType([1,2,3,2,3,2,2,2]) == Oscar.CycleType([1 => 1, 2 => 5, 3 => 2])
#   @test_throws ArgumentError Oscar.CycleType([1 => 1, 1 => 1])
#T should the input be checked?
  end

  @testset "degree, order, sign" begin
    g = cperm(1:3, 4:5, 6:7, 8:10, 11:15)
    c = cycle_structure(g)
    @test degree(g) == degree(c)
    @test order(g) == order(c)
    @test order(Int, g) == order(Int, c)
    @test sign(g) == sign(c)
    @test iseven(g) == iseven(c)
    @test isodd(g) == isodd(c)

    for g in symmetric_group(6)
      c = cycle_structure(g)
      @test degree(g) == degree(c)
      @test order(g) == order(c)
      @test order(Int, g) == order(Int, c)
      @test sign(g) == sign(c)
      @test iseven(g) == iseven(c)
      @test isodd(g) == isodd(c)
    end
  end

  @testset "cycles" begin
    g = cperm(1:3, 4:5, 6:7, 8:10, 11:15)
    c = cycle_structure(g)
    cc = cycles(g)
    for i in 1:length(c)
      @test length(filter(x -> length(x) == c[i][1], cc)) == c[i][2]
    end
    @test cc == [[1, 2, 3], [4, 5], [6, 7], [8, 9, 10], [11, 12, 13, 14, 15]]
  end
end

@testset "smaller_degree_permuation_group" begin
  c = cperm(1:3,4:6,7:9)
  G,_ = sub(symmetric_group(9), [c])
  H,iso = smaller_degree_permutation_representation(G)
  @test degree(H)<degree(G)
end

@testset "fixed_points" begin
  # Trivial group (only identity) -- all points fixed
  G = group((), 1:5)
  @test sort(fixed_points(G)) == [1, 2, 3, 4, 5]

  # Full symmetric group -- no points fixed
  G = symmetric_group(5)
  @test fixed_points(G) == Int[]

  # Subgroup generated by (1,2); 3, 4, 5 are fixed
  G = group(@perm((1,2)), 1:5)
  @test sort(fixed_points(G)) == [3, 4, 5]
end