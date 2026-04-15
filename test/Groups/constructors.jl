@testset "The groups Sym(n) and Alt(n)" begin

  for n = 5:8
    G = @inferred symmetric_group(n)
    A = @inferred alternating_group(n)

    @test degree(G) isa Integer
    @test degree(G) == n
    @test degree(A) == n

    @test length(moved_points(G)) == n
    nmp = number_of_moved_points(G)
    @test nmp == n
    @test nmp isa Int

    @test is_finite(G)

    @test order(G) isa ZZRingElem
    @test order(G) == factorial(n)
    @test 2 * order(A) == factorial(n)
    @testset "order with type" begin
      for T in [Int, BigInt, ZZRingElem]
        @test order(T, G) == factorial(n)
        @test order(T, G) isa T
      end
      if n < 13
        @test order(Int32, G) isa Int32
      end
    end

    @test exponent(G) isa ZZRingElem
    for T in [Int, BigInt, ZZRingElem]
      @test exponent(T, A) isa T
    end
    @test exponent(G) == lcm(1:n)

    @test gens(G) isa Vector{PermGroupElem}
    @test gens(A) isa Vector{PermGroupElem}
    @test ngens(G) == ngens(G)
  end

  @test_throws ArgumentError symmetric_group(0)
  @test_throws ArgumentError alternating_group(-1)

  @test is_natural_alternating_group(alternating_group(4))
  @test !is_natural_alternating_group(omega_group(3,3))
  @test is_isomorphic_to_alternating_group(alternating_group(4))
  @test is_isomorphic_to_alternating_group(omega_group(3,3))
  @test !is_isomorphic_to_alternating_group(symmetric_group(4))

  @test is_natural_symmetric_group(symmetric_group(4))
  @test !is_natural_symmetric_group(PcGroup(symmetric_group(4)))
  @test is_isomorphic_to_symmetric_group(symmetric_group(4))
  @test is_isomorphic_to_symmetric_group(PcGroup(symmetric_group(4)))
  @test !is_isomorphic_to_symmetric_group(alternating_group(4))
end

@testset "Special Constructors" begin

  @test isa(symmetric_group(5), PermGroup)
  
  @test isa(alternating_group(5), PermGroup)
    
  @test isa(dihedral_group(6), PcGroup)
  @test isa(dihedral_group(PermGroup, 6), PermGroup)

  @test is_quaternion_group(small_group(8, 4))
  @test ! is_quaternion_group(small_group(12, 3))
  @test is_dicyclic_group(small_group(8, 4))
  @test ! is_dicyclic_group(small_group(13, 1))

  @test small_group_identification(small_group(8, 4)) == (8, 4)
  @test isa(small_group(8, 4), PcGroup)
  @test isa(small_group(60, 5), PermGroup)
  
  @test isa(transitive_group(5, 5), PermGroup)
  
  @test isa(cyclic_group(5), PcGroup)
  @test isa(cyclic_group(PermGroup, 5), PermGroup)
  @test_throws ArgumentError cyclic_group(-1)
  @test_throws ArgumentError cyclic_group(PermGroup, -1)

  @test_throws ArgumentError cyclic_generator(symmetric_group(3))

  @test isa(elementary_abelian_group(27), PcGroup)
  @test isa(elementary_abelian_group(PermGroup, 27), PermGroup)
  @test isa(elementary_abelian_group(FinGenAbGroup, 27), FinGenAbGroup)
  @test_throws ArgumentError elementary_abelian_group(6)
  @test_throws ArgumentError elementary_abelian_group(PermGroup, 6)

  for p in [1, next_prime(2^62), next_prime(ZZRingElem(2)^66)]
    g = cyclic_group(p)
    @test is_cyclic(g)
    @test is_elementary_abelian(g)
    @test order(cyclic_generator(g)) == order(g)
    @test !is_dihedral_group(g)
    @test is_finite(g)
    @test order(g) == p
  end

  for p in [next_prime(2^62), next_prime(ZZRingElem(2)^66)]
    n = 2*ZZRingElem(p)
    g = dihedral_group(n)
    @test !is_cyclic(g)
    #@test is_dihedral_group(g)
    @test is_finite(g)
    @test order(g) == n

    n = p^2
    g = elementary_abelian_group(n)
    @test is_abelian(g)
    @test is_finite(g)
    @test !is_cyclic(g)
    @test exponent(g) == p
    @test order(g) == n
  end

  g = cyclic_group(PosInf())
  @test is_cyclic(g)
  @test !is_finite(g)
  @test_throws InfiniteOrderError{PcGroup} order(g)

  g = dihedral_group(PosInf())
  @test !is_cyclic(g)
  @test !is_finite(g)
  @test_throws InfiniteOrderError{PcGroup} order(g)

  G = abelian_group(PcGroup,[2, 3])
  @test isa(G, PcGroup)
  @test is_cyclic(G)
  G1 = abelian_group(PermGroup, [2, 3])
  @test is_isomorphic(G, G1)
# G = abelian_group(PcGroup, [ZZ(2)^70])
  G = abelian_group(SubPcGroup, [ZZ(2)^70])

# FIXME: a function `free_abelian_group` is not defined in GAPGroups, since it is already defined in Hecke
#=
  H = free_abelian_group(2)
  @test !is_finite(H)
  @test is_abelian(H)
=#
  
  @test mathieu_group(10) isa PermGroup
  @test order(mathieu_group(10))==720
  @test_throws ArgumentError mathieu_group(-1)
  @test_throws ArgumentError mathieu_group(8)
  @test_throws ArgumentError mathieu_group(13)
  @test_throws ArgumentError mathieu_group(20)
  @test_throws ArgumentError mathieu_group(25)

  @testset "free_group($args)" for args in [
          ("x","y"), (:x,:y), ('x','y'),
          (["x","y"],), ([:x,:y],), (['x','y'],),
          (2, ), (2, "x"), (2, :x), (2, 'x'),
      ]
    F = free_group(args...)
    @test F isa FPGroup
    @test_throws InfiniteOrderError{FPGroup} order(F)
    @test_throws ArgumentError index(F, trivial_subgroup(F)[1])
    @test_throws MethodError degree(F)
    @test !is_finite(F)
    @test !is_abelian(F)
    @test ngens(F) == 2
    @test length(gens(F)) == 2
  end

  F = free_group(3,"y")
  @test F isa FPGroup
  
  F = free_group(3,:y)
  @test F isa FPGroup

  @test_throws ArgumentError free_group(-1)

  Q8 = quaternion_group(8)
  @test isa(Q8, PcGroup)

  Dic12 = dicyclic_group(12)
  @test isa(Dic12, PcGroup)
  
  gl = GL(2, 3)
  @test isa(gl, MatrixGroup)
  
  sl = SL(2, 3)
  @test isa(sl, MatrixGroup)
  
end

@testset "Extraspecial groups" begin
   @testset for p in [2, 3, 5], n in [1, 2], type in [:+, :-]
      G = extraspecial_group(p, n, type)
      @test is_extraspecial_group(G)
      @test G isa PcGroup

      for T in [PcGroup, SubPcGroup, PermGroup, FPGroup, SubFPGroup]
         G = extraspecial_group(T, p, n, type)
         @test is_extraspecial_group(G)
         @test G isa T
         @test (p == 2) ? (exponent(G) == 4) :
                          ((type == :+) ? exponent(G) == p : exponent(G) == p^2)
      end
   end
   @test_throws ArgumentError extraspecial_group(2, 0, :+)
   @test_throws ArgumentError extraspecial_group(4, 1, :+)
   @test_throws MethodError extraspecial_group(2, 1, "+")
end

@testset "Classical groups" begin
   @testset for n in [2,5], q in [4,9]
      G = GL(n,q)
      S = SL(n,q)
      @test G==general_linear_group(n,q)
      @test S==special_linear_group(n,q)
      @test order(S)==prod(BigInt[q^n-q^i for i in 0:(n-1)])÷(q-1)
      @test index(G,S)==q-1
   end

   @testset for n in 1:3, q in [2,3]
      @test unitary_group(n,q)==GU(n,q)
      @test special_unitary_group(n,q)==SU(n,q)
      @test index(GU(n,q),SU(n,q))==q+1
   end

   @testset for n in [2,4,6], q in [4,9]
      @test symplectic_group(n,q)==Sp(n,q)
   end

   @testset for q in [3,4]
      @testset for n in [4,6], e in [+1,-1]
         @test GO(e,n,q)==orthogonal_group(e,n,q)
         @test SO(e,n,q)==special_orthogonal_group(e,n,q)
         @test index(GO(e,n,q), SO(e,n,q)) == gcd(2, q-1)
         @test index(SO(e,n,q), omega_group(e,n,q)) == 2
         @test index(GO(e,n,q), omega_group(e,n,q)) == 2 * gcd(2, q-1)
      end
      @testset for n in [3,5]
         @test GO(n,q)==orthogonal_group(n,q)
         @test SO(n,q)==special_orthogonal_group(n,q)
         @test index(GO(n, q), SO(n, q)) == gcd(2, q-1)
         @test index(SO(n, q), omega_group(n, q)) == gcd(2, q-1)
      end
   end

   @test order(omega_group(+1,4,3))==288
   @test order(omega_group(-1,4,3))==360
   @test order(omega_group(3,3))==12
end
