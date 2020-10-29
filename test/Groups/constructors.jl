@testset "The groups Sym(n) and Alt(n)" begin

  for n = 5:8
    G = @inferred symmetric_group(n)
    A = @inferred alternating_group(n)
    @test degree(G) isa Integer
    @test degree(G) == n
    @test degree(A) == n
    @test isfinite(G)
    @test order(G) isa Int64
    @test exponent(G) isa Integer
    @test exponent(A) isa Integer
    @test exponent(G) == lcm(1:n)
    if n < 13 
      @test order(Int32, G) isa Int32 
    end
    @test order(G) == factorial(n)
    @test order(A) == factorial(n)/2
    @test order(BigInt, G) == factorial(n)
    @test gens(G) isa Vector{PermGroupElem}
    @test gens(A) isa Vector{PermGroupElem}
    @test ngens(G) == length(gens(G))
  end
  @test_throws ArgumentError symmetric_group(0)
  @test_throws ArgumentError alternating_group(-1)
  @test isalternating_group(omega_group(3,3))
  @test issymmetric_group(symmetric_group(PcGroup,4))
end

@testset "Special Constructors" begin

  @test isa(symmetric_group(5), PermGroup)
  
  @test isa(alternating_group(5), PermGroup)
    
  @test isa(dihedral_group(6), PcGroup)
  @test isa(dihedral_group(PermGroup, 6), PermGroup)
  
  @test isa(alternating_group(PcGroup,3), PcGroup)
  @test isa(symmetric_group(PcGroup,3), PcGroup)
  @test isisomorphic(symmetric_group(4), symmetric_group(PcGroup,4))[1]
  
  
  @test isquaternion_group(small_group(8, 4))
  @test small_group_identification(small_group(8, 4)) == (8, 4)
  @test isa(small_group(8, 4), PcGroup)
  @test isa(small_group(60, 5), PermGroup)
  
  @test isa(transitive_group(5, 5), PermGroup)
  
  @test isa(cyclic_group(5), PcGroup)
  @test isa(cyclic_group(PermGroup, 5), PermGroup)
  
  G = abelian_group(PcGroup,[2, 3])
  @test isa(G, PcGroup)
  @test iscyclic(G)
  G1 = abelian_group(PermGroup, [2, 3])
  @test isisomorphic(G, G1)[1]

# FIXME: a function `free_abelian_group` is not defined in GAPGroups, since it is already defined in Hecke
#=
  H = free_abelian_group(2)
  @test !isfinite(H)
  @test isabelian(H)
=#
  
  @test mathieu_group(10) isa PermGroup
  @test order(mathieu_group(10))==720

  F = free_group("x","y")
  @test F isa FPGroup
  @test !isfinite(F)
  @test !isabelian(F)

  F = free_group(2)
  @test F isa FPGroup
  
  F = free_group(["x","y"])
  @test F isa FPGroup
  
  F = free_group(3,"y")
  @test F isa FPGroup
  
  Q8 = quaternion_group(8)
  @test isa(Q8, PcGroup)
  
  gl = GL(2, 3)
  @test isa(gl, MatrixGroup)
  
  sl = SL(2, 3)
  @test isa(sl, MatrixGroup)
  
end
