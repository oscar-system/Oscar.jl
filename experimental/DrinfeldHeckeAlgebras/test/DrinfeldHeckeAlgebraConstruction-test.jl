@testset "DrinfeldHeckeAlgebras.DrinfeldHeckeAlgebraConstruction" begin
  @testset "create trivial drinfeld-hecke algebra for C2" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = drinfeld_hecke_algebra(G)
    
    @test is_trivial(A)
  end

  @testset "create drinfeld-hecke algebra from drinfeld-hecke form for C2" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    kappa_1 = matrix(QQ, [0 1; -1 0])
    kappa_g = matrix(QQ, [0 2//3; -2//3 0])
    forms = Dict(one(G) => kappa_1, G[1] => kappa_g)
    A = drinfeld_hecke_algebra(forms)
    
    @test !is_trivial(A)
  end

  @testset "can't create from empty forms" begin
    @test_throws ArgumentError drinfeld_hecke_algebra(Dict())
  end

  @testset "can't create from non-alternating forms" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    kappa_1 = matrix(QQ, [1 1; -1 0])
    forms = Dict(one(G) => kappa_1)
    
    @test_throws ArgumentError drinfeld_hecke_algebra(forms)
  end

  @testset "create generic drinfeld-hecke algebra for C2 over QQ" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    @test !is_trivial(A)
    @test length(parameters(A)) == 2
  end

  @testset "create generic drinfeld-hecke algebra for C2 over F5" begin
    R, _ = residue_field(ZZ, 5)
    G = matrix_group(matrix(R, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    @test !is_trivial(A)
    @test length(parameters(A)) == 2
  end

  @testset "create generic drinfeld-hecke algebra for C2 over F2" begin
    R, _ = residue_field(ZZ, 2)
    G = matrix_group(matrix(R, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    @test !is_trivial(A)
    @test length(parameters(A)) == 1
  end

  @testset "create generic drinfeld-hecke algebra for C2 over F3" begin
    R, _ = residue_field(ZZ, 3)
    G = matrix_group(matrix(R, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    @test !is_trivial(A)
    @test length(parameters(A)) == 2
  end

  @testset "create generic drinfeld-hecke algebra for C2 over algebraic closure of QQ" begin
    F = algebraic_closure(QQ)
    G = matrix_group(matrix(F, [-1 0;0 -1]))
    A = generic_drinfeld_hecke_algebra(G)
    
    @test !is_trivial(A)
    @test length(parameters(A)) == 2
  end

  @testset "create generic drinfeld-hecke algebra for symplectic doubling of Weyl group B2 over QQ" begin
    w1 = matrix(QQ,[-1 2; 0 1])
    w2 = matrix(QQ,[1 0; 1 -1])
    W = matrix_group([w1,w2])
    G = symplectic_doubling(W)
    A = generic_drinfeld_hecke_algebra(G)
    
    @test !is_trivial(A)
    @test length(parameters(A)) == 3
  end
end
