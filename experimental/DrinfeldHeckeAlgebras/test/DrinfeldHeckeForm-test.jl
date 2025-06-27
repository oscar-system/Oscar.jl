@testset "DrinfeldHeckeAlgebras.DrinfeldHeckeForm" begin

  @testset "create zero drinfeld-hecke form" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    κ = drinfeld_hecke_form(G)
    @test is_zero(κ)
  end

  @testset "create parametrized drinfeld-hecke form" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    κ = generic_drinfeld_hecke_form(G)
    
    @test nparams(S) == 2
    
    S = base_ring(κ)
    t1 = S[1]
    t2 = S[2]
    
    g = one(G)
    h = G[1]
    
    κ_g = alternating_bilinear_form(matrix(S, [0 t1; -t1 0]))
    κ_h = alternating_bilinear_form(matrix(S, [0 t2; -t2 0]))
    
    @test κ[g] == κ_g
    @test κ[h] == κ_h
    @test nforms(κ) == 2
  end

  @testset "evaluate parameter" begin
    G = matrix_group(matrix(QQ, [-1 0;0 -1]))
    κ = generic_drinfeld_hecke_form(G)
    
    κ2 = evaluate_parameters(κ, [2,3])
    
    S = base_ring(κ2)
    g = one(G)
    h = G[1]

    κ_g = alternating_bilinear_form(matrix(S, [0 2; -2 0]))
    κ_h = alternating_bilinear_form(matrix(S, [0 3; -3 0]))
    
    @test κ2[g] == κ_g
    @test κ2[h] == κ_h
  end

  @testset "evaluate parameter not enough values" begin
    κ = generic_drinfeld_hecke_form(G)
    @test_throws ArgumentError evaluate_parameters(κ, [2])
  end

  @testset "evaluate parameter too many values" begin
    κ = generic_drinfeld_hecke_form(G)
    @test_throws ArgumentError evaluate_parameters(κ, [2,3,4])
  end

  @testset "evaluate parameter wrong element type" begin
    κ = generic_drinfeld_hecke_form(G)
    @test_throws ArgumentError evaluate_parameters(κ, [1, "hi"])
  end

  @testset "evaluate parameter not parametrized" begin
    κ = drinfeld_hecke_form(G)
    @test_throws ArgumentError evaluate_parameters(κ, [1, 2])
  end
      
  @testset "create drinfeld-hecke form from forms input" begin
    T = elem_type(typeof(QQ))
    forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
    g = one(G)
    forms[g] = MS([0 1;-1 0])
    
    κ = drinfeld_hecke_form(forms)
    
    @test !is_zero(κ)
    @test nforms(κ) == 1
    
    κ_g = alternating_bilinear_form(MS([0 1;-1 0]))
    
    @test κ[g] == κ_g
  end

  @testset "change drinfeld-hecke form" begin
    κ = drinfeld_hecke_form(G)
    
    @test is_zero(κ)
    
    T = elem_type(typeof(QQ))
    forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
    g = one(G)
    h = G[1]
    forms[g] = MS([0 3;-3 0])
    forms[h] = MS([0 -2;2 0])
    set_forms(κ, forms)
    
    @test !is_zero(κ)
    @test nforms(κ) == 2
    
    κ_g = alternating_bilinear_form(MS([0 3;-3 0]))
    κ_h = alternating_bilinear_form(MS([0 -2;2 0]))
    
    @test κ[g] == κ_g
    @test κ[h] == κ_h
  end

  @testset "can't create from empty forms" begin
    T = elem_type(typeof(QQ))
    forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
    
    @test_throws ArgumentError drinfeld_hecke_form(forms)
  end

  @testset "can't create from non-alternating forms" begin
    g = one(G)
    T = elem_type(typeof(QQ))
    forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
    forms[g] = MS([1 0;-1 0])
    
    @test_throws ArgumentError drinfeld_hecke_form(forms)
  end

  @testset "apply form" begin
    g = one(G)
    T = elem_type(typeof(QQ))
    forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
    forms[g] = MS([0 1;-1 0])
    κ = drinfeld_hecke_form(forms)
    RG = κ.group_algebra
    v1 = [QQ(1), QQ()]
    v2 = [QQ(), QQ(1)]
    @test κ(v1,v2) == RG(1)
  end

  @testset "parametrized drinfeld-hecke forms over various fields" begin
    @testset "algebraic closure of Q" begin
      F = algebraic_closure(QQ)
      G = matrix_group(matrix(F, [-1 0;0 -1]))
      κ = generic_drinfeld_hecke_form(G)
      @test ngens(base_ring(κ)) == 2
    end
  
    @testset "Z/5Z" begin
      F, _ = residue_field(ZZ, 5)
      G = matrix_group(matrix(F, [-1 0;0 -1]))
      κ = generic_drinfeld_hecke_form(G)
      @test ngens(base_ring(κ)) == 2
    end
  end

  @testset "drinfeld-hecke form over symmetric group S_3 acting on 3-dimensional QQ-VS" begin
    S3 = symmetric_group(3)
    mat_gens = [permutation_matrix(QQ, g) for g in gens(S3)]
    G = matrix_group(mat_gens)
    
    @testset "create zero drinfeld-hecke form" begin
      κ = drinfeld_hecke_form(G)
      @test is_zero(κ)
    end
  
    @testset "create parametrized drinfeld-hecke form" begin
      κ = generic_drinfeld_hecke_form(G)
      S = base_ring(κ)
  
      @test ngens(S) == 1
  
      t = S[1]
      MS = matrix_space(S,3,3)
  
      g = G([0 0 1; 1 0 0; 0 1 0])
      h = G([0 1 0; 0 0 1; 1 0 0])
  
      κ_g = alternating_bilinear_form(MS([0 t -t; -t 0 t; t -t 0]))
      κ_h = alternating_bilinear_form(MS([0 -t t; t 0 -t; -t t 0]))
  
      @test κ[g] == κ_g
      @test κ[h] == κ_h
      @test nforms(κ) == 2
    end
  
    @testset "evaluate parameter" begin
      κ = generic_drinfeld_hecke_form(G)
      κ = evaluate_parameters(κ, [-2])
  
      g = G(matrix(QQ, [0 0 1; 1 0 0; 0 1 0]))
      h = G(matrix(QQ, [0 1 0; 0 0 1; 1 0 0]))
  
      κ_g = alternating_bilinear_form(matrix(QQ, [0 -2 2; 2 0 -2; -2 2 0]))
      κ_h = alternating_bilinear_form(matrix(QQ, [0 2 -2; -2 0 2; 2 -2 0]))
  
      @test κ[g] == κ_g
      @test κ[h] == κ_h
      @test nforms(κ) == 2
    end
  
    @testset "create drinfeld-hecke form from input" begin
      T = elem_type(typeof(QQ))
      forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
      
      MS = matrix_space(QQ,3,3)
      
      g = G([0 0 1; 1 0 0; 0 1 0])
      h = G([0 1 0; 0 0 1; 1 0 0])
      
      forms[g] = MS([0 -1 1; 1 0 -1; -1 1 0])
      forms[h] = MS([0 1 -1; -1 0 1; 1 -1 0])
  
      κ = drinfeld_hecke_form(forms)
      
      @test !is_zero(κ)
      @test nforms(κ) == 2
      
      κ_g = alternating_bilinear_form(forms[g])
      κ_h = alternating_bilinear_form(forms[h])
            
      @test κ[g] == κ_g
      @test κ[h] == κ_h
    end
  end
end
