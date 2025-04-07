@testset "DrinfeldHeckeAlgebras.DrinfeldHeckeForm" begin
  @testset "drinfeld-hecke form over cyclic group C_2 acting on 2-dimensional QQ-VS" begin
    MS = matrix_space(QQ,2,2)
    G = matrix_group(matrix(MS([-1 0;0 -1])))
    
    @testset "create zero drinfeld-hecke form" begin
      κ = drinfeld_hecke_form(G)
      @test is_zero(κ)
    end
  
    @testset "create parametrized drinfeld-hecke form" begin
      κ = parametrized_drinfeld_hecke_form(G)
      
      S = base_ring(κ)
      G_S = change_base_ring(S, G)
      
      @test ngens(S) == 2
      
      t1 = S[1]
      t2 = S[2]
      
      g = G_S(matrix(S, [-1 0; 0 -1]))
      h = G_S(matrix(S, [1 0; 0 1]))
      
      κ_g = alternating_bilinear_form(matrix(S, [0 t1; -t1 0]))
      κ_h = alternating_bilinear_form(matrix(S, [0 t2; -t2 0]))
      
      @test κ[g] == κ_g
      @test κ[h] == κ_h
      @test nforms(κ) == 2
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
  
    @testset "change drinfeld-hecke form locally" begin
      κ = drinfeld_hecke_form(G)
      
      @test is_zero(κ)
      
      g = one(G)
      κ[g] = MS([0 1;-1 0])
      
      @test !is_zero(κ)
      @test nforms(κ) == 1
      
      κ_g = alternating_bilinear_form(MS([0 1;-1 0]))
      
      @test κ[g] == κ_g
    end
  
    @testset "change drinfeld-hecke form globally" begin
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
      forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
      
      @test_throws ErrorException drinfeld_hecke_form(forms)
    end
  
    @testset "can't create from non-alternating forms" begin
      κ = drinfeld_hecke_form(G)
      
      @test_throws ErrorException κ[g] = MS([1 0;-1 0])
    end
  
    @testset "apply form" begin
      κ = drinfeld_hecke_form(G)
      κ[g] = MS([0 1;-1 0])
      RG = κ.group_algebra
      v1 = [QQ(1), QQ()]
      v2 = [QQ(), QQ(1)]
      @test κ(v1,v2) == RG(1)
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
      κ = parametrized_drinfeld_hecke_form(G)
      S = base_ring(κ)
      G_S = change_base_ring(S, G)
  
      @test ngens(S) == 1
  
      t1 = S[1]
      MS = matrix_space(S,3,3)
  
      g = G_S([0 0 1; 1 0 0; 0 1 0])
      h = G_S([0 1 0; 0 0 1; 1 0 0])
  
      κ_g = alternating_bilinear_form(MS([0 -t1 t1; t1 0 -t1; -t1 t1 0]))
      κ_h = alternating_bilinear_form(MS([0 t1 -t1; -t1 0 t1; t1 -t1 0]))
  
      @test κ[g] == κ_g
      @test κ[h] == κ_h
      @test nforms(κ) == 2
    end
  
    @testset "substitute parameter" begin
      κ = parametrized_drinfeld_hecke_form(G)
      κ = substitute_parameter(κ, 1, 2)
  
      g = G(matrix(QQ, [0 0 1; 1 0 0; 0 1 0]))
      h = G(matrix(QQ, [0 1 0; 0 0 1; 1 0 0]))
  
      κ_g = alternating_bilinear_form(matrix(QQ, [0 -2 2; 2 0 -2; -2 2 0]))
      κ_h = alternating_bilinear_form(matrix(QQ, [0 2 -2; -2 0 2; 2 -2 0]))
  
      @test κ[g] == κ_g
      @test κ[h] == κ_h
      @test nforms(κ) == 2
    end
  
    @testset "substitute parameter index out of bounds" begin
      κ = parametrized_drinfeld_hecke_form(G)
      @test_throws ErrorException substitute_parameter(κ, 2, 2)
    end
  
    @testset "substitute parameter wrong element type" begin
      κ = parametrized_drinfeld_hecke_form(G)
      @test_throws ErrorException substitute_parameter(κ, 1, "hi")
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
  
    @testset "change drinfeld-hecke form locally not possible" begin
      κ = drinfeld_hecke_form(G)
      @test is_zero(κ)
  
      g = G([0 0 1; 1 0 0; 0 1 0])
      MS = matrix_space(QQ,3,3)
  
      @test_throws ErrorException κ[g] = MS([0 -1 1; 1 0 -1; -1 1 0])
    end
  end
end
