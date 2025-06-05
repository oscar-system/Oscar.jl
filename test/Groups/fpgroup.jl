@testset "fp group elements letters" begin
  F = @free_group(:F1, :F2)
  x = F1^5 * F2^-3
  l_1 = [1, 1, 1, 1, 1, -2, -2, -2]
  
  # test trivial cases
  @test letters(one(F)) == []
  @test F(Vector{Int}()) == one(F)
  
  # test general letters usage
  @test letters(x) == l_1
  @test F(l_1) == x

  # test for generator number check
  l_2 = [l_1 ; [3, 4]]
  @test_throws ArgumentError F(l_2)

  G, epi = quo(F, [F1^10, F2^10])
  
  # test trivial cases for quotient
  @test letters(one(G)) == []
  @test G(Vector{Int}()) == one(G)

  # test general letters usage for quotient
  @test letters(epi(x)) == l_1
  @test G(l_1) == epi(x)
end

@testset "fp group elements syllables" begin
  F = @free_group(:F1, :F2)
  x = F1^5 * F2^-3
  s_1 = [1 => 5, 2 => -3]
  
  # test trivial cases
  @test syllables(one(F)) == []
  @test F(Vector{Pair{Int, Int}}()) == one(F)

  # check general syllable usage
  @test syllables(x) == s_1
  @test F(s_1) == x

  # test for generator number check
  s_2 = [s_1 ; [3 => 1, 4 => 2]]
  @test_throws ArgumentError F(s_2)

  G, epi = quo(F, [F1^10, F2^10])
  
  # test trivial cases for quotient
  @test syllables(one(G)) == []
  @test G(Vector{Pair{Int, Int}}()) == one(G)

  # check general syllable usage for quotient
  @test syllables(epi(x)) == s_1
  @test G(s_1) == epi(x)
end
