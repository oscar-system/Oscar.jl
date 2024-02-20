@testset "Algebraic Shifting" begin
  K = simplicial_complex([[1, 2, 3], [1, 4]])
  k_ext = Oscar.exterior_shift(GF(2), K)
  k_sym = Oscar.symmetric_shift(QQ, K)
end
