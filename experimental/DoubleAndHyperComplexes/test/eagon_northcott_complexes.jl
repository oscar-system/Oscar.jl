@testset "Eagon-Northcott complexes and their restrictions" begin
  R, (x, y, z, w) = QQ[:x, :y, :z, :w]

  A = R[x y z; y z w]
  enc = Oscar.EagonNorthcottComplex(A)
  @test ngens(enc[0]) == 1
  @test ngens(enc[1]) == 3
  @test ngens(enc[2]) == 2
  @test is_zero(compose(map(enc, 2), map(enc, 1)))
  @test is_zero(homology(enc, 1)[1])
  
  F0 = Oscar.original_module(enc)
  I, _ = sub(F0, [F0[1], F0[2], x*F0[1]])
  ind = Oscar.InducedENC(enc, I)
  @test is_zero(homology(ind, 1)[1])

  m = 3
  n = 7
  r = 4
  A = matrix_space(R, m, n)([rand(R, 0:r, 0:r, 0:r) for i in 1:m, j in 1:n])
  enc = Oscar.EagonNorthcottComplex(A)

  for i in 2:n-m+1
    @test is_zero(compose(map(enc, i), map(enc, i-1)))
  end
  
  F0 = Oscar.original_module(enc)
  I, _ = sub(F0, [F0[1], F0[2], x*F0[3], y*F0[5], F0[7]])
  ind = Oscar.InducedENC(enc, I)
  for i in 2:m-n+1
    @test is_zero(compose(map(ind, i), map(ind, i-1)))
  end

  m = 2
  n = 4
  R, x = QQ[["x_$(i)_$j" for i in 1:m for j in 1:n]...]
  A = matrix_space(R, m, n)(x)
  enc = Oscar.EagonNorthcottComplex(A)
  
  for i in 2:n-m+1
    @test is_zero(compose(map(enc, i), map(enc, i-1)))
    @test is_zero(homology(enc, i-1)[1])
  end
  
  F0 = Oscar.original_module(enc)
  I, _ = sub(F0, [7*x[5]*F0[1], x[7]*F0[2] + F0[4], F0[n-1], F0[n]])
  ind = Oscar.InducedENC(enc, I)
  for i in 2:m-n+1
    @test is_zero(compose(map(ind, i), map(ind, i-1)))
    @test is_zero(homology(ind, i-1)[1])
  end

end
