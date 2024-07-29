@testset "heap based sparse linear algebra" begin
  n = Oscar.HeapNode(1, ZZ(5))
  n.left = Oscar.HeapNode(2, ZZ(7))
  n.right = Oscar.HeapNode(2, ZZ(-3))
  n.right.right = Oscar.HeapNode(3, ZZ(8))
  v = Oscar.HeapSRow(ZZ, n)
  w = Oscar.HeapSRow(ZZ, [1 => ZZ(5), 2 => ZZ(4), 3 => ZZ(11), 3=>ZZ(-3)])
  @test !isdefined(v, :index_cache)
  @test v[3] == 8
  @test isdefined(v, :index_cache) && v.index_cache[3] == 8
  @test Oscar.tree_size(v.top) == 4
  @test v[2] == 4
  @test v.index_cache[2] == 4
  @test Oscar.tree_size(v.top) == 3

  u = v + w
  @test u == 2*v
end

@testset "Gauss algorithm on different matrix types" begin
  s = 4
  m = 10*s
  n = 20*s

  A = Oscar.HeapSMat(ZZ, 0, n)
  for i in 1:m
    l = [(rand(1:20*s), ZZ(rand(-10:10))) for i in 1:s];
    v = Oscar.HeapSRow(ZZ, l)
    push!(A, v)
  end
  AA = Oscar.sparse_matrix(A)

  B0, S0, T0 = Oscar.upper_triangular_form!(copy(A))
  B1, S1, T1 = Oscar.upper_triangular_form!(copy(AA))

  @test S0*A == B0
  @test T0*B0 == A

  S1 = Oscar.dense_matrix(Oscar.HeapSMat(S1))
  B1 = Oscar.dense_matrix(Oscar.HeapSMat(B1))
  T1 = Oscar.dense_matrix(Oscar.HeapSMat(T1))
  A1 = Oscar.dense_matrix(A)
  @test S1*A1 == B1
  @test T1*B1 == A1

  @test Oscar.HeapSMat(B1) == B0
  @test Oscar.HeapSMat(T1) == T0
  @test Oscar.HeapSMat(S1) == S0
end

@testset "Gauss algorithm on different matrix types over univariate polynomial rings" begin
  s = 1
  m = 10*s
  n = 20*s

  P, (x,) = polynomial_ring(QQ, [:x])
  QQt, t = QQ[:t]
  A = Oscar.HeapSMat(QQt, 0, n)
  for i in 1:m
    l = [(rand(1:20*s), evaluate(rand(P, 1:10, 1:10, 1:10), [t])) for i in 1:s];
    v = Oscar.HeapSRow(QQt, l)
    push!(A, v)
  end
  AA = Oscar.sparse_matrix(A)

  B0, S0, T0 = Oscar.upper_triangular_form!(copy(A))
  B1, S1, T1 = Oscar.upper_triangular_form!(copy(AA))

  @test S0*A == B0
  @test T0*B0 == A

  S1 = Oscar.dense_matrix(Oscar.HeapSMat(S1))
  B1 = Oscar.dense_matrix(Oscar.HeapSMat(B1))
  T1 = Oscar.dense_matrix(Oscar.HeapSMat(T1))
  A1 = Oscar.dense_matrix(A)
  @test S1*A1 == B1
  @test T1*B1 == A1

  @test Oscar.HeapSMat(B1) == B0
  @test Oscar.HeapSMat(T1) == T0
  @test Oscar.HeapSMat(S1) == S0
end

