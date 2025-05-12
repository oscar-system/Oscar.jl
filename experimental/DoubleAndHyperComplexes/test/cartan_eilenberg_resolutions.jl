@testset "Cartan-Eilenberg resolutions" begin
  R, (x, y, z, w) = QQ[:x, :y, :z, :w]

  A = R[x y z; y z w]

  R1 = free_module(R, 1)
  I, inc = sub(R1, [a*R1[1] for a in minors(A, 2)])
  M = cokernel(inc)
  R4 = free_module(R, 4)
  theta = sum(a*g for (a, g) in zip(gens(R), gens(R4)); init=zero(R4))
  K = koszul_complex(Oscar.KoszulComplex, theta)

  comp = tensor_product(K, Oscar.ZeroDimensionalComplex(M))
  res = Oscar.CartanEilenbergResolution(comp);
  tot = total_complex(res);
  tot_simp = simplify(tot);

  res_M, _ = free_resolution(Oscar.SimpleFreeResolution, M)
  comp2 = tensor_product(K, res_M)
  tot2 = total_complex(comp2)
  tot_simp2 = simplify(tot2);

  @test [ngens(tot_simp[i]) for i in 0:5] == [ngens(tot_simp2[i]) for i in 0:5]
end

