@testset "double complexes of vector spaces" begin
  U = vector_space(QQ, 1)
  V = vector_space(QQ, 3)
  W = vector_space(QQ, 1)
  Z = vector_space(QQ, 0)

  phi = hom(V, U, [U[1], zero(U), zero(U)])
  psi = hom(W, V, [V[3]])
  inc_Z = hom(Z, W, elem_type(W)[])
  pr_Z = hom(U, Z, [zero(Z)])

  C = chain_complex([inc_Z, psi, phi, pr_Z])


  D = tensor_product(C, C)

  for i in left_bound(D)+1:right_bound(D)
    for j in lower_bound(D)+2:upper_bound(D)
      @test iszero(compose(vertical_map(D, i, j), vertical_map(D, i, j-1)))
    end
  end

  for j in lower_bound(D)+1:upper_bound(D)
    for i in left_bound(D)+2:right_bound(D)
      @test iszero(compose(horizontal_map(D, i, j), horizontal_map(D, i-1, j)))
    end
  end

  for j in lower_bound(D)+1:upper_bound(D)
    for i in left_bound(D)+1:right_bound(D)
      @test compose(horizontal_map(D, i, j), vertical_map(D, i-1, j)) == compose(vertical_map(D, i, j), horizontal_map(D, i, j-1))
    end
  end

  tot = total_complex(D)
  for i in map_range(tot)
    i == last(map_range(tot)) && continue
    if i != first(map_range(tot)) && i != last(map_range(tot)) + 1
      @test !iszero(map(tot, i)) || !iszero(map(tot, i-1))
    end
    @test iszero(compose(map(tot, i), map(tot, i-1)))
  end

end
