@testset "vinberg_algorithm" begin
  # example of Vinberg's original paper "Some arithmetical discrete groups in Lobacevskii spaces"
  Q = diagonal_matrix(ZZ, [1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1])
  v0 = ZZ[-1 0 0 0 0 0 0 0 0 0 0]
  alpha = 9
  roots = vinberg_algorithm(Q, ZZ(alpha); v0)
  @test length(roots) == 12
  roots[12] = abs.(roots[12])
  @test roots[12] == ZZ[3 1 1 1 1 1 1 1 1 1 1]

  # rank 3 matrix generating a lot of roots with big distance to v0. Jacob Grosek's paper "Fundamental Domains for Hyperbolic Tesselations" 
  v0 = ZZ[-10 0 2]
  Q = diagonal_matrix(ZZ, [1, -1, -13])
  d = ZZ.([-1, -2, -13])
  alpha = 6000
  roots = vinberg_algorithm(Q, ZZ(alpha), v0 = v0, root_lengths = d)
  l = length(roots)
  @test l == 8
  cusp_counter = 0
  for i in 1:l
    for j in i+1:l
      vi = roots[i]
      vj = roots[j]
      k = (vi*Q*transpose(vj))[1, 1]
      if k > 0 && (vi*Q*transpose(vi))[1, 1] * (vj*Q*transpose(vj))[1, 1] == k*k  # formula for Coxeter diagrams (in dimension 2 parallel hyperplanes are cusps)
        cusp_counter = cusp_counter + 1
      end
    end
  end
  @test cusp_counter == 2

  # non diagonal example of Rudolf Scharlau's paper "On the classification of arithmetic reflection groups on hyperbolic-3-spaces"
  Q = matrix(ZZ, 4, 4, [0 -1 0 0; -1 0 0 0; 0 0 -2 -1; 0 0 -1 -2])
  alpha = 500
  roots = vinberg_algorithm(Q, ZZ(alpha))
  l = length(roots)
  @test l == 4
  edge_counter = 0
  for i in 1:l
    for j in i+1:l
      vi = roots[i]
      vj = roots[j]
      k = (vi*Q*transpose(vj))[1, 1]
      if 1 >= (k*k)//((vi*Q*transpose(vi))[1, 1] * (vj*Q*transpose(vj))[1, 1]) # formula for Coxeter diagrams
        edge_counter = edge_counter + 1
      end
    end
  end
  @test edge_counter == 6

  # basic example of Jacob Grosek's paper "Fundamental Domains for Hyperbolic Tesselations" with reflection on one of the roots
  Q = diagonal_matrix(ZZ, [1, -1, -3])
  v0 = ZZ[1 0 0]
  alpha = 4
  roots = vinberg_algorithm(Q, ZZ(alpha), v0 = v0, direction_vector = ZZ[0 9 6]) 
  @test roots[1] == ZZ[0 -1 0]
  @test roots[2] == ZZ[0 0 -1]
  @test roots[3] == ZZ[1 0 1]
  @test roots[4] == ZZ[3 3 1]
  e = ZZ.(roots[3]) 
  v0_ref = v0 - (numerator((2*(v0*Q*transpose(e))[1, 1])//((e*Q*transpose(e))[1, 1])) * e)  # we want to reflect v0 on e
  roots_after_reflection = vinberg_algorithm(Q, ZZ(alpha), v0 = v0_ref, direction_vector = ZZ[-3 1 -2]) # roots we get by replacing v0 with the reflection of v0
  @test roots_after_reflection[3] == ((-1) * roots[3])
  @test roots_after_reflection[4] == roots[4]
end
