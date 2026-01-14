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


@testset "Isometry Positivity Checker" begin
  Qb = algebraic_closure(QQ);
  R, x = QQ[:x]
  # Gram matrix
  G = QQ[0 1 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1; 1 -2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 -2 0 0 0 0 0 0 0 0 1 1 1 0 1 0 1; 0 0 0 -2 0 0 0 0 0 0 1 0 1 0 1 1 1 0; 0 0 0 0 -2 0 0 0 0 0 0 0 1 1 1 0 0 1; 0 0 0 0 0 -2 0 0 0 0 1 1 0 1 0 1 0 0; 0 0 0 0 0 0 -2 0 0 0 1 0 0 1 1 0 1 1; 0 0 0 0 0 0 0 -2 0 0 1 1 1 0 0 0 1 0; 0 0 0 0 0 0 0 0 -2 0 0 1 0 0 0 0 1 1; 0 0 0 0 0 0 0 0 0 -2 0 0 0 0 1 1 0 0; 1 0 0 1 0 1 1 1 0 0 -2 2 2 2 2 0 0 2; 1 0 1 0 0 1 0 1 1 0 2 -2 1 0 2 2 1 1; 1 0 1 1 1 0 0 1 0 0 2 1 -2 0 0 1 2 1; 1 0 1 0 1 1 1 0 0 0 2 0 0 -2 1 2 2 1; 1 0 0 1 1 0 1 0 0 1 2 2 0 1 -2 1 1 1; 1 0 1 1 0 1 0 0 0 1 0 2 1 2 1 -2 2 2; 1 0 0 1 0 0 1 1 1 0 0 1 2 2 1 2 -2 1; 1 0 1 0 1 0 1 0 1 0 2 1 1 1 1 2 1 -2]

  # isometry
  f = QQ[3 1 0 -1 0 1 1 -1 -1 -1 1 0 -1 1 0 0 0 0; 2 1 0 -1 0 0 0 -1 0 0 0 0 -1 1 0 0 0 0; 2 1 0 -1 0 1 1 -1 0 -1 1 0 -1 1 0 0 0 0; 2 1 0 -1 0 1 1 -1 -1 0 1 0 -1 1 0 0 0 0; 3 1 0 -1 0 1 0 -1 -1 -1 1 0 -1 1 0 0 0 0; 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 1 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 3 1 0 -1 0 0 1 -1 -1 -1 1 0 -1 1 0 0 0 0; 1 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 11 5 -3//2 -5//2 -1//2 1//2 1//2 -7//2 -3//2 -3//2 0 -1 -3 2 0 0 0 0; 6 2 0 -3 0 3 2 -2 -2 -2 2 1 -3 3 0 0 0 0; 0 -1 1//2 -1//2 1//2 3//2 3//2 -1//2 -1//2 -1//2 1 0 -1 2 0 0 0 0; 2 0 0 -1 0 2 2 -1 -1 -1 2 1 -2 2 0 0 0 0; 3 1 0 -1 0 1 2 -1 -1 -2 2 0 -1 2 0 -1 0 0; 6 2 -1 -2 0 1 2 -3 -2 -1 1 -1 -3 3 1 0 0 0; 14 7 -3//2 -5//2 -3//2 -1//2 -3//2 -3//2 -3//2 -5//2 0 0 -1 0 -1 -1 0 -1; 2 1 1 -1 0 1 1 -1 -1 -1 0 0 -1 2 -1 1 1 0]

  # check that f is an isometry
  @assert f*G*transpose(f)  == G

  # define the lattice
  B = QQ[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1//2 1//2 1//2 1//2 1//2 1//2 1//2 1//2 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]

  V = quadratic_space(QQ, G)
  Vf = quadratic_space_with_isometry(V, f)
  L = lattice(Vf, B)
  @assert lattice(Vf,basis_matrix(L)*f) == L
  tau = Oscar._get_tau(change_base_ring(ZZ, isometry(L)), Qb)

  C0 = x^14 + 2*x^13 + 5*x^12 + 8*x^11 + 11*x^10 + 14*x^9 + 15*x^8 + 16*x^7 + 15*x^6 + 14*x^5 + 11*x^4 + 8*x^3 + 5*x^2 + 2*x + 1

  bi_form = Oscar._get_bilinearform(L, Qb)

  v = eigenspace(isometry(L), tau)
  w = eigenspace(isometry(L), tau^(-1))
  h = Oscar._get_h(lattice(L), v, w, bi_form)

  @test bi_form(v,v) == 0
  @test bi_form(w,w) == 0
  @test iterate(Oscar._get_Cfancy(L, C0)) === nothing
  @test bi_form(h,h)>0
  Rh = Oscar._get_R(lattice(L), h)
  if Rh != []
    @test Oscar._check_R(Rh[1], v, w, bi_form) == false
  end
  @test isometry_is_positive(L, h)[1] == true
end
