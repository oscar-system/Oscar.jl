using Random

RNG = Random.MersenneTwister(42)

"""
    randpoly(R::Ring,coeffs=0:9,max_exp=4,max_terms=8)

Return a random Polynomial from the Polynomial Ring `R` with coefficients in `coeffs`
with exponents between `0` and `max_exp` und between `0` and `max_terms` terms.
"""
function randpoly(R::Oscar.Ring,coeffs=0:2,max_exp=2,max_terms=3)
  n = nvars(R)
  K = base_ring(R)
  E = [[Random.rand(RNG,0:max_exp) for i=1:n] for j=1:max_terms]
  C = [K(Random.rand(RNG,coeffs)) for i=1:max_terms]
  M = MPolyBuildCtx(R)
  for i=1:max_terms
    push_term!(M,C[i],E[i])
  end
  return finish(M)
end

@testset "Modules: Constructors" begin
  Oscar.set_seed!(235)
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  F = FreeMod(R,3)
  v = [x, x^2*y+z^3, R(-1)]
  @test v == Vector(F(v))

  M = sub_object(F, [F(v), F([z, R(1), R(0)])])
  N = quo_object(M, [SubquoModuleElem([x+y^2, y^3*z^2+1], M)])
  AN, ai = ambient_module(N, :with_morphism)
  @test AN.quo === N.quo
  for i=1:ngens(N)
    @test AN(repres(N[i])) == ai(N[i])
  end

  G = FreeMod(R,2)
  @test F(v) in F
  @test !(F(v) in G)
  @test (F(v) + F([z, R(1), R(0)])) in M
  @test !(F([R(1), R(0), R(0)]) in M)
  @test N[1] in M

  @test gen(F, 1) == deepcopy(gen(F, 1))
  @test gen(M, 1) == deepcopy(gen(M, 1))
  @test gen(N, 1) == deepcopy(gen(N, 1))

  M = SubquoModule(F, [x*F[1]])
  N = SubquoModule(F, [y*F[1]])
  G = FreeMod(R,3,"f")
  M_2 = SubquoModule(G, [x*G[1]])
  @test !is_canonically_isomorphic(M,N)
  is_iso, phi = is_canonically_isomorphic_with_map(M,M_2)
  @test is_iso
  @test is_welldefined(phi)
  @test is_bijective(phi)
end

@testset "Submodule Membership" begin
    R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
    F = FreeMod(R, 3)
    gens_submodule = [x*F[1], 3*y*F[2]]
    S, _ = sub(F, gens_submodule)
    x = x*y*F[1]+3*y^2*F[2]
    @test in(x, S)
    coord = coordinates(x, S)
    @test coord == sparse_row(R, [1, 2], [y, y])
end

@testset "Modules: Simplify elements of subquotients" begin
  Oscar.set_seed!(235)
    R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
    F1 = free_module(R, 3)
    F2 = free_module(R, 1)
    G = free_module(R, 2)
    V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
    V2 = [z*G[2]+y*G[1]]
    a1 = hom(F1, G, V1)
    a2 = hom(F2, G, V2)
    M = subquotient(a1,a2)
    m3 = x*M[1]+M[2]+x*M[3]
    @test repres(simplify(m3)) == x*G[1] + (y - z)*G[2]
end

@testset "Intersection of modules" begin
  Oscar.set_seed!(235)
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])

  A1 = R[x y;
        2*x^2 3*y^2]
  A2 = R[x^3 x^2*y;
        (2*x^2+x) (2*y^3+y)]
  B = R[4*x*y^3 (2*x+y)]
  F2 = FreeMod(R,2)

  M1 = SubquoModule(F2, A1, B)
  M2 = SubquoModule(F2, A2, B)

  A1 = R[x R(1); y x^2]
  A2 = R[x y]
  B1 = R[x^2 y^2]
  res = R[(x*y^2) (y^3-x*y+y^2); (-x*y) (-x^2*y^2+x*y^3-x^2*y+y^3-x*y)]

  P,i = intersect(M1,M1)
  @test M1 == P
  @test image(i)[1] == M1

  P,i = intersect(SubquoModule(F2, A1, B1), SubquoModule(F2, A2, B1))
  @test is_welldefined(i)
  @test is_injective(i)
  @test SubquoModule(F2, res, B1) == P #macaulay2

  A1 = R[x y; x^2 y^2]
  A2 = R[x+x^2 y+y^2]
  P,i = intersect(SubquoModule(F2, A1,B1), SubquoModule(F2, A2,B1))
  @test is_welldefined(i)
  @test is_injective(i)
  @test SubquoModule(F2, A2, B1) == P

  #Test that no obvious zeros are in the generator set
  F = free_module(R,1)
  AM = R[x;]
  BM = R[x^2; y^3; z^4]
  M = SubquoModule(F, AM, BM)
  AN = R[y;]
  BN = R[x^2; y^3; z^4]
  N = SubquoModule(F, AN, BN)
  P,_ = intersect(M, N)
  for g in gens(P)
    @test !iszero(ambient_representative(g))
  end


  F = FreeMod(R, 1)
  I, _ = sub(F, [F[1]])
  K, _ = sub(F, [zero(R)*F[1]])
  I, _ = intersect(I, K)
  @test iszero(I)
end

@testset "module quotients" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);
  F = free_module(R, 1);
  V = [x^2*F[1]; y^3*F[1]; z^4*F[1]];
  N, _ = sub(F, V);
  @test iszero(annihilator(N))
  
  R, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);
  F = graded_free_module(R, 1);
  B = R[x^2; y^3; z^4];
  AM = R[x;];
  M = SubquoModule(F, AM, B)
  AN = R[y;];
  N = SubquoModule(F, AN, B)
  L = quotient(M, N)
  @test ngens(L) == 3
end

@testset "Presentation" begin
  Oscar.set_seed!(235)

  # over Integers
  R, (x,y,z) = polynomial_ring(ZZ, [:x, :y, :z])
  generator_matrices = [R[x x^2*y; y y^2*x^2], R[x x^2*y; y y^2*x^2], R[x y; x^2*y y^2*x], R[x+R(1) x^10; x^2+y^2 y^4-x^4], R[42*x*y 7*x^2; 6*x 9*y^2]]
  relation_matrices  = [R[x^3 y^4], R[x^3 y^4; x^2*y x*y^2], R[x*y^2 x; y^3 x*y^3], R[x x*y], R[3*x*y 7*x*y; 42 7]]
  true_pres_matrices = [R[x^5*y-y^4 -x^5+x*y^3], R[-x^2*y-x*y-y x^2+x; x*y^4-x^3*y+y^4-x^2*y -x*y^3+x^3; -y^4+x^2*y x*y^3-x^3; x^2*y^3+x*y^3+y^3 -x^2*y^2-x*y^2], R[-x*y R(1); -x^2*y^5+x*y^3 R(0)], R[x^5-x*y^4+x^3*y+x*y^3 x^11-x^2*y-x*y; -x^5*y^7+x*y^11+2*x^5*y^6-x^3*y^8-3*x*y^10+2*x^5*y^5+2*x^3*y^7-3*x^5*y^4+2*x^3*y^6+5*x*y^8+2*x^5*y^3-3*x^3*y^5-5*x*y^7+2*x^3*y^4+2*x*y^6 -x^11*y^7+2*x^11*y^6+2*x^11*y^5-3*x^11*y^4+2*x^11*y^3+x^2*y^8-2*x^2*y^7+x*y^8-2*x^2*y^6-2*x*y^7+3*x^2*y^5-2*x*y^6-2*x^2*y^4+3*x*y^5-2*x*y^4], R[-13*y R(0); -377 -2639*x; -39 -273*x; -13*y-39 -273*x; -13 -91*x; y^2-42*x -294*x^2+21*x*y; 9*y^2-x+26*y+78 -7*x^2+189*x*y+546*x; -y^2+3*x 21*x^2-21*x*y]]
  for (A,B,true_pres_mat) in zip(generator_matrices, relation_matrices, true_pres_matrices)
    SQ = SubquoModule(A,B)
    pres_mat = generator_matrix(present_as_cokernel(SQ).quo)
    F = FreeMod(R,ncols(pres_mat))
    @test cokernel(F,pres_mat) == cokernel(F,true_pres_mat)

    pres_SQ, i = present_as_cokernel(SQ, :both)
    p = i.inverse_isomorphism
    @test is_welldefined(i)
    @test is_welldefined(p)
    @test is_bijective(i)
    @test is_bijective(p)
  end

  # issue 3797
  R, (x,y) = graded_polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [x, y])
  FIm = free_resolution(I, algorithm = :mres)
  @test is_graded(FIm)

  # over Rationals
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  generator_matrices = [R[x x^2*y; y y^2*x^2], R[x x^2*y; y y^2*x^2], R[x y; x^2*y y^2*x], R[x+R(1) x^10; x^2+y^2 y^4-x^4], R[x+R(1) x^10; x^2+y^2 y^4-x^4]]
  relation_matrices  = [R[x^3 y^4], R[x^3 y^4; x^2*y x*y^2], R[x*y^2 x; y^3 x*y^3], R[x+y x+y; x*y x*y], R[x+y x+y; y x; x y]]
  true_pres_matrices = [R[x^5*y-y^4 -x^5+x*y^3], R[-x^2*y-x*y-y x^2+x; -x^2*y^4+x^4*y x^2*y^3-x^4], R[-x*y R(1); -x^2*y^5+x*y^3 R(0)], R[-x^4+y^4-x^2-y^2 -x^10+x+R(1)], R[R(0) -x^2+y^2; -2*x^2 -x^9*y+x+1; -2*x*y^2 -x^8*y^3+y^2+x; -2*x*y^2+2*y^2 -x^8*y^3+x^7*y^3+y^2-1]]
  for (A,B,true_pres_mat) in zip(generator_matrices, relation_matrices, true_pres_matrices)
    SQ = SubquoModule(A,B)
    pres_mat = generator_matrix(present_as_cokernel(SQ).quo)
    F = FreeMod(R,ncols(pres_mat))
    @test cokernel(F,pres_mat) == cokernel(F,true_pres_mat)

    pres_SQ, i = present_as_cokernel(SQ, :both)
    p = i.inverse_isomorphism
    @test is_welldefined(i)
    @test is_welldefined(p)
    @test is_bijective(i)
    @test is_bijective(p)
  end

	R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
	A = R[x; y]
	B = R[x^2; x*y; y^2; z^4]
	M = SubquoModule(A, B)
	free_res = free_resolution(M, length=1)
  @test is_complete(free_res) == false
  free_res[2]
  @test length(free_res.C.maps) == 4
	@test free_res[3] == free_module(R, 2)
	@test free_res[4] == free_module(R, 0)
	@test free_res[100] == free_module(R, 0)
  @test is_complete(free_res) == true
	free_res = free_resolution(M)
  @test get_attribute(free_res.C, :algorithm) == :fres
	@test all(iszero, homology(free_res.C))
	free_res = free_resolution_via_kernels(M)
	@test all(iszero, homology(free_res))
	free_res = free_resolution(M, algorithm = :mres)
  @test get_attribute(free_res.C, :algorithm) == :mres
	@test all(iszero, homology(free_res.C))
	free_res = free_resolution(M, algorithm = :nres)
	@test all(iszero, homology(free_res.C))

  N = SubquoModule(R[x+2*x^2; x+y], R[z^4;])
  tensor_resolution = tensor_product(N,free_res)
  @test chain_range(tensor_resolution) == chain_range(free_res)
  for i in Hecke.map_range(tensor_resolution)
    f = map(free_res,i)
    M_i = domain(f)
    tensored_f = map(tensor_resolution,i)
    to_pure_tensors_i = get_attribute(domain(tensored_f),:tensor_pure_function)
    to_pure_tensors_i_plus_1 = get_attribute(codomain(tensored_f), :tensor_pure_function)
    for (n,mi) in zip(gens(N),gens(M_i))
      @test tensored_f(to_pure_tensors_i((n,mi))) == to_pure_tensors_i_plus_1(n,f(mi))
    end
  end

  N = SubquoModule(R[x+2*x^2*z; x+y-z], R[z^4;])
  tensor_resolution = tensor_product(free_res,N)
  @test chain_range(tensor_resolution) == chain_range(free_res)
  for i in Hecke.map_range(tensor_resolution)
    f = map(free_res,i)
    M_i = domain(f)
    tensored_f = map(tensor_resolution,i)
    to_pure_tensors_i = get_attribute(domain(tensored_f),:tensor_pure_function)
    to_pure_tensors_i_plus_1 = get_attribute(codomain(tensored_f), :tensor_pure_function)
    for (mi,n) in zip(gens(M_i),gens(N))
      @test tensored_f(to_pure_tensors_i((mi,n))) == to_pure_tensors_i_plus_1(f(mi),n)
    end
  end

  N = SubquoModule(R[x+2*x^2; x+y], R[z^4;])
  hom_resolution = hom(N,free_res)
  @test chain_range(hom_resolution) == chain_range(free_res)
  for i in Hecke.map_range(hom_resolution)
    f = map(free_res,i)
    hom_f = map(hom_resolution,i)
    hom_N_M_i = domain(hom_f)
    for v in gens(hom_N_M_i)
      @test element_to_homomorphism(hom_f(v)) == element_to_homomorphism(v)*f
    end
  end

  N = SubquoModule(R[x+2*x^2; x+y], R[z^4; x^2-y*z])
  hom_resolution = hom(free_res,N)
  @test last(chain_range(hom_resolution)) == first(chain_range(free_res))
  @test first(chain_range(hom_resolution)) == last(chain_range(free_res))
  for i in Hecke.map_range(hom_resolution)
                f = map(free_res,i+1)         #f[i]: M[i] -> M[i-1]
                hom_f = map(hom_resolution,i) #f[i]: M[i] -> M[i+1]
    hom_M_i_N = domain(hom_f)
    for v in gens(hom_M_i_N)
      @test element_to_homomorphism(hom_f(v)) == f*element_to_homomorphism(v)
    end
  end

  hom_hom_resolution = hom(hom_resolution,N)
  @test chain_range(hom_hom_resolution) == chain_range(free_res)

  hom_resolution = hom_without_reversing_direction(free_res,N)
  @test last(chain_range(hom_resolution)) == -first(chain_range(free_res))
  @test first(chain_range(hom_resolution)) == -last(chain_range(free_res))
  for i in Hecke.map_range(hom_resolution)
    f = map(free_res,-i+1)
    hom_f = map(hom_resolution,i)
    hom_M_i_N = domain(hom_f)
    for v in gens(hom_M_i_N)
      @test element_to_homomorphism(hom_f(v)) == f*element_to_homomorphism(v)
    end
  end
  hom_hom_resolution = hom_without_reversing_direction(hom_resolution,N)
  @test chain_range(hom_hom_resolution) == chain_range(free_res)

  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  F = free_module(R, 2)
  C, isom = present_as_cokernel(F, :cache_morphism)
  @test ambient_free_module(C) == F
  @test relations(C) == [zero(F)]
  @test domain(isom) == F
  @test codomain(isom) == C

  R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);
  A, p = quo(R, ideal(R, x^5))
  M1 = identity_matrix(A, 2)
  M2 = A[-x-y 2*x^2+x; z^4 0; 0 z^4; 8*x^3*y - 4*x^3 - 4*x^2*y + 2*x^2 + 2*x*y - x - y x; x^4 0]
  M = SubquoModule(M1, M2)
  fr = free_resolution(M, length = 9)
  @test all(iszero, homology(fr)[2:end])

  R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);
  Z = R(0)
  O = R(1)
  B = [Z Z Z O; w*y w*z-x*y x*z-y^2 Z];
  A = transpose(matrix(B));
  M = graded_cokernel(A)
  FM1 = free_resolution(M, length = 2)
  FM1[4]
  @test all(iszero, homology(FM1))
  FM2 = free_resolution(M, length = 2, algorithm = :mres)
  FM2[4]
  @test all(iszero, homology(FM2))
  FM3 = free_resolution(M, length = 2, algorithm = :nres)
  FM3[4]
  @test all(iszero, homology(FM3))

  R, (v, w, x, y, z) = polynomial_ring(QQ, [:v, :w, :x, :y, :z]);
  U = complement_of_point_ideal(R, [0, 0, 0, 0, 0]);
  RL, _ = localization(R, U);
  I = ideal(RL, [v, w,x, y, z])
  MI = ideal_as_module(I)
  FMI = free_resolution_via_kernels(MI)
  @test is_complete(FMI)
  @test length(FMI) == 4

  R, (v, w, x, y, z) = graded_polynomial_ring(QQ, [:v, :w, :x, :y, :z]);
  I = ideal(R, [v, w, x, y, z])
  A, _ = quo(R, I)
  M = quotient_ring_as_module(A)
  FM = free_resolution(M, length = 1, algorithm = :mres)
  @test length(FM) == 1
  FM[3]
  @test length(FM) == 3
  @test rank(FM[3]) == 10
end

@testset "Homology of ComplexOfMorphisms{FPModule{FqFieldElem}}" begin
    F = GF(2)
    function _rep_mat(n)
        I = [i for i in 1:n-1 for d in (0,1)]
        J = [(i+d-1) % n + 1 for i in 1:n-1 for d in (0,1)]
        matrix(F, sparse(I, J, trues(2*(n-1)), n-1, n))
    end
    for n in 2:50
        d = hom(free_module(F,n-1), free_module(F,n), _rep_mat(n))
        C = chain_complex([d])
        H = homology(C)
        ker_d = kernel(d)[1]
        @test dim(ker_d) == 0
        @test H[1] == ker_d
        coker_d, coker_proj = cokernel(d)
        @test dim(coker_d) == 1
        @test H[2][1] == coker_d
        @test all(g -> iszero(coker_proj(d(g))), gens(domain(d)))
    end
end

@testset "Prune With Map" begin
  Oscar.set_seed!(235)
  # ungraded
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  M = SubquoModule(identity_matrix(R, 3), R[1 x x])
  N, phi = prune_with_map(M)
  @test rank(ambient_free_module(N)) == 2
  @test (phi).(gens(N)) == gens(M)[2:3]

  M = SubquoModule(identity_matrix(R, 3), R[x 1 x])
  N, phi = prune_with_map(M)
  @test rank(ambient_free_module(N)) == 2
  @test (phi).(gens(N)) == [gens(M)[1], gens(M)[3]]

  # graded
  R, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
  F = graded_free_module(R, [2, 1])
  M = SubquoModule(F, identity_matrix(R, 2), R[1 x])
  N, phi = prune_with_map(M)
  @test rank(ambient_free_module(N)) == 1
  @test (phi).(gens(N)) == [gens(M)[2]]
  @test degrees_of_generators(N) == [degrees_of_generators(M)[2]]
end

@testset "Ext, Tor" begin
  Oscar.set_seed!(235)
  # These tests are only meant to check that the ext and tor function don't throw any error
  # These tests don't check the correctness of ext and tor

  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  A = R[x; y]
  B = R[x^2; x*y; y^2; z^4]
  M = SubquoModule(A, B)
  F = free_module(R, 1)
  Q, _ = quo(F, [x*F[1]])
  G = free_module(R, 2)
  M_coker = present_as_cokernel(M)

  T0 = tor(Q, M, 0)
  T1 = tor(Q, M, 1)
  T2 =  tor(Q, M, 2)
  @test_broken is_canonically_isomorphic(T0, M) # Should probably no longer be tested for after #4809.
  @test is_canonically_isomorphic(present_as_cokernel(T1), M_coker)
  @test iszero(T2)
  T0 = tor(M, Q, 0)
  T1 = tor(M, Q, 1)
  T2 = tor(M, Q, 2)
  @test is_canonically_isomorphic(present_as_cokernel(T0), M_coker)
  @test is_canonically_isomorphic(Oscar._old_simplify(present_as_cokernel(T1))[1], M_coker)
  @test iszero(T2)

  E0 = ext(Q, M, 0)
  E1 = ext(Q, M, 1)
  E2 = ext(Q, M, 2)
  @test is_canonically_isomorphic(present_as_cokernel(E0), M_coker)
  # test was dependend on internal design which has changed
  #@test is_canonically_isomorphic(E1, M_coker)
  @test iszero(E2)
  E0 = ext(M, Q, 0)
  E1 = ext(M, Q, 1)
  E2 = ext(M, Q, 2)
  E3 = ext(M, Q, 3)
  E4 = ext(M, Q, 4)
  @test iszero(E0)
  @test iszero(E1)
  @test is_canonically_isomorphic(present_as_cokernel(Oscar._old_simplify(E2)[1]), M_coker)
  @test is_canonically_isomorphic(E3, M_coker)
  @test iszero(E4)

  # Test that tor, ext don't crash outside of "sensible" arguments
  T3 = tor(Q, M, 20)
  E5 = ext(Q, M, 20)
end

@testset "Gröbner bases" begin
  Oscar.set_seed!(235)
  R, (x,y) = polynomial_ring(QQ, [:x, :y])
  F = FreeMod(R, 1)

  J = SubquoModule(F, [x*F[1], (x^2)*F[1], (x+y)*F[1]])
  @test leading_module(J) == SubquoModule(F, [x*F[1], y*F[1]])

  J = SubquoModule(F, [(x*y^2+x*y)*F[1], (x^2*y+x^2-y)*F[1]])
  @test leading_module(J) == SubquoModule(F, [x^2*F[1], y^2*F[1], x*y*F[1]]) # Example 1.5.7 in Singular book

  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  F = FreeMod(R, 2)
  lp = lex(gens(base_ring(F)))*lex(gens(F))

  M = SubquoModule(F, [(x^2*y^2*F[1]+y*z*F[2]), x*z*F[1]+z^2*F[2]])
  @test leading_module(M,lp) == SubquoModule(F, [x*z*F[1], x*y^2*z^2*F[2], x^2*y^2*F[1]])

  R, x = polynomial_ring(QQ, :x => 1:4)
  F = FreeMod(R, 1)
  lp = lex(gens(base_ring(F)))*lex(gens(F))

  J = SubquoModule(F, [(x[1]+x[2]+R(1))*F[1], (x[1]+x[2]+2*x[3]+2*x[4]+1)*F[1],(x[1]+x[2]+x[3]+x[4]+1)*F[1]])
  @test Oscar.oscar_generators(reduced_groebner_basis(J, lp)) == Oscar.oscar_generators(Oscar.ModuleGens([(x[3]+x[4])*F[1], (x[1]+x[2]+1)*F[1]], F))
  @test haskey(J.groebner_basis, lp)

  R, (x,y) = polynomial_ring(QQ, [:x, :y])
  F = FreeMod(R, 1)
  lp = lex(gens(base_ring(F)))*lex(gens(F))
  I = SubquoModule(F, [(x-1)*F[1], (y^2-1)*F[1]])
  f = (x*y^2+y)*F[1]
  @test Oscar.reduce(f, I) == (y+1)*F[1]

  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  F = FreeMod(R, 2)

  A = R[x+1 y*z+x^2; (y+2*z) z^3]
  B = R[2*z*(x*z+y^2) (x*z)^5]
  M = SubquoModule(F, A, B)
  gb = groebner_basis(M)
  P = sum(Oscar.SubModuleOfFreeModule(F, gb), Oscar.SubModuleOfFreeModule(F, gb.quo_GB))
  Q = Oscar.SubModuleOfFreeModule(F, groebner_basis(M.sum))
  @test P == Q
  v = x*((x+1)*F[1] + (y*z+x^2)*F[2]) + (y-z)*((y+2*z)*F[1] + z^3*F[2]) + 2*z*(x*z+y^2)*F[1] + (x*z)^5*F[2]
  @test represents_element(v,M)
  R, (x,y,z) = polynomial_ring(AcbField(64), [:x, :y, :z])
  F = FreeMod(R, 2)

  A = R[x+1 y*z+x^2; (y+2*z) z^3]
  B = R[2*z*(x*z+y^2) (x*z)^5]
  M = SubquoModule(F, A, B)
  @test_throws ArgumentError groebner_basis(M)
end

@testset "Singular Ordering Test" begin
  R, x = polynomial_ring(QQ, :x => 1:4)
  F = FreeMod(R, 1)
  lp = lex(gens(base_ring(F))) * lex(gens(F))
  J = SubquoModule(F, [
        (x[1] + x[2] + R(1)) * F[1],
        (x[1] + x[2] + 2*x[3] + 2*x[4] + 1) * F[1],
        (x[1] + x[2] + x[3] + x[4] + 1) * F[1]
  ])
  mg = reduced_groebner_basis(J, lp)
  odr = Oscar.singular_ordering(mg)
  @test Singular.ordering_as_symbol(odr) == :lex
end

@testset "Test kernel" begin
  Oscar.set_seed!(235)

  # over Integers
  R, (x,y,z) = polynomial_ring(ZZ, [:x, :y, :z])
  matrices = [R[x^2+x y^2+y; x^2+y y^2; x y], R[5*x^5+x*y^2 4*x*y+y^2+R(1); 4*x^2*y 2*x^2+3*y^2-R(5)]]
  kernels = [R[x^2-x*y+y -x^2+x*y -x^2+x*y-y^2-y], R[R(0) R(0)]]
  for (A,Ker) in zip(matrices, kernels)
    F1 = FreeMod(R, nrows(A))
    F2 = FreeMod(R, ncols(A))
    K,emb = kernel(FreeModuleHom(F1,F2,A))
    @test K == image(emb)[1]
    @test image(emb)[1] == image(map(F1,Ker))[1]
  end
  for k=1:3
    A = matrix([randpoly(R,0:2,2,2) for i=1:3,j=1:2])
    A = map(A)
    K,emb = kernel(A)
    @test iszero(emb*A)
  end

  # over Rationals
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  matrices = [R[x^2+x y^2+y; x^2+y y^2; x y], R[5*x^5+x*y^2 4*x*y+y^2+R(1); 4*x^2*y 2*x^2+3*y^2-R(5)],
        R[8*x^2*y^2*z^2+13*x*y*z^2  12*x^2+7*y^2*z;
        13*x*y^2+12*y*z^2  4*x^2*y^2*z+8*x*y*z;
        9*x*y^2+4*z  12*x^2*y*z^2+9*x*y^2*z]]
  kernels = [R[x^2-x*y+y -x^2+x*y -x^2+x*y-y^2-y], R[R(0) R(0)],
         R[-36*x^3*y^4*z+156*x^3*y^3*z^2+144*x^2*y^2*z^4+117*x^2*y^4*z+108*x*y^3*z^3-72*x^2*y^3*z-16*x^2*y^2*z^2-32*x*y*z^2 -96*x^4*y^3*z^4-72*x^3*y^4*z^3-156*x^3*y^2*z^4-117*x^2*y^3*z^3+63*x*y^4*z+108*x^3*y^2+28*y^2*z^2+48*x^2*z 32*x^4*y^4*z^3+116*x^3*y^3*z^3+104*x^2*y^2*z^3-91*x*y^4*z-84*y^3*z^3-156*x^3*y^2-144*x^2*y*z^2]]
  for (A,Ker) in zip(matrices, kernels)
    F1 = FreeMod(R, nrows(A))
    F2 = FreeMod(R, ncols(A))
    K,emb = kernel(FreeModuleHom(F1,F2,A))
    @test K == image(emb)[1]
    @test image(emb)[1] == image(map(F1,Ker))[1]
  end
  for k=1:3
    A = matrix([randpoly(R,0:2,2,2) for i=1:3,j=1:2])
    A = map(A)
    K,emb = kernel(A)
    @test image(emb)[1] == K
    @test iszero(emb*A)
  end
end

@testset "iszero(SubquoModule)" begin
  Oscar.set_seed!(235)
  R, (x,y) = polynomial_ring(QQ, [:x, :y])
  A = R[x^2+2*x*y y^2*x-2*x^2*y;-y x*y]
  B = R[x^2 y^2*x;-y x*y]
  @test iszero(SubquoModule(A,B))
  for k=1:3
    A = matrix([randpoly(R,0:2,2,2) for i=1:3,j=1:2])
    B = matrix([randpoly(R,0:2,2,2) for i=1:2,j=1:2])
    @test iszero(SubquoModule(A,A))
    @test !iszero(SubquoModule(A,B)) # could go wrong
  end
end

@testset "simplify subquotient" begin
  Oscar.set_seed!(235)
  R, (x,y) = polynomial_ring(QQ, [:x, :y])
  A1 = R[x*y R(0)]
  B1 = R[R(0) R(1)]
  M1 = SubquoModule(A1,B1)
  M2,i2,p2 = Oscar._old_simplify(M1)
  for k=1:5
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:1])), M1)
    @test elem == i2(p2(elem))
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
    @test elem == p2(i2(elem))
  end

  @test is_welldefined(i2)
  @test is_welldefined(p2)
  @test is_bijective(i2)
  @test is_bijective(p2)

  M2,i2,p2 = simplify_with_same_ambient_free_module(M1)
  @test ambient_free_module(M2) === ambient_free_module(M1)
  @test is_welldefined(i2)
  @test is_welldefined(p2)
  @test is_bijective(i2)
  @test is_bijective(p2)
  @test i2*p2 == id_hom(M2)
  @test p2*i2 == id_hom(M1)

  A1 = matrix([randpoly(R,0:2,2,1) for i=1:3,j=1:2])
  B1 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:2])
  M1 = SubquoModule(A1,B1)
  M2,i2,p2 = Oscar._old_simplify(M1)
  for k=1:5
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:3])), M1)
    @test elem == i2(p2(elem))
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:ngens(M2)])), M2)
    @test elem == p2(i2(elem))
  end

  @test is_welldefined(i2)
  @test is_welldefined(p2)
  @test is_bijective(i2)
  @test is_bijective(p2)

  A1 = matrix([randpoly(R,0:2,2,1) for i=1:3,j=1:3])
  B1 = matrix([randpoly(R,0:2,2,1) for i=1:2,j=1:3])
  M1 = SubquoModule(A1,B1)
  M2,i2,p2 = Oscar._old_simplify(M1)

  for k=1:5
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:3])), M1)
    @test elem == i2(p2(elem))
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
    @test elem == p2(i2(elem))
  end

  @test is_welldefined(i2)
  @test is_welldefined(p2)
  @test is_bijective(i2)
  @test is_bijective(p2)

  M1 = SubquoModule(B1,A1)
  M2,i2,p2 = Oscar._old_simplify(M1)
  for k=1:5
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:2])), M1)
    @test elem == i2(p2(elem))
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M2)])), M2)
    @test elem == p2(i2(elem))
  end

  @test is_welldefined(i2)
  @test is_welldefined(p2)
  @test is_bijective(i2)
  @test is_bijective(p2)
end

@testset "quotient modules" begin
  Oscar.set_seed!(235)
  R, (x,y) = polynomial_ring(QQ, [:x, :y])

  F3 = FreeMod(R,3)
  M1 = SubquoModule(F3,R[x^2*y^3-x*y y^3 x^2*y; 2*x^2 3*y^2*x 4],R[x^4*y^5 x*y y^4])
  N1 = SubquoModule(F3,R[x^4*y^5-4*x^2 x*y-6*y^2*x y^4-8],R[x^4*y^5 x*y y^4])
  Q1,p1 = quo(M1,N1, cache_morphism=true)

  @test Q1 == SubquoModule(F3,R[x^2*y^3-x*y y^3 x^2*y],R[x^4*y^5 x*y y^4; x^4*y^5-4*x^2  -6*x*y^2+x*y  y^4-8])
  @test p1 == find_morphism(M1, Q1)
  for k=1:5
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:2])), M1)
    @test p1(elem) == transport(Q1, elem)
  end

  F2 = FreeMod(R,2)
  M2 = SubquoModule(F2,R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],R[x^3-y^2 y^4-x-y])
  elems = [SubquoModuleElem(sparse_row(R[x*y -x*y^2 x*y]), M2), SubquoModuleElem(sparse_row(R[x R(0) R(-1)]), M2)]
  Q2,p2 = quo(M2,elems, cache_morphism=true)

  @test Q2 == SubquoModule(F2,R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],
            R[x^3-y^2 y^4-x-y; x^2*y^3+x^2*y^2+x^3*y^3+x*y^3-x^5*y^2 x^4*y+2*x*y^2-x*y^5+x^2*y^2; x^2*y-y^2 x^4+x*y])
  for k=1:5
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:3])), M2)
    @test p2(elem) == transport(Q2, elem)
  end

  M3 = SubquoModule(F3,R[x^2*y+13*x*y+2x-1 x^4 2*x*y; y^4 3*x -1],R[y^2 x^3 y^2])
  #N3 = SubquoModule(F3,R[x^2*y+13*x*y+2x-1-x*y^2 0 x^4-x*y^2; y^4-x*y^2 3*x-x^4 -1-x*y^2],R[2*y^2 2*x^3 2*y^2])
  N3 = SubquoModule(F3,R[x^2*y+13*x*y+2x-1-x*y^2 0 2*x*y-x*y^2; y^4-x*y^2 3*x-x^4 -1-x*y^2],R[2*y^2 2*x^3 2*y^2])
  Q3,p3 = quo(M3,N3, cache_morphism=true)

  @test iszero(quo_object(M3,M3))
  @test iszero(Q3)
  for k=1:5
    elem = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1,i=1:1])), M3)
    @test p3(elem) == transport(Q3, elem)
    @test iszero(p3(elem))
  end
end

@testset "submodules" begin
  Oscar.set_seed!(235)
  R, (x,y) = polynomial_ring(QQ, [:x, :y])

  F2 = FreeMod(R,2)
  M1 = SubquoModule(F2,R[x^2*y+x*y x*y^2-x; x+x*y^2 y^3],R[x^2 y^3-x])
  S1,i1 = sub(M1, [M1(sparse_row(R[1 1])),M1(sparse_row(R[y -x]))], cache_morphism=true)

  @test S1 == SubquoModule(F2,R[x*y^2+x^3-x^2 x*y^3-x*y-x^2; x^2*y+x*y^2+x*y-x^2+x x*y^2],R[x^2 y^3-x])
  @test i1 == find_morphism(S1, M1)
  for k=1:5
      elem = S1(sparse_row(matrix([randpoly(R) for _=1:1,i=1:2])))
      @test i1(elem) == transport(M1, elem)
  end

  M2 = SubquoModule(F2,R[x*y^2+x*y x^3+2*y; x^4 y^3; x^2*y^2+y^2 x*y],R[x^3-y^2 y^4-x-y])
  S2,i2 = sub(M2,[M2(sparse_row(R[x*y -x*y^2 x*y])),M2(sparse_row(R[x 0 -1]))], cache_morphism=true)

  @test S2 == SubquoModule(F2,R[x^2*y^3+x^2*y^2+x^3*y^3+x*y^3-x^5*y^2 x^4*y+2*x*y^2-x*y^5+x^2*y^2; x^2*y-y^2 x^4+x*y],R[x^3-y^2 y^4-x-y])
  @test i2 == find_morphisms(S2, M2)[1]
  for k=1:5
      elem = S2(sparse_row(matrix([randpoly(R) for _=1:1,i=1:2])))
      @test i2(elem) == transport(M2, elem)
  end

  M3 = SubquoModule(F2,R[x*y^2 x^3+2*y; x^4 y^3; x*y+y^2 x*y],R[x^3-y^2 y^4-x-y])
  elems = [M3(sparse_row(R[0 6 0])),M3(sparse_row(R[9 0 -x])),M3(sparse_row(R[0 0 -42]))]
  S3,i3 = sub(M3,elems; cache_morphism=true)

  @test S3 == M3
  for k=1:5
      elem = S3(sparse_row(matrix([randpoly(R) for _=1:1, i=1:3])))
      @test i3(elem) == transport(M3, elem)
  end
end

@testset "Hom module" begin
  Oscar.set_seed!(235)
  R, (x0,x1,x2,x3,x4,x5) = polynomial_ring(QQ, [:x0, :x1, :x2, :x3, :x4, :x5])
  f1= transpose(R[-x2*x3 -x4*x5 0; x0*x1 0 -x4*x5; 0 x0*x1 -x2*x3])
  g1 = transpose(R[x0*x1 x2*x3 x4*x5])
  M = cokernel(f1)
  N = cokernel(g1)
  SQ = hom(M,N)[1]

  f2 = R[0 0]
  M = cokernel(f1)
  N = cokernel(f2)
  SQ = hom(M,N)[1]


  function free_module_SQ(ring,n)
    F = FreeMod(ring,n)
    return SubquoModule(F,gens(F))
  end

  function free_module_SQ(F::FreeMod)
    return SubquoModule(F,gens(F))
  end

  # test if Hom(R^n,R^m) gives R^(m*n): (there is a dedicated function for free modules but this is a simple test for the function for subquos)
  for n=1:5
    for m=1:5
      M=free_module_SQ(R,n)
      N=free_module_SQ(R,m)
      SQ = hom(M,N)[1]
      #SQ = simplify(hom(M,N)[1])[1]
      @test SQ == free_module_SQ(ambient_free_module(SQ))
    end
  end
  R, (x,y) = polynomial_ring(QQ, [:x, :y])
  A1 = R[x^2+1 x*y; x^2+y^3 x*y]
  B1 = R[x+x^4+y^2+1 x^5; y^4-3 x*y^2-1]
  M1 = SubquoModule(A1, B1)
  A2 = R[x;]
  B2 = R[y^3;]
  M2 = SubquoModule(A2,B2)
  SQ = hom(M1,M2)[1]
  for v in gens(SQ)
    @test v == homomorphism_to_element(SQ, element_to_homomorphism(v))
  end

  R, (x,y) = polynomial_ring(QQ, [:x, :y])
  A1 = R[x^2+1 x*y; x^2+y^3 x*y]
  B1 = R[x+x^4+y^2+1 x^5; y^4-3 x*y^2-1]
  M1 = SubquoModule(A1, B1)
  A2 = R[x;]
  B2 = R[y^3;]
  M2 = SubquoModule(A2,B2)
  SQ = hom(M1,M2, algorithm=:matrices)[1]
  for v in gens(SQ)
    @test v == homomorphism_to_element(SQ, element_to_homomorphism(v))
  end

  End_M = hom(M1,M1)[1]
  R_as_module = FreeMod(R,1)
  phi = multiplication_induced_morphism(R_as_module, End_M)
  @test element_to_homomorphism(phi(R_as_module[1])) == id_hom(M1)
  @test image(element_to_homomorphism(phi((x+y)*R_as_module[1])))[1] == (ideal(R,x+y)*M1)[1]

  # test if hom(zero-module, ...) is zero
  Z = FreeMod(R,0)
  @test iszero(hom(Z,Z)[1])
  for k=1:10
    A = matrix([randpoly(R,0:2,2,1) for i=1:3,j=1:2])
    B = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:2])
    N = SubquoModule(A,B)

    @test iszero(hom(N,Z)[1])
    @test iszero(hom(Z,N)[1])
  end

  # test welldefinedness of randomly generated homomorphisms (using hom() and element_to_homomorphism())
  for k=1:10
    A1 = matrix([randpoly(R,0:2,2,1) for i=1:3,j=1:2])
    A2 = matrix([randpoly(R,0:2,2,1) for i=1:2,j=1:2])
    B1 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:2])
    B2 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:2])
    N = SubquoModule(A1,B1)
    M = SubquoModule(A2,B2)
    HomNM = k <= 5 ? hom(N,M)[1] : hom(N,M, algorithm=:matrices)[1]
    for l=1:10
      v = sparse_row(matrix([randpoly(R,0:2,2,1) for _=1:1, j=1:AbstractAlgebra.ngens(HomNM)]))
      H = HomNM(v)
      H = element_to_homomorphism(H)
      @test is_welldefined(H)
    end
  end
end

@testset "tensoring morphisms" begin
  Oscar.set_seed!(235)
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])

  F2 = FreeMod(R,2)
  F3 = FreeMod(R,3)
  F4 = FreeMod(R,4)

  for _=1:3 # At some later point we run into non-finishing computations with this seed. 
    A1 = matrix([randpoly(R,0:2,4,3) for i=1:3,j=1:2])
    B1 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:2])


    A2 = matrix([randpoly(R,0:2,2,1) for i=1:3,j=1:3])
    B2 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:3])

    M1 = SubquoModule(F2,A1,B1)
    M2 = SubquoModule(F3,A2,B2)

    M,pure_M = tensor_product(M1,M2, task=:map)
    phi = hom_tensor(M,M,[id_hom(M1),id_hom(M2)])

    for _=1:3
      v = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M)])), M)
      @test phi(v) == v
    end

    A3 = matrix([randpoly(R,0:2,2,1) for i=1:2,j=1:2])
    #B3 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:2])
    M3 = SubquoModule(Oscar.SubModuleOfFreeModule(F2,A3))

    N,pure_N = tensor_product(M3,F4, task=:map)

    M3_to_M1 = SubQuoHom(M3,M1, matrix([randpoly(R,0:2,2,2) for i=1:ngens(M3), j=1:ngens(M1)]))
    is_welldefined(M3_to_M1) || continue
    F4_to_M2 = FreeModuleHom(F4,M2, matrix([randpoly(R,0:2,2,2) for i=1:ngens(F4), j=1:ngens(M2)]))

    phi = hom_tensor(N,M,[M3_to_M1,F4_to_M2])
    u1 = SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(M3)])), M3)
    u2 = F4(sparse_row(matrix([randpoly(R) for _=1:1, i=1:ngens(F4)])))
    @test phi(pure_N((u1,u2))) == pure_M((M3_to_M1(u1),F4_to_M2(u2)))
  end
end

@testset "direct product" begin
  Oscar.set_seed!(235)
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])

  F2 = FreeMod(R,2)
  F3 = FreeMod(R,3)

  A1 = matrix([randpoly(R,0:2,2,2) for i=1:3,j=1:2])
  B1 = matrix([randpoly(R,0:2,2,2) for i=1:1,j=1:2])
  M1 = SubquoModule(F2,A1,B1)

  A2 = matrix([randpoly(R,0:2,2,1) for i=1:2,j=1:3])
  B2 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:3])
  M2 = SubquoModule(F3,A2,B2)

  sum_M, emb = direct_sum(M1,M2)

  @test domain(emb[1]) === M1
  @test domain(emb[2]) === M2
  @test codomain(emb[1]) === sum_M
  @test codomain(emb[2]) === sum_M

  sum_M, pr = direct_sum(M1,M2, task=:prod)
  @test codomain(pr[1]) === M1
  @test codomain(pr[2]) === M2
  @test domain(pr[1]) === sum_M
  @test domain(pr[2]) === sum_M

  prod_M, emb, pr = direct_sum(M1,M2,task=:both)
  @test length(pr) == length(emb) == 2
  @test ngens(prod_M) == ngens(M1) + ngens(M2)

  for g in gens(prod_M)
    @test g == sum([emb[i](pr[i](g)) for i=1:length(pr)])
  end
  for g in gens(M1)
    @test g == pr[1](emb[1](g))
  end
  for g in gens(M2)
    @test g == pr[2](emb[2](g))
  end

  A1 = matrix([randpoly(R,0:2,2,2) for i=1:3,j=1:2])
  B1 = matrix([randpoly(R,0:2,2,2) for i=1:1,j=1:2])
  N1 = SubquoModule(F2,A1,B1)

  A2 = matrix([randpoly(R,0:2,2,1) for i=1:2,j=1:3])
  B2 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:3])
  N2 = SubquoModule(F3,A2,B2)

  prod_N = direct_product(N1,N2,task=:none)
  @test ngens(prod_M) == ngens(M1) + ngens(M2)

  for g in gens(prod_N)
    @test g == sum([canonical_injection(prod_N,i)(canonical_projection(prod_N,i)(g)) for i=1:2])
  end
  for g in gens(N1)
    @test g == canonical_projection(prod_N,1)(canonical_injection(prod_N,1)(g))
  end
  for g in gens(N2)
    @test g == canonical_projection(prod_N,2)(canonical_injection(prod_N,2)(g))
  end

  # testing hom_product

  M1_to_N1 = SubQuoHom(M1,N1,zero_matrix(R,3,3))
  H12 = hom(M1,N2)[1]
  H21 = hom(M2,N1)[1]
  M1_to_N2 = iszero(H12) ? SubQuoHom(M1,N2,zero_matrix(R,3,2)) : element_to_homomorphism(H12[1])
  M2_to_N1 = iszero(H21) ? SubQuoHom(M2,N1,zero_matrix(R,2,3)) : element_to_homomorphism(H21[1])
  M2_to_N2 = SubQuoHom(M2,N2,R[0 0; 0 0])
  @assert is_welldefined(M1_to_N1)
  @assert is_welldefined(M1_to_N2)
  @assert is_welldefined(M2_to_N1)
  @assert is_welldefined(M2_to_N2)

  phi = hom_product(prod_M,prod_N,[M1_to_N1 M1_to_N2; M2_to_N1 M2_to_N2])
  for g in gens(M1)
    @test M1_to_N1(g) == canonical_projection(prod_N,1)(phi(emb[1](g)))
    @test M1_to_N2(g) == canonical_projection(prod_N,2)(phi(emb[1](g)))
  end
  for g in gens(M2)
    @test M2_to_N1(g) == canonical_projection(prod_N,1)(phi(emb[2](g)))
    @test M2_to_N2(g) == canonical_projection(prod_N,2)(phi(emb[2](g)))
  end

  # testing mixed typed modules

  prod_FN,prod,emb = direct_product(F2,N2,task=:both)
  @test ngens(prod_FN) == ngens(F2) + ngens(N2)
  for g in gens(prod_FN)
    @test g == sum([emb[i](prod[i](g)) for i=1:2])
  end
  for g in gens(F2)
    @test g == prod[1](emb[1](g))
  end
  for g in gens(N2)
    @test g == prod[2](emb[2](g))
  end
end

@testset "Coordinates (lift)" begin
  Oscar.set_seed!(235)
  Z3, a = finite_field(3,1,"a")
  R, (x,y) = polynomial_ring(Z3, [:x, :y])
  coeffs = [Z3(i) for i=0:1]

  A = R[x*y x^2+y^2; y^2 x*y;x^2+1 1]
  B = R[2*x^2 x+y; x^2+y^2-1 x^2+2*x*y]
  F = FreeMod(R,2)
  M = SubquoModule(F,A,B)

  monomials = [x,y]
  coeff_generator = ([c1,c2] for c1 in coeffs for c2 in coeffs)
  for coefficients in ([c1,c2,c3] for c1 in coeff_generator for c2 in coeff_generator for c3 in coeff_generator)
    v = sparse_row(R, [(i,sum(coefficients[i][j]*monomials[j] for j=1:2)) for i=1:3])
    v_as_FreeModElem = sum([v[i]*repres(M[i]) for i=1:ngens(M)])

    elem1 = SubquoModuleElem(v_as_FreeModElem,M)
    elem2 = SubquoModuleElem(v,M)

    @test elem1 == elem2
  end
end

@testset "module homomorphisms" begin
  Oscar.set_seed!(235)
  R, (x,y) = polynomial_ring(QQ, [:x, :y])

  F3 = FreeMod(R,3)
  F4 = FreeMod(R,4)
  phi = hom(F3,F4, [F4[1],F4[3]+x*F4[4],(x+y)*F4[4]] )
  z = hom(F3,F4, [zero(F4) for _ in gens(F3)])
  for v in gens(F3)
      @test (z-phi)(v) == (-phi)(v)
  end
  @test iszero(preimage(phi,zero(F4)))
  @test phi(preimage(phi,y*(F4[3]+x*F4[4]+(x+y)*F4[3]))) == y*(F4[3]+x*F4[4]+(x+y)*F4[3])

  A1 = R[9*y 11*x; 0 8; 14*x*y^2 x]
  A2 = R[4*x*y 15; 4*y^2 6*x*y]
  B1 = R[2*x*y^2 6*x^2*y^2]
  B2 = R[15*x*y 3*y]


  #1) H: N --> M where N is a cokernel, H should be an isomorphism
  F2 = FreeMod(R,2)
  M = SubquoModule(F2,A1,B1)
  N, H = present_as_cokernel(M, :cache_morphism)
  Hinv = H.inverse_isomorphism
  @test is_welldefined(H)

  ## testing the homomorphism theorem: #################################
  KerH,iKerH = kernel(H)
  ImH,iImH = image(H)

  NmodKerH, pNmodKerH = quo(N,KerH, cache_morphism=true)
  Hbar = SubQuoHom(NmodKerH,M,matrix(H))
  Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

  @test is_welldefined(Hbar)
  @test is_bijective(Hbar)

  Hbar_inv = inv(Hbar)

  # test, if caching of inverse maps works:
  @test Hbar_inv === Hbar.inverse_isomorphism
  @test Hbar === Hbar_inv.inverse_isomorphism

  # test if Hbar and Hbar_inv are inverse to each other:
  @test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
  @test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
  #######################################################################

  # test, if H is bijective with inverse Hinv:
  @test is_bijective(H)
  @test all([inv(H)(H(g))==g for g in gens(N)])
  @test all([H(inv(H)(g))==g for g in gens(M)])
  @test inv(H) === Hinv
  @test ImH == M
  @test iszero(KerH)


  #2) H: N --> M = N/(submodule of N) canonical projection
  M,H = quo(N,[N(sparse_row(R[1 x^2-1 x*y^2])),N(sparse_row(R[y^3 y*x^2 x^3]))], cache_morphism=true)
  @test is_welldefined(H)

  ## test addition/subtraction of morphisms
  H_1 = H+H-H
  for v in gens(N)
      @test H_1(v) == H(v)
  end
  for v in gens(N)
    @test ((x+y+R(1))*H_1)(v) == H_1((x+y+R(1))*v)
  end

  ## testing the homomorphism theorem: #################################
  KerH,iKerH = kernel(H)
  ImH,iImH = image(H)

  NmodKerH, pNmodKerH = quo(N,KerH, cache_morphism=true)
  Hbar = SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
  Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

  @test is_welldefined(Hbar)
  @test is_bijective(Hbar)

  Hbar_inv = inv(Hbar)

  # test, if caching of inverse maps works:
  @test Hbar_inv === Hbar.inverse_isomorphism
  @test Hbar === Hbar_inv.inverse_isomorphism

  # test if Hbar and Hbar_inv are inverse to each other:
  @test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
  @test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
  #######################################################################

  # test, if H is surjective:
  @test ImH == M

  #3) H:N --> M neither injective nor surjective, also: tests 'restrict_domain()'
  MM = SubquoModule(F2,A1,B1)
  M, iM, iSQ = sum(MM, SubquoModule(F2,A2,B1))
  NN, p, i = direct_product(MM,SubquoModule(A2,B2), task = :both)
  i1,i2 = i[1],i[2]
  p1,p2 = p[1],p[2]
  nn = ngens(NN)
  u1 = R[3*y 14*y^2 6*x*y^2 x^2*y 3*x^2*y^2]
  u2 = R[5*x*y^2 10*y^2 4*x*y 7*x^2*y^2 7*x^2]
  u3 = R[13*x^2*y 4*x*y 2*x 7*x^2 9*x^2]
  N,iN = sub(NN,[NN(sparse_row(u1)), NN(sparse_row(u2)), NN(sparse_row(u3))], cache_morphism=true)

  H = restrict_domain(p1*iM,N)
  @test is_welldefined(H)

  ## testing the homomorphism theorem: #################################
  KerH,iKerH = kernel(H)
  ImH,iImH = image(H)

  NmodKerH, pNmodKerH = quo(N,KerH, cache_morphism=true)
  Hbar = SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
  Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

  @test is_welldefined(Hbar)
  @test is_bijective(Hbar)

  Hbar_inv = inv(Hbar)

  # test, if caching of inverse maps works:
  @test Hbar_inv === Hbar.inverse_isomorphism
  @test Hbar === Hbar_inv.inverse_isomorphism

  # test if Hbar and Hbar_inv are inverse to each other:
  @test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
  @test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
  #######################################################################



  #4) H: M --> N random map created via the hom() function
  N = SubquoModule(F2,A1,B1)
  M = SubquoModule(F2,A2,B2)
  HomNM = hom(N,M)[1]
  u1 = R[x^2*y^2 4*x^2*y^2] # 0 5*x*y^2] # original vector was longer but current output is too small.
  H = HomNM(sparse_row(u1))
  H = element_to_homomorphism(H)
  @test is_welldefined(H)

  ## testing the homomorphism theorem: #################################
  KerH,iKerH = kernel(H)
  ImH,iImH = image(H)

  NmodKerH, pNmodKerH = quo(N,KerH, cache_morphism=true)
  Hbar = SubQuoHom(NmodKerH,M,matrix(H)) # induced map N/KerH --> M
  Hbar = restrict_codomain(Hbar,ImH) # induced map N/KerH --> ImH

  @test is_welldefined(Hbar)
  @test is_bijective(Hbar)

  Hbar_inv = inv(Hbar)

  # test, if caching of inverse maps works:
  @test Hbar_inv === Hbar.inverse_isomorphism
  @test Hbar === Hbar_inv.inverse_isomorphism

  # test if Hbar and Hbar_inv are inverse to each other:
  @test all([Hbar_inv(Hbar(g))==g for g in gens(NmodKerH)])
  @test all([Hbar(Hbar_inv(g))==g for g in gens(ImH)])
end


@testset "lift of homomorphisms" begin
  K = GF(3)
  S, (x0, x1, x2, x3, x4) = graded_polynomial_ring(K, [:x0, :x1, :x2, :x3, :x4]);
  m = ideal(S, [x1^2+(-x1+x2+x3-x4)*x0, x1*x2+(x1-x3+x4)*x0, x1*x3+(-x1+x4+x0)*x0, x1*x4+(-x1+x3+x4-x0)*x0, x2^2+(x1-x2-x4-x0)*x0, x2*x3+(x1-x2+x3+x4-x0)*x0, x2*x4+(x1+x2-x3-x4-x0)*x0, x3^2+(x3+x4-x0)*x0,x3*x4+(-x3-x4+x0)*x0, x4^2+(x1+x3-x4-x0)*x0]);
  A, _ = quo(S, m);
  FA = free_resolution(A, algorithm = :mres, length = 2)
  phi1 = map(FA, 1)
  e1 = gen(codomain(phi1),1)
  phi2 = map(FA, 2)
  #phi3 = map(FA, 3)
  L = monomial_basis(A, 2);
  a = [[i == ((j-1) ÷ 5 + 1) ? S(L[((j-1) % 5) + 1]) : S(0) for i in 1:10] for j in 1:50]
  A1 = hom(domain(phi1), codomain(phi1), [p*e1 for p in a[1]])
  A1phi2 = phi2 * A1
  A2 = lift(A1phi2, phi1)
  @test matrix(A2)[1,1] == 2*x1 + 2*x2 + x3 + 2*x4
  #A2phi3 = phi3 * A2
  #A3 = lift2(A2phi3, phi2)
  #A3m = matrix(A3)
  #@test A3m[2,17]==S(2)
end


@testset "preimage" begin
  Oscar.set_seed!(235)
  R, (x,y) = polynomial_ring(QQ, [:x, :y])

  for _=1:10
    A1 = matrix([randpoly(R,0:2,2,1) for i=1:3,j=1:1])
    A2 = matrix([randpoly(R,0:2,2,1) for i=1:2,j=1:2])
    B1 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:1])
    B2 = matrix([randpoly(R,0:2,2,1) for i=1:1,j=1:2])

    N = SubquoModule(A1,B1)
    M = SubquoModule(A2,B2)
    HomNM = hom(N,M)[1]
    if iszero(HomNM)
      continue
    end
    H = HomNM(sparse_row(matrix([randpoly(R,0:2,2,1) for _=1:1,j=1:ngens(HomNM)])))
    H = element_to_homomorphism(H)

    u = [SubquoModuleElem(sparse_row(matrix([randpoly(R) for _=1:1, _=1:ngens(N)])), N) for _=1:3]
    image_of_u = sub_object(M, map(H, u))
    preimage_test_module = image_of_u + sub_object(M,[M[1]])
    _,emb = preimage(H,preimage_test_module,:with_morphism)
    @test issubset(sub_object(N,u), image(emb)[1])
  end
end

@testset "change of base rings" begin
  Oscar.set_seed!(235)
  R, (x,y) = QQ[:x, :y]
  U = Oscar.MPolyPowersOfElement(x)
  S = Oscar.MPolyLocRing(R, U)
  F = FreeMod(R, 2)
  FS, mapF = change_base_ring(S, F)
  @test 1//x*mapF(x*F[1]) == FS[1]

  shift = hom(R, R, [x-1, y-2])
  FSshift, mapFSshift = change_base_ring(shift, F)
  @test mapFSshift(x*F[1]) == (x-1)*FSshift[1]

  A = R[x y]
  B = R[x^2 x*y]
  M = SubquoModule(F, A, B)
  MS, mapM = change_base_ring(S, M)
  @test iszero(mapM(M[1]))

  f = MapFromFunc(R, S, S)
  MS, mapM = change_base_ring(f, M)
  @test iszero(mapM(M[1]))
end

@testset "duals" begin
  Oscar.set_seed!(235)
  R, (x,y,z) = QQ[:x, :y, :z]
  F1 = FreeMod(R, 1)
  F2 = FreeMod(R, 2)
  F2v, ev = Oscar.dual(F2, codomain=F1)
  @test ev(F2v[1])(F2[1]) == F1[1] # the first generator

  FF, psi = Oscar.double_dual(F2)
  @test is_injective(psi)
  @test is_surjective(psi)

  M, inc = sub(F2, [x*F2[1], y*F2[1]])
  F1 = FreeMod(R, 1)
  Mv, ev = dual(M, codomain=F1)
  @test ev(Mv[1])(M[1]) == x*F1[1]

  Mvv, psi = Oscar.double_dual(M, codomain=F1)
  @test matrix(psi) == R[x; y]

  ### Quotient rings

  A = R[x y; z x-1]
  Q, _ = quo(R, ideal(R, det(A)))

  F1 = FreeMod(Q, 1)
  F2 = FreeMod(Q, 2)
  F2v, ev = Oscar.dual(F2, codomain=F1)
  @test ev(F2v[1])(F2[1]) == F1[1] # the first generator

  FF, psi = Oscar.double_dual(F2)
  @test is_injective(psi)
  @test is_surjective(psi)

  M, pr = quo(F2, [sum(A[i, j]*F2[j] for j in 1:ngens(F2)) for i in 1:nrows(A)])
  Mv, ev = Oscar.dual(M, codomain=F1)
  Mvv, psi = Oscar.double_dual(M, codomain=F1)
  @test is_injective(psi)
  @test is_surjective(psi) # works correctly!
end

@testset "free resolution in case of no relations" begin
  Oscar.set_seed!(235)
  R, (x,y,z) = polynomial_ring(QQ, [:x, :y, :z])
  Z = abelian_group(0)
  F = free_module(R, 3)
  G = free_module(R, 2)
  V = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
  a = hom(F, G, V)
  Mk = kernel(a)
  Mk1 = Mk[1]
  fr = free_resolution(Mk1)
  @test is_zero(map(fr.C,1))
end

@testset "length of free resolution" begin
  Oscar.set_seed!(235)
  S, (x0, x1, x2, x3, x4) = graded_polynomial_ring(GF(3), [:x0, :x1, :x2, :x3, :x4]);
  m = ideal(S, [x1^2+(-x1+x2+x3-x4)*x0, x1*x2+(x1-x3+x4)*x0,
            x1*x3+(-x1+x4+x0)*x0,
            x1*x4+(-x1+x3+x4-x0)*x0, x2^2+(x1-x2-x4-x0)*x0,
            x2*x3+(x1-x2+x3+x4-x0)*x0, x2*x4+(x1+x2-x3-x4-x0)*x0,
            x3^2+(x3+x4-x0)*x0,x3*x4+(-x3-x4+x0)*x0,
            x4^2+(x1+x3-x4-x0)*x0]);
  A, _ = quo(S, m);
  FA = free_resolution(A)
  @test length(FA.C.maps) == 9
end

@testset "vector_space_dim and vector_space_basis" begin
  Oscar.set_seed!(235)
  R,(x,y,z,w) = QQ[:x, :y, :z, :w]
  U=complement_of_point_ideal(R,[0,0,0,0])
  RL, loc_map = localization(R,U)
  Floc = free_module(RL,2)
  v = gens(Floc)
  Mloc, _=quo(Floc,[x*v[1],y*v[1],z*v[1],w^2*(w+1)^3*v[1],v[2]])
  @test vector_space_dim(Mloc) == 2
  @test vector_space_dim(Mloc,1) == 1
  @test length(vector_space_basis(Mloc)) == 2
  @test length(vector_space_basis(Mloc,0)) == 1
end

@testset "canonical maps and garbage collection" begin
  Oscar.set_seed!(235)
  R, (x, y) = QQ[:x, :y]

  F = FreeMod(R, 1)

  # we need to wrap the creation of I in a scope of its own so
  # that gc works within the test suite.
  function dummy(F::FreeMod)
    I, inc = sub(F, [x*F[1]], cache_morphism=true)

    @test haskey(F.incoming, I)
    @test Oscar._recreate_morphism(I, F, F.incoming[I]) == inc

    @test haskey(I.outgoing, F)
    @test Oscar._recreate_morphism(I, F, I.outgoing[F]) == inc

    I = 5
    inc = "a"
  end
  dummy(F)

  GC.gc()
  GC.gc()

  @test length(keys(F.incoming)) == 0
  # The other way around it will not work, because I has a reference to its ambient_free_module f.

  function dummy2(F::FreeMod)
    I, inc_I = sub(F, [x*F[1]], cache_morphism=true)
    J, inc_J = sub(I, [x^2*I[1]], cache_morphism=true)

    @test haskey(J.outgoing, I)
    @test haskey(I.incoming, J)
    @test Oscar._recreate_morphism(J, I, J.outgoing[I]) == inc_J
    @test Oscar._recreate_morphism(J, I, I.incoming[J]) == inc_J

    I = 5
    inc_I = 6
  end
  dummy2(F)
  
  GC.gc()
  GC.gc()

  @test length(keys(F.incoming)) == 0
end


@testset "issue 3107" begin
  Oscar.set_seed!(235)
  X = veronese();
  I = defining_ideal(X);
  Pn = base_ring(I)
  FI = free_resolution(I)
  F = graded_free_module(Pn, 1)
  dualFIC = hom(FI.C, F)
end

@testset "issue 4203" begin
  R, (x,y) = polynomial_ring(GF(2), ["x","y"]);
  
  g = [x+1, y+1]
  s = syzygy_generators(g)
  F = parent(first(s))
  s2 = syzygy_generators(g; parent=F)
  @test s == s2
  s3 = syzygy_generators(g)
  @test s2 != s3
  
  A, _ = quo(R, ideal(R, [x^2+1, y^2+1]));
  
  g = A.([x+1, y+1])
  s = syzygy_generators(g)
  F = parent(first(s))
  s2 = syzygy_generators(g; parent=F)
  @test s == s2
  s3 = syzygy_generators(g)
  @test s2 != s3
  
  F = free_module(A,2);
  M,_ = sub(F, [F([A(x+1),A(y+1)])]);
  FF = free_module(A, ngens(M)); 
  phi = hom(FF, M, gens(M)); 
  K, = kernel(phi); 
  @test ambient_representatives_generators(K) == [(x*y + x + y + 1)*FF[1]]
  @test_throws ErrorException standard_basis(K)
  s1 = syzygy_generators(ambient_representatives_generators(M))
  G = parent(first(s1))
  s2 = syzygy_generators(ambient_representatives_generators(M); parent=G)
  @test s1 == s2
  s3 = syzygy_generators(ambient_representatives_generators(M))
  @test s3 != s2
  @test_throws ArgumentError syzygy_generators(ambient_representatives_generators(M); parent=F)
end

@testset "Issue #4809" begin
  R, (x,) = polynomial_ring(QQ, [:x])
  F = free_module(R, 1)
  M = SubquoModule(F, [x*F[1]], [x^2*F[1]])
  MoM = tensor_product(M, M)
  @test !is_zero(MoM)
  decomp = Oscar.tensor_generator_decompose_function(MoM)
  m = Oscar.tensor_pure_function(MoM)
  b0 = (M[1], M[1])
  a = m(b0...)
  b = decomp(a)
  @test b == b0
end

@testset "composition of morphisms" begin
  R, (x, y) = QQ[:x, :y]
  P, t = QQ[:t]
  kk, _ = extension_field(t^2 + 1)
  S, (u, v) = kk[:u, :v]
  R1 = free_module(R, 1)
  S1 = free_module(S, 1)
  id_R1 = hom(R1, R1, [R1[1]])
  id_S1 = hom(S1, S1, [S1[1]])
  f = hom(R1, S1, [S1[1]], hom(R, S, [u, v]))
  i = first(gens(kk))
  conj = hom(kk, kk, -i)
  conj_S1 = hom(S1, S1, [S1[1]], hom(S, S, conj, [u, v]))
  # recreate the same morphism, but with an anonymous function as ring_map
  conj_S1_alt = hom(S1, S1, [S1[1]], f->evaluate(map_coefficients(conj, f), [u, v]))
  @test_throws ErrorException conj_S1 == conj_S1_alt

  a = compose(compose(id_R1, f), compose(conj_S1, id_S1))
  b = compose(compose(id_R1, compose(f, conj_S1)), id_S1)
  c = compose(compose(compose(id_R1, f), conj_S1), id_S1)
  d = compose(id_R1, compose(f, compose(conj_S1, id_S1)))

  @test a == b
  a_alt = compose(compose(id_R1, f), compose(conj_S1_alt, id_S1))
  @test_throws ErrorException f == a # non-comparable ring_maps
  @test_throws ErrorException a_alt == a # same
  @test_throws MethodError !(conj_S1 == id_S1) # ring map vs. no ring map
  @test conj_S1 == conj_S1 # identical ring maps are OK
  @test all(a(x) == b(x) for x in gens(R1))
  @test all(a(x) == c(x) for x in gens(R1))
  @test all(a(x) == d(x) for x in gens(R1))

  @test all(i*conj_S1(x) == conj_S1(-i*x) for x in gens(S1))
  @test conj_S1 == compose(conj_S1, id_S1) # identical ring_map
  @test conj_S1 == compose(id_S1, conj_S1)
end

@testset "check exactness for mres" begin
  R, (w, x, y, z) = graded_polynomial_ring(QQ, [:w, :x, :y, :z]);
  Z = R(0)
  O = R(1)
  B = [Z Z Z O; w*y w*z-x*y x*z-y^2 Z];
  A = transpose(matrix(B));
  M = graded_cokernel(A)
  FM2 = free_resolution(M, algorithm = :mres)
  @test all(iszero, homology(FM2))
end

@testset "Modules over ZZ and QQ" begin

  @testset "Module Constructors over ZZ" begin
    F0 = free_module(FreeMod, ZZ,3)
    @test rank(F0) == 3
    @test base_ring(F0) == ZZ

    F = FreeMod(ZZ, 3)
    M, inc = sub(F, [2*F[1], 3*F[2]])
    N, proj = quo(F, [F[1] + 2*F[3]])

    @test rank(F) == 3
    @test is_welldefined(inc)
    @test is_welldefined(proj)
  end

  @testset "Submodule Membership over ZZ" begin
    R = ZZ
    F = FreeMod(R, 3)
    gens_submodule = [2*F[1], 3*F[2]]
    S, _ = sub(F, gens_submodule)
    x = 4*F[1]+3*F[2]
    @test in(x, S)
    coord = coordinates(x, S)
    @test coord == sparse_row(R, [1, 2], [2, 1])
  end

  @testset "Module Homomorphisms over ZZ" begin
    F1 = FreeMod(ZZ, 2)
    F2 = FreeMod(ZZ, 2)
    phi = hom(F1, F2, [F2[1] + 2*F2[2], 3*F2[1]+ 6*F2[2]])

    @test is_welldefined(phi)

    K, emb = kernel(phi)
    @test ngens(K) == 1
    @test is_welldefined(emb)

    I, im = image(phi)
    @test ngens(I) == 2
    @test is_welldefined(im)
  end

  @testset "Presentations over ZZ" begin
    F = FreeMod(ZZ, 2)
    M, inc = sub(F, [2*F[1] + 3*F[2]])
    pres_M, iso = present_as_cokernel(M, :both)

    @test is_welldefined(iso)
    @test is_bijective(iso)
    @test pres_M isa SubquoModule
    @test ngens(pres_M) == 1

    F = FreeMod(ZZ, 3)
    M, inc = sub(F, [5*F[1],3*F[2],F[3]])
    N, inc = sub(F, [10*F[1], 3*F[2]])
    W, _ = quo(M, N)
    pres_W, iso = present_as_cokernel(W, :both)
    @test is_welldefined(iso)
    @test is_bijective(iso)
    @test pres_W isa SubquoModule
    @test ngens(pres_W) == 3
    presentation(W)
    ppW, iso2 = prune_with_map(pres_W)
    @test is_welldefined(iso2)
    @test is_bijective(iso2)
    @test ppW isa SubquoModule
    @test ngens(ppW) == 2
  end

  @testset "Ext and Tor over ZZ" begin
    F = FreeMod(ZZ, 1)
    M, _ = sub(F, [2*F[1]])
    N, _ = quo(F, [3*F[1]])

    T0 = tor(M, N, 0)
    T1 = tor(M, N, 1)

    @test T0 isa SubquoModule
    @test T1 isa SubquoModule

    E0 = ext(M, N, 0)
    E1 = ext(M, N, 1)

    @test E0 isa SubquoModule
    @test E1 isa SubquoModule
  end

  @testset "Free resolutions over ZZ" begin
    F = FreeMod(ZZ, 1)
    M, _ = quo(F, [4*F[1]])

    fr = free_resolution(M, length=3)

    @test length(fr.C.maps) == 3
    @test iszero(homology(fr.C)[1])
  end

  @testset "Tensoring morphisms over ZZ" begin
    R = ZZ
    F2 = FreeMod(R, 2)
    F4 = FreeMod(R, 4)
    A1 = matrix(ZZ, [1 2; 3 4; 5 6])
    B1 = matrix(ZZ, [7 8])
    A2 = matrix(ZZ, [2 0 1; 0 1 3; 1 1 1])
    B2 = matrix(ZZ, [0 1 1])
    M1 = SubquoModule(F2, A1, B1)
    F3 = FreeMod(R, 3)
    M2 = SubquoModule(F3, A2, B2)
    M, pure_M = tensor_product(M1, M2, task=:map)
    A3 = matrix(ZZ, [1 2; 3 4])
    M3 = SubquoModule(Oscar.SubModuleOfFreeModule(F2, A3))
    N, pure_N = tensor_product(M3, F4, task=:map)
    M3_to_M1 = SubQuoHom(M3, M1, matrix(ZZ, [1 0 0; 0 1 0]))
    @test is_welldefined(M3_to_M1)
    F4_to_M2 = FreeModuleHom(F4, M2, matrix(ZZ, [1 0 0; 0 1 0; 0 0 1; 0 0 0]))
    phi2 = hom_tensor(N, M, [M3_to_M1, F4_to_M2])
    u1 = M3[1]
    u2 = F4[1]
    @test phi2(pure_N((u1, u2))) == pure_M((M3_to_M1(u1), F4_to_M2(u2)))
    F3 = FreeMod(R, 3)
    A1 = matrix(ZZ, [1 2; 3 4; 5 6])
    B1 = matrix(ZZ, [7 8])
    A2 = matrix(ZZ, [2 0 1; 0 1 3; 1 1 1])
    B2 = matrix(ZZ, [0 1 1])
    M1 = SubquoModule(F2, A1, B1)
    M2 = SubquoModule(F3, A2, B2)
    M, pure_M = tensor_product(M1, M2, task=:map)
    phi_M1 = hom(M1, M1, identity_matrix(ZZ, 3))
    phi_M2 = hom(M2, M2, identity_matrix(ZZ, 3))
    phi = hom_tensor(M, M, [phi_M1, phi_M2])
    @test is_welldefined(phi)
    @test is_bijective(phi)
    for u1 in gens(M1), u2 in gens(M2)
      input_elem = pure_M((u1, u2))
      output_elem = pure_M((phi_M1(u1), phi_M2(u2)))
      @test phi(input_elem) == output_elem
    end
  end

  @testset "Direct product over ZZ" begin
    R = ZZ
    F2 = FreeMod(R, 2)
    F3 = FreeMod(R, 3)
    A1 = matrix(ZZ, [1 0; 0 1; 2 2])
    B1 = matrix(ZZ, [0 1])
    M1 = SubquoModule(F2, A1, B1)
    A2 = matrix(ZZ, [1 2 3; 0 0 1])
    B2 = matrix(ZZ, [1 0 1])
    M2 = SubquoModule(F3, A2, B2)
    sum_M, emb = direct_sum(M1, M2)
    @test domain(emb[1]) === M1
    @test domain(emb[2]) === M2
    @test codomain(emb[1]) === sum_M
    @test codomain(emb[2]) === sum_M
    sum_M2, pr = direct_sum(M1, M2, task=:prod)
    @test codomain(pr[1]) === M1
    @test codomain(pr[2]) === M2
    @test domain(pr[1]) === sum_M2
    @test domain(pr[2]) === sum_M2
    prod_M, emb, pr = direct_sum(M1, M2, task=:both)
    @test length(pr) == length(emb) == 2
    @test ngens(prod_M) == ngens(M1) + ngens(M2)
    for g in gens(prod_M)
      @test g == sum([emb[i](pr[i](g)) for i in 1:length(pr)])
    end
    prod_N = direct_product(M1, M2, task=:none)
    @test ngens(prod_N) == ngens(M1) + ngens(M2)
    for g in gens(prod_N)
      @test g == sum([canonical_injection(prod_N, i)(canonical_projection(prod_N, i)(g)) for i in 1:2])
    end
  end



  @testset "Module Constructors over QQ" begin
    F0 = free_module(FreeMod, QQ,3)
    @test rank(F0) == 3
    @test base_ring(F0) == QQ

    F = FreeMod(QQ, 3)
    M, inc = sub(F, [QQ(1//2)*F[1], QQ(2//3)*F[2]])
    N, proj = quo(F, [F[1] + QQ(1//3)*F[3]])

    @test ngens(M) == 2
    @test ngens(N) == 3
    F = FreeMod(QQ, 3)
    M, inc = sub(F, [QQ(1//2)*F[1], QQ(2//3)*F[2]])
    N, proj = quo(F, [F[1] + QQ(1//3)*F[3]])

    @test is_welldefined(inc)
    @test is_welldefined(proj)
    @test is_welldefined(inc)
    @test is_welldefined(proj)
  end

  @testset "Module Homomorphisms over QQ" begin
    F1 = FreeMod(QQ, 2)
    F2 = FreeMod(QQ, 2)
    phi = hom(F1, F2, [F2[1] + QQ(1//2)*F2[2], QQ(3//4)*F2[1]])

    @test is_welldefined(phi)

    K, emb = kernel(phi)
    @test is_welldefined(emb)

    I, im = image(phi)
    @test is_welldefined(im)
  end

  @testset "Presentations over QQ" begin
    F = FreeMod(QQ, 2)
    M, inc = sub(F, [QQ(1//2)*F[1] + QQ(1//3)*F[2]])
    pres_M, iso = present_as_cokernel(M, :both)

    @test is_welldefined(iso)
    @test is_bijective(iso)
    @test pres_M isa SubquoModule
    @test ngens(pres_M) == 1

    F = FreeMod(QQ, 3)
    M, inc = sub(F, [5*F[1], 3*F[2], F[3]])
    N, inc = sub(F, [10*F[1], 3*F[2]])
    W, _ = quo(M, N)
    pres_W, iso = present_as_cokernel(W, :both)
    MP, iso2 = prune_with_map(pres_W)
    @test is_welldefined(iso2)
    @test is_bijective(iso2)
    @test MP isa SubquoModule
    @test ngens(MP) == 1
  end

  @testset "Ext and Tor over QQ" begin
    F = FreeMod(QQ, 1)
    M, _ = sub(F, [QQ(1//2)*F[1]])
    N, _ = quo(F, [QQ(1//3)*F[1]])

    T0 = tor(M, N, 0)
    T1 = tor(M, N, 1)

    @test T0 isa SubquoModule
    @test T1 isa SubquoModule

    E0 = ext(M, N, 0)
    E1 = ext(M, N, 1)

    @test E0 isa SubquoModule
    @test E1 isa SubquoModule
  end

  @testset "Free resolutions over QQ" begin
    F = FreeMod(QQ, 1)
    M, _ = quo(F, [QQ(1//4)*F[1]])

    fr = free_resolution(M, length=3)

    @test length(fr.C.maps) == 3
    @test iszero(homology(fr.C)[1])
    MP, _ = prune_with_map(M)
    @test is_zero(MP)
  end

  @testset "Tensoring morphisms over QQ" begin
    R = QQ
    F2 = FreeMod(R, 2)
    F3 = FreeMod(R, 3)
    A1 = matrix(R, [1 2; 3 4; 5 6])
    B1 = matrix(R, [7 8])
    A2 = matrix(R, [2 0 1; 0 1 3; 1 1 1])
    B2 = matrix(R, [0 1 1])
    M1 = SubquoModule(F2, A1, B1)
    M2 = SubquoModule(F3, A2, B2)
    M, pure_M = tensor_product(M1, M2, task=:map)
    phi_M1 = SubQuoHom(M1, M1, matrix(R, [1 1 0; 0 2 1; 1 0 1]))
    phi_M2 = SubQuoHom(M2, M2, matrix(R, [0 1 0; 1 0 1; 1 1 0]))
    phi = hom_tensor(M, M, [phi_M1, phi_M2])
    v = M[1] + 2*M[2]
    @test phi(v) == pure_M((phi_M1(M1[1]), phi_M2(M2[1]))) + 2*pure_M((phi_M1(M1[1]), phi_M2(M2[2])))
    F4 = FreeMod(R, 4)
    A3 = matrix(R, [1 2; 3 4])
    M3 = SubquoModule(Oscar.SubModuleOfFreeModule(F2, A3))
    N, pure_N = tensor_product(M3, F4, task=:map)
    M3_to_M1 = SubQuoHom(M3, M1, matrix(R, [1 2 0; 0 1 3]))
    @test is_welldefined(M3_to_M1)
    F4_to_M2 = FreeModuleHom(F4, M2, matrix(R, [1 0 2; 0 1 0; 1 1 1; 2 0 0]))
    phi2 = hom_tensor(N, M, [M3_to_M1, F4_to_M2])
    u1 = M3[1]
    u2 = F4[1]
    @test phi2(pure_N((u1, u2))) == pure_M((M3_to_M1(u1), F4_to_M2(u2)))
  end

  @testset "Direct product over QQ" begin
    R = QQ
    F2 = FreeMod(R, 2)
    F3 = FreeMod(R, 3)
    A1 = matrix(R, [1 0; 0 1; 2 2])
    B1 = matrix(R, [0 1])
    M1 = SubquoModule(F2, A1, B1)
    A2 = matrix(R, [1 2 3; 0 0 1])
    B2 = matrix(R, [1 0 1])
    M2 = SubquoModule(F3, A2, B2)
    sum_M, emb = direct_sum(M1, M2)
    @test domain(emb[1]) === M1
    @test domain(emb[2]) === M2
    @test codomain(emb[1]) === sum_M
    @test codomain(emb[2]) === sum_M
    sum_M2, pr = direct_sum(M1, M2, task=:prod)
    @test codomain(pr[1]) === M1
    @test codomain(pr[2]) === M2
    @test domain(pr[1]) === sum_M2
    @test domain(pr[2]) === sum_M2
    prod_M, emb, pr = direct_sum(M1, M2, task=:both)
    @test length(pr) == length(emb) == 2
    @test ngens(prod_M) == ngens(M1) + ngens(M2)
    for g in gens(prod_M)
      @test g == sum([emb[i](pr[i](g)) for i in 1:length(pr)])
    end
    prod_N = direct_product(M1, M2, task=:none)
    @test ngens(prod_N) == ngens(M1) + ngens(M2)
    for g in gens(prod_N)
      @test g == sum([canonical_injection(prod_N, i)(canonical_projection(prod_N, i)(g)) for i in 1:2])
    end
  end

  @testset "gradings over non-graded rings" begin
    @test_throws AssertionError graded_free_module(ZZ, [1])
  end

  @testset "size of modules" begin
    R = ZZ
    F = free_module(FreeMod, R, 2)
    M = cokernel(hom(F, F, matrix(ZZ, [2 0; 0 3])))
    @test size(M) == 6

    N = cokernel(hom(F, F, matrix(ZZ, [7 0; 0 0])))
    @test size(N) == PosInf()

    K, a = finite_field(7, "a")
    G = free_module(FreeMod, K, 3)
    H = free_module(FreeMod, K, 1)
    L = cokernel(hom(H, G, matrix(K, [1 0 0])))
    @test size(L) == 49

    K, a = finite_field(7, "a")
    G = free_module(FreeMod, K, 3)
    L = cokernel(hom(G, G, matrix(K, [1 0 0; 0 1 0; 0 0 1])))
    @test size(L) == 1

    G = free_module(FreeMod, QQ, 3)
    L = cokernel(hom(G, G, matrix(QQ, [1 0 0; 0 1 0; 0 0 1])))
    @test size(L) == 1
  end
end

@testset "minimization of resolutions" begin
  R, (x, y, z) = QQ[:x, :y, :z]

  @test AbstractAlgebra.is_known(is_local, R)
  @test !is_local(R)

  A = R[x y z; y-1 z-1 x]
  I = ideal(R, minors(A, 2))
  Q, _ = quo(R, I)
  @test !AbstractAlgebra.is_known(is_local, Q)
  M = quotient_ring_as_module(Q)
  res = free_resolution(M)
  @test_throws ErrorException minimize(res)

  U1 = complement_of_point_ideal(R, [0, 0, 0])
  L1, _ = localization(R, U1)
  @test AbstractAlgebra.is_known(is_local, L1)
  M1, _ = change_base_ring(L1, M)
  res1 = free_resolution(M1)
  res1[3]
  res1min = minimize(res1)
  @test [ngens(res1min[i]) for i in 0:2] == [1, 2, 1]

  U2 = complement_of_prime_ideal(ideal(R, gens(R)))
  L2, _ = localization(R, U2)
  @test AbstractAlgebra.is_known(is_local, L2)
  M2, _ = change_base_ring(L2, M)
  res2 = free_resolution(M2)
  res2min = minimize(res2)
  @test [ngens(res2min[i]) for i in 0:2] == [1, 2, 1]

  f = x^2 + y^2 + z^2
  I = ideal(f)
  Q, _ = quo(R, f)
  M0, _ = change_base_ring(Q, M)

  res = free_resolution(M0; length=5)
  @test_throws ErrorException minimize(res)

  L1, _ = localization(Q, U1)
  @test !AbstractAlgebra.is_known(is_local, L1)
  M1, _ = change_base_ring(L1, M)
  res1 = free_resolution(M1; length=5)
  res1[3]
  res1min = minimize(res1)
  @test [ngens(res1min[i]) for i in 0:2] == [1, 2, 1]

  U2 = complement_of_prime_ideal(ideal(R, gens(R)))
  L2, _ = localization(Q, U2)
  @test !AbstractAlgebra.is_known(is_local, L2)
  M2, _ = change_base_ring(L2, M)
  res2 = free_resolution(M2; length=5)
  res2min = minimize(res2)
  @test [ngens(res2min[i]) for i in 0:2] == [1, 2, 1]
end

@testset "Krull dimension for subquos" begin
    R, (a, b, c, d) = graded_polynomial_ring(QQ, [:a, :b, :c, :d])
    I = ideal(R, [a*d - b*c])
    M = quotient_ring_as_module(I)

    E2 = ext(M, graded_free_module(R, 1), 0 + 2)
    @test is_zero(E2)
    @test krull_dim(E2) == -inf

    E3 = ext(M, graded_free_module(R, 1), 1 + 2)
    @test is_zero(E3)
    @test krull_dim(E3) == -inf

    @test typeof(E2) == typeof(E3)
end

@testset "New minimization of free resolutions" begin
  R, (x, y, z, w) = QQ[:x, :y, :z, :w]
  U = complement_of_point_ideal(R, [0, 0, 0, 0])
  L, _ = localization(R, U)
  A = L[x y z; y-1 z-1 w]
  I = ideal(L, minors(A, 2))
  II = ideal_as_module(I)
  M = quotient_ring_as_module(I)
  res = free_resolution(M; length=1)
  @test ngens(res[0]) == 1
  @test ngens(res[1]) == 3
  @test ngens(res[2]) == 2

  min_res = minimize(res)
  @test min_res isa FreeResolution
  @test ngens(min_res[0]) == 1
  @test ngens(min_res[1]) == 2
  @test ngens(min_res[2]) == 1
end

@testset "caching of kernels" begin
  R, (x, y) = QQ[:x, :y]
  F = free_module(R, 1)
  phi = hom(F, F, [x*F[1]])
  I, _ = kernel(phi)
  @test I === kernel(phi)[1]
  phi = hom(F, F, [x*F[1]])
  I, _ = kernel(phi; cached=false)
  @test I !== kernel(phi)[1]
end

@testset "Issue #5490" begin
  # This is about saturation not working properly over quotient rings
  P, (a, b, x, y, z) = QQ[:a, :b, :x, :y, :z]
  I = ideal(P, [x^3 + x*y^2 - z^2, a*y - b*z, b*x - y, a*x - z, -a*z + x^2 + y^2, -a^2 + b*y + x, -a^3 + b^2*z + z, -a^2*b + b^2*y + y])
  A, _ = quo(P, I)

  FA = free_module(A, 8)
  rel_mat = [0 0 0 3*a*z - 2*y^2 2*x*y -2*z 0 0; 0 0 0 0 3*a*z - 2*y^2 0 2*x*y -2*z; -6*a*z^2 + 4*y^2*z -4*x*y*z 4*z^2 -6*x*z -4*y*z 3*a*z - 2*y^2 -2*x*z 2*x*y; 0 0 0 6*x*y^2 - 9*z^2 2*y^3 6*x*z 4*x*y^2 -4*y*z; 0 0 0 -3*a*z + 2*y^2 -2*x*y 2*z 0 0; 0 0 0 0 -2*x*y^2 + 3*z^2 0 2*b*z^2 - 2*y^3 -2*x*z; 4*x*y^2*z - 6*z^3 -4*b*z^3 + 4*y^3*z 4*x*z^2 -6*a*z^2 + 6*y^2*z -4*x*y*z -2*x*y^2 + 3*z^2 -2*a*z^2 + 2*y^2*z 2*b*z^2 - 2*y^3; 0 0 0 -6*x*y^2 + 9*z^2 -2*y^3 -6*x*z -4*x*y^2 4*y*z]
  rel_mat = matrix_space(A, 8, 8)(rel_mat)
  U, inc_U = sub(FA, rel_mat)

  J = ideal(A, x)
  U_sat = Oscar._saturation(U.sub, J)
  U_sat, _ = sub(FA, gens(U_sat))
  @test all(x in U_sat for x in gens(U))
end

