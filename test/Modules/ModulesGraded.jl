using Random

RNG = Random.MersenneTwister(42)

##################################################################
# Tests graded modules
##################################################################

@testset "Graded free modules constructors" begin
    Rg, _ = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    M1 = graded_free_module(Rg, 2)
    @test degrees(M1) == [Z[0], Z[0]]
    M2 = graded_free_module(Rg, [2,2,4,3,3,5])
    @test degrees(M2) == [2*Z[1], 2*Z[1], 4*Z[1], 3*Z[1], 3*Z[1], 5*Z[1]]
    M3 = grade(M2,[2,2,4,3,3,6])
    @test degrees(M3) == [2*Z[1], 2*Z[1], 4*Z[1], 3*Z[1], 3*Z[1], 6*Z[1]]
    @test degrees_of_generators(M3) == [2*Z[1], 2*Z[1], 4*Z[1], 3*Z[1], 3*Z[1], 6*Z[1]]
    @test is_graded(M1)
    @test is_graded(M2)
    @test is_graded(M3)
end

@testset "Graded free modules hom constructor" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 3)
    G = graded_free_module(Rg, 2)
    V = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
    a = hom(F, G, V)
    @test is_graded(a)
    @test domain(a) == F
    @test codomain(a) == G
    @test a(F[1]) == V[1]
    W = Rg[y 0; x y; 0 z]
    a = hom(F, G, W)
    @test is_graded(a)
    @test domain(a) == F
    @test codomain(a) == G
    @test a(F[1]) == y*G[1]
end

@testset "Kernel and image free modules hom" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 3)
    G = graded_free_module(Rg, 2)
    V = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
    a = hom(F, G, V)
    Mk = kernel(a)
    Mi = image(a)
    @test is_graded(Mk[1])
    @test is_graded(Mi[1])
    @test ngens(Mk[1]) == 1
    @test ngens(Mi[1]) == 3
    @test codomain(Mk[2]) == F
    @test codomain(Mi[2]) == G
end

@testset "Presentations 1" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, 3)
    G = graded_free_module(Rg, 2)
    V = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
    a = hom(F, G, V)
    Mk = kernel(a)
    Mk1 = Mk[1]
    p = presentation(Mk1)
    @test all(is_graded(p[i]) for i in -2:1)
    @test degrees(p[0]) == [2 * Z[1]]
    Mp, f = present_as_cokernel(Mk1, :all)
    @test is_graded(Mp)
    @test is_graded(f)
    @test degrees_of_generators(Mp) == [2 * Z[1]]
end

@testset "Degrees of free module homomorphisms" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, 3)
    G = graded_free_module(Rg, 2)
    V = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
    a = hom(F, G, V)
    @test degree(a) == Z[1]
    @test is_homogeneous(a) == false
    F = graded_free_module(Rg, [1,1,1])
    G = graded_free_module(Rg, [0,0])
    V = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
    a = hom(F, G, V)
    @test is_graded(a)
    @test degree(a) == Z[0]
    @test is_homogeneous(a)
end

@testset "Subquotient constructors" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, [1,1,1])
    G = graded_free_module(Rg, [0,0])
    B = Rg[y 0; x y; 0 z]
    b = hom(F, G, B)
    F1 = graded_free_module(Rg, [2,2,2])
    F2 = graded_free_module(Rg, [2])
    G = graded_free_module(Rg, [1,1])
    V1 = [y*G[1], (x+y)*G[1]+y*G[2], z*G[2]]
    V2 = [z*G[2]+y*G[1]]
    a1 = hom(F1, G, V1)
    a2 = hom(F2, G, V2)
    V = subquotient(a1,a2)
    @test is_graded(V)
    @test degrees_of_generators(V) == [2*Z[1], 2*Z[1], 2*Z[1]]
    A1 = Rg[x y;
        2*x^2 3*y^2]
    A2 = Rg[x^3 x^2*y;
        (2*x^2+x*y)*x (2*y^3+y*x^2)]
    B = Rg[4*x*y^3 (2*x+y)^4]
    F2 = graded_free_module(Rg, [0,0])
    M1 = SubquoModule(F2, A1, B)
    @test degrees_of_generators(M1) == [Z[1], 2*Z[1]]
    F = graded_free_module(Rg, 2)
    O = [x*F[1]+y*F[2],y*F[2]]
    M = SubquoModule(F, O)
    @test degrees_of_generators(M) == [Z[1], Z[1]]
    O1 = [x*F[1]+y*F[2],y*F[2]]
    O2 = [x^2*F[1]+y^2*F[2],y^2*F[2]]
    M = SubquoModule(F, O1, O2)
    @test degrees_of_generators(M) == [Z[1], Z[1]]
    @test degrees_of_generators(M.quo) == [2*Z[1], 2*Z[1]]
end

@testset "Graded modules from matrices constructors" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F2 = graded_free_module(Rg,[8,8])
    A1 = Rg[x y;
        2*x^2 3*y^2]
    MM = graded_cokernel(F2, A1)
    @test degrees_of_generators(MM) == [8*Z[1], 8*Z[1]]
    MM = cokernel(F2, A1)
    @test degrees_of_generators(MM) == [8*Z[1], 8*Z[1]]
    MM = graded_cokernel(A1)
    @test degrees_of_generators(MM) == [Z[0], Z[0]]
    MM = graded_image(F2, A1)
    @test degrees_of_generators(MM) == [9*Z[1], 10*Z[1]]
    MM = image(F2, A1)
    @test degrees_of_generators(MM) == [9*Z[1], 10*Z[1]]
    MM = graded_image(A1)
    @test degrees_of_generators(MM) == [Z[1], 2*Z[1]]
end

@testset "Graded morphisms from matrix constructor" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F = graded_free_module(Rg,[1,2, 3])
    A2 = Rg[x^3 x^2 x;
       (2*x^2+x*y)*x^2 (2*y^2+x^2)*x x^2]
    f = graded_map(F, A2)
    @test degrees(domain(f)) == [4*Z[1], 5*Z[1]]
    @test degrees(codomain(f)) == [Z[1], 2*Z[1], 3*Z[1]]
end

@testset "Graded morphisms of subquotients constructor" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, 1);
    A = Rg[x; y]
    B = Rg[x^2; y^3; z^4]
    M = SubquoModule(F, A, B)
    N = M;
    # Problem with the previous test: V[2] is zero and
    # the homomorphism is hence not graded.
    # V = [y^2*N[1], x^2*N[2]]
    V = [y^2*N[1], x*y*N[2]]
    a = hom(M, N, V);
    @test is_graded(a)
    @test degree(a) == 2*Z[1]
end

@testset "Graded morphisms of subquotients and kernel" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1);
    A = Rg[x; y]
    B = Rg[x^2; y^3; z^4]
    M = SubquoModule(F, A, B)
    N = M;
    #V = [y^2*N[1], x^2*N[2]]
    V = [y^2*N[1], x*y*N[2]]
    a = hom(M, N, V);
    K, incl = kernel(a);
    @test ngens(K) == 2
    @test domain(incl) == K
end

@testset "Is isomorphic for graded free modules" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, [1,1,3,2])
    G = graded_free_module(Rg, [1,1,2,3])
    @test is_isomorphic(F, G)
    F = graded_free_module(Rg, [1,1,5,2])
    G = graded_free_module(Rg, [1,1,2,3])
    @test !is_isomorphic(F, G)
end

@testset "Canonical isomorphism between graded subquotients" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F1 = graded_free_module(Rg,[2,3, 4])
    A2 = Rg[x^3 x^2 x;
       (2*x^2+x*y)*x^2 (2*y^2+x^2)*x x^2]
    MM1 = image(F1, A2)
    F2 = graded_free_module(Rg,[2,4, 3])
    A2a = Rg[x^3 x x^2;
       (2*x^2+x*y)*x^2 x^2 (2*y^2+x^2)*x]
    MM2 = image(F2, A2a)
    @test is_canonically_isomorphic(MM1.sub, MM2.sub)
    _ , f = is_canonically_isomorphic_with_map(MM1, MM2)
    @test f(MM1[1]) == MM2[1]
end

@testset "Being a subset and equality for graded subquotients" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg,2)
    O1 = [x*F[1]+y*F[2],y*F[2]]
    O1a = [x*F[1],y*F[2]]
    O2 = [x^2*F[1]+y^2*F[2],y^2*F[2]]
    M1 = SubquoModule(F, O1, O2)
    M2 = SubquoModule(F, O1a, O2)
    @test M1 == M2
    @test is_subset(M1,M2)
    @test is_subset(M2,M1)
    @test M1.sub == M2.sub
    @test is_subset(M1.sub,M2.sub)
    @test is_subset(M2.sub,M1.sub)
end

@testset "Sum of graded subquotients" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    AM = Rg[x;]
    BM = Rg[x^2; y^3; z^4]
    M = SubquoModule(F, AM, BM)
    AN = Rg[y;]
    BN = Rg[x^2; y^3; z^4]
    N = SubquoModule(F, AN, BN)
    O, f1, f2 = sum(M, N)
    @test ngens(O) == 2
    @test codomain(f1) == O
    @test codomain(f2) == O
    @test domain(f1) == M
    @test domain(f2) == N
end


@testset "Being well defined and equality for graded morphisms of subquotients" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    A = Rg[x; y]
    B = Rg[x^2; y^3; z^4]
    M = SubquoModule(F, A, B)
    N = M
    # V = [y^2*N[1], x^2*N[2]]
    V = [y^2*N[1], x*y*N[2]]
    a = hom(M, N, V)
    @test is_welldefined(a)
    #W = Rg[y^2 0; 0 x^2]
    W = Rg[y^2 0; 0 x*y]
    b = hom(M, N, W)
    @test a == b
    @test nrows(matrix(a)) == 2
    @test ncols(matrix(a)) == 2
    F = graded_free_module(Rg, 1)
    W = [y*N[1], x*N[2]]
    c = hom(M, N, W);
    @test !is_welldefined(c)
end

@testset "Being zero for elements of graded subquotients" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = free_module(Rg, 1)
    A = Rg[x; y]
    B = Rg[x^2; y^3; z^4]
    M = SubquoModule(F, A, B)
    @test !is_zero(M[1])
    @test is_zero(x*M[1])
end

@testset "Hom module of graded free modules" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F1 = graded_free_module(Rg, [1,2,2])
    F2 = graded_free_module(Rg, [3,5])
    V, f = hom(F1, F2)
    @test is_graded(V)
    @test degrees(V) == [2*Z[1], 4*Z[1], Z[1], 3*Z[1], Z[1], 3*Z[1]]
    F1 = graded_free_module(Rg, 3)
    F2 = graded_free_module(Rg, 5)
    V = hom(F1, F2)
    g = V[2](V[1][1])
    @test is_graded(g)
    @test g(F1[1])==F2[1]
end

@testset "Presentations 2" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, [1,2,2]);
    p = presentation(F)
    @test rank(p[-2]) == 0
    @test degrees(p[-1]) == [Z[1], 2*Z[1], 2*Z[1]]
    @test degrees(p[0]) == [Z[1], 2*Z[1], 2*Z[1]]
    @test rank(p[1]) == 0
    @test is_zero(map(p,-1))
    @test degree(map(p,0)) == Z[0]
    @test is_zero(map(p,1))
    F = graded_free_module(Rg, 1);
    A = Rg[x; y];
    B = Rg[x^2; y^3; z^4];
    M = SubquoModule(F, A, B);
    P = presentation(M);
    @test rank(P[-2]) == 0
    @test is_graded(P[-2])
    @test P[-1] == M
    @test degrees(P[0]) == [Z[1], Z[1]]
    @test degrees(P[1]) == [2*Z[1], 2*Z[1], 3*Z[1], 5*Z[1], 5*Z[1]]
    @test is_homogeneous(map(P, -1))
    @test is_homogeneous(map(P, 0))
    @test is_homogeneous(map(P, 1))
    @test domain(map(P, -1)) == M
    @test domain(map(P, 0)) == P[0]
    @test codomain(map(P, 0)) == M
    @test codomain(map(P, 1)) == P[0]
    @test domain(map(P, 1)) == P[1]
    N, phi = present_as_cokernel(M, :all)
    @test degrees_of_generators(N) == [Z[1], Z[1]]
    @test is_homogeneous(phi)
    @test domain(phi) == N
    @test codomain(phi) == M
    psi = phi.inverse_isomorphism
    @test is_homogeneous(psi)
    @test domain(psi) == M
    @test codomain(psi) == N
    @test is_bijective(phi)
    @test is_bijective(psi)
  I = ideal(Rg, [x^2, x*y, y])
  A, _ = quo(Rg, I)
  M = quotient_ring_as_module(A)
  mp = presentation(M, minimal = true)
  @test rank(mp[0]) == 1
  @test rank(mp[1]) == 2
end


@testset "Equality morphism for subquotients" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1);
    A = Rg[x-y; x+y];
    B = Rg[x^2; y^3; z^4];
    M = SubquoModule(F, A, B);
    A1 = Rg[x; y; x^2*y+y^3];
    B1 = Rg[x^2; y^3; z^4];
    M1 = SubquoModule(F, A1, B1);
    f = is_equal_with_morphism(M, M1)
    @test is_homogeneous(f)
    @test domain(f) == M
    @test codomain(f) == M1
    @test is_bijective(f)
end

@testset "Free resolutions" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, 3)
    rs = free_resolution(F)
    @test is_graded(rs)
    @test rank(rs[0]) == 3
    @test rank(rs[-1]) == 3
    A = Rg[x; y]
    B = Rg[x^2; x*y; y^2; z^4]
    F = graded_free_module(Rg, 1)
    M = SubquoModule(F, A, B)
    fr = free_resolution_via_kernels(M)
    @test is_graded(fr)
    @test degrees(fr[0]) == [Z[1], Z[1]]
    @test degrees(fr[1]) == vcat(fill(2*Z[1], 4), fill(5*Z[1], 2))
    @test degrees(fr[2]) == vcat(fill(3*Z[1], 2), fill(6*Z[1], 4))
    @test degrees(fr[3]) == fill(7*Z[1], 2)
end

@testset "Hom module of graded subquotient modules" begin
    Rg, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, 2);
    V = [x*F[1], y^2*F[2]];
    M = quo(F, V)[1]
    H, f = hom(M, M)
    @test is_graded(H)
    @test degrees_of_generators(H) == [Z[0], Z[0]]
    @test degrees_of_generators(H.quo) == [Z[1], 2*Z[1], Z[1], 2*Z[1]]
    @test is_homogeneous(f(H[1]))
    a = element_to_homomorphism(x*H[1] + y*H[2])
    @test matrix(a) == Rg[x 0; 0 y]
    W =  [x*M[1], y*M[2]];
    a = hom(M, M, W);
    @test degree(a) == Z[1]
    @test is_welldefined(a)
    m = homomorphism_to_element(H, a)
    @test m == y*H[2]
end

@testset "Dual and double dual" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F1 = graded_free_module(Rg, 1)
    F2 = graded_free_module(Rg, 2)
    F2v, ev = dual(F2, codomain=F1)
    @test ev(F2v[1])(F2[1]) == F1[1]
    FF, psi = double_dual(F2)
    @test degrees_of_generators(FF) == [Z[0], Z[0]]
    @test is_injective(psi)
    M, inc = sub(F2, [x*F2[1], y*F2[1]])
    F1 = graded_free_module(Rg, 1)
    Mv, ev = dual(M, codomain=F1)
    @test degrees_of_generators(Mv) == [Z[0]]
    @test ev(Mv[1])(M[1]) == x*F1[1]
    Mvv, psi = double_dual(M, codomain=F1)
    @test matrix(psi) == Rg[x; y]
    F = graded_free_module(Rg, 2);
    V = [x*F[1], y^2*F[2]];
    M = quo(F, V)[1]
    dM, f = dual(M)
    @test is_graded(dM)
    @test is_zero(dM)
    @test domain(f)== dM
    ddM, g = double_dual(M)
    @test is_graded(ddM)
    @test is_zero(ddM)
    @test codomain(g)== ddM
    F = graded_free_module(Rg, 1);
    M = SubquoModule(F, [(x^2-y^2)*F[1], (x^2+y^2)*F[1]]);
    dM, f = dual(M)
    @test is_graded(dM)
    @test degrees_of_generators(M) == [2*Z[1], 2*Z[1]]
    @test domain(f)== dM
    ddM, g = double_dual(M)
    @test is_graded(ddM)
    @test degrees_of_generators(ddM) == [Z[0]]
    @test codomain(g)== ddM
    @test is_injective(g)
end

@testset "Hom module element morphism conversion" begin
    Rg, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
    A1 = Rg[x^2+y^2 x*y; x^2+y^2 x*y]
    B1 = Rg[x*y^3+x^4+y^2*x^2 x^4; y^4-3*x^4 x*y^3-x^4]
    F = graded_free_module(Rg, 2)
    M1 = SubquoModule(F, A1, B1)
    A2 = Rg[x;]
    B2 = Rg[y^3;]
    F1 = graded_free_module(Rg, 1)
    M2 = SubquoModule(F1, A2,B2)
    SQ = hom(M1,M2)[1] # This is the zero module
end

@testset "Multiplication morphism" begin
    Rg, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
    A1 = Rg[x^2+y^2 x*y; x^2+y^2 x*y]
    B1 = Rg[x*y^3+x^4+y^2*x^2 x^4; y^4-3*x^4 x*y^3-x^4]
    F = graded_free_module(Rg, 2)
    M1 = SubquoModule(F, A1, B1)
    End_M = hom(M1,M1)[1]
    R_as_module = graded_free_module(Rg,1)
    phi = multiplication_induced_morphism(R_as_module, End_M)
    @test is_homogeneous(phi)
    @test matrix(phi) == Rg[1 0]
end

@testset "Graded tensor product of graded free modules" begin
    Rg, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
    Z = grading_group(Rg)
    F1 = graded_free_module(Rg, [1,2]);
    F2 = graded_free_module(Rg, [3]);
    F3 = graded_free_module(Rg, [5,6]);
    tensor_product(F1, F2, F3)
    M,t = tensor_product(F1, F2, F3 ,task=:map)
    v = t((F1[1],F2[1],F3[1]))
    @test degree(v) == 9*Z[1]
    @test degrees_of_generators(M) == [9*Z[1], 10*Z[1], 10*Z[1], 11*Z[1]]
end

@testset "Graded tensor product of graded subquotients" begin
    Rg, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
    Z = grading_group(Rg)
    A1 = Rg[x^2+y^2 x*y; x^2+y^2 x*y]
    B1 = Rg[x*y^3+x^4+y^2*x^2 x^4; y^4-3*x^4 x*y^3-x^4]
    F = graded_free_module(Rg, 2)
    M1 = SubquoModule(F, A1, B1)
    A2 = Rg[x;]
    B2 = Rg[y^3;]
    F1 = graded_free_module(Rg, 1)
    M2 = SubquoModule(F1, A2, B2)
    M, pure_M = tensor_product(M1, M2, task=:map)
    @test is_graded(M)
    v = pure_M((M1[1],M2[1]))
    @test degree(v) == 3*Z[1]
    @test degrees_of_generators(M) == [3*Z[1], 3*Z[1]]
end

@testset "Betti tables 1" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    A = Rg[x; y]
    B = Rg[x^2; x*y; y^2; z^4]
    M = SubquoModule(F, A, B)
    free_res1 = free_resolution(M, length=1)
    free_res2 = free_resolution(M)
    free_res3 = free_resolution_via_kernels(M)
    @test first(chain_range(free_res1)) == 1
    @test first(chain_range(free_res2)) == 4
    @test first(chain_range(free_res3)) == 4
end

@testset "Betti tables 2" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    A = Rg[x; y]
    B = Rg[x^2; x*y; y^2; z^4]
    F = graded_free_module(Rg, 1)
    M = SubquoModule(F, A, B)
    fr = free_resolution(M)
    T1 = Dict{Tuple{Int64, Any}, Int64}(
    (0, Z[1]) => 2,
    (3, 7*Z[1]) => 2,
    (2, 3*Z[1]) => 2,
    (1, 2*Z[1]) => 4,
    (1, 5*Z[1]) => 2,
    (2, 6*Z[1]) => 4)
    T = betti_table(fr)
    @test as_dictionary(T) == T1
    T = betti(fr)
    @test as_dictionary(T) == T1
    reverse_direction!(T)
    @test as_dictionary(T) == T1
    fr2 = free_resolution_via_kernels(M)
    @test as_dictionary(betti(fr2)) == T1
end

@testset "Resolution and Betti tables 2" begin
    K = GF(31991)
    Rg, (x,y,z,u,v) = graded_polynomial_ring(K, [:x, :y, :z, :u, :v])
    Z = grading_group(Rg)
    I = ideal([10000*x^3-8565*x^2*y-7937*x*y^2+14060*x^2*z+1416*x*y*z-3771*x*z^2-15712*x^2*u-5990*x*y*u+3271*x*z*u-6141*x*u^2-14457*x^2*v-194*x*y*v-15529*y^2*v+12735*x*z*v+11618*y*z*v+9223*z^2*v+14387*x*u*v+7102*y*u*v+13097*z*u*v+3485*u^2*v+12211*x*v^2-4061*y*v^2+2878*z*v^2-6247*u*v^2-1256*v^3,-10000*x^2*y+8565*x*y^2+7937*y^3-14060*x*y*z-1416*y^2*z+3771*y*z^2+15712*x*y*u+5990*y^2*u-3271*y*z*u+6141*y*u^2-x^2*v+3340*x*y*v-797*y^2*v-7781*x*z*v+14337*y*z*v+7675*z^2*v+9192*x*u*v+1747*y*u*v-13195*z*u*v+9126*u^2*v-5013*x*v^2+15010*y*v^2-6899*z*v^2-6373*u*v^2-6580*v^3,-10000*x^2*z+8565*x*y*z+7937*y^2*z-14060*x*z^2-1416*y*z^2+3771*z^3+15712*x*z*u+5990*y*z*u-3271*z^2*u+6141*z*u^2+12149*x^2*v+669*x*y*v+9621*y^2*v+1543*x*z*v+8614*y*z*v-649*z^2*v-12278*x*u*v+14655*y*u*v+10594*z*u*v+10462*u^2*v-2455*x*v^2-10875*y*v^2-7761*z*v^2-10318*u*v^2-6636*v^3,-10000*x^2*u+8565*x*y*u+7937*y^2*u-14060*x*z*u-1416*y*z*u+3771*z^2*u+15712*x*u^2+5990*y*u^2-3271*z*u^2+6141*u^3-9289*x^2*v+13297*x*y*v+13024*y^2*v-12477*x*z*v+9258*y*z*v+10372*z^2*v-4128*x*u*v-3551*y*u*v+4570*z*u*v+15473*u^2*v+4446*x*v^2-7043*y*v^2-14100*z*v^2+5002*u*v^2+9799*v^3,-9289*x^3+13297*x^2*y+13024*x*y^2-12477*x^2*z+9258*x*y*z+10372*x*z^2+13406*x^2*u-3745*x*y*u-15529*y^2*u-14686*x*z*u+11618*y*z*u+9223*z^2*u-2131*x*u^2+7102*y*u^2+13097*z*u^2+3485*u^3+4446*x^2*v-7043*x*y*v-14100*x*z*v-14778*x*u*v-4061*y*u*v+2878*z*u*v-6247*u^2*v+9799*x*v^2-1256*u*v^2,9289*x^2*y-13297*x*y^2-13024*y^3+12477*x*y*z-9258*y^2*z-10372*y*z^2-x^2*u+7468*x*y*u+2754*y^2*u-7781*x*z*u+9767*y*z*u+7675*z^2*u+9192*x*u^2-13726*y*u^2-13195*z*u^2+9126*u^3-4446*x*y*v+7043*y^2*v+14100*y*z*v-5013*x*u*v+10008*y*u*v-6899*z*u*v-6373*u^2*v-9799*y*v^2-6580*u*v^2,9289*x^2*z-13297*x*y*z-13024*y^2*z+12477*x*z^2-9258*y*z^2-10372*z^3+12149*x^2*u+669*x*y*u+9621*y^2*u+5671*x*z*u+12165*y*z*u-5219*z^2*u-12278*x*u^2+14655*y*u^2-4879*z*u^2+10462*u^3-4446*x*z*v+7043*y*z*v+14100*z^2*v-2455*x*u*v-10875*y*u*v-12763*z*u*v-10318*u^2*v-9799*z*v^2-6636*u*v^2,12149*x^3+669*x^2*y+9621*x*y^2-12914*x^2*z+8420*x*y*z-15529*y^2*z+12086*x*z^2+11618*y*z^2+9223*z^3-12278*x^2*u+14655*x*y*u-7010*x*z*u+7102*y*z*u+13097*z^2*u+10462*x*u^2+3485*z*u^2-2455*x^2*v-10875*x*y*v+4450*x*z*v-4061*y*z*v+2878*z^2*v-10318*x*u*v-6247*z*u*v-6636*x*v^2-1256*z*v^2,-12149*x^2*y-669*x*y^2-9621*y^3-x^2*z+1797*x*y*z-9411*y^2*z-7781*x*z^2+14986*y*z^2+7675*z^3+12278*x*y*u-14655*y^2*u+9192*x*z*u-8847*y*z*u-13195*z^2*u-10462*y*u^2+9126*z*u^2+2455*x*y*v+10875*y^2*v-5013*x*z*v-9220*y*z*v-6899*z^2*v+10318*y*u*v-6373*z*u*v+6636*y*v^2-6580*z*v^2,x^3+11117*x^2*y+991*x*y^2+15529*y^3+7781*x^2*z+4919*x*y*z-11618*y^2*z-7675*x*z^2-9223*y*z^2-9192*x^2*u+15857*x*y*u-7102*y^2*u+13195*x*z*u-13097*y*z*u-9126*x*u^2-3485*y*u^2+5013*x^2*v+4770*x*y*v+4061*y^2*v+6899*x*z*v-2878*y*z*v+6373*x*u*v+6247*y*u*v+6580*x*v^2+1256*y*v^2])
    V = gens(I)
    F = graded_free_module(Rg,1)
    V1 = [p * F[1] for p in V]
    M = quo(F, V1)[1]
    free_res0 = free_resolution(M)
    @test is_graded(free_res0)
    @test all(iszero, homology(free_res0.C))
    free_res0a = free_resolution_via_kernels(M)
    @test is_graded(free_res0a)
    @test all(iszero, homology(free_res0a.C))
    # a minimal generating set of I as produced by the old `minimal_generating_set` function
    V = [x^2*u + 25341*x^2*v + 18282*x*y*u + 28685*x*y*v + 29881*x*z*u + 6569*x*z*v + 6499*x*u^2 + 31326*x*u*v + 24012*x*v^2 + 12818*y^2*u + 5322*y^2*v + 28229*y*z*u + 30832*y*z*v + 15643*y*u^2 + 21879*y*u*v + 3926*y*v^2 + 24252*z^2*u + 29725*z^2*v + 25334*z*u^2 + 23257*z*u*v + 15677*z*v^2 + 2082*u^3 + 14404*u^2*v + 10307*u*v^2 + 27796*v^3, x^2*z + 30515*x^2*v + 18282*x*y*z + 16824*x*y*v + 29881*x*z^2 + 6499*x*z*u + 26466*x*z*v + 16451*x*u*v + 7982*x*v^2 + 12818*y^2*z + 22172*y^2*v + 28229*y*z^2 + 15643*y*z*u + 6890*y*z*v + 5453*y*u*v + 25194*y*v^2 + 24252*z^3 + 25334*z^2*u + 22269*z^2*v + 2082*z*u^2 + 6186*z*u*v + 30485*z*v^2 + 29693*u^2*v + 8645*u*v^2 + 30085*v^3, x^2*z + 18213*x^2*u + 3291*x*y*z + 9316*x*y*u + 12069*x*z^2 + 11593*x*z*u + 1353*x*z*v + 20876*x*u^2 + 10094*x*u*v + 9101*y^2*z + 22652*y^2*u + 7221*y*z^2 + 5429*y*z*u + 3059*y*z*v + 31352*y*u^2 + 19564*y*u*v + 1559*z^3 + 11151*z^2*u + 26923*z^2*v + 27410*z*u^2 + 24499*z*u*v + 23597*z*v^2 + 22504*u^3 + 22395*u^2*v + 23649*u*v^2, x^2*y + 2844*x^2*v + 18282*x*y^2 + 29881*x*y*z + 6499*x*y*u + 2367*x*y*v + 23383*x*z*v + 26590*x*u*v + 20977*x*v^2 + 12818*y^3 + 28229*y^2*z + 15643*y^2*u + 27298*y^2*v + 24252*y*z^2 + 25334*y*z*u + 14097*y*z*v + 2082*y*u^2 + 22128*y*u*v + 19545*y*v^2 + 22153*z^2*v + 1137*z*u*v + 10273*z*v^2 + 22348*u^2*v + 17906*u*v^2 + 30776*v^3, x^2*y + 799*x^2*u + 3291*x*y^2 + 12069*x*y*z + 15385*x*y*u + 1353*x*y*v + 10765*x*z*u + 13522*x*u^2 + 6512*x*u*v + 9101*y^3 + 7221*y^2*z + 6933*y^2*u + 3059*y^2*v + 1559*y*z^2 + 1971*y*z*u + 26923*y*z*v + 26152*y*u^2 + 1358*y*u*v + 23597*y*v^2 + 9947*z^2*u + 17766*z*u^2 + 9849*z*u*v + 2274*u^3 + 5458*u^2*v + 10896*u*v^2, x^2*y + 3123*x^2*z + 9872*x*y^2 + 18385*x*y*z + 13015*x*y*u + 10875*x*y*v + 18894*x*z^2 + 21302*x*z*u + 12000*x*z*v + 6834*y^3 + 22815*y^2*z + 20435*y^2*u + 11817*y^2*v + 1555*y*z^2 + 20948*y*z*u + 2160*y*z*v + 10015*y*u^2 + 23814*y*u*v + 5940*y*v^2 + 24225*z^3 + 3577*z^2*u + 15634*z^2*v + 3483*z*u^2 + 4477*z*u*v + 11118*z*v^2, x^3 + 18282*x^2*y + 29881*x^2*z + 6499*x^2*u + 24718*x^2*v + 12818*x*y^2 + 28229*x*y*z + 15643*x*y*u + 24102*x*y*v + 24252*x*z^2 + 25334*x*z*u + 4528*x*z*v + 2082*x*u^2 + 139*x*u*v + 17849*x*v^2 + 15095*y^2*v + 26880*y*z*v + 11767*y*u*v + 31258*y*v^2 + 29583*z^2*v + 10344*z*u*v + 27327*z*v^2 + 26121*u^2*v + 20528*u*v^2 + 10928*v^3]
    F = graded_free_module(Rg,1)
    V1 = [p * F[1] for p in V]
    M = quo(F, V1)[1]
    # Free resolution via Singular gives wrong result
    free_res1 = free_resolution(M)
    @test is_graded(free_res1)
    @test all(iszero, homology(free_res1.C))
    free_res3 = free_resolution_via_kernels(M)
    @test is_graded(free_res3)
    B = betti(free_res3)
    B1 = Dict{Tuple{Int64, Any}, Int64}(
    (0, Z[0]) => 1,
    (1, 3 * Z[1]) => 7,
    (2, 4 * Z[1]) => 10,
    (3, 5 * Z[1]) => 5,
    (4, 6 * Z[1]) => 1)
    @test as_dictionary(B) == B1
    @test all(iszero, homology(free_res3.C))
end

@testset "Tensor module resolution" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    A = Rg[x; y]
    B = Rg[x^2; x*y; y^2; z^4]
    M = SubquoModule(F, A, B)
    free_res = free_resolution_via_kernels(M)
    F1 = graded_free_module(Rg, 1)
    N = SubquoModule(F1, Rg[x*y+2*x^2; x+y], Rg[z^4;])
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
end

@testset "Tensor resolution module" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    A = Rg[x; y]
    B = Rg[x^2; x*y; y^2; z^4]
    M = SubquoModule(F, A, B)
    free_res = free_resolution_via_kernels(M)
    F1 = graded_free_module(Rg, 1)
    N = SubquoModule(F1, Rg[x+2*x^2*z; x+y-z], Rg[z^4;])
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
end

@testset "Hom module resolution" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    A = Rg[x; y]
    B = Rg[x^2; x*y; y^2; z^4]
    M = SubquoModule(F, A, B)
    free_res = free_resolution_via_kernels(M)
    F1 = graded_free_module(Rg, 1)
    N = SubquoModule(F1, Rg[x*y+2*x^2; x+y], Rg[z^4;])
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
    N = SubquoModule(F1, Rg[x*y+2*x^2; x+y], Rg[z^4; x^2-y*z])
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
end

@testset "Hom resolution module" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    A = Rg[x; y]
    B = Rg[x^2; x*y; y^2; z^4]
    M = SubquoModule(F, A, B)
    free_res = free_resolution_via_kernels(M)
    F1 = graded_free_module(Rg, 1)
    N = SubquoModule(F1, Rg[x*y+2*x^2; x+y], Rg[z^4;])
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
end

@testset "Tor and Ext" begin
  Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
  A = Rg[x; y]
  B = Rg[x^2; x*y; y^2; z^4]
  F = graded_free_module(Rg, 1)
  M = SubquoModule(F, A, B)
  Q, _ = quo(F, [x*F[1]])
  G = free_module(Rg, 2)
  M_coker = present_as_cokernel(M)

  T0 = tor(Q, M, 0)
  T1 = tor(Q, M, 1)
  T2 =  tor(Q, M, 2)
  @test is_canonically_isomorphic(T0, M)
  @test ngens(present_as_cokernel(T1)) == ngens(M_coker)
  # Todo twist
  @test iszero(T2)
  T0 = tor(M, Q, 0)
  T1 = tor(M, Q, 1)
  T2 = tor(M, Q, 2)
  @test is_canonically_isomorphic(present_as_cokernel(T0), M_coker)
  # todo simplify
  @test iszero(T2)

  E0 = ext(Q, M, 0)
  E1 = ext(Q, M, 1)
  E2 = ext(Q, M, 2)
  @test is_canonically_isomorphic(present_as_cokernel(E0), M_coker)
  @test ngens(present_as_cokernel(E1)) == ngens(M_coker)
  @test iszero(E2)
  E0 = ext(M, Q, 0)
  E1 = ext(M, Q, 1)
  E2 = ext(M, Q, 2)
  E3 = ext(M, Q, 3)
  E4 = ext(M, Q, 4)
  @test iszero(E0)
  @test iszero(E1)
  # Todo simplify
  @test ngens(E3) == ngens(M_coker)
  # Todo twist
  @test iszero(E4)
end

@testset "Groebner bases graded" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    J = SubquoModule(F, [x*F[1], (x^2)*F[1], (x+y)*F[1]])
    @test leading_module(J) == SubquoModule(F, [x*F[1], y*F[1]])
    J = SubquoModule(F, [(x*y^2+x*y^2)*F[1], (x^2*y+x^2*z-y^3)*F[1]])
    @test leading_module(J) == SubquoModule(F, [x*y^2*F[1], x^2*y*F[1], y^4*F[1], x^3*z^2*F[1]])
    F = graded_free_module(Rg, 2)
    lp = lex(gens(base_ring(F)))*lex(gens(F))
    M = SubquoModule(F, [(x^2*y^2*F[1]+y*z^3*F[2]), x*z*F[1]+z^2*F[2]])
    @test leading_module(M,lp) == SubquoModule(F, [x*z*F[1], x^2*y^2*F[1], x*y^2*z^2*F[2]])
    Rg, x = graded_polynomial_ring(QQ, "x#" => 1:4)
    F = graded_free_module(Rg, 1)
    lp = lex(gens(base_ring(F)))*lex(gens(F))
    J = SubquoModule(F, [(x[1]+x[2])*F[1], (x[1]+x[2]+2*x[3]+2*x[4])*F[1],(x[1]+x[2]+x[3]+x[4])*F[1]])
    @test reduced_groebner_basis(J, lp).O == Oscar.ModuleGens([(x[3]+x[4])*F[1], (x[1]+x[2])*F[1]], F).O
    @test haskey(J.groebner_basis, lp)
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F = graded_free_module(Rg, 1)
    lp = lex(gens(base_ring(F)))*lex(gens(F))
    I = SubquoModule(F, [(x-y)*F[1], (y^2-z^2)*F[1]])
    f = (x*y^2+y*z^2)*F[1]
    @test Oscar.reduce(f, I) == 2*y*z^2*F[1]
    F = graded_free_module(Rg, 2)
    A = Rg[(x+z)^2 y*z+x^2; (y*x+2*z^2) z^2]
    B = Rg[2*z*(x*z+y^2) x*z^2]
    M = SubquoModule(F, A, B)
    gb = groebner_basis(M)
    P = sum(Oscar.SubModuleOfFreeModule(F, gb), Oscar.SubModuleOfFreeModule(F, gb.quo_GB))
    Q = Oscar.SubModuleOfFreeModule(F, groebner_basis(M.sum))
    @test P == Q
    v = x*((x^2 + 2*x*z + z^2)*F[1] + (x^2 + y*z)*F[2]) + (2*x*z^2 + 2*y^2*z)*F[1] + x*z^2*F[2]
    @test represents_element(v, M)
end

@testset "Tensor product of morphisms graded" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    F2 = graded_free_module(Rg,2)
    F3 = graded_free_module(Rg,3)
    F4 = graded_free_module(Rg,4)
    A1 = Rg[6*x*y^2*z^2 + 4*y*z^4 + 12*y^2*z^3 6*x^2*y^2*z + 11*x^2*z^3; x y;  y*z   3*x*z + 3*x*y + 5*y*z]
    B1 = Rg[12*z^2   x*y]
    A2 = Rg[8*x^2*y*z         6*x^2*z^2       4*y^4;
            9*x*z^5   3*x^2*y^2*z^2    11*y*z^5;
            8*x*y*z            15*y^3   9*x*z^2]
    B2 = Rg[0   14*x*y   6*x^2]
    M1 = SubquoModule(F2,A1,B1)
    M2 = SubquoModule(F3,A2,B2)
    M, pure_M = tensor_product(M1,M2, task=:map)
    @test is_graded(M)
    phi = hom_tensor(M, M, [identity_map(M1), identity_map(M2)])
    @test is_homogeneous(phi)
    v = M[1]
    @test phi(v) == v
    A3 = Rg[7*x^2*y^2   11*x^2*y*z;
            6*x^2*y        4*x^2*z]
    M3 = SubquoModule(Oscar.SubModuleOfFreeModule(F2, A3))
    N,pure_N = tensor_product(M3,F4, task=:map)
    @test is_graded(N)
    M3_to_M1 = SubQuoHom(M3, M1, [M1[1], x^3*M1[2]])
    @test is_welldefined(M3_to_M1)
    F4_to_M2 = FreeModuleHom(F4,M2, [x^2*M2[1], M2[2], y^3*M2[3], z^3*M2[3]])
    phi = hom_tensor(N,M,[M3_to_M1,F4_to_M2])
    u1 = M3[1]
    u2 = F4[1]
    @test phi(pure_N((u1,u2))) == pure_M((M3_to_M1(u1),F4_to_M2(u2)))
end

@testset "direct product graded" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])
    Z = grading_group(Rg)
    F2 = graded_free_module(Rg, 2)
    F3 = graded_free_module(Rg, 3)
    A1 = Rg[6*x*y^2*z^2 + 4*y*z^4 + 12*y^2*z^3 6*x^2*y^2*z + 11*x^2*z^3; x y;  y*z   3*x*z + 3*x*y + 5*y*z]
    B1 = Rg[12*z^2   x*y]
    A2 = Rg[8*x^2*y*z         6*x^2*z^2       4*y^4;
            9*x*z^5   3*x^2*y^2*z^2    11*y*z^5;
            8*x*y*z            15*y^3   9*x*z^2]
    B2 = Rg[0   14*x*y   6*x^2]
    M1 = SubquoModule(F2, A1, B1)
    M2 = SubquoModule(F3, A2, B2)
    sum_M, emb = direct_sum(M1, M2)
    @test is_graded(sum_M)
    @test is_homogeneous(emb[1])
    @test is_homogeneous(emb[2])
    @test domain(emb[1]) === M1
    @test domain(emb[2]) === M2
    @test codomain(emb[1]) === sum_M
    @test codomain(emb[2]) === sum_M
    sum_M, pr = direct_sum(M1,M2, task=:prod)
    @test is_graded(sum_M)
    @test is_homogeneous(pr[1])
    @test is_homogeneous(pr[2])
    @test codomain(pr[1]) === M1
    @test codomain(pr[2]) === M2
    @test domain(pr[1]) === sum_M
    @test domain(pr[2]) === sum_M
    prod_M, emb, pr = direct_sum(M1,M2,task=:both)
    @test is_graded(sum_M)
    @test is_homogeneous(emb[1])
    @test is_homogeneous(emb[2])
    @test is_homogeneous(pr[1])
    @test is_homogeneous(pr[2])
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
    A1 = Rg[4*x*y^2*z^2 + 6*x*z^4               9*x^2*z^3;
              5*x^2*y*z + 12*x*y^3   8*x^2*y^2 + z^4;
              11*x^2*z^2          8*y^4 + 4*y*z^3]
    B1 = Rg[10*x^2*y*z    15*y^2*z^2 + 3*y*z^3]
    N1 = SubquoModule(F2,A1,B1)
    A2 = Rg[   14*x^4   4*x*y^2*z   8*x^2*y^2;
       2*x*y^2*z   4*x^2*z^2        15*y^4]
    B2 = Rg[12*x*y*z^2   13*x^2*z^2   11*x*z^3]
    N2 = SubquoModule(F3,A2,B2)
    prod_N = direct_product(N1,N2,task=:none)
    @test is_graded(prod_N)
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

    M1_to_N1 = SubQuoHom(M1,N1,zero_matrix(Rg,3,3))
    @test is_homogeneous(M1_to_N1)
    H12 = hom(M1,N2)[1]
    @test is_graded(H12)
    H21 = hom(M2,N1)[1]
    @test is_graded(H21)
    M1_to_N2 = iszero(H12) ? SubQuoHom(M1,N2,zero_matrix(R,3,2)) : element_to_homomorphism(H12[1])
    M2_to_N1 = iszero(H21) ? SubQuoHom(M2,N1,zero_matrix(R,2,3)) : element_to_homomorphism(x^3*H21[1])
    M2_to_N2 = SubQuoHom(M2, N2, [0*N2[1],0*N2[1],0*N2[1]])
    @test degree(M1_to_N1) == Z[0]
    @test degree(M1_to_N2) == 6*Z[1]
    @test degree(M2_to_N1) == 6*Z[1]
    @test degree(M2_to_N2) == Z[0]
    @test is_welldefined(M1_to_N1)
    @test is_welldefined(M1_to_N2)
    @test is_welldefined(M2_to_N1)
    @test is_welldefined(M2_to_N2)

    phi = hom_product(prod_M,prod_N,[M1_to_N1 M1_to_N2; M2_to_N1 M2_to_N2])
    @test degree(phi) == 6*Z[1]
    for g in gens(M1)
        @test M1_to_N1(g) == canonical_projection(prod_N,1)(phi(emb[1](g)))
        @test M1_to_N2(g) == canonical_projection(prod_N,2)(phi(emb[1](g)))
    end
    for g in gens(M2)
        @test M2_to_N1(g) == canonical_projection(prod_N,1)(phi(emb[2](g)))
        @test M2_to_N2(g) == canonical_projection(prod_N,2)(phi(emb[2](g)))
    end
    prod_FN,prod,emb = direct_product(F2,N2,task=:both)
    @test is_graded(prod_FN)
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

@testset "Coordinates" begin
    Z3, a = finite_field(3,1,"a")
    Rg, (x,y) = graded_polynomial_ring(Z3, [:x, :y])
    coeffs = [Z3(i) for i=0:1]
    A = Rg[x*y x^2+y^2; y^2 x*y; x^2+y^2 y^2]
    B = Rg[2*x^2 (x+y)^2; x^2+y^2 x^2+2*x*y]
    F = graded_free_module(Rg, 2)
    M = SubquoModule(F, A, B)

    monomials = [x,y]
    coeff_generator = ([c1,c2] for c1 in coeffs for c2 in coeffs)
    for coefficients in ([c1,c2,c3] for c1 in coeff_generator for c2 in coeff_generator for c3 in coeff_generator)
      v = sparse_row(Rg, [(i,sum(coefficients[i][j]*monomials[j] for j=1:2)) for i=1:3])
      v_as_FreeModElem = sum([v[i]*repres(M[i]) for i=1:ngens(M)])
      elem1 = SubquoModuleElem(v_as_FreeModElem,M)
      elem2 = SubquoModuleElem(v,M)
      @test elem1 == elem2
    end
end


@testset "Presentations 3" begin
    Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights=[-1,1,-1])
    F = graded_free_module(Rg, 1)
    B = Rg[x^2; y^3; z^4]
    M, inc = sub(F,B)
    M = cokernel(inc)
    fr = free_resolution(M)
    phi = map(fr,4)
    N = cokernel(phi)
    f = present_as_cokernel(N)
    @test ngens(f) == 1
end

######################################################################
# Tests filtrated modules
######################################################################

"""
    randpoly(R::Ring,coeffs=0:9,max_exp=4,max_terms=8)

> Return a random Polynomial from the Polynomial Ring `R` with coefficients in `coeffs`
> with exponents between `0` and `max_exp` und between `0` and `max_terms` terms
"""
function randpoly(R::Oscar.Ring,coeffs=0:9,max_exp=4,max_terms=8)
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

@testset "Basic degree and homogeneity tests for free modules" begin
    R, x = graded_polynomial_ring(QQ, 3)
    F = FreeMod_dec(R, 3) #standard grading

    @test degree(F[1]) == decoration(R)[0]
    @test degree(x[2]*x[3]^2*F[1]+x[1]^3*F[3]) == 3*decoration(R)[1]
    @test !is_homogeneous((x[1]+x[2]^2)*F[1])
    @test !is_homogeneous(F[1]+x[3]*F[3])
    @test is_homogeneous(x[1]*x[2]*F[1]+x[3]^2*F[2])
    @test homogeneous_component(x[1]*F[1]+x[2]^2*F[2]+x[3]*F[3], decoration(R)[1]) == x[1]*F[1]+x[3]*F[3]
    D = homogeneous_components(x[1]*F[1]+x[2]^2*F[2]+x[3]*F[3])
    @test length(D) == 2
    @test D[2*decoration(R)[1]] == x[2]^2*F[2]

    g = abelian_group(2,0)
    R, x = graded_polynomial_ring(QQ, 2; weights=[g[1], g[2]])
    F = FreeMod_dec(R, [g[2], g[1], g[1]-2*g[2]])

    @test degree(F[1]) == g[2]
    @test degree(x[1]*F[1]+x[2]^3*F[3]) == g[1]+g[2]
    @test !is_homogeneous((x[1]+x[2]^2)*F[1])
    @test !is_homogeneous(F[1]+x[2]*F[3])
    @test is_homogeneous(x[1]*x[2]*F[1]+x[2]^2*F[2])
    @test homogeneous_component(x[1]*F[1]+x[2]^2*F[2]+x[2]*F[3], g[1]+g[2]) == x[1]*F[1]
    D = homogeneous_components(x[1]^2*F[1]+x[2]^2*F[2]+x[1]*x[2]^3*F[3])
    @test length(D) == 2
    @test D[g[1]+2*g[2]] == x[2]^2*F[2]
    @test_throws ErrorException degree(x[1]*F[1]+x[2]^2*F[2]+x[2]*F[3])

    # Filtered case
    Z = abelian_group(0)
    Qx, x = polynomial_ring(QQ, 2)
    R, x = filtrate(Qx, [Z[1], Z[1]], (x,y) -> x[1] > y[1])
    F = FreeMod_dec(R, 2)
    @test !is_homogeneous(x[1]*F[1] + x[2]^3*F[2])
    @test degree(x[1]*F[1] + x[2]^3*F[2]) == Z[1]
end

@testset "Tensor product of decorated free modules" begin
    Z = abelian_group(0)
    R, x = graded_polynomial_ring(QQ, 3; weights=[Z[1], 5*Z[1], -Z[1]])
    F2 = FreeMod_dec(R, [Z[0], Z[1]])
    F3 = FreeMod_dec(R, [Z[1], -2*Z[1], Z[0]])

    F6, pure = tensor_product(F2,F3, task=:map)
    @test is_homogeneous(F6[1])
    v = pure((F2[1], F3[2]))
    @test degree(v) == -2*Z[1]
    v = pure(((x[1]^5+x[2])*F2[2], x[3]*F3[3]))
    @test degree(v) == 5*Z[1]
    for v in gens(F6)
        a,b = inv(pure)(v)
        @test degree(v) == degree(a) + degree(b)
    end
end

@testset "Hom module for decorated free modules" begin
    R, x = graded_polynomial_ring(QQ, 3, :x)
    Z = grading_group(R)
    F = FreeMod_dec(R, [Z[1], 2*Z[1]])
    G = FreeMod_dec(R, [Z[1], 2*Z[1], 3*Z[1]])

    H,_ = hom(G,F)
    for v in gens(H)
        @test v == homomorphism_to_element(H, element_to_homomorphism(v))
    end

    phi = H[1]
    v = x[1]*F[1] + F[2]
    @test degree(phi) == degree(element_to_homomorphism(phi)(v)) - degree(v)
end

Oscar.@_AuxDocTest "minimal Betti tables", (fix = false),
raw"""
# The Veronese surface in IP^4 due to Wolfram

```jldoctest; setup = :(using Oscar)
julia> F = GF(31991)
Prime field of characteristic 31991

julia> R, (x,y,z,u,v) = graded_polynomial_ring(F, [:x, :y, :z, :u, :v]);

julia> I = ideal([F(45)/16*x^3-8565*x^2*y-7937*x*y^2+14060*x^2*z+F(53)/113*x*y*z-F(125)/17*x*z^2-15712*x^2*u-5990*x*y*u-F(71)/88*x*z*u-6141*x*u^2+F(49)/104*x^2*v-194*x*y*v+F(63)/103*y^2*v+F(74)/103*x*z*v+11618*y*z*v+F(41)/111*z^2*v+14387*x*u*v-F(64)/9*y*u*v+13097*z*u*v+F(84)/101*u^2*v+12211*x*v^2+F(85)/63*y*v^2-F(119)/100*z*v^2-6247*u*v^2-F(74)/51*v^3,-F(45)/16*x^2*y+8565*x*y^2+7937*y^3-14060*x*y*z-F(53)/113*y^2*z+F(125)/17*y*z^2+15712*x*y*u+5990*y^2*u+F(71)/88*y*z*u+6141*y*u^2-x^2*v+3340*x*y*v+F(111)/40*y^2*v+F(22)/37*x*z*v-F(110)/29*y*z*v-F(71)/25*z^2*v-F(71)/87*x*u*v+F(112)/55*y*u*v+F(103)/80*z*u*v-F(100)/7*u^2*v-5013*x*v^2+15010*y*v^2+F(52)/51*z*v^2+F(126)/5*u*v^2-6580*v^3,-F(45)/16*x^2*z+8565*x*y*z+7937*y^2*z-14060*x*z^2-F(53)/113*y*z^2+F(125)/17*z^3+15712*x*z*u+5990*y*z*u+F(71)/88*z^2*u+6141*z*u^2+F(41)/79*x^2*v+F(121)/48*x*y*v+9621*y^2*v+F(105)/83*x*z*v+F(27)/26*y*z*v-649*z^2*v-12278*x*u*v+14655*y*u*v+10594*z*u*v+10462*u^2*v+F(76)/13*x*v^2+F(97)/50*y*v^2-7761*z*v^2+F(52)/31*u*v^2-6636*v^3,-F(45)/16*x^2*u+8565*x*y*u+7937*y^2*u-14060*x*z*u-F(53)/113*y*z*u+F(125)/17*z^2*u+15712*x*u^2+5990*y*u^2+F(71)/88*z*u^2+6141*u^3-F(40)/31*x^2*v+13297*x*y*v+F(126)/113*y^2*v-F(51)/100*x*z*v-F(97)/38*y*z*v+10372*z^2*v-F(4)/31*x*u*v+F(32)/9*y*u*v-F(1)/7*z*u*v+15473*u^2*v+F(101)/36*x*v^2+F(97)/109*y*v^2-14100*z*v^2+F(109)/32*u*v^2-F(5)/111*v^3,-F(40)/31*x^3+13297*x^2*y+F(126)/113*x*y^2-F(51)/100*x^2*z-F(97)/38*x*y*z+10372*x*z^2+F(26)/105*x^2*u-3745*x*y*u+F(63)/103*y^2*u-F(98)/61*x*z*u+11618*y*z*u+F(41)/111*z^2*u+F(26)/15*x*u^2-F(64)/9*y*u^2+13097*z*u^2+F(84)/101*u^3+F(101)/36*x^2*v+F(97)/109*x*y*v-14100*x*z*v-14778*x*u*v+F(85)/63*y*u*v-F(119)/100*z*u*v-6247*u^2*v-F(5)/111*x*v^2-F(74)/51*u*v^2,F(40)/31*x^2*y-13297*x*y^2-F(126)/113*y^3+F(51)/100*x*y*z+F(97)/38*y^2*z-10372*y*z^2-x^2*u+F(103)/30*x*y*u+2754*y^2*u+F(22)/37*x*z*u+F(126)/95*y*z*u-F(71)/25*z^2*u-F(71)/87*x*u^2-F(109)/7*y*u^2+F(103)/80*z*u^2-F(100)/7*u^3-F(101)/36*x*y*v-F(97)/109*y^2*v+14100*y*z*v-5013*x*u*v+10008*y*u*v+F(52)/51*z*u*v+F(126)/5*u^2*v+F(5)/111*y*v^2-6580*u*v^2,F(40)/31*x^2*z-13297*x*y*z-F(126)/113*y^2*z+F(51)/100*x*z^2+F(97)/38*y*z^2-10372*z^3+F(41)/79*x^2*u+F(121)/48*x*y*u+9621*y^2*u+5671*x*z*u-F(42)/71*y*z*u-5219*z^2*u-12278*x*u^2+14655*y*u^2+F(58)/59*z*u^2+10462*u^3-F(101)/36*x*z*v-F(97)/109*y*z*v+14100*z^2*v+F(76)/13*x*u*v+F(97)/50*y*u*v-12763*z*u*v+F(52)/31*u^2*v+F(5)/111*z*v^2-6636*u*v^2,F(41)/79*x^3+F(121)/48*x^2*y+9621*x*y^2-F(22)/109*x^2*z+F(25)/19*x*y*z+F(63)/103*y^2*z+F(23)/45*x*z^2+11618*y*z^2+F(41)/111*z^3-12278*x^2*u+14655*x*y*u+F(126)/73*x*z*u-F(64)/9*y*z*u+13097*z^2*u+10462*x*u^2+F(84)/101*z*u^2+F(76)/13*x^2*v+F(97)/50*x*y*v-F(106)/115*x*z*v+F(85)/63*y*z*v-F(119)/100*z^2*v+F(52)/31*x*u*v-6247*z*u*v-6636*x*v^2-F(74)/51*z*v^2,-F(41)/79*x^2*y-F(121)/48*x*y^2-9621*y^3-x^2*z-F(22)/89*x*y*z-F(32)/17*y^2*z+F(22)/37*x*z^2-F(86)/111*y*z^2-F(71)/25*z^3+12278*x*y*u-14655*y^2*u-F(71)/87*x*z*u+F(74)/47*y*z*u+F(103)/80*z^2*u-10462*y*u^2-F(100)/7*z*u^2-F(76)/13*x*y*v-F(97)/50*y^2*v-5013*x*z*v-9220*y*z*v+F(52)/51*z^2*v-F(52)/31*y*u*v+F(126)/5*z*u*v+6636*y*v^2-6580*z*v^2,x^3+11117*x^2*y+991*x*y^2-F(63)/103*y^3-F(22)/37*x^2*z-F(35)/13*x*y*z-11618*y^2*z+F(71)/25*x*z^2-F(41)/111*y*z^2+F(71)/87*x^2*u+F(68)/115*x*y*u+F(64)/9*y^2*u-F(103)/80*x*z*u-13097*y*z*u+F(100)/7*x*u^2-F(84)/101*y*u^2+5013*x^2*v-F(67)/114*x*y*v-F(85)/63*y^2*v-F(52)/51*x*z*v+F(119)/100*y*z*v-F(126)/5*x*u*v+6247*y*u*v+6580*x*v^2+F(74)/51*y*v^2]);

julia> I = ideal(R, minimal_generating_set(I));

julia> F = graded_free_module(R, 1);

julia> sub_F, inc = sub(F, [g*F[1] for g in gens(I)]);

julia> M = cokernel(inc);

julia> A, _ = quo(R, I);

julia> minimal_betti_table(A)
degree: 0  1   2  3  4
----------------------
     0: 1  -   -  -  -
     1: -  -   -  -  -
     2: -  7  10  5  1
----------------------
 total: 1  7  10  5  1

julia> minimal_betti_table(M)
degree: 0  1   2  3  4
----------------------
     0: 1  -   -  -  -
     1: -  -   -  -  -
     2: -  7  10  5  1
----------------------
 total: 1  7  10  5  1

julia> minimal_betti_table(I)
degree: 0   1  2  3
-------------------
     3: 7  10  5  1
-------------------
 total: 7  10  5  1
```

# small example due to Janko

```jldoctest; setup = :(using Oscar)
julia> R, x = graded_polynomial_ring(QQ, :x => 1:7);

julia> I = ideal(R, [x[1]*x[2]*x[5], x[1]*x[2]*x[6], x[3]*x[4]*x[6], x[3]*x[4]*x[7], x[5]*x[7]]);

julia> A, _ = quo(R, I);

julia> minimal_betti_table(A)
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  1  -  -
     2: -  4  4  -
     3: -  -  1  -
     4: -  -  -  1
------------------
 total: 1  5  5  1
```

# another example due to Wolfram

```jldoctest; setup = :(using Oscar)
julia> R, (x, y, z, w) = graded_polynomial_ring(QQ, [:x, :y, :z, :w]);

julia> I = ideal(R, [w^2 - x*z, w*x - y*z, x^2 - w*y, x*y - z^2, y^2 - w*z]);

julia> A, _ = quo(R, I);

julia> minimal_betti_table(free_resolution(A))
degree: 0  1  2  3
------------------
     0: 1  -  -  -
     1: -  5  5  -
     2: -  -  -  1
------------------
 total: 1  5  5  1
```
"""

@testset "sheaf cohomology" begin
  S, _ = graded_polynomial_ring(QQ, :x => 1:4)
  I = ideal(S, gens(S))
  FI = free_resolution(I)
  M = cokernel(map(FI, 2))
  tbl = Oscar._sheaf_cohomology_bgg(M, -6, 2)
  lbt = sheaf_cohomology(M, -6, 2, algorithm = :bgg)
  @test tbl.values == lbt.values
  @test tbl[3, -6] == 70
  @test tbl[1, 0] == 1
  @test iszero(tbl[2, -2])

  F = free_module(S, 1)
  @test_throws AssertionError Oscar._sheaf_cohomology_bgg(F, -6, 2)

  S, _ = graded_polynomial_ring(QQ, :x => 1:4; weights=[1,2,3,4])
  F = graded_free_module(S, 1)
  @test_throws AssertionError Oscar._sheaf_cohomology_bgg(F, -6, 2)

  S, _ = graded_polynomial_ring(QQ, :x => 1:5)
  F = graded_free_module(S, 1)
  tbl = sheaf_cohomology(F, -7, 2)
  a = tbl.values
  b = transpose(a) * a
  @test is_symmetric(b)
end

@testset "twist" begin
  R, (x, y) = graded_polynomial_ring(QQ, [:x, :y], [1 0; 0 1]);
  I = ideal(R, [x, y])
  M = quotient_ring_as_module(I)
  N = twist(M, [1, 2])
  @test Int(degree(gen(N, 1))[1]) == -1
  @test Int(degree(gen(N, 1))[2]) == -2
end

@testset "all_monomials for graded free modules" begin
  R, (x, y, u, v, w) = QQ[:x, :y, :u, :v, :w]
  S, (x, y, u, v, w) = grade(R)

  F = graded_free_module(S, [-1, 2])

  for d in -4:4
    amm = Oscar.all_monomials(F, d)
    length(amm)
    eltype(amm)
    @test d < -1 || !isempty(amm)
    @test all(x->degree(x) == grading_group(F)([d]), amm)
  end
end

##################################################################
# Tests random elements
##################################################################

@testset "random free module hom" begin
    Rg, (x, y, z) = graded_polynomial_ring(GF(101), [:x, :y, :z])
    Z = grading_group(Rg)
    F1 = graded_free_module(Rg, [1,2,2])
    F2 = graded_free_module(Rg, [3,5])
    V, f = hom(F1, F2)
    ff = rand_homogeneous(Rg,8)
    @test is_homogeneous(ff)
    @test degree(ff) == 8*Z[1]
    v = rand_homogeneous(V, 8)
    @test is_homogeneous(v)
    @test degree(v) == 8*Z[1]
end

@testset "random subquo module hom" begin
    Rg, (x, y) = graded_polynomial_ring(GF(101), [:x, :y])
    Z = grading_group(Rg)
    F = graded_free_module(Rg, [3,5]);
    V = [x*F[1], y^2*F[2]];
    M = quo(F, V)[1]
    F2 = graded_free_module(Rg, [2,2]);
    H, f = hom(F2, M)
    v = rand_homogeneous(H, 8)
    @test is_homogeneous(v)
    @test degree(v) == 8*Z[1]
end

#= Disabled for the moment but continued soon.
@testset "monomials of subquos" begin
  S, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z])

  S1 = graded_free_module(S, [0])
  I = ideal(S, [u^2 for u in gens(S)])
  IS1, inc = I*S1
  M = cokernel(inc)

  a = Oscar.AllSubquoMonomials(M, 3)
  b = Oscar.all_exponents(M, 3)

  v = collect(a)
  @test length(v) == length(a) == 1 == length(collect(b))
  @test !(M(x^3 * S1[1]) in v)
  @test M(x*y*z * S1[1]) in v

  I = ideal(S, [u^3 for u in gens(S)])
  IS1, inc = I*S1
  M = cokernel(inc)

  a = Oscar.AllSubquoMonomials(M, 3)
  b = Oscar.all_exponents(M, 3)

  v = collect(a)
  @test length(v) == length(a) == length(collect(b))
  @test !(M(x^3 * S1[1]) in v)
  @test M(x*y*z * S1[1]) in v
  @test M(x^2*y * S1[1]) in v

  J, _ = sub(S1, [x*y*z*S1[1]])
  I = ideal(S, [u^4 for u in gens(S)])
  IS1, inc = I*S1
  M, _ = quo(J, IS1)
  a = Oscar.AllSubquoMonomials(M, 4)
  b = Oscar.all_exponents(M, 4)
  @test length(collect(a)) == length(a) == 3 == length(collect(b))

  a = Oscar.AllSubquoMonomials(M, 6)
  b = Oscar.all_exponents(M, 6)
  @test length(collect(a)) == length(a) == 7 == length(collect(b))
end
=#

@testset "simplification of graded subquos, issue #3108" begin
  X = rational_d10_pi9_quart_1();
  I = defining_ideal(X);
  Pn = base_ring(I)
  n = ngens(Pn)-1
  c = codim(I)
  FI = free_resolution(I)
  FIC = FI.C;
  r = range(FIC)
  C = shift(FIC[first(r):-1:1], -c)
  F = free_module(Pn, 1)
  OmegaPn = grade(F, [n+1])
  D = hom(C, OmegaPn)
  Omega = homology(D, 0);
  is_graded(Omega)
  SOmega, b = simplify(Omega)
  a = get_attribute(b, :inverse)
  @test is_graded(SOmega)
  @test is_isomorphism(a)
  @test is_isomorphism(b)
  M, iso = forget_grading(Omega)
  @test is_isomorphism(iso)
  inv_iso = get_attribute(iso, :inverse)
  @test is_isomorphism(inv_iso)
  M, _ = forget_grading(Omega)
  prune_with_map(M)
end

@testset "kernels FreeMod -> Subquo with gradings" begin
  S, (x, y) = graded_polynomial_ring(QQ, [:x, :y])
  G = grading_group(S)
  F0 = graded_free_module(S, [zero(G)])
  F1 = graded_free_module(S, [-1])
  M, _ = sub(F1, [F1[1]])
  phi = hom(F0, M, [x*M[1]])
  @test is_graded(phi)
  K, inc = kernel(phi)
end

@testset "ideal and quotient ring as module" begin
  R, (x, y) = polynomial_ring(QQ,[:x, :y])
  I = ideal(R, [y^2-x^3-x^2])
  II = ideal_as_module(I)
  M = quotient_ring_as_module(I)
  Q, pr = quo(R, I)
  J = ideal(Q, gens(Q))
  ideal_as_module(J)
  quotient_ring_as_module(J)
  U = complement_of_point_ideal(R, [0, 0])
  Rloc, _ = localization(R, U)
  Iloc = ideal(Rloc, [y^2-x^3-x^2])
  Mloc = ideal_as_module(Iloc)
  @test is_empty(relations(Mloc))
  QI = quotient_ring_as_module(I)
  @test rank(ambient_free_module(QI)) == 1
  Q, pr = quo(Rloc, Iloc)
  J = ideal(Q, gens(Q))
  ideal_as_module(J)
  quotient_ring_as_module(J)
end
