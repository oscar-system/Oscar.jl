@testset "Modules over MPolyQuos" begin
  P, (x, y) = QQ["x", "y"]
  I = ideal(P, x)
  R, _ = quo(P, I)
  FR = FreeMod(R, 2, cached=false)
  M = oscar._as_poly_module(FR)
  @test M isa SubQuo{<:MPolyElem}
  @test FR isa FreeMod{<:MPolyQuoElem}
  FP = oscar._poly_module(FR)
  @test FP isa FreeMod{<:MPolyElem}

  M = FreeMod(R, 1, cached=false)
  psi = hom(FR, M, [2*M[1], M[1]])

  K, inc = kernel(psi)

  v = FR[1] - 2*FR[2] 
  @test v in K

  N = oscar.SubModuleOfFreeModule(FR, [R(y-1)*FR[1], R(y)*FR[2]])
  c = coordinates(y*(y-1-x)*FR[1] + y^5*(1+x)*FR[2], N)
  @test c[1] == R(y)
  @test c[2] == R(y^4)
end

@testset "Issues in #1806 part 1" begin
    R, (x, y, z) = PolynomialRing(QQ, ["x$i" for i in 1:3])

    M = R[x y; y-1 z]
    f = det(M)
    I = ideal(R, f)
    A, _ = quo(R, I)
    A2 = FreeMod(A, 2) 
    MA = change_base_ring(A, M)
    phi = FreeModuleHom(A2, A2, MA)
    K, inc = kernel(phi)
    @test iszero(kernel(inc)[1])

    psi = FreeModuleHom(A2, A2, ambient_representatives_generators(K))
    K2, inc2 = kernel(psi)
    @test iszero(kernel(inc2)[1])

    phi2 = hom(A2, A2, ambient_representatives_generators(K2))

    K3, inc3 = kernel(phi2)

    @test !(K3 == K2)
end

@testset "Issues in #1806 part 2" begin
    R, (x,y) = QQ["x", "y"]
    Q, phi = quo(R, ideal(R, x))
    F = FreeMod(R, 1)
    FF = FreeMod(Q, 1)
    f = hom(F, FF, gens(FF), phi)
    @test matrix(f) isa MatrixElem
end

@testset "Issues in #1806 part 3" begin
    R, (x,y) = QQ["x", "y"]
    Q, phi = quo(R, ideal(R, x^2))
    F = FreeMod(R, 1)
    FQ = FreeMod(Q, 1)
    f = hom(F, FQ, gens(FQ), phi)
    W, (t,) = PolynomialRing(QQ, ["t"])
    psi = hom(Q, W, [zero(W), W[1]])
    FW = FreeMod(W, 1)
    g = hom(FQ, FW, gens(FW), psi)
    h = compose(f, g) 
    @test h(F[1]) == FW[1]
end
