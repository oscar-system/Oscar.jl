@testset "algebra homomorphisms" begin
  r, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"])
  s, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"])
  S = quo(s, ideal(s, [c-b^3]))[1]
  T, (X, Y, Z) = grade(r)
  U = [X, Y, Z]
  V = S.([2*a+b^6, 7*b-a^2, c^2])
  f = hom(r, S, V)
  g1 = hom(r, T, U)
  g2 = hom(T, T, U)
  K = kernel(f)
  R = quo(r, K)[1]
  phi = hom(R, S, V)
  psi = inverse(phi)
  #W = psi.image
  W = Oscar._images(psi)
  prphi = preimage(phi, S(a))
  prpsi = preimage(psi, R(x))

  @test is_surjective(f) == true
  @test is_injective(f) == false
  @test is_bijective(f) == false
  @test is_surjective(phi) == true
  @test is_injective(phi) == true
  @test is_bijective(phi) == true
  @test prphi == R(1//2 * x-1//2*z)
  @test prpsi == S(2*a+c^2)
  @test W[1] == R(1//2 * x-1//2*z)
  @test W[2] == R(1//28*x^2 - 1//14*x*z + 1//7*y + 1//28*z^2)
  @test W[3] == R(1//21952*x^6 - 3//10976*x^5*z + 3//5488*x^4*y + 15//21952*x^4*z^2 - 3//1372*x^3*y*z - 5//5488*x^3*z^3 + 3//1372*x^2*y^2 + 9//2744*x^2*y*z^2 + 15//21952*x^2*z^4 - 3//686*x*y^2*z - 3//1372*x*y*z^3 - 3//10976*x*z^5 + 1//343*y^3 + 3//1372*y^2*z^2 + 3//5488*y*z^4 + 1//21952*z^6)
end

@testset "finiteness tests for alghoms" begin
  r, (X, Y, Z) = polynomial_ring(QQ, ["X", "Y", "Z"])
  s1, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"])
  S1 = quo(s1, ideal(s1, [c-b^3]))[1]
  V = S1.([2*a+b^6, 7*b-a^2, c^2])
  f1 = hom(r, S1, V)

  S2 = quo(r, ideal(r, [X*Y]))[1]
  a, b, c = S2(X), S2(Y), S2(Z)
  f2 = hom(r, S2, [(a*b)^3+a^2+c, b^2-1, c^3])

  r2 , (x, y) = polynomial_ring(QQ, ["x", "y"])
  s3 = quo(r, ideal(r, [X^2-Y^2*Z]))[1]
  f3 = hom(r2, s3, [gen(s3, 1), gen(s3, 2)])

  @test isfinite(f1) == true
  @test isfinite(f2) == true
  @test isfinite(f3) == false
end

@testset "subalgebra membership" begin
  s, (a, b, c) = polynomial_ring(QQ, ["a", "b", "c"])
  S = quo(s, ideal(s, [c-b^3]))[1]
  t = subalgebra_membership(S(c+a^2-b^6), S.([a,b^3]))
  T = parent(t[2])
  @test t[1] == true
  @test t[2] == gen(T, 1)^2-gen(T, 2)^2+gen(T, 2)

  R, (x, y, z) = GradedPolynomialRing(QQ, [ "x", "y", "z" ], [ 3, 1, 3 ])
  f = x^2 - y^6 + z^2
  v = [ x, y^3, z - y^3 ]
  fl, t = subalgebra_membership_homogeneous(f, v)
  @test fl
  @test t(v...) == f

  I = ideal(R, [ z - y^3 ])
  Q, RtoQ = quo(R, I)
  f = Q(x)^2 + Q(y)^6 + Q(z)^2
  v = [ Q(x), Q(y)^3 ]
  fl, t = subalgebra_membership_homogeneous(f, v)
  @test fl
  @test t(v...) == f
end

@testset "#655" begin
  R, vars = QQ["x","y"]
  x = vars[1]
  y = vars[2]
  f = hom(R, R, vars)
  push!(vars, x)
  @test f(x) == x
end
