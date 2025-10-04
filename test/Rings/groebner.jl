@testset "groebner" begin
    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    I = ideal(R,[x*y^2 - x, x^3 - 2*y^5])
    groebner_basis(I, ordering=lex(R), algorithm=:fglm)
    @test gens(I.gb[lex(R)]) == QQMPolyRingElem[y^7 - y^5, x*y^2 - x, x^3 - 2*y^5]
    groebner_basis(I, ordering=degrevlex(R))
    @test gens(I.gb[degrevlex(R)]) == QQMPolyRingElem[x*y^2 - x, x^4 - 2*x*y, -x^3 + 2*y^5]
    @test_throws ErrorException groebner_basis(I, ordering=neglex(R))
    standard_basis(I, ordering=neglex(R))
    @test gens(I.gb[neglex(R)]) == QQMPolyRingElem[-x^3 + 2*y^5, x]
    @test leading_ideal(I, ordering=degrevlex(gens(R))) == ideal(R,[x*y^2, x^4, y^5])
    @test leading_ideal(I) == ideal(R,[x*y^2, x^4, y^5])
    @test leading_ideal(I, ordering=lex(gens(R))) == ideal(R,[y^7, x*y^2, x^3])

    # algorithm = :modular option
    R = @polynomial_ring(QQ, [:x, :y])
    I = ideal(R, [x^2+y,y*x-1])
    # uses f4 in msolve/AlgebraicSolving
    groebner_basis(I, ordering=degrevlex(R), algorithm=:modular)
    @test gens(I.gb[degrevlex(R)]) == QQMPolyRingElem[x + y^2, x*y - 1, x^2 + y]
    # uses multi-modular implementation in Oscar applying Singular finite field
    # computations
    groebner_basis(I, ordering=lex(R), algorithm=:modular)
    @test gens(I.gb[lex(R)]) == QQMPolyRingElem[y^3 + 1, x + y^2]
    T = @polynomial_ring(QQ, :t)
    R = @polynomial_ring(T, [:x, :y])
    I = ideal(R, [x^2+y,y*x-1])
    @test_throws ErrorException groebner_basis(I, ordering=lex(R), algorithm=:modular)

    # issue 3665
    kt,t = polynomial_ring(GF(2),:t)
    Ft = fraction_field(kt)
    P,(x,y) = polynomial_ring(Ft,[:x,:y])
    I = ideal(P,[x,y])
    @test normal_form([x], I) == [0]

    R, (x, y) = polynomial_ring(GF(5), [:x, :y])
    I = ideal(R, [x])
    gb = groebner_basis(I)
    @test normal_form(y, I) == y
    @test Oscar._normal_form_singular([y], I, degrevlex(R)) == [y]
    @test Oscar._normal_form_f4([y], I) == [y]
    @test normal_form(y, I) == y

    G = groebner_basis(I)
    J = ideal(R, y)
    @test reduce(J.gens, G) == [y]
    @test reduce(y, G) == y
    J = ideal(R, x*y^2)
    matrix, res = reduce_with_quotients(J.gens, G)
    @test matrix * gens(G) + res == gens(J)
    matrix, res = reduce_with_quotients(x*y^2, G)
    @test matrix * gens(G) + [res] == [x*y^2]
    units, quots, res = reduce_with_quotients_and_unit(J.gens, G)
    @test matrix * gens(G) + res == units * gens(J)
    unit, quots, res = reduce_with_quotients_and_unit(x*y^2, G)
    @test matrix * gens(G) + [res] == units * [x*y^2]
    @test reduce(y^3, [y^2 - x, x^3 - 2*y^2]) == x*y
    @test reduce(y^3, elem_type(R)[]) == y^3
    @test reduce(y^3, Oscar.IdealGens(R, elem_type(R)[])) == y^3
    @test reduce([y^3], [y^2 - x, x^3 - 2*y^2]) == [x*y]
    @test reduce([y^3], elem_type(R)[]) == [y^3]
    @test reduce([y^3], Oscar.IdealGens(R, elem_type(R)[])) == [y^3]
    @test reduce([y^3], [y^2 - x, x^3 - 2*y^2], ordering=lex(R)) == [y^3]
    f = x+y^3
    g = x
    @test reduce(f, [g], ordering=degrevlex(R)) == x+y^3
    @test reduce(f, [g], ordering=degrevlex(R), complete_reduction=true) == y^3
    f = x*y^3-y^4
    F = [x*y^2-x,x^3-2*x*y^2]
    q, r = reduce_with_quotients(f, F)
    @test q * F + [r] == [f]
    q, r = reduce_with_quotients(f, elem_type(R)[])
    @test q * elem_type(R)[] + [r] == [f]
    q, r = reduce_with_quotients(f, Oscar.IdealGens(R, elem_type(R)[]))
    @test q * elem_type(R)[] + [r] == [f]
    q, r = reduce_with_quotients([f], F)
    @test q * F + r == [f]
    q, r = reduce_with_quotients([f], elem_type(R)[])
    @test q * elem_type(R)[] + r == [f]
    q, r = reduce_with_quotients([f], Oscar.IdealGens(R, elem_type(R)[]))
    @test q * elem_type(R)[] + r == [f]
    u, q, r = reduce_with_quotients_and_unit(f, F)
    @test q * F + [r] == u * [f]
    u, q, r = reduce_with_quotients_and_unit(f, elem_type(R)[])
    @test q * elem_type(R)[] + [r] == u * [f]
    u, q, r = reduce_with_quotients_and_unit(f, Oscar.IdealGens(R, elem_type(R)[]))
    @test q * elem_type(R)[] + [r] == u * [f]
    u, q, r = reduce_with_quotients_and_unit([f], F)
    @test q * F + r == u * [f]
    u, q, r = reduce_with_quotients_and_unit([f], elem_type(R)[])
    @test q * elem_type(R)[] + r == u * [f]
    u, q, r = reduce_with_quotients_and_unit([f], Oscar.IdealGens(R, elem_type(R)[]))
    @test q * elem_type(R)[] + r == u * [f]
    f = x
    F = [1-x]
    q, r = reduce_with_quotients(f, F, ordering=neglex(R))
    @test q * F + [r] != [f]
    u, q, r = reduce_with_quotients_and_unit(f, F, ordering=neglex(R))
    @test q * F + [r] == u * [f]
    # Issue 3105
    R, (x,y,z) = QQ[:x, :y, :z]
    f = x^3 - x^2*y - x^2*z + x
    f1 = x^2*y - z
    f2 = x*y - 1
    _,r = reduce_with_quotients(f, [f1, f2], ordering = deglex(R))
    @test r == x^3-x^2*z
    I = ideal(R,[y^2 - x, x^3 - 2*y^2])
    @test is_groebner_basis(I.gens, ordering=degrevlex(R)) == true
    @test is_groebner_basis(I.gens, ordering=lex(R)) == false
    @test_throws ErrorException is_groebner_basis(I.gens, ordering=neglex(R))
    @test is_standard_basis(I.gens, ordering=degrevlex(R)) == true
    @test is_standard_basis(I.gens, ordering=neglex(R)) == false
    R, (x1, x2, x3, x4) = polynomial_ring(QQ, [:x1, :x2, :x3, :x4])
    J = ideal(R, [x1+2*x2+2*x3+2*x4-1,
       x1^2+2*x2^2+2*x3^2+2*x4^2-x1,
       2*x1*x2+2*x2*x3+2*x3*x4-x2,
       x2^2+2*x1*x3+2*x2*x4-x3])
    standard_basis(J, ordering=neglex(R))
    Oscar.groebner_assure(J, true, true)
    @test J.gb[degrevlex(R)].isReduced == true 
    @test_throws ErrorException standard_basis(J, ordering=negdeglex(R), algorithm=:fglm)
    standard_basis(J, ordering=lex(R), algorithm=:fglm)
    H = QQMPolyRingElem[128304*x4^8 - 93312*x4^7 + 15552*x4^6 + 3144*x4^5 - 1120*x4^4 + 36*x4^3 + 15*x4^2 - x4,
        5913075*x3 + 371438283744*x4^7 - 237550027104*x4^6 + 22645939824*x4^5 + 11520686172*x4^4 - 2024910556*x4^3 - 132524276*x4^2 + 30947828*x4,
        1971025*x2 - 97197721632*x4^7 + 73975630752*x4^6 - 12121915032*x4^5 - 2760941496*x4^4 + 814792828*x4^3 - 1678512*x4^2 - 9158924*x4,
        5913075*x1 - 159690237696*x4^7 + 31246269696*x4^6 + 27439610544*x4^5 - 6475723368*x4^4 - 838935856*x4^3 + 275119624*x4^2 + 4884038*x4 - 5913075]
    @test elements(J.gb[lex(R)]) == H
    R, (x,y) = polynomial_ring(AcbField(64),[:x,:y])
    I = ideal([x^2+y^2+1//3,x^2+x*y+1//3*x])
    @test_throws ArgumentError groebner_basis(I)

    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    I = ideal(R, [x+y^2, x^2+x*y+3*y])
    H = QQMPolyRingElem[x,y]
    standard_basis_highest_corner(I, ordering=negdeglex(R))
    @test elements(I.gb[negdeglex(R)]) == H
    @test_throws ArgumentError standard_basis_highest_corner(I)

    R, (x, y) = polynomial_ring(QQ, [:x, :y])
    I = ideal(R,[x^2+x*y+y,x+y^2])
    G = standard_basis(I, ordering=lex(R), complete_reduction=true)
    @test G.ord == lex(R)
end

@testset "normal form graded" begin
    R, (a, b, c) = graded_polynomial_ring(QQ, [:a, :b, :c])
    I = ideal([a^2+b^2, c^3])
    # verify normal form of inhomogeneous polynomial work
    @test normal_form(a^2+b, I) == -b^2+b
end

@testset "groebner leading ideal" begin
   R, (t, x, y, z) = polynomial_ring(QQ, [:t, :x, :y, :z])

   I = ideal(R, [x + y + t + z, x^2 + y^2 + t^3])

   @test leading_ideal(groebner_basis(I, ordering=lex([t, x, y, z]))) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, ordering=degrevlex([t, x, y, z]))) == ideal(R, [t, x^3])
   @test leading_ideal(groebner_basis(I, ordering=lex([t, x, y, z]))) == ideal(R, [x^3, t])
   @test leading_ideal(groebner_basis(I, ordering=invlex([t, x, y, z]))) == ideal(R, [y^2, z])
   @test leading_ideal(groebner_basis(I, ordering=wdeglex([t, x, y, z], [2, 3, 1, 4]))) == ideal(R, [z, t^3])
   @test leading_ideal(groebner_basis(I, ordering=wdegrevlex([t, x, y, z], [2, 1, 1, 1]))) == ideal(R, [t, x^3])
end

@testset "groebner orderings" begin
   R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])

   I = ideal(R, [x + y + z, x^2 + y^2 + z^3])

   @test groebner_basis(I, ordering=lex([x])*lex([y,z])) == groebner_basis(I, ordering=lex([x, y, z]))
   @test groebner_basis(I, ordering=lex([z])*lex([y])*lex([x])) == groebner_basis(I, ordering=invlex([x, y, z]))
   @test groebner_basis(I, ordering=degrevlex([x, y, z])*invlex([y])) == groebner_basis(I, ordering=degrevlex([x, y, z]))
   @test groebner_basis(I, ordering=deglex([z])*deglex([x])*deglex([y])) == groebner_basis(I, ordering=lex([z])*lex([x, y]))
   @test groebner_basis(I, ordering=deglex([x, y, z])) == groebner_basis(I, ordering=wdeglex([x, y, z], [1, 1, 1]))
   M = matrix_ordering([x, y, z], [1 1 1; 0 1 0; 1 0 0])
   @test gens(groebner_basis(I, ordering = M)) == [ x + y + z, 2*x^2 + 2*x*z + z^3 + z^2 ]
   @test_throws ErrorException groebner_basis(I, ordering = negdeglex([x, y, z]))
   @test gens(standard_basis(I, ordering=negdeglex([x, y, z]))) == [ x + y + z, 2*y^2 + 2*y*z + z^3 + z^2 ]

   G, U = groebner_basis_with_transformation_matrix(I, ordering=lex([x])*lex([y,z]))
   H, V = groebner_basis_with_transformation_matrix(I, ordering=lex([x, y, z]))
   @test gens(G) == gens(H) && U == V
   G, U = groebner_basis_with_transformation_matrix(I, ordering=lex([z])*lex([y])*lex([x]))
   H, V = groebner_basis_with_transformation_matrix(I, ordering=invlex([x, y, z]))
   @test gens(G) == gens(H) && U == V
   G, U = groebner_basis_with_transformation_matrix(I, ordering=degrevlex([x, y, z])*invlex([y]))
   H, V = groebner_basis_with_transformation_matrix(I, ordering=degrevlex([x, y, z]))
   @test gens(G) == gens(H) && U == V
   G, U = groebner_basis_with_transformation_matrix(I, ordering=deglex([z])*deglex([x])*deglex([y]))
   H, V = groebner_basis_with_transformation_matrix(I, ordering=lex([z])*lex([x, y]))
   @test gens(G) == gens(H) && U == V
   G, U = groebner_basis_with_transformation_matrix(I, ordering=deglex([x, y, z]))
   H ,V = groebner_basis_with_transformation_matrix(I, ordering=wdeglex([x, y, z], [1, 1, 1]))
   @test gens(G) == gens(H) && U == V
   @test_throws ErrorException groebner_basis_with_transformation_matrix(I, ordering = negdeglex([x, y, z]))
   G, M = standard_basis_with_transformation_matrix(I, ordering = negdeglex([x, y, z]))
   @test gens(I) * M == gens(G)
end

@testset "non-global orderings" begin
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, [x^2*(x-1)-y^2, y^3*(x+y-6)])
  o = negdegrevlex(gens(R))
  G = standard_basis(I, ordering=o)
  @test normal_form(x^5-5, I, ordering=o) == -5
  J = ideal(R, [x^5-5])
  units, matr, res = reduce_with_quotients_and_unit(J.gens, G)
  @test matr * gens(G) + res == units * gens(J)
  u = negdegrevlex([x])*negdegrevlex([y])
  @test ideal_membership(x^4, I, ordering=u)
end

@testset "f4" begin
  R, (x1,x2,x3,x4) = polynomial_ring(GF(next_prime(2^28)), [:x1, :x2, :x3, :x4])
  I = ideal(R,[x1+2*x2+2*x3+2*x4-1,
          x1^2+2*x2^2+2*x3^2+2*x4^2-x1,
          2*x1*x2+2*x2*x3+2*x3*x4-x2,
          x2^2+2*x1*x3+2*x2*x4-x3])
  H = groebner_basis_f4(I);
  G = [x1 + 2*x2 + 2*x3 + 2*x4 + 268435458
                x3^2 + 2*x2*x4 + 76695850*x3*x4 + 115043772*x4^2 + 115043768*x2 + 191739613*x3 + 230087535*x4
                x2*x3 + 268435457*x2*x4 + 230087533*x3*x4 + 76695842*x4^2 + 210913575*x2 + 38347923*x3 + 153391692*x4
                x2^2 + 2*x2*x4 + 153391692*x3*x4 + 230087538*x4^2 + 230087536*x2 + 115043768*x3 + 191739613*x4
                x3*x4^2 + 238609298*x4^3 + 14913081*x2*x4 + 56338306*x3*x4 + 129246702*x4^2 + 263464432*x2 + 260150414*x3 + 258493405*x4
                x2*x4^2 + 89478486*x4^3 + 29826162*x2*x4 + 263464432*x3*x4 + 238609297*x4^2 + 141674270*x2 + 9942054*x3
                x4^4 + 35851649*x4^3 + 232885084*x2*x4 + 251179133*x3*x4 + 231479137*x4^2 + 191484965*x2 + 85955250*x3 + 167408122*x4]
  @test elements(H) == G
  @test isdefined(I, :gb)
  @test Oscar.oscar_generators(I.gb[degrevlex(gens(base_ring(I)))]) == G
  @test length(I.gb) == 1
  H = groebner_basis_f4(I, eliminate=2);
  G = [x3^2*x4 + 73209671*x3*x4^2 + 260301051*x4^3 + 188447115*x3^2 + 167207272*x3*x4 + 120660383*x4^2 + 210590781*x3 + 109814506*x4
                x3^3 + 156877866*x3*x4^2 + 59264971*x4^3 + 224858274*x3^2 + 183605206*x3*x4 + 130731555*x4^2 + 110395535*x3 + 158620953*x4
                x4^4 + 167618101*x3*x4^2 + 102789335*x4^3 + 193931678*x3^2 + 156155981*x3*x4 + 60823186*x4^2 + 239040667*x3 + 127377432*x4
                x3*x4^3 + 99215126*x3*x4^2 + 261328123*x4^3 + 132228634*x3^2 + 93598185*x3*x4 + 85654356*x4^2 + 3613010*x3 + 240673711*x4]
  @test elements(H) == G
  @test length(I.gb) == 1
  # issue 5216
  R, (a,b,c,d,x,y,z,w) = polynomial_ring(QQ, ["a", "b", "c", "d", "x", "y", "z", "w"])
  I = ideal(R, [x - a*c, y - a*c*d, z - a*c^2 - b, w - a*c^2*d - b*d])
  H = groebner_basis_f4(I; eliminate=4);
  G= [-x*w + y*z]
  @test elements(H) == G
end

@testset "fglm" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  A = Oscar.IdealGens(R, [x*y-1, x^2+y^2])
  @test_throws ErrorException Oscar._fglm(A, lex(R))
  I = ideal(R, gens(A))
  groebner_basis(I, algorithm=:buchberger)
  @test_throws ErrorException Oscar._fglm(I.gb[degrevlex(R)], lex(R))
  groebner_basis(I, complete_reduction=true)
  G = Oscar._fglm(I.gb[degrevlex([x, y])], lex(R))
  @test gens(G) == QQMPolyRingElem[y^4 + 1, x + y^3]
  J = ideal(R, [x])
  @test_throws ErrorException groebner_basis(J, algorithm=:fglm)
  J.gb = Dict()
  fglm(I, destination_ordering=lex(R))
  @test haskey(I.gb, degrevlex(R))
  @test I.gb[degrevlex(R)].isReduced
  @test haskey(I.gb, lex(R))
  @test gens(I.gb[lex(R)]) == QQMPolyRingElem[y^4 + 1, x + y^3]
  R, (x, y) = polynomial_ring(ZZ, [:x, :y])
  I = ideal(R, [x*y-1, x^2+y^2])
  @test_throws ErrorException groebner_basis(I, algorithm=:fglm)
  @test_throws ErrorException fglm(I, destination_ordering=lex(R))

  # Issue #4026
  R, (x, y) = QQ[:x, :y]
  I = ideal([x, y])
  G = groebner_basis(I, ordering = lex([y, x]), algorithm = :fglm)
  @test gens(G) == elem_type(R)[x, y]
end

@testset "groebner hilbert driven" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal([x*(x+1), x^2 - y^2 + (x-2) * y])
  gb = standard_basis(I, ordering = lex(R), algorithm = :hilbert)
  @test is_groebner_basis(gb, ordering = lex(R))
  S, (x, y) = grade(R)
  I = ideal([x*(x+y), x^2 - y^2 + (x-2*y) * y])
  gb = standard_basis(I, ordering = deglex(R), algorithm = :hilbert)
  @test is_groebner_basis(gb, ordering = deglex(R))
  R, (x, y) = polynomial_ring(GF(65521), [:x, :y])
  I = ideal([x*(x+1), x^2 - y^2 + (x-2) * y])
  gb = standard_basis(I, ordering = lex(R), algorithm = :hilbert)
  @test is_groebner_basis(gb, ordering = lex(R))
  @test haskey(I.gb, lex(R))
end

@testset "groebner basis modular" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R, [x^2, x*y + 32771*y^2])
  gb = Oscar.groebner_basis_modular(I, ordering = degrevlex(R))
  @test is_groebner_basis(I.gb[degrevlex(R)], ordering = degrevlex(R))
  @test all(iszero, Oscar.reduce(groebner_basis(I), gb))
  @test all(iszero, Oscar.reduce(gb, groebner_basis(I)))
  J = ideal(R, [x+y^2, x*y+y^3])
  J.gb[degrevlex(R)] = Oscar.IdealGens(R, [x^3])
  @test Oscar._certify_modular_groebner_basis(J, degrevlex(R)) == false
  groebner_basis_modular(J, ordering=wdegrevlex(R,[2,1]), certify=true)
  @test gens(J.gb[wdegrevlex([x, y], [2, 1])]) == QQMPolyRingElem[x+y^2]
end

@testset "known groebner bases" begin
  R, (x, y) = QQ[:x, :y]
  I = ideal(R, x)
  # A priori we expect nothing to be known
  @test !Oscar.is_known(is_one, I)
  @test !Oscar.is_known(krull_dim, I)
  @test !Oscar.is_known(groebner_basis, I)
  # ...except things which are really easy to check.
  @test Oscar.is_known(is_zero, I)
  # The following triggers a groebner basis computation
  is_one(I)
  # and then stuff is known.
  @test Oscar.is_known(is_one, I)
  # Is any groebner basis known?
  @test Oscar.is_known(groebner_basis, I)
  # Is a groebner basis for lex known?
  @test !Oscar.is_known(groebner_basis, I; ordering=lex(gens(R)))
end

@testset "factoring standard resp. groebner bases" begin
  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R,[x*y])
  L = factoring_standard_basis(I)
  @test gens(L[1]) == QQMPolyRingElem[y]
  @test gens(L[2]) == QQMPolyRingElem[x]
  @test L[1].isGB == true
  @test L[2].isGB == true

  R, (x, y, z) = polynomial_ring(GF(32003), [:x, :y, :z])
  I = ideal(R,[x*y*z,y*z+x])

  L = factoring_standard_basis(I)

  @test length(L) == 1
  @test gens(L[1]) == FqMPolyRingElem[x, y*z+x]

  L = factoring_standard_basis(I, ordering=lex(R))

  @test length(L) == 2
  @test gens(L[1]) == FqMPolyRingElem[z, x + y*z]
  @test gens(L[2]) == FqMPolyRingElem[y, x + y*z]

  R, (x, y) = polynomial_ring(QQ, [:x, :y])
  I = ideal(R,[x*y])
  L = factoring_groebner_basis(I)
  @test gens(L[1]) == QQMPolyRingElem[y]
  @test gens(L[2]) == QQMPolyRingElem[x]
  @test L[1].isGB == true
  @test L[2].isGB == true

  @test_throws ErrorException L = factoring_groebner_basis(I, ordering=neglex(R))
end
