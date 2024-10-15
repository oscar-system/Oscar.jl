@testset "localizations of polynomial rings" begin
  mktempdir() do path
    R, (x, y) = QQ[:x, :y]
    U = powers_of_element(x*y)
    L, _ = localization(R, U)
    test_save_load_roundtrip(path, L) do loaded
      @test is_unit(L(x^2*y))
      @test !is_unit(L(x-1))
    end
    a = L(x)
    b = L(y*(x-1))
    test_save_load_roundtrip(path, a) do loaded
      @test is_unit(a)
    end
    test_save_load_roundtrip(path, b) do loaded
      @test !is_unit(b)
    end
    J = ideal(L, gens(L))
    test_save_load_roundtrip(path, J) do loaded
      @test is_one(J)
    end

    U = complement_of_prime_ideal(ideal(R, x+y))
    L, _ = localization(R, U)
    test_save_load_roundtrip(path, L) do loaded
      @test is_unit(L(x^2*y))
      @test !is_unit(L(x+y))
    end
    a = L(x)
    b = L(y*(x+y))
    test_save_load_roundtrip(path, a) do loaded
      @test is_unit(a)
    end
    test_save_load_roundtrip(path, b) do loaded
      @test !is_unit(b)
    end
    J = ideal(L, gens(L))
    test_save_load_roundtrip(path, J) do loaded
      @test is_one(J)
    end


    U = complement_of_point_ideal(R, QQ.([0, 0]))
    L, _ = localization(R, U)
    test_save_load_roundtrip(path, L) do loaded
      @test !is_unit(L(x^2*y))
      @test is_unit(L(x+y+1))
    end
    a = L(x-1)
    b = L(y*(x+y))
    test_save_load_roundtrip(path, a) do loaded
      @test is_unit(a)
    end
    test_save_load_roundtrip(path, b) do loaded
      @test !is_unit(b)
    end

    I = ideal(R, [x+y])
    Q, _ = quo(R, I)
    U = powers_of_element(x*y)
    L, _ = localization(Q, U)
    test_save_load_roundtrip(path, L) do loaded
      @test is_unit(L(x^2*y))
      @test !is_unit(L(x-1))
    end
    a = L(x)
    b = L(y*(x-1))
    test_save_load_roundtrip(path, a) do loaded
      @test is_unit(a)
    end
    test_save_load_roundtrip(path, b) do loaded
      @test !is_unit(b)
    end
    J = ideal(L, gens(L))
    test_save_load_roundtrip(path, J) do loaded
      @test is_one(J)
    end
  end
end

@testset "affine schemes" begin
  mktempdir() do path
    A = affine_space(QQ, 3)
    test_save_load_roundtrip(path, A) do loaded
      @test OO(A) isa MPolyRing{<:FieldElem}
      @test ngens(OO(A)) == 3
    end

    (x, y, z) = gens(OO(A))
    U = PrincipalOpenSubset(A, [x, y])
    test_save_load_roundtrip(path, U) do loaded
      @test OO(U) isa Oscar.MPolyLocRing
      @test ngens(OO(U)) == 3
      @test is_unit(OO(U)(x))
    end
  end
end

@testset "maps from polynomial rings" begin
  mktempdir() do path
    R, (x, y) = GF(29^2)[:x, :y]
    phi = hom(R, R, x->x^29, gens(R))
    @test_throws ArgumentError save(path, phi)
    # the following should in principle be supported for some serializations 
    # which are not for long term storage (to come)
    # test_save_load_roundtrip(path, phi) do loaded
    #   @test phi(one(R)) == one(R)
    # end
  end
end

@testset "covered schemes" begin
  mktempdir() do path
    IP = projective_space(GF(29), 3)
    X = covered_scheme(IP)
    test_save_load_roundtrip(path, X) do loaded
      @test length(patches(default_covering(X))) == 4
      @test length(gluings(default_covering(X))) == 16
    end
    cov = default_covering(X)
    for U in affine_charts(X)
      for V in affine_charts(X)
        glue = cov[U, V]
        glue isa SimpleGluing || continue
        test_save_load_roundtrip(path, glue) do loaded
          f, g = gluing_morphisms(glue)
          @test is_isomorphism(f)
          @test is_isomorphism(g)
        end
      end
    end
  end
end


@testset "covered schemes II" begin
  mktempdir() do path
    IA3 = affine_space(QQ, 3)
    I = ideal(OO(IA3), gens(OO(IA3)))
    pr = blow_up(IA3, I)
    Y = domain(pr)
    und = Oscar.underlying_morphism(pr)
    test_save_load_roundtrip(path, und) do loaded
      @test length(patches(default_covering(domain(und)))) == 3
    end
  end
end
