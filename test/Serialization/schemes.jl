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

  end
end

