@testset "Serialization of sparse modules" begin
  mktempdir() do path
    R, (x, y) = QQ[:x, :y]

    F = FreeMod(R, 2)
    test_save_load_roundtrip(path, F) do loaded
      #@test symbols(F) == symbols(loaded)
      @test F === loaded
    end
    v = x*F[1]
    test_save_load_roundtrip(path, v) do loaded
      @test v == loaded
    end

    I = ideal(R, [x, y])
    II, inc = I*F
    test_save_load_roundtrip(path, II) do loaded
      #@test symbols(F) == symbols(loaded)
      @test II === loaded
    end
    v = x*II[1]
    test_save_load_roundtrip(path, v) do loaded
      @test v == loaded
    end
  end
end

