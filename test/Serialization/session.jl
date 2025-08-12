@testset "Session Functionality" begin
  mktempdir() do path
    R, x = polynomial_ring(QQ, :x)

    test_save_load_roundtrip(path, R) do loaded
      @test loaded == R
    end

    poly_path = joinpath(path, "poly")
    R_path = joinpath(path, "R")
    save(poly_path, x^2 + 1)
    save(R_path, R)

    # simulates loading a new session
    Oscar.reset_global_serializer_state()

    loaded_R = load(R_path)
    loaded_poly = load(poly_path)

    @test loaded_R == parent(loaded_poly)

    # simulates loading a new session
    Oscar.reset_global_serializer_state()

    loaded_poly = load(poly_path)
    loaded_R = load(R_path)

    @test loaded_R == parent(loaded_poly)

  end
end
