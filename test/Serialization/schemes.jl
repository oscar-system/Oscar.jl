@testset "localizations of polynomial rings" begin
  mktempdir() do path
    R, (x, y) = QQ[:x, :y]
    U = powers_of_element(x*y)
    L, _ = localization(R, U)
    test_save_load_roundtrip(path, L) do loaded
      is_unit(L(x^2*y))
      !is_unit(L(x-1))
    end
  end
end

