@testset "MPolyAnyMap" begin
  mktempdir() do path
    R, (x, y) = QQ[:x, :y]
    S, (s, t, m) = QQ[:s, :t, :m]
    phi = hom(R, S, [s * t, m])

    test_save_load_roundtrip(path, phi) do loaded
      loaded == phi
    end
  end
end
