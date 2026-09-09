@testset "MPolyAnyMap" begin
  mktempdir() do path
    R, (x, y) = QQ[:x, :y]
    S, (s, t, m) = QQ[:s, :t, :m]
    phi = hom(R, S, [s * t, m])

    test_save_load_roundtrip(path, phi) do loaded
      loaded == phi
    end

    Q, (z, w) = R[:z, :w]
    T, (u, v) = S[:u, :v]
    psi = hom(Q, T, phi, [u * v, v - u])

    test_save_load_roundtrip(path, psi) do loaded
      loaded == psi
    end
  end
end
