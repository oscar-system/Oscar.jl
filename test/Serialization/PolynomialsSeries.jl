# Setup for coefficient rings
R, x = polynomial_ring(QQ, :x)
q = x^2 + 3//4
L, (e, f) = number_field([x^2 + 5, x^3 - 6])
K, a = number_field(q)
Ky, y = K[:y]
Tow, b = number_field(y^2 + 1, "b")
NonSimRel, c = number_field([y^2 - 5 * a, y^2 - 7 * a])
Qu, u = rational_function_field(QQ, "u")
Zt, t = polynomial_ring(residue_ring(ZZ, 2)[1], :t)
Fin, d = Nemo.Native.finite_field(t^2 + t + 1)
Frac = fraction_field(R)
P7 = PadicField(7, 30)
T = tropical_semiring()
PF = GF(7)
F  = GF(2, 2)
Fs, s = F[:s]
FF, r = finite_field(s^2 + gen(F) * s + 1, "r")

cases = [
  (QQ, QQFieldElem(3, 4), QQFieldElem(1, 2), "Rationals"),
  (R, x^2, x + 1, "Iterated Multivariate PolyRing"),
  (residue_ring(ZZ, 6)[1], 3, 5, "Integers Modulo 6"),
  (L, e, f, "Non Simple Extension"),
  (K, a, a + 1, "Simple Extension"),
  (Tow, a^2 * b, a + b, "Tower Extension"),
  (NonSimRel, c[1], c[2] * a, "Non Simple Rel Extension"),
  (Fin, d, 1, "Finite Field"),
  (Qu, u, 1 // u, "RationalFunctionField"),
  (Frac, 1 // x, x^2, "Fraction Field"),
  (T, T(1), T(3)^2, "Tropical Semiring"),
  (PF, PF(1), PF(6), "Default Prime Field"),
  (F, F(1), F(1), "Default Prime Field Extension"),
  (FF, FF(1), r, "Default Finite Field"),
  (QQBarField(), sqrt(QQBarField()(-7)), QQBarField()(5)^(QQ(4//5)), "QQBar"),
  (P7, 7 + 3*7^2, 7^5, "Padic Field"),
]

@testset "Serialization.Polynomials.and.Series" begin
  mktempdir() do path
    @testset "Empty Ideal" begin
      i = Oscar.ideal(QQ[:x, :y][1], [])
      test_save_load_roundtrip(path, i) do loaded
        @test loaded == i
      end
    end

    @testset "Graded Ring" begin
      R, (x, y) = QQ[:x, :y]
      A = [1 3; 2 1]
      M, (m1, m2) = grade(R, A)

      test_save_load_roundtrip(path, m1 * m2) do loaded
        @test loaded == m1 * m2
        @test grading_group(parent(loaded)) == grading_group(M)
      end

      GM, _ = grade(M, A)
      test_save_load_roundtrip(path, GM) do loaded
        @test loaded == GM
        @test forget_grading(loaded) == forget_grading(GM)
      end
    end

    for case in cases
      @testset "Univariate Polynomial over $(case[4])" begin
        R, z = polynomial_ring(case[1], :z)
        p = z^2 + case[2] * z + case[3]
        test_save_load_roundtrip(path, p) do loaded
          @test loaded == p
        end
        
        @testset "Load with params" begin
          test_save_load_roundtrip(path, p; params=R) do loaded
            @test loaded == z^2 + case[2] * z + case[3]
          end
        end
      end

      @testset "Multivariate Polynomial over $(case[4])"  begin
        R, (z, w) = polynomial_ring(case[1], [:z, :w])
        p = z^2 + case[2] * z * w + case[3] * w^3
        test_save_load_roundtrip(path, p) do loaded
          @test loaded == z^2 + case[2] * z * w + case[3] * w^3
        end

        @testset "Load with params" begin
          test_save_load_roundtrip(path, p; params=R) do loaded
            @test loaded == z^2 + case[2] * z * w + case[3] * w^3
          end
        end

        if R isa MPolyRing{T} where T <: Union{QQFieldElem, ZZRingElem, zzModRingElem}
          @testset "MPoly Ideals over $(case[4])" begin
            q = z
            i = Oscar.ideal(R, [p, q])
            test_save_load_roundtrip(path, i) do loaded_i
              S = parent(loaded_i[1])
              h = hom(R, S, gens(S))
              @test h(i) == loaded_i
            end

            S = parent(i[1])
            test_save_load_roundtrip(path, i; params=S) do loaded_i
              @test i == loaded_i
            end

            gb = groebner_basis(i)
            test_save_load_roundtrip(path, gb;) do loaded_gb
              @test gens(gb) == gens(loaded_gb)
              @test ordering(gb) == ordering(loaded_gb)
            end

            # need a cleaner way to setup type params in general
            params = Dict(Oscar.params(Oscar.type_params(gb))...)
            params[:ordering_type] = Oscar.type(params[:ordering_type])
            test_save_load_roundtrip(path, gb; params=params) do loaded_gb
              @test gens(gb) == gens(loaded_gb)
              @test ordering(gb) == ordering(loaded_gb)
            end
          end
        end
      end

      @testset "Universal Polynomial over $(case[4])" begin
        R = universal_polynomial_ring(case[1])
        z, w = gens(R, ["z", "w"])
        p = z^2 + case[2] * z * w + case[3] * w^3
        test_save_load_roundtrip(path, p) do loaded
          test_p = z^2 + case[2] * z * w + case[3] * w^3
          @test loaded.p == test_p.p
        end

        @testset "Load with params" begin
          test_save_load_roundtrip(path, p; params=R) do loaded
            @test p == loaded
          end
        end
      end

      # Tropical Semirings currently can't have formal power series
      filter!(case-> case[4] != "Tropical Semiring", cases)
      
      @testset "Multivariate Laurent Polynomial over $(case[4])" begin
        R, (z, w) = laurent_polynomial_ring(case[1], [:z, :w])
        p = z^2 + case[2] * z * w^(-4) + case[3] * w^(-3)
        test_save_load_roundtrip(path, p) do loaded
          @test loaded == z^2 + case[2] * z * w^(-4) + case[3] * w^(-3)
        end

        @testset "Load with params" begin
          test_save_load_roundtrip(path, p; params=R) do loaded
            @test p == loaded
          end
        end

        if R isa AbstractAlgebra.Generic.LaurentMPolyWrapRing{T} where T <: Union{QQFieldElem, ZZRingElem, zzModRingElem}
          @testset "Laurent MPoly Ideals over $(case[4])" begin
            q = w^2 + z
            i = Oscar.ideal(R, [p, q])
            test_save_load_roundtrip(path, i) do loaded_i
              S = parent(gens(loaded_i)[1])
              h = hom(R, S, gens(S))
              @test gens(i) == gens(loaded_i)
            end

            S = parent(gens(i)[1])
            test_save_load_roundtrip(path, i; params=S) do loaded_i
              @test gens(i) == gens(loaded_i)
            end
          end
        end
      end

      @testset "Series" begin
        @testset "Power Series over $(case[4])" begin
          rel_R, rel_z = power_series_ring(case[1], 10, "z")
          rel_p = rel_z^2 + case[2] * rel_z + case[3] * rel_z^3
          test_save_load_roundtrip(path, rel_p) do loaded
            @test rel_p == loaded
          end

          test_save_load_roundtrip(path, rel_p; params=rel_R) do loaded
            @test rel_p == loaded
          end

          abs_R, abs_z = power_series_ring(case[1], 10, "z"; model=:capped_absolute)
          abs_p = abs_z^2 + case[2] * abs_z + case[3]
          test_save_load_roundtrip(path, abs_p) do loaded
            @test abs_p == loaded
          end

          test_save_load_roundtrip(path, abs_p; params=abs_R) do loaded
            @test abs_p == loaded
          end
        end

        @testset "Laurent Series over $(case[4])" begin
          L, z = laurent_series_ring(case[1], 10, "z")
          p = z^(-1) + case[2] * z + case[3] * z^2
          test_save_load_roundtrip(path, p) do loaded
            @test p == loaded
          end

          test_save_load_roundtrip(path, p; params=L) do loaded
            @test p == loaded
          end
        end
      end
    end
  end
end

@testset "localizations and quotients" begin
  mktempdir() do path
    R, (x, y, z) = GF(103)[:x, :y, :z]
    for U in [powers_of_element(x), 
              complement_of_point_ideal(R, GF(103).([1, 2, 3])),
              complement_of_prime_ideal(ideal(R, x))
             ]
      L, _ = localization(R, U)
      Q, _ = quo(L, ideal(L, y))
      Qz = Q(z)
      # Test MPolyQuoLocRings and their elements
      test_save_load_roundtrip(path, Qz; params=Q) do loaded
        @test Qz == loaded
      end
      # Test ideals in these rings. 
      J = ideal(Q, z)
      test_save_load_roundtrip(path, J; params=Q) do loaded
        @test J == loaded
      end
      # Test MPolyLocalizedRingHom
      phi = hom(L, L, gens(L))
      test_save_load_roundtrip(path, phi) do loaded
        @test phi == loaded
      end   
      # Test MPolyQuoLocalizedRingHom
      phi = hom(Q, Q, gens(Q))
      test_save_load_roundtrip(path, phi) do loaded
        @test phi == loaded
      end
    end
  end
end
