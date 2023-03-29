# Setup for coefficient rings
R, x = polynomial_ring(QQ, "x")
q = x^2 + 3//4
L, (e, f) = number_field([x^2 + 5, x^3 - 6])
K, a = number_field(q)
Ky, y = K["y"]
Tow, b = number_field(y^2 + 1, "b")
NonSimRel, c = number_field([y^2 - 5 * a, y^2 - 7 * a])
Qu, u = RationalFunctionField(QQ, "u")
Zt, t = polynomial_ring(residue_ring(ZZ, 2), "t")
Fin, d = FiniteField(t^2 + t + 1)
Frac = fraction_field(R)
P7 = PadicField(7, 30)
T = TropicalSemiring()

cases = [
    (QQ, QQFieldElem(3, 4), QQFieldElem(1, 2), "Rationals"),
    (R, x^2, x + 1, "Iterated Multivariate PolyRing"),
    (residue_ring(ZZ, 6), 3, 5, "Integers Modulo 6"),
    (L, e, f, "Non Simple Extension"),
    (K, a, a + 1, "Simple Extension"),
    (Tow, a^2 * b, a + b, "Tower Extension"),
    (NonSimRel, c[1], c[2] * a, "Non Simple Rel Extension"),
    (Fin, d, 1, "Finite Field"),
    (Qu, u, 1 // u, "RationalFunctionField"),
    (Frac, 1 // x, x^2, "Fraction Field"),
    (P7, 7 + 3*7^2, 7^5, "Padic Field"),
    (T, T(1), T(3)^2, "Tropical Semiring")
]

function get_hom(R1::T, R2::T) where T <: Union{
    MPolyRing{NfAbsNSElem}, PolyRing{NfAbsNSElem}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)
    return hom(D, I, gens(I))
end

function get_hom(R1::T, R2::T) where T <: Union{
    SeriesRing{S}, Generic.LaurentSeriesField{S}} where S <: NfAbsNSElem
    D = base_ring(R1)
    I = base_ring(R2)
    return hom(D, I, gens(I))
end

function get_hom(R1::T, R2::T) where T <: Union{
    MPolyRing{Hecke.NfRelElem{nf_elem}}, PolyRing{Hecke.NfRelElem{nf_elem}},
    AbstractAlgebra.PolyRing{AbstractAlgebra.Generic.Frac{QQPolyRingElem}},
    AbstractAlgebra.PolyRing{AbstractAlgebra.Generic.RationalFunctionField{QQFieldElem, QQPolyRingElem}}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)
    D_base_field = base_field(D)
    I_base_field = base_field(I)
    h_1 = hom(D_base_field, I_base_field, gen(I_base_field))
    return hom(D, I, h_1, gen(I))
end

function get_hom(R1::T, R2::T) where T <: (
    Union{SeriesRing{S}, Generic.LaurentSeriesField{S}} where S <: Union{
        Hecke.NfRelElem{nf_elem},
        AbstractAlgebra.Generic.Frac{QQPolyRingElem}}
    )
    D = base_ring(R1)
    I = base_ring(R2)
    D_base_field = base_field(D)
    I_base_field = base_field(I)
    h_1 = hom(D_base_field, I_base_field, gen(I_base_field))
    return hom(D, I, h_1, gen(I))
end

function get_hom(R1::T, R2::T) where T <: Union{
    MPolyRing{Hecke.NfRelNSElem{nf_elem}},
    PolyRing{Hecke.NfRelNSElem{nf_elem}}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)
    D_base_field = base_field(D)
    I_base_field = base_field(I)
    h_1 = hom(D_base_field, I_base_field, gen(I_base_field))
    return hom(D, I, h_1, gens(I))
end

function get_hom(R1::T, R2::T) where T <: Union{
    Generic.LaurentSeriesField{S}, SeriesRing{S}} where S <: Hecke.NfRelNSElem{nf_elem}
    D = base_ring(R1)
    I = base_ring(R2)
    D_base_field = base_field(D)
    I_base_field = base_field(I)
    h_1 = hom(D_base_field, I_base_field, gen(I_base_field))
    return hom(D, I, h_1, gens(I))
end

function get_hom(R1::T, R2::T) where {
    T <: Union{MPolyRing{S}, PolyRing{S}} where S <: Union{
        nf_elem, zzModRingElem, ZZRingElem, QQFieldElem, fqPolyRepFieldElem}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)
    return hom(D, I, gen(I))
end

function get_hom(R1::T, R2::T) where T <: Union{
    SeriesRing{S}, Generic.LaurentSeriesField{S}} where S <: Union{
        nf_elem, zzModRingElem, ZZRingElem, QQFieldElem, fqPolyRepFieldElem}
    D = base_ring(R1)
    I = base_ring(R2)
    return hom(D, I, gen(I))
end

function test_equality(p::T, l::T) where T <: (
    MPolyRingElem{S} where S <:Union{
        QQFieldElem, ZZRingElem, zzModRingElem, padic, Oscar.TropicalSemiringElem
    })
    P = parent(p)
    L = parent(l)
    h = hom(P, L, gens(L))
    return h(p) == l
end

function test_equality(p::T, l::T) where T <: (
    PolyRingElem{S} where S <: Union{
        QQFieldElem, ZZRingElem, zzModRingElem, padic, Oscar.TropicalSemiringElem})
    P = parent(p)
    L = parent(l)
    return L(collect(coefficients(p))) == l
end

function test_equality(p::T, l::T) where T <: (
    RelPowerSeriesRingElem{S} where S <: Union{QQFieldElem, ZZRingElem, zzModRingElem, padic})
    L = parent(l)
    v = valuation(p)
    pl = pol_length(p)
    coeffs = map(o -> coeff(p, o), v:v + pl)
    return L(coeffs, pl, precision(p), v) == l
end

function test_equality(p::T, l::T) where T <: (
    AbsPowerSeriesRingElem{S} where S <: Union{QQFieldElem, ZZRingElem, zzModRingElem, padic})
    L = parent(l)
    coeffs = map(o -> coeff(p, o), 0:pol_length(p))
    
    return L(coeffs, pol_length(p), precision(p)) == l
end

function test_equality(p::ZZLaurentSeriesRingElem, l::ZZLaurentSeriesRingElem)
    L = parent(l)
    v = valuation(p)
    coeffs = map(o -> coeff(p, o), v : v + pol_length(p))
    return L(coeffs, pol_length(p), precision(p), v, Nemo.scale(p)) == l
end

function test_equality(p::T, l::T) where T <: (
    Generic.LaurentSeriesElem{S} where S <: Union{QQFieldElem, padic})
    L = parent(l)
    v = valuation(p)
    coeffs = map(o -> coeff(p, o), v : v + pol_length(p))
    return L(coeffs, pol_length(p), precision(p), v, Generic.scale(p)) == l
end

function test_equality(p::T, l::T) where T <: (
    Generic.LaurentSeriesRingElem{S} where S <: zzModRingElem)
    L = parent(l)
    v = valuation(p)
    coeffs = map(o -> coeff(p, o), v : v + pol_length(p))
    return L(coeffs, pol_length(p), precision(p), v, Generic.scale(p)) == l
end

function test_equality(p::T, l:: T) where T  <: Union{
    MPolyRingElem{S}, PolyRingElem{S}} where S <: Union{
        AbstractAlgebra.Generic.Frac{QQPolyRingElem},
        AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem, QQPolyRingElem},
        QQPolyRingElem}
    g = gen(base_ring(parent(l)))
    mapped_coeffs = map(i -> evaluate(i, g), coefficients(p))
    return mapped_coeffs == collect(coefficients(l))
end

function test_equality(p::T, l:: T) where T  <: Union{
    Generic.LaurentSeriesElem{S}, SeriesElem{S}} where S <: Union{
        AbstractAlgebra.Generic.Frac{QQPolyRingElem},
        AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem, QQPolyRingElem},
        QQPolyRingElem}
    dom = base_ring(parent(p))
    codom = base_ring(parent(l))
    g = gen(base_ring(parent(l)))
    evaluate_on_gen = map_from_func(y -> evaluate(y, g), dom, codom)
    return compare_series_coeffs(p, l, evaluate_on_gen)
end


function test_equality(p::T, l::T) where T <: (
    MPolyRingElem{S} where S <: Union{
        Hecke.NfRelNSElem{nf_elem},
        Hecke.NfRelElem{nf_elem},
        NfAbsNSElem,
        fqPolyRepFieldElem,
        nf_elem})
    P = parent(p)
    L = parent(l)
    h = get_hom(P, L)
    return [h(c) for c in coefficients(p)] == collect(coefficients(l))
end

function test_equality(p::T, l::T) where T <: (
    PolyRingElem{S} where S <: Union{
        Hecke.NfRelNSElem{nf_elem},
        Hecke.NfRelElem{nf_elem},
        NfAbsNSElem,
        fqPolyRepFieldElem,
        nf_elem})
    P = parent(p)
    L = parent(l)
    h = get_hom(P, L)
    return [h(c) for c in coefficients(p)] == collect(coefficients(l))
end

function test_equality(p::T, l::T) where T <: (
    Union{SeriesElem{S}, Generic.LaurentSeriesElem} where S <: Union{
        Hecke.NfRelNSElem{nf_elem},
        Hecke.NfRelElem{nf_elem},
        NfAbsNSElem,
        fqPolyRepFieldElem,
        nf_elem})
    L = parent(l)
    P = parent(p)
    h = get_hom(P, L)
    return compare_series_coeffs(p, l, h)
end

function compare_series_coeffs(p::T, l::T,
                               h::Union{Map, typeof(identity)}) where T <: RelPowerSeriesRingElem
    v = valuation(p)
    pl = pol_length(p)
    coeffs_p = map(o -> h(coeff(p, o)), v:v + pl)
    L = parent(l)
    return L(coeffs_p, pl, precision(p), v) == l
end

function compare_series_coeffs(p::T, l::T,
                               h::Union{Map, typeof(identity)}) where T <: AbsPowerSeriesRingElem
    coeffs_p = map(o -> h(coeff(p, o)), 0:pol_length(p))
    L = parent(l)
    return L(coeffs_p, pol_length(p), precision(p)) == l
end

function compare_series_coeffs(p::T, l::T,
                               h::Union{Map, typeof(identity)}
                               ) where T <: Generic.LaurentSeriesElem
    v = valuation(p)
    coeffs_p = map(o -> h(coeff(p, o)), v:v + pol_length(p))
    L = parent(l)
    return L(coeffs_p, pol_length(p), precision(p), v, Generic.scale(p)) == l
end

@testset "Serialization.Polynomials.and.Series" begin
    mktempdir() do path
        for case in cases
            @testset "Univariate Polynomial over $(case[4])" begin
                R, z = polynomial_ring(case[1], "z")
                p = z^2 + case[2] * z + case[3]
                test_save_load_roundtrip(path, p) do loaded
                    @test test_equality(p, loaded)
                end

                @testset "Load with parent" begin
                    test_save_load_roundtrip(path, p; parent=R) do loaded
                        @test p == loaded
                    end
                end
            end
            
            @testset "Multivariate Polynomial over $(case[4])" begin
                R, (z, w) = polynomial_ring(case[1], ["z", "w"])
                p = z^2 + case[2] * z * w + case[3] * w^3
                test_save_load_roundtrip(path, p) do loaded
                  @test test_equality(p, loaded)
                end

                @testset "Load with parent" begin
                    test_save_load_roundtrip(path, p; parent=R) do loaded
                        @test p == loaded
                    end
                end

                if R isa MPolyRing{T} where T <: Union{QQFieldElem, ZZRingElem, zzModRingElem}
                    @testset "MPoly Ideals over $(case[4])" begin
                        q = w^2 + z
                        i = Oscar.ideal(R, [p, q])
                        test_save_load_roundtrip(path, i) do loaded_i
                            S = parent(loaded_i[1])
                            h = hom(R, S, gens(S))
                            @test h(i) == loaded_i
                        end
                    end
                end
            end

            # Tropical Semirings currently can't have formal power series
            filter!(case-> case[4] != "Tropical Semiring", cases)

            @testset "Series" begin
                @testset "Power Series over $(case[4])" begin
                    rel_R, rel_z = power_series_ring(case[1], 10, "z")
                    rel_p = rel_z^2 + case[2] * rel_z + case[3] * rel_z^3
                    test_save_load_roundtrip(path, rel_p) do loaded
                        @test test_equality(rel_p, loaded)
                    end

                    test_save_load_roundtrip(path, rel_p; parent=rel_R) do loaded
                        @test rel_p == loaded
                    end

                    abs_R, abs_z = power_series_ring(case[1], 10, "z"; model=:capped_absolute)
                    abs_p = abs_z^2 + case[2] * abs_z + case[3]
                    test_save_load_roundtrip(path, abs_p) do loaded
                        @test test_equality(abs_p, loaded)
                    end

                    test_save_load_roundtrip(path, abs_p; parent=abs_R) do loaded
                        @test abs_p == loaded
                    end
                end

                @testset "Laurent Series over $(case[4])" begin
                    L, z = laurent_series_ring(case[1], 10, "z")
                    p = z^(-1) + case[2] * z + case[3] * z^2
                    test_save_load_roundtrip(path, p) do loaded
                        @test test_equality(p, loaded)
                    end

                    test_save_load_roundtrip(path, p; parent=L) do loaded
                        @test p == loaded
                    end
                end
            end
        end
    end
end
