# Setup for coefficient rings
R, x = PolynomialRing(QQ, "x")
q = x^2 + 3//4
L, (e, f) = NumberField([x^2 + 5, x^3 - 6])
K, a = NumberField(q)
Ky, y = K["y"]
Tow, b = NumberField(y^2 + 1, "b")
NonSimRel, c = NumberField([y^2 - 5 * a, y^2 - 7 * a])
Zt, t = PolynomialRing(ResidueRing(ZZ, 2), "t")
Fin, d = FiniteField(t^2 + t + 1)
Frac = FractionField(R)
P7 = PadicField(7, 30)
T = TropicalSemiring()

cases = [
    (QQ, fmpq(3, 4), fmpq(1, 2), "Rationals"),
    (ZZ, 3, 4, "Integers"),
    (R, x^2, x + 1, "Iterated Multivariate PolyRing"),
    (ResidueRing(ZZ, 6), 3, 5, "Integers Modulo 6"),
    (L, e, f, "Non Simple Extension"),
    (K, a, a + 1, "Simple Extension"),
    (Tow, a^2 * b, a + b, "Tower Extension"),
    (NonSimRel, c[1], c[2] * a, "Non Simple Rel Extension"),
    (Fin, d, 1, "Finite Field"),
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
    AbstractAlgebra.PolyRing{AbstractAlgebra.Generic.Frac{fmpq_poly}}}
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
        AbstractAlgebra.Generic.Frac{fmpq_poly}}
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
        nf_elem, nmod, fmpz, fmpq, fq_nmod}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)
    return hom(D, I, gen(I))
end

function get_hom(R1::T, R2::T) where T <: Union{
    SeriesRing{S}, Generic.LaurentSeriesField{S}} where S <: Union{
        nf_elem, nmod, fmpz, fmpq, fq_nmod}
    D = base_ring(R1)
    I = base_ring(R2)
    return hom(D, I, gen(I))
end

function test_equality(p::T, l::T) where T <: (
    MPolyElem{S} where S <:Union{
        fmpq, fmpz, nmod, padic, Oscar.TropicalSemiringElem
    })
    P = parent(p)
    L = parent(l)
    h = hom(P, L, gens(L))
    return h(p) == l
end

function test_equality(p::T, l::T) where T <: (
    PolyElem{S} where S <: Union{
        fmpq, fmpz, nmod, padic, Oscar.TropicalSemiringElem})
    P = parent(p)
    L = parent(l)
    return L(collect(coefficients(p))) == l
end

function test_equality(p::T, l::T) where T <: (
    RelSeriesElem{S} where S <: Union{fmpq, fmpz, nmod, padic})
    L = parent(l)
    coeffs = map(o -> coeff(p, o), 0:pol_length(p))
    return L(coeffs, pol_length(p), precision(p), valuation(p)) == l
end

function test_equality(p::T, l::T) where T <: (
    AbsSeriesElem{S} where S <: Union{fmpq, fmpz, nmod, padic})
    L = parent(l)
    coeffs = map(o -> coeff(p, o), 0:pol_length(p))
    
    return L(coeffs, pol_length(p), precision(p)) == l
end

function test_equality(p::fmpz_laurent_series, l::fmpz_laurent_series)
    L = parent(l)
    v = valuation(p)
    coeffs = map(o -> coeff(p, o), v : v + pol_length(p))
    return L(coeffs, pol_length(p), precision(p), v, Nemo.scale(p)) == l
end

function test_equality(p::T, l::T) where T <: (
    Generic.LaurentSeriesElem{S} where S <: Union{fmpq, padic})
    L = parent(l)
    v = valuation(p)
    coeffs = map(o -> coeff(p, o), v : v + pol_length(p))
    return L(coeffs, pol_length(p), precision(p), v, Generic.scale(p)) == l
end

function test_equality(p::T, l::T) where T <: (
    Generic.LaurentSeriesRingElem{S} where S <: nmod)
    L = parent(l)
    v = valuation(p)
    coeffs = map(o -> coeff(p, o), v : v + pol_length(p))
    return L(coeffs, pol_length(p), precision(p), v, Generic.scale(p)) == l
end

function test_equality(p::T, l:: T) where T  <: Union{
    MPolyElem{S}, PolyElem{S}} where S <: Union{
    AbstractAlgebra.Generic.Frac{fmpq_poly}, fmpq_poly}
    g = gen(base_ring(parent(l)))
    mapped_coeffs = map(i -> evaluate(i, g), coefficients(p))
    return mapped_coeffs == collect(coefficients(l))
end

function test_equality(p::T, l:: T) where T  <: Union{
    Generic.LaurentSeriesElem{S},SeriesElem{S}} where S <: Union{
    AbstractAlgebra.Generic.Frac{fmpq_poly}, fmpq_poly}
    dom = base_ring(parent(p))
    codom = base_ring(parent(l))
    g = gen(base_ring(parent(l)))
    evaluate_on_gen = map_from_func(y -> evaluate(y, g), dom, codom)
    return compare_series_coeffs(p, l, evaluate_on_gen)
end


function test_equality(p::T, l::T) where T <: (
    MPolyElem{S} where S <: Union{
        Hecke.NfRelNSElem{nf_elem},
        Hecke.NfRelElem{nf_elem},
        NfAbsNSElem,
        fq_nmod,
        nf_elem})
    P = parent(p)
    L = parent(l)
    h = get_hom(P, L)
    return [h(c) for c in coefficients(p)] == collect(coefficients(l))
end

function test_equality(p::T, l::T) where T <: (
    PolyElem{S} where S <: Union{
        Hecke.NfRelNSElem{nf_elem},
        Hecke.NfRelElem{nf_elem},
        NfAbsNSElem,
        fq_nmod,
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
        fq_nmod,
        nf_elem})
    L = parent(l)
    P = parent(p)
    h = get_hom(P, L)
    return compare_series_coeffs(p, l, h)
end

function compare_series_coeffs(p::T, l::T,
                               h::Union{Map, typeof(identity)}) where T <: RelSeriesElem
    coeffs_p = map(o -> h(coeff(p, o)), 0:pol_length(p))
    L = parent(l)
    return L(coeffs_p, pol_length(p), precision(p), valuation(p)) == l
end

function compare_series_coeffs(p::T, l::T,
                               h::Union{Map, typeof(identity)}) where T <: AbsSeriesElem
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

@testset "Polynomials and Series" begin
    mktempdir() do path
        for case in cases
            @testset "Univariate Polynomial over $(case[4])" begin
                R, z = PolynomialRing(case[1], "z")
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
                R, (z, w) = PolynomialRing(case[1], ["z", "w"])
                p = z^2 + case[2] * z * w + case[3] * w^3
                test_save_load_roundtrip(path, p) do loaded
                  @test test_equality(p, loaded)
                end

                @testset "Load with parent" begin
                    test_save_load_roundtrip(path, p; parent=R) do loaded
                        @test p == loaded
                    end
                end

                @testset "MPoly Ideals over $(case[4])" begin
                    q = w^2 + z
                    i = Oscar.ideal(R, [p, q])
                    test_save_load_roundtrip(path, i) do loaded_i
                        if R isa MPolyRing{T} where T <: Union{fmpq, fmpz, nmod}
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
                    rel_R, rel_z = PowerSeriesRing(case[1], 10, "z")
                    rel_p = rel_z^2 + case[2] * rel_z + case[3]
                    test_save_load_roundtrip(path, rel_p) do loaded
                        @test test_equality(rel_p, loaded)
                    end

                    test_save_load_roundtrip(path, rel_p; parent=rel_R) do loaded
                        @test rel_p == loaded
                    end

                    abs_R, abs_z = PowerSeriesRing(case[1], 10, "z"; model=:capped_absolute)
                    abs_p = abs_z^2 + case[2] * abs_z + case[3]
                    test_save_load_roundtrip(path, abs_p) do loaded
                        @test test_equality(abs_p, loaded)
                    end

                    test_save_load_roundtrip(path, abs_p; parent=abs_R) do loaded
                        @test abs_p == loaded
                    end
                end

                @testset "Laurent Series over $(case[4])" begin
                    L, z = LaurentSeriesRing(case[1], 10, "z")
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
