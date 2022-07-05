# Setup for fields
R, x = PolynomialRing(QQ, "x")
q = x^2 + 3//4
K, a = NumberField(q)
Ky, y = K["y"]
Tow, b = NumberField(y^2 + 1, "b")
NonSimRel, c = NumberField([y^2 - 5 * a, y^2 - 7 * a])
Zt, t = PolynomialRing(ResidueRing(ZZ, 2), "t")
Fin, d = FiniteField(t^2 + t + 1)
Frac = FractionField(R)

cases = [
    (QQ, fmpq(3, 4), fmpq(1, 2), "Rationals"),
    (ZZ, 3, 4, "Integers"),
    (ResidueRing(ZZ, 6), 3, 5, "Integers Modulo 6"),
    (K, a, a + 1, "Simple Extension"),
    (Tow, a^2 * b, a + b, "Tower Extension"),
    (NonSimRel, c[1], c[2] * a, "Non Simple Extension"),
    (Fin, d, 0, "Finite Field"),
    (Frac, 1 // x, x^2, "Fraction Field")
]


function get_hom(R1::T, R2::T) where T <: Union{
    MPolyRing{NfAbsNSElem}, PolyRing{NfAbsNSElem}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)
    return hom(D, I, gens(I))
end

function get_hom(R1::T, R2::T) where T <: Union{
    MPolyRing{Hecke.NfRelElem{nf_elem}}, PolyRing{Hecke.NfRelElem{nf_elem}},
    AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.Frac{fmpq_poly}}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)
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


function get_hom(R1::T, R2::T) where {
    T <: Union{MPolyRing{S}, PolyRing{S}} where S <: Union{
        nf_elem, nmod, fmpz, fmpq, fq_nmod, fmpq_poly}}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)

    return hom(D, I, gen(I))
end

function test_equality(p::MPolyElem{T}, l::MPolyElem{T}) where T <: Union{fmpq, fmpz, nmod}
    P = parent(p)
    L = parent(l)
    h = hom(P, L, gens(L))
    return h(p) == l
end

function test_equality(p::AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_poly}},
                       l::AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.Frac{fmpq_poly}}) 
    P = parent(p)
    L = parent(l)
 
    mapped_coeffs = map(i -> evaluate(i, x), coefficients(l))
    return mapped_coeffs == collect(coefficients(p))
end

function test_equality(p::MPolyElem{T}, l::MPolyElem{T}) where T <: Union{
    Hecke.NfRelNSElem{nf_elem}, Hecke.NfRelElem{nf_elem}, fq_nmod, nf_elem}
    P = parent(p)
    L = parent(l)
    h = get_hom(P, L)
    return [h(c) for c in coefficients(p)] == collect(coefficients(l))
end

function test_equality(p::PolyElem, l::PolyElem)
    P = parent(p)
    L = parent(l)

    if p isa PolyElem{T} where T <: Union{fmpq, fmpz, nmod}
        return L(collect(coefficients(p))) == l
        
    elseif P isa AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.Frac{fmpq_poly}}
        mapped_coeffs = map(i -> evaluate(i, x), coefficients(l))
        return mapped_coeffs == collect(coefficients(p))
    else 
        h = get_hom(P, L)
        return [h(c) for c in coefficients(p)] == collect(coefficients(l))
    end
end

@testset "Polynomials" begin
    mktempdir() do path
        for case in cases
            @testset "Univariate Polynomial over $(case[4])" begin
                R, z = PolynomialRing(case[1], "z")
                p = z^2 + case[2] * z + case[3]
                test_save_load_roundtrip(path, p) do loaded
                  S = parent(loaded)
                  @test test_equality(p, loaded)
                end
            end
            
            @testset "Multivariate Polynomial over $(case[4])" begin
                R, (z, w) = PolynomialRing(case[1], ["z", "w"])
                p = z^2 + case[2] * z * w + case[3] * w^3
                test_save_load_roundtrip(path, p) do loaded
                  @test test_equality(p, loaded)
                end

                @testset "MPoly Ideals over $(case[4])" begin
                    q = w^2 - z
                    i = ideal(R, [p, q])
                    test_save_load_roundtrip(path, i) do loaded_i
                        if R isa MPolyRing{T} where T <: Union{fmpq, fmpz, nmod}
                            S = parent(loaded_i[1])
                            h = hom(R, S, gens(S))
                            @test h(i) == loaded_i
                        end
                    end
                end
            end
        end
    end
end
