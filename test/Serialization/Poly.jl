# Setup for fields
R, x = PolynomialRing(QQ, "x")
q = x^2 + 3//4
K, a = NumberField(q)
Ky, y = K["y"]
Tow, b = NumberField(y^2 + 1, "b")
NonSimRel, c = NumberField([y^2 - 5 * a, y^2 - 7 * a])
Zt, t = PolynomialRing(ResidueRing(ZZ, 2), "t")
Fin, d = FiniteField(t^2 + t + 1)

cases = [
    [QQ, fmpq(3, 4), fmpq(1, 2)],
    [ZZ, 3, 4],
    [ResidueRing(ZZ, 6), 3, 5],
    [K, a, a + 1],
    [Tow, a^2 * b, a + b],
    [NonSimRel, c[1], c[2] * a],
    [Fin, d, 0]
]

function get_hom(R1::T, R2::T) where T <: Union{MPolyRing, PolyRing}
    D = coefficient_ring(R1)
    I = coefficient_ring(R2)

    if D isa NfAbsNS
        return hom(D, I, gens(I))

    elseif D isa Hecke.NfRel{nf_elem}
        D_base_field = base_field(D)
        I_base_field = base_field(I)
        h_1 = hom(D_base_field, I_base_field, gen(I_base_field))
        return hom(D, I, h_1, gen(I))

    elseif D isa NfRelNS{nf_elem}
        D_base_field = base_field(D)
        I_base_field = base_field(I)
        h_1 = hom(D_base_field, I_base_field, gen(I_base_field))
        return hom(D, I, h_1, gens(I))
    end

    return hom(D, I, gen(I))
end

function test_equality(p::MPolyElem, l::MPolyElem)
    P = parent(p)
    L = parent(l)

    if p isa MPolyElem{T} where T <: Union{fmpq, fmpz, nmod}
        h = hom(P, L, gens(L))
        return h(p) == l
    else
        h = get_hom(P, L)
        return [h(c) for c in coefficients(p)] == collect(coefficients(l))
    end
end

function test_equality(p::PolyElem, l::PolyElem)
    P = parent(p)
    L = parent(l)
    
    if p isa PolyElem{T} where T <: Union{fmpq, fmpz, nmod}
        return L(collect(coefficients(p))) == l
    else 
        h = get_hom(P, L)
        return [h(c) for c in coefficients(p)] == collect(coefficients(l))
    end
end

@testset "Polynomials" begin
    mktempdir() do path
        for case in cases
            @testset "Univariate Polynomial over $(cases[1])" begin
                R, x = PolynomialRing(case[1], "x")
                p = x^2 + case[2] * x + case[3]
                filename = joinpath(path, "polynomial.uv")
                save(p, filename)
                loaded = load(filename)
                S = parent(loaded)
                @test test_equality(p, loaded)
            end
            
            @testset "Multivariate Polynomial over $(case[1])" begin
                R, (x, y) = PolynomialRing(case[1], ["x", "y"])
                p = x^2 + case[2] * x * y + case[3] * y^3
                filename = joinpath(path, "polynomial_.mv")
                save(p, filename)
                loaded = load(filename)
                @test test_equality(p, loaded)

                #@testset "MPoly Ideals over $(case[1])" begin
                    #q = y^2 - x
                    #i = ideal(R, [p, q])
                    #filename = joinpath(path, "ideal.mv")
                    #save(i, filename)
                    #loaded_i = load(filename)
                    #S = parent(loaded_i[1])
                    #h = hom(R, S, gens(S))
                    #@test h(i) == loaded_i
                #end
            end
        end
    end
end
