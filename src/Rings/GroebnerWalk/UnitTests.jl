using Test
include("GroebnerWalk.jl")
include("Examples.jl")

@testset "UnitTests" begin
    @testset "Testing GroebnerwalkUtilitys" begin

        R, (x, y, z) = polynomial_ring(
            QQ,
            ["x", "y", "z"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 3)),
        )

        f1 = 3 * x^2 + 16 * x^2 * z + 14 * x^2 * y^3
        f2 = y^3 * z + 17 * x^2 * z^2 + 7 * x^2 * y^2 * z^2 + 13 * x^3 * z^2
        I = Singular.Ideal(R, [f1, f2])
        sol = [14 * x^2 * y^3, y^3 * z + 7 * x^2 * y^2 * z^2]
        @test initials(R, gens(I), [1, 3, 1]) == sol

        @test difference_lead_tail(I) ==
              [[0, 3, -1], [0, 3, 0], [-1, 2, 0], [2, -1, 1], [0, 2, 0]]

        F = [
            13 * x^3 * z^2,
            14 * x^2 * y^3,
            98 * x * y^5 * z^2,
            y^7 * z + x^2 * z^3,
            14 * x * y^6 * z,
        ]
        g = y^7 * z + x^2 * z^3 + 28 * x^2 * y^4
        q = Array{Singular.elem_type(R),1}(undef, 5)
        q[1] = R(0)
        q[2] = R(2 * y)
        q[3] = R(0)
        q[4] = R(1)
        q[5] = R(0)
        @test division_algorithm(g, F, R) == q

        J = Singular.Ideal(R, [f2, f1])
        f1 = 4 * x^2 + 16 * x^2 * z + 14 * x^2 * y^3
        f2 = y^3 * z + 17 * x^2 * z^2 + 7 * x^2 * y^2 * z^2 + 13 * x^3 * z^2
        K = Singular.Ideal(R, [f1, f2])
        @test equalitytest(I, J) == true
        @test equalitytest(I, K) == false

        @test deg(f1, 3) == 5

        id = trinks1()
        R = base_ring(id)
        dim = nvars(R)
        ve = ones(Int, dim)
        StartOrd = ordering_as_matrix(:degrevlex, dim)
        TarOrd = ordering_as_matrix(:lex, dim)
        I = Singular.std(id, complete_reduction = true)
        @test pertubed_vector(I, StartOrd, 1) == [1, 1, 1, 1, 1, 1]
        @test pertubed_vector(I, StartOrd, 2) == [4, 4, 4, 4, 4, 3]
        @test pertubed_vector(I, StartOrd, 3) == [49, 49, 49, 49, 48, 42]
        @test pertubed_vector(I, StartOrd, 4) ==
              [1000, 1000, 1000, 999, 990, 900]
        @test pertubed_vector(I, TarOrd, 1) == [1, 0, 0, 0, 0, 0]
        @test pertubed_vector(I, TarOrd, 2) == [4, 1, 0, 0, 0, 0]
        @test pertubed_vector(I, TarOrd, 3) == [49, 7, 1, 0, 0, 0]
        @test pertubed_vector(I, TarOrd, 4) == [1000, 100, 10, 1, 0, 0]

        @test inCone(I, [1000, 1000, 1000, 999, 990, 900]) == true
        @test inCone(I, [100, 1000, 1000, 999, 990, 900]) == false

        @test dot([1, 2, 3, 4], [2, 2, 2, 2]) == 20

    end

    @testset "Testing FraktalWalkUtilitys" begin
        R, (x, y, z) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y", "z"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 3)),
        )

        f1 = 3 * x^2
        f2 = y^3 * z + 17 * x^2 * z^2
        I = Singular.Ideal(R, [f1, f2])
        f1 = 3 * x^2
        f2 = y^3 * z
        J = Singular.Ideal(R, [f1, f2])
        f1 = x^3 * y^2 + z^2 + y^2
        f2 = y^3
        K = Singular.Ideal(R, [f1, f2])

        @test ismonomial(gens(I)) == false
        @test ismonomial(gens(J)) == true
        @test isbinomial(gens(I)) == true
        @test isbinomial(gens(J)) == true
        @test isbinomial(gens(K)) == false
        @test ismonomial(gens(K)) == false
    end

    @testset "Testing GenericWalkUtilitys" begin
        R, (x, y) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y"],
            ordering = Singular.ordering_M(ordering_as_matrix(:degrevlex, 2)),
        )

        f1 = x^2 - y^3
        f2 = x^3 - y^2 - x
        I = Singular.Ideal(R, [f1, f2])
        G = Singular.std(I, complete_reduction = true)

        @test next_gamma(
            gens(G),
            [Singular.leading_term(g) for g in gens(G)],
            [0],
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
        ) == [-2, 3]

        Rn, (x, y) = Singular.PolynomialRing(
            Singular.QQ,
            ["x", "y"],
            ordering = Singular.ordering_M(ordering_as_matrix(:lex, 2)),
        )
        g = [Rn(y^3 - x^2), Rn(x^3)]
        @test facet_initials(
            [change_ring(x, Rn) for x in gens(G)],
            [change_ring(Singular.leading_term(g), Rn) for g in gens(G)],
            [-2, 3],
        ) == g

        @test difference_lead_tail(
            gens(G),
            [Singular.leading_term(g) for g in gens(G)],
        ) == [[-2, 3], [3, -2], [2, 0]]

        @test isparallel([1, 2], [1, 4]) == false
        @test isparallel([1, 2], [2, 4]) == true
        @test isparallel([-1, 0], [-2, 1]) == false
        @test isparallel([-1, 0], [2, 0]) == true

        @test less_facet(
            [-2, 3],
            [-1, 4],
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
        ) == true
        @test less_facet(
            [-1, 7],
            [-1, 4],
            ordering_as_matrix(:degrevlex, 2),
            ordering_as_matrix(:lex, 2),
        ) == false

        dim = 5
        ve = [1, 1, 1, 1, 1]
        StartOrd = ordering_as_matrix(:degrevlex, dim)
        TarOrd = ordering_as_matrix(:lex, dim)
        R, (a, b, c, d, e) = Singular.PolynomialRing(
            Singular.QQ,
            ["a", "b", "c", "d", "e"],
            ordering = Singular.ordering_M(StartOrd),
        )
        S = change_order(R, TarOrd)
        J = Singular.Ideal(
            R,
            [
                b + 3 * b^3 + 2 * b * c * e + 5 * b * d * e,
                4 + b^2 + 4 * b * c + 5 * b^3 + c * d * e,
                d * e + 5 * b^2 * e,
            ],
        )
        I = Singular.std(J, complete_reduction = true)
        f1 = R(
            a^3 +
            a^2 +
            b^5 * a^3 * c^9 +
            e^3 +
            b^2 * a^2 * c^4 +
            d^3 +
            e^3 +
            b^2 * d^2 * e^7,
        )
        f2 = R(a^3 * b^3)
        f3 = R(a^2 * b^4 + c^2 + a^3 * 4 + a * e^3)
        f4 = R(0)

        @test (reduce(f3, I)) ==
              reduce_walk(f3, gens(I), [Singular.leading_term(g) for g in gens(I)])
        @test (reduce(f1, I)) ==
              reduce_walk(f1, gens(I), [Singular.leading_term(g) for g in gens(I)])
        @test (reduce(f2, I)) ==
              reduce_walk(f2, gens(I), [Singular.leading_term(g) for g in gens(I)])
        @test (reduce(f4, I)) ==
              reduce_walk(f4, gens(I), [Singular.leading_term(g) for g in gens(I)])
        J = Singular.std(J)
        
        @test equalitytest(
            Singular.std(J, complete_reduction = true),
            Singular.Ideal(
                R,
                interreduce(
                    collect(gens(J)),
                    [Singular.leading_term(g) for g in gens(J)],
                ),
            ),
        )

    end
end
