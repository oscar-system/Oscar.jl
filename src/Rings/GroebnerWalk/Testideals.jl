include("GroebnerWalk.jl")
include("Examples.jl")
using Test
 
@testset "Groebnerwalks" begin
    @testset "Testing Groebnerwalks" begin
        let id = katsura4()
            R = base_ring(id)
            dim = nvars(R)
            TarOrd = ordering_as_matrix(:lex, dim)
            S = change_order(R, TarOrd)
            I = groebner_basis(id, complete_reduction = true)

            ideals = []
            for i = 2:nvars(S)-2
                push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i,2))
            end
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :standard, 2,1))
            #push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :tran))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_start_order))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_look_ahead))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_lex))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

            s = Singular.std(
                Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                complete_reduction = true,
            )

            for id in ideals
                @test equalitytest(
                    Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                    s,
                )
            end
        end

        let id = katsura4()
            R = base_ring(id)
            dim = nvars(R)
            StartOrd = ordering_as_matrix([1,3,1,3,1],:lex)
            TarOrd = ordering_as_matrix([1,0,0,0,2],:lex)
            S = change_order(R, TarOrd)

            I = standard_basis(ideal(change_order(R,StartOrd), [change_ring(x, change_order(R,StartOrd)) for x in gens(id)]), complete_reduction = true)

            ideals = []
            for i = 1:nvars(S)-1
                push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :pertubed, i))
            end
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :standard))
            #push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :tran))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :fractal_start_order))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :fractal_look_ahead))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :fractal_lex))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :fractal_combined))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :generic))

            s = Singular.std(
                Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                complete_reduction = true,
            )

            for id in ideals
                @test equalitytest(
                    Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                    s,
                )
            end
        end

        let id = cyclic4()
            R = base_ring(id)
            dim = nvars(R)
            TarOrd = ordering_as_matrix(:lex, dim)
            S = change_order(R, TarOrd)

            I = Singular.std(id, complete_reduction = true)

            ideals = []
            for i = 2:nvars(S)-1
                push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i))
            end
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :standard))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :tran))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_start_order))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_look_ahead))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_lex))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

            s = Singular.std(
                Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                complete_reduction = true,
            )


            for id in ideals
                @test equalitytest(
                    Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                    s,
                )
            end
        end

        let id = ex2()
            R = base_ring(id)
            dim = nvars(R)
            TarOrd = ordering_as_matrix(:lex, dim)
            S = change_order(R, TarOrd)
            I = Singular.std(id, complete_reduction = true)

            ideals = []
            for i = 2:nvars(S)
                push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i))
            end
            push!(ideals, groebnerwalk(I, :degrevlex, :lex, :standard))
            #push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :tran))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_look_ahead))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_start_order))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_lex))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
            push!(ideals, groebnerwalk(I, :degrevlex, :lex, :generic))

            s = Singular.std(
                Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                complete_reduction = true,
            )

            for id in ideals
                @test equalitytest(
                    Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                    s,
                )
            end
        end


        id = redeco7()
        R = base_ring(id)
        dim = nvars(R)
        TarOrd = ordering_as_matrix(:lex, dim)
        S = change_order(R, TarOrd)
        I = Singular.std(id, complete_reduction = true)

        ideals = []
        for i = 2:nvars(S)-3
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i))
        end
        push!(ideals, groebnerwalk(I, :degrevlex, :lex, :standard))
        push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
        push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

        s = Singular.std(
            Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
            complete_reduction = true,
        )

        for id in ideals
            @test equalitytest(
                Singular.Ideal(S, [change_ring(x, S) for x in Singular.gens(id)]),
                s,
            )
        end
    end
end
