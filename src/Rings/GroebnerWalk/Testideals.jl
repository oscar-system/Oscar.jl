include("GroebnerWalk.jl")
include("Examples.jl")
using Test
 
@testset "Groebnerwalks" begin
    @testset "Testing Groebnerwalks" begin
        let id = katsura5()
            R = base_ring(id)
            dim = nvars(R)
            TarOrd = ordering_as_matrix(:lex, dim)
            I = groebner_basis(id, complete_reduction = true)

            ideals = []
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined,2))

            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

            for i = 2:nvars(R)-2
                push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i))
            end
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :standard, 2))
            #push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :tran))

            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_start_order))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

            S = matrix_ordering(R,TarOrd)

        
            s = groebner_basis(id, ordering=S,complete_reduction = true)


            for i in ideals
                @test equalitytest(
                    Oscar.IdealGens(R,gens(i), S),
                    s,
                )
            end
        end

        let id = katsura4()
            R = base_ring(id)
            dim = nvars(R)
            StartOrd = ordering_as_matrix([1,3,1,3,1],:lex)
            TarOrd = ordering_as_matrix([1,0,0,0,2],:lex)
            I = groebner_basis(id, ordering=matrix_ordering(R, StartOrd), complete_reduction = true)


            ideals = []
            for i = 1:nvars(R)-1
                push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :pertubed, i))
            end
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :standard))
            #push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :tran))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :fractal_combined))
            push!(ideals, groebnerwalk(I, ordering_as_matrix([1,3,1,3,1],:lex), ordering_as_matrix([1,0,0,0,2],:lex), :generic))

            S = matrix_ordering(R,TarOrd)

        
            s = groebner_basis(id, ordering=S,complete_reduction = true)


            for i in ideals
                @test equalitytest(
                    Oscar.IdealGens(R,gens(i), S),
                    s,
                )
            end
        end

        let id = cyclic4()
            R = base_ring(id)
            dim = nvars(R)
            TarOrd = ordering_as_matrix(:lex, dim)
            I = groebner_basis(id, complete_reduction = true)

            ideals = []
            for i = 2:nvars(R)-1
                push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i))
            end
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :standard))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :tran))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

            S = matrix_ordering(R,TarOrd)

        
            s = groebner_basis(id, ordering=S,complete_reduction = true)


            for i in ideals
                @test equalitytest(
                    Oscar.IdealGens(R,gens(i), S),
                    s,
                )
            end
        end

        let id = ex2()
            R = base_ring(id)
            dim = nvars(R)
            TarOrd = ordering_as_matrix(:lex, dim)
            I = groebner_basis(id, complete_reduction = true)

            ideals = []
            for i = 2:nvars(R)
                push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i))
            end
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :standard))
            #push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :tran))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_start_order))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

            S = matrix_ordering(R,TarOrd)

        
            s = groebner_basis(id, ordering=S,complete_reduction = true)


            for i in ideals
                @test equalitytest(
                    Oscar.IdealGens(R,gens(i), S),
                    s,
                )
            end
        end


        id = redeco7()
        R = base_ring(id)
        dim = nvars(R)
        TarOrd = ordering_as_matrix(:lex, dim)
        I = groebner_basis(id, complete_reduction = true)

        ideals = []
        for i = 2:nvars(R)-3
            push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :pertubed, i))
        end
        push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :standard))
        push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :fractal_combined))
        push!(ideals, groebnerwalk(I, ordering_as_matrix(:degrevlex, dim), ordering_as_matrix(:lex, dim), :generic))

        S = matrix_ordering(R,TarOrd)

        
        s = groebner_basis(id, ordering=S,complete_reduction = true)


        for i in ideals
            @test equalitytest(
                Oscar.IdealGens(R,gens(i), S),
                s,
            )
        end
    end
end
