using Oscar
using Test
using TestSetExtensions
#using SparseArrays

include("MBOld.jl")

G = Oscar.GAP.Globals
forGap = Oscar.GAP.julia_to_gap
fromGap = Oscar.GAP.gap_to_julia

"""
We are testing our code in multiple ways. First, we calculated two small examples per hand and compare those. Then we check basic properties of the result.
For example we know the size of our monomial basis. These properties get partially used in the algorithm and could therefore be true for false results. We
have another basic algorithm that solves the problem without the recursion, weightspaces and saving of computations. The third test compares the results we 
can compute with the weaker version.
"""

function compare_algorithms(dynkin::Char, n::Int64, lambda::Vector{Int64})
    #print("TESTETSETSET", dynkin, n, lambda)
    dim, m, v = MBOld.basisLieHighestWeight(string(dynkin), n, lambda) # basic algorithm
    w = BasisLieHighestWeight.basisLieHighestWeight2(string(dynkin), n, lambda) # algorithm that needs to be tested
    L = G.SimpleLieAlgebra(forGap(string(dynkin)), n, G.Rationals)
    gapDim = G.DimensionOfHighestWeightModule(L, forGap(lambda)) # dimension
    @test Set(m) == w # compare if result of basic and sophisticated algorithm match
    @test gapDim == length(w) # check if dimension is correct
end

function check_dimension(dynkin::Char, n::Int64, lambda::Vector{Int64}, monomial_order::String)
    w = BasisLieHighestWeight.basisLieHighestWeight2(string(dynkin), n, lambda, monomial_order=monomial_order) # algorithm that needs to be tested
    L = G.SimpleLieAlgebra(forGap(string(dynkin)), n, G.Rationals)
    gapDim = G.DimensionOfHighestWeightModule(L, forGap(lambda)) # dimension
    @test gapDim == length(w) # check if dimension is correct
end


@testset ExtendedTestSet "Test basisLieHighestWeight" begin
    # TODO: add test for basis (not just dimension)
    @testset "Known examples" begin
        @test BasisLieHighestWeight.basisLieHighestWeight2("A", 2, [1,0]) == Set([[0,0,0], [0,0,1], [1,0,0]])
        @test BasisLieHighestWeight.basisLieHighestWeight2("A", 2, [1,0], ops=[1,2,1]) == Set([[0,0,0], [0,1,1], [1,0,0]])
    end
    @testset "Compare with simple algorithm and check dimension" begin
        @testset "Dynkin type $dynkin" for dynkin in ('A', 'B', 'C', 'D')
            @testset "n = $n" for n in 1:4
                if (!(dynkin == 'B' && n < 2) && !(dynkin == 'C' && n < 2) && !(dynkin == 'D' && n < 4))
                    for i in 1:n                                # w_i
                       lambda = zeros(Int64,n)
                        lambda[i] = 1
                       compare_algorithms(dynkin, n, lambda)
                    end
                    
                    if (n > 1)
                        lambda = [1, (0 for i in 1:n-2)..., 1]  # w_1 + w_n
                        compare_algorithms(dynkin, n, lambda)
                    end
    
                    if (n < 4)
                        lambda = ones(Int64,n)                  # w_1 + ... + w_n
                        compare_algorithms(dynkin, n, lambda)
                    end
                end 
            end
        end
    end
    @testset "Check dimension" begin
        @testset "Monomial order $monomial_order" for monomial_order in ("Lex", "GLex", "GRevLex", "GRevLex")
           #@testset "Operators $ops" for ops in ("regular", "longest-word")
            check_dimension('A', 3, [1,1,1], monomial_order)
    #            #check_dimension('B', 3, [2,1,0], monomial_order, ops)
    #            #check_dimension('C', 3, [1,1,1], monomial_order, ops)
    #            #check_dimension('D', 4, [3,0,1,1], monomial_order, ops)
    #            #check_dimension('F', 4, [2,0,1,0], monomial_order, ops)
    #            #check_dimension('G', 2, [1,0], monomial_order, ops)
    #            #check_dimension('G', 2, [2,2], monomial_order, ops)
            #end
        end
    end
end
