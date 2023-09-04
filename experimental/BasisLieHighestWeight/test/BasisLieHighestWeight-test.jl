using Oscar
using Test
# using TestSetExtensions
include("MBOld.jl")
forGap = Oscar.GAP.julia_to_gap

"""
We are testing our code in multiple ways. First, we calculated two small examples per hand and compare those. Then we 
check basic properties of the result. For example we know the size of our monomial basis. These properties get partially
used in the algorithm and could therefore be true for false results. We have another basic algorithm that solves the 
problem without the recursion, weightspaces and saving of computations. The third test compares the results we can 
compute with the weaker version.
"""

function compare_algorithms(dynkin::Char, n::Int64, lambda::Vector{Int64})
    # old algorithm
    mons_old = MBOld.basisLieHighestWeight(string(dynkin), n, lambda) # basic algorithm

    # new algorithm
    base = BasisLieHighestWeight.basis_lie_highest_weight(string(dynkin), n, lambda)
    mons_new = base.monomial_basis.set_mon
    L = Oscar.GAP.Globals.SimpleLieAlgebra(forGap(string(dynkin)), n, Oscar.GAP.Globals.Rationals)
    gap_dim = Oscar.GAP.Globals.DimensionOfHighestWeightModule(L, forGap(lambda)) # dimension

    # comparison
    # convert set of monomials over different ring objects to string representation to compare for equality
    @test issetequal(string.(mons_old), string.(mons_new)) # compare if result of old and new algorithm match
    @test gap_dim == length(mons_new) # check if dimension is correct
    print(".")
end

function check_dimension(dynkin::Char, n::Int64, lambda::Vector{Int64}, monomial_order::String)
    base = BasisLieHighestWeight.basis_lie_highest_weight(string(dynkin), n, lambda, monomial_order=monomial_order) 
    mons_new = base.monomial_basis.set_mon
    L = Oscar.GAP.Globals.SimpleLieAlgebra(forGap(string(dynkin)), n, Oscar.GAP.Globals.Rationals)
    gap_dim = Oscar.GAP.Globals.DimensionOfHighestWeightModule(L, forGap(lambda)) # dimension
    @test gap_dim == length(mons_new) # check if dimension is correct
end

@testset "Test BasisLieHighestWeight" begin
    @testset "is_fundamental" begin
        @test BasisLieHighestWeight.is_fundamental([ZZ(0), ZZ(1), ZZ(0)])
        @test !BasisLieHighestWeight.is_fundamental([ZZ(0), ZZ(1), ZZ(1)])
    end

    @testset "compute_sub_weights" begin
        @test isequal(BasisLieHighestWeight.compute_sub_weights([ZZ(0), ZZ(0), ZZ(0)]), [])
        sub_weights =  convert(Vector{Vector{ZZRingElem}}, [[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1], [2, 0, 0],
                        [0, 2, 0], [2, 1, 0], [1, 2, 0], [2, 0, 1], [0, 2, 1], [2, 1, 1], [1, 2, 1], [2, 2, 0],
                        [0, 3, 0], [2, 2, 1], [1, 3, 0], [0, 3, 1], [1, 3, 1], [2, 3, 0]])
        @test isequal(BasisLieHighestWeight.compute_sub_weights([ZZ(2), ZZ(3), ZZ(1)]), sub_weights)
    end

    @testset "Known examples basis_lie_highest_weight" begin
        base = BasisLieHighestWeight.basis_lie_highest_weight("A", 2, [1,0])
        mons = base.monomial_basis.set_mon
        @test issetequal(string.(mons), Set(["1", "x3", "x1"]))
        base = BasisLieHighestWeight.basis_lie_highest_weight("A", 2, [1,0], operators=[1,2,1])
        mons = base.monomial_basis.set_mon
        @test issetequal(string.(mons), Set(["1", "x2*x3", "x3"]))
    end
    @testset "Compare basis_lie_highest_weight with algorithm of Johannes and check dimension" begin
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
        @testset "Monomial order $monomial_order" for monomial_order in ("lex", "revlex", "degrevlex")
            # the functionality longest-word was temporarily removed because it required coxeter groups from 
            # https://github.com/jmichel7/Gapjm.jl
            #@testset "Operators $ops" for ops in ("regular", "longest-word") 
            check_dimension('A', 3, [1,1,1], monomial_order)
                #check_dimension('B', 3, [2,1,0], monomial_order, ops)
                #check_dimension('C', 3, [1,1,1], monomial_order, ops)
                #check_dimension('D', 4, [3,0,1,1], monomial_order, ops)
                #check_dimension('F', 4, [2,0,1,0], monomial_order, ops)
                #check_dimension('G', 2, [1,0], monomial_order, ops)
                #check_dimension('G', 2, [2,2], monomial_order, ops)
            #end
        end
    end
end
