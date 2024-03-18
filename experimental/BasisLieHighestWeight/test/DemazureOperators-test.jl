using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion, _is_weighted
using Test

include("../src/LieAlgebras.jl")
include("../src/BirationalSequence.jl")
include("../src/NewMonomial.jl")
include("../src/RootConversion.jl")

include("../src/DemazureOperatorsDictionary.jl")

function calc_base_dicts(base::Oscar.BasisLieHighestWeight.MonomialBasis)::Vector{Dict{Vector{Int}, Int}}
  L = base.lie_algebra
  weight_w = base.birational_sequence.weights_w
  
  base_dicts = []
  n = length(base.birational_sequence.operators)
  monomials_list = collect(base.monomials)

  for i in reverse(0:n)
    indices = [idx for (idx, vec) in enumerate([leading_exponent_vector(mon) for mon in monomials_list]) 
                    if length(vec) == n && ((i == 0) || all(vec[1:i] .== 0))]
    filtered_monomials = monomials_list[indices]

    dict = Dict{Vector{Int}, Int}()
    for mon in filtered_monomials
      mon_weight = [Int(x) for x in weight(mon, weight_w)]  # Convert to Vector{Int}
      dict[mon_weight] = get(dict, mon_weight, 0) + 1  # Increment the count for this weight
    end
    push!(base_dicts, dict)
  end
  return base_dicts
end

# TODO Test with B3, highest_weight = [2, 2, 2], reduced_expression = [3, 2, 3, 1, 2, 3, 1, 2, 1]
# Oder kleiner test


@testset "Test DemazureOperatorsDimensions" begin
  type = :A
  rank = 4
  highest_weight = [1, 1, 1, 1]
  reduced_expression = [4, 3, 2, 1, 2, 3, 4, 3, 2, 3]
  #type = :A
  #rank = 2
  #highest_weight = [1, 1]
  #reduced_expression = [1, 2, 1]

  base = BasisLieHighestWeight.basis_lie_highest_weight(type, rank, highest_weight, reduced_expression)
  demazure_dimensions =  demazure_operators_summary_dictionary(type, rank, highest_weight, reduced_expression)
  base_dicts = calc_base_dicts(base)

  sub_word = []
  for i in (length(reduced_expression)+1):-1:1
      demazure_dict_old = demazure_dimensions[i]
      # subtract highest_weight from each key
      demazure_dict = Dict{Vector{Int}, Int}()
      for (key, value) in demazure_dict_old
          new_key = - key .+ highest_weight
          demazure_dict[new_key] = value
      end
      # demazure_dict = demazure_dict_old

      base_dict = base_dicts[i]
      @testset "sub_word" begin
        println("sub_word: ", sub_word)
        println("total dimensions demazure_dict: ", sum(values(demazure_dict)))
        println("total dimensions base_dict: ", sum(values(base_dict)))
        @test isequal(sum(values(demazure_dict)), sum(values(base_dict)))  
        @test isequal(demazure_dict, base_dict)
      end

      if i <= length(reduced_expression)
        alpha_i = reduced_expression[i]
        push!(sub_word, alpha_i)
      end
  end

  println("")

  # Convert keys of last base_dict to ZZRingElem
  demazure_dict = get_dim_weightspace_demazure(lie_algebra(type, rank), ZZ.(highest_weight), reduced_expression)
  base_dict = Dict{Vector{ZZRingElem}, Int}()
  for (key, value) in base_dicts[end]
      transformed_key = ZZ.(key)
      base_dict[transformed_key] = value
  end


  println("sub_word: ", sub_word)
  println("total dimensions demazure_dict: ", sum(values(demazure_dict)))
  println("total dimensions base_dict: ", sum(values(base_dict)))
  @test isequal(sum(values(demazure_dict)), sum(values(base_dict)))  
  @test isequal(demazure_dict, base_dict)
end
