include("MBOld.jl")

# We are testing our code in multiple ways. First, we calculated two small examples per hand and compare those. Then we 
# check basic properties of the result. For example we know the size of our monomial basis. These properties get partially
# used in the algorithm and could therefore be true for false results. We have another basic algorithm that solves the 
# problem without the recursion, weightspaces and saving of computations. The third test compares the results we can 
# compute with the weaker version.

function compare_algorithms(dynkin::Symbol, n::Int64, lambda::Vector{Int64})
  # old algorithm
  mons_old = MBOld.basisLieHighestWeight(string(dynkin), n, lambda) # basic algorithm

  # new algorithm
  basis = basis_lie_highest_weight(dynkin, n, lambda)
  mons_new = monomials(basis)
  L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(dynkin), n, GAP.Globals.Rationals)
  gap_dim = GAP.Globals.DimensionOfHighestWeightModule(L, GAP.Obj(lambda)) # dimension

  # comparison
  # convert set of monomials over different ring objects to string representation to compare for equality
  @test issetequal(string.(mons_old), string.(mons_new)) # compare if result of old and new algorithm match
  @test gap_dim == length(mons_new) # check if dimension is correct
end

function check_dimension(
  dynkin::Symbol, n::Int64, lambda::Vector{Int64}, monomial_ordering::Symbol
)
  basis = basis_lie_highest_weight(dynkin, n, lambda; monomial_ordering)
  L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(dynkin), n, GAP.Globals.Rationals)
  gap_dim = GAP.Globals.DimensionOfHighestWeightModule(L, GAP.Obj(lambda)) # dimension
  @test gap_dim == dim(basis) == length(monomials(basis)) # check if dimension is correct
end

@testset "sub_weights(_proper)" begin
  sub_weights = BasisLieHighestWeight.sub_weights
  sub_weights_proper = BasisLieHighestWeight.sub_weights_proper
  R = root_system(:B, 3)

  w_zero = zero(weight_lattice(R))
  @test issetequal(sub_weights(w_zero), [w_zero])
  @test isempty(sub_weights_proper(w_zero))

  w_231 = WeightLatticeElem(R, [2, 3, 1])
  sub_weights_proper_231 = [
    WeightLatticeElem(R, coeffs) for coeffs in [
      [0, 0, 1],
      [0, 1, 0],
      [0, 1, 1],
      [0, 2, 0],
      [0, 2, 1],
      [0, 3, 0],
      [0, 3, 1],
      [1, 0, 0],
      [1, 0, 1],
      [1, 1, 0],
      [1, 1, 1],
      [1, 2, 0],
      [1, 2, 1],
      [1, 3, 0],
      [1, 3, 1],
      [2, 0, 0],
      [2, 0, 1],
      [2, 1, 0],
      [2, 1, 1],
      [2, 2, 0],
      [2, 2, 1],
      [2, 3, 0],
    ]
  ]
  @test issetequal(sub_weights(w_231), [w_zero, w_231, sub_weights_proper_231...])
  @test issetequal(sub_weights_proper(w_231), sub_weights_proper_231)
end

@testset "Test BasisLieHighestWeight" begin
  @testset "Known examples basis_lie_highest_weight" begin
    base = basis_lie_highest_weight(:A, 2, [1, 0])
    mons = monomials(base)
    @test issetequal(string.(mons), Set(["1", "x3", "x1"]))
    base = basis_lie_highest_weight(:A, 2, [1, 0], [1, 2, 1])
    mons = monomials(base)
    @test issetequal(string.(mons), Set(["1", "x2*x3", "x3"]))
  end

  @testset "Compare basis_lie_highest_weight with algorithm of Johannes and check dimension" begin
    @testset "Dynkin type $dynkin" for dynkin in (:A, :B, :C, :D)
      @testset "n = $n" for n in 1:4
        if (
          !(dynkin == :B && n < 2) && !(dynkin == :C && n < 2) && !(dynkin == :D && n < 4)
        )
          for i in 1:n                                # w_i
            lambda = zeros(Int64, n)
            lambda[i] = 1
            compare_algorithms(dynkin, n, lambda)
          end

          if (n > 1)
            lambda = [1, (0 for i in 1:(n - 2))..., 1]  # w_1 + w_n
            compare_algorithms(dynkin, n, lambda)
          end

          if (n < 4)
            lambda = ones(Int64, n)                  # w_1 + ... + w_n
            compare_algorithms(dynkin, n, lambda)
          end
        end
      end
    end
  end

  @testset "Compare against GAP algorithm of Xin on some examples" begin
    basis_lusztig = basis_lie_highest_weight_lusztig(:A, 3, [2, 1, 1], [2, 3, 1, 2, 3, 1])

    @test issetequal(
      [only(exponents(m)) for m in monomials(basis_lusztig)],
      [
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1],
        [1, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 0, 1],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 2, 0],
        [0, 0, 0, 0, 1, 1],
        [1, 0, 0, 0, 1, 1],
        [0, 0, 1, 0, 0, 1],
        [0, 1, 0, 0, 1, 0],
        [0, 0, 0, 1, 0, 0],
        [1, 0, 0, 0, 2, 0],
        [0, 0, 1, 0, 1, 0],
        [2, 0, 0, 0, 1, 0],
        [2, 0, 0, 0, 0, 1],
        [0, 1, 0, 0, 0, 1],
        [0, 0, 0, 0, 2, 1],
        [1, 0, 0, 0, 2, 1],
        [0, 0, 1, 0, 1, 1],
        [0, 1, 0, 0, 2, 0],
        [0, 0, 0, 1, 1, 0],
        [2, 0, 0, 0, 1, 1],
        [1, 0, 1, 0, 0, 1],
        [1, 1, 0, 0, 1, 0],
        [1, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 1, 1],
        [0, 0, 0, 1, 0, 1],
        [2, 0, 0, 0, 2, 0],
        [1, 0, 1, 0, 1, 0],
        [1, 1, 0, 0, 0, 1],
        [0, 0, 1, 0, 2, 0],
        [2, 0, 0, 0, 2, 1],
        [1, 0, 1, 0, 1, 1],
        [0, 0, 2, 0, 0, 1],
        [1, 1, 0, 0, 2, 0],
        [1, 0, 0, 1, 1, 0],
        [0, 1, 1, 0, 1, 0],
        [1, 1, 0, 0, 1, 1],
        [1, 0, 0, 1, 0, 1],
        [0, 1, 1, 0, 0, 1],
        [0, 2, 0, 0, 1, 0],
        [0, 0, 1, 0, 2, 1],
        [0, 0, 0, 1, 2, 0],
        [0, 1, 0, 0, 2, 1],
        [0, 0, 0, 1, 1, 1],
        [1, 0, 1, 0, 2, 0],
        [3, 0, 0, 0, 2, 0],
        [3, 0, 0, 0, 1, 1],
        [1, 1, 0, 0, 2, 1],
        [1, 0, 0, 1, 1, 1],
        [0, 1, 1, 0, 1, 1],
        [0, 0, 1, 1, 0, 1],
        [0, 2, 0, 0, 2, 0],
        [0, 1, 0, 1, 1, 0],
        [1, 0, 1, 0, 2, 1],
        [0, 0, 2, 0, 1, 1],
        [1, 0, 0, 1, 2, 0],
        [0, 1, 1, 0, 2, 0],
        [3, 0, 0, 0, 2, 1],
        [2, 0, 1, 0, 1, 1],
        [2, 1, 0, 0, 2, 0],
        [2, 0, 0, 1, 1, 0],
        [2, 1, 0, 0, 1, 1],
        [2, 0, 0, 1, 0, 1],
        [0, 2, 0, 0, 1, 1],
        [0, 0, 0, 1, 2, 1],
        [2, 0, 1, 0, 2, 0],
        [1, 0, 0, 1, 2, 1],
        [0, 1, 1, 0, 2, 1],
        [0, 0, 1, 1, 1, 1],
        [0, 1, 0, 1, 2, 0],
        [2, 1, 0, 0, 2, 1],
        [2, 0, 0, 1, 1, 1],
        [1, 1, 1, 0, 1, 1],
        [1, 0, 1, 1, 0, 1],
        [1, 2, 0, 0, 2, 0],
        [1, 1, 0, 1, 1, 0],
        [0, 2, 0, 0, 2, 1],
        [0, 1, 0, 1, 1, 1],
        [2, 0, 1, 0, 2, 1],
        [1, 0, 2, 0, 1, 1],
        [2, 0, 0, 1, 2, 0],
        [1, 1, 1, 0, 2, 0],
        [1, 2, 0, 0, 1, 1],
        [0, 0, 2, 0, 2, 1],
        [4, 0, 0, 0, 2, 1],
        [2, 0, 0, 1, 2, 1],
        [1, 1, 1, 0, 2, 1],
        [1, 0, 1, 1, 1, 1],
        [0, 1, 2, 0, 1, 1],
        [1, 1, 0, 1, 2, 0],
        [0, 2, 1, 0, 2, 0],
        [1, 2, 0, 0, 2, 1],
        [1, 1, 0, 1, 1, 1],
        [0, 2, 1, 0, 1, 1],
        [0, 3, 0, 0, 2, 0],
        [0, 0, 1, 1, 2, 1],
        [0, 1, 0, 1, 2, 1],
        [1, 0, 2, 0, 2, 1],
        [3, 0, 1, 0, 2, 1],
        [3, 0, 0, 1, 2, 0],
        [3, 1, 0, 0, 2, 1],
        [3, 0, 0, 1, 1, 1],
        [1, 1, 0, 1, 2, 1],
        [0, 2, 1, 0, 2, 1],
        [0, 1, 1, 1, 1, 1],
        [0, 2, 0, 1, 2, 0],
        [1, 0, 1, 1, 2, 1],
        [0, 1, 2, 0, 2, 1],
        [3, 0, 0, 1, 2, 1],
        [2, 1, 1, 0, 2, 1],
        [2, 0, 1, 1, 1, 1],
        [2, 1, 0, 1, 2, 0],
        [2, 2, 0, 0, 2, 1],
        [2, 1, 0, 1, 1, 1],
        [0, 3, 0, 0, 2, 1],
        [2, 0, 2, 0, 2, 1],
        [0, 1, 1, 1, 2, 1],
        [2, 1, 0, 1, 2, 1],
        [1, 2, 1, 0, 2, 1],
        [1, 1, 1, 1, 1, 1],
        [1, 2, 0, 1, 2, 0],
        [0, 2, 0, 1, 2, 1],
        [2, 0, 1, 1, 2, 1],
        [1, 1, 2, 0, 2, 1],
        [1, 3, 0, 0, 2, 1],
        [4, 0, 0, 1, 2, 1],
        [1, 1, 1, 1, 2, 1],
        [0, 2, 2, 0, 2, 1],
        [1, 2, 0, 1, 2, 1],
        [0, 3, 1, 0, 2, 1],
        [3, 0, 1, 1, 2, 1],
        [3, 1, 0, 1, 2, 1],
        [0, 2, 1, 1, 2, 1],
        [2, 1, 1, 1, 2, 1],
        [2, 2, 0, 1, 2, 1],
        [1, 2, 1, 1, 2, 1],
      ],
    )

    basis_string = basis_lie_highest_weight_string(:A, 3, [2, 1, 1], [2, 3, 1, 2, 3, 1])

    @test issetequal(
      [only(exponents(m)) for m in monomials(basis_string)],
      [
        [0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [1, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0],
        [0, 0, 1, 1, 0, 0],
        [1, 1, 0, 0, 0, 0],
        [0, 1, 0, 1, 0, 0],
        [0, 0, 2, 0, 0, 0],
        [0, 1, 1, 0, 0, 0],
        [1, 1, 1, 0, 0, 0],
        [0, 1, 1, 1, 0, 0],
        [0, 1, 0, 1, 0, 1],
        [0, 0, 1, 1, 1, 0],
        [1, 0, 2, 0, 0, 0],
        [0, 0, 2, 1, 0, 0],
        [2, 0, 1, 0, 0, 0],
        [2, 1, 0, 0, 0, 0],
        [0, 2, 0, 1, 0, 0],
        [0, 1, 2, 0, 0, 0],
        [1, 1, 2, 0, 0, 0],
        [0, 1, 2, 1, 0, 0],
        [0, 1, 1, 1, 0, 1],
        [0, 0, 2, 1, 1, 0],
        [2, 1, 1, 0, 0, 0],
        [1, 1, 1, 1, 0, 0],
        [1, 1, 0, 1, 0, 1],
        [1, 0, 1, 1, 1, 0],
        [0, 2, 1, 1, 0, 0],
        [0, 2, 0, 1, 0, 1],
        [2, 0, 2, 0, 0, 0],
        [1, 0, 2, 1, 0, 0],
        [1, 2, 0, 1, 0, 0],
        [0, 0, 3, 1, 0, 0],
        [2, 1, 2, 0, 0, 0],
        [1, 1, 2, 1, 0, 0],
        [1, 1, 1, 1, 0, 1],
        [1, 0, 2, 1, 1, 0],
        [0, 1, 1, 2, 0, 1],
        [0, 0, 2, 2, 1, 0],
        [1, 2, 1, 1, 0, 0],
        [1, 2, 0, 1, 0, 1],
        [0, 2, 0, 2, 0, 1],
        [0, 1, 1, 2, 1, 0],
        [0, 1, 3, 1, 0, 0],
        [0, 0, 3, 1, 1, 0],
        [0, 2, 2, 1, 0, 0],
        [0, 2, 1, 1, 0, 1],
        [1, 0, 3, 1, 0, 0],
        [3, 0, 2, 0, 0, 0],
        [3, 1, 1, 0, 0, 0],
        [1, 2, 2, 1, 0, 0],
        [1, 2, 1, 1, 0, 1],
        [0, 2, 1, 2, 0, 1],
        [0, 2, 0, 2, 0, 2],
        [0, 1, 2, 2, 1, 0],
        [0, 1, 1, 2, 1, 1],
        [1, 1, 3, 1, 0, 0],
        [1, 0, 3, 1, 1, 0],
        [0, 1, 2, 2, 0, 1],
        [0, 0, 3, 2, 1, 0],
        [3, 1, 2, 0, 0, 0],
        [2, 1, 2, 1, 0, 0],
        [2, 1, 1, 1, 0, 1],
        [2, 0, 2, 1, 1, 0],
        [2, 2, 1, 1, 0, 0],
        [2, 2, 0, 1, 0, 1],
        [0, 3, 0, 2, 0, 1],
        [0, 2, 3, 1, 0, 0],
        [2, 0, 3, 1, 0, 0],
        [1, 2, 3, 1, 0, 0],
        [0, 2, 2, 2, 0, 1],
        [0, 1, 3, 2, 1, 0],
        [0, 1, 2, 2, 1, 1],
        [2, 2, 2, 1, 0, 0],
        [2, 2, 1, 1, 0, 1],
        [1, 2, 1, 2, 0, 1],
        [1, 2, 0, 2, 0, 2],
        [1, 1, 2, 2, 1, 0],
        [1, 1, 1, 2, 1, 1],
        [0, 3, 1, 2, 0, 1],
        [0, 3, 0, 2, 0, 2],
        [2, 1, 3, 1, 0, 0],
        [2, 0, 3, 1, 1, 0],
        [1, 1, 2, 2, 0, 1],
        [1, 0, 3, 2, 1, 0],
        [1, 3, 0, 2, 0, 1],
        [0, 0, 4, 2, 1, 0],
        [4, 1, 2, 0, 0, 0],
        [2, 2, 3, 1, 0, 0],
        [1, 2, 2, 2, 0, 1],
        [1, 1, 3, 2, 1, 0],
        [1, 1, 2, 2, 1, 1],
        [0, 2, 1, 3, 0, 2],
        [0, 1, 2, 3, 1, 1],
        [1, 3, 1, 2, 0, 1],
        [1, 3, 0, 2, 0, 2],
        [0, 3, 0, 3, 0, 2],
        [0, 2, 1, 3, 1, 1],
        [0, 1, 4, 2, 1, 0],
        [0, 3, 2, 2, 0, 1],
        [1, 0, 4, 2, 1, 0],
        [3, 1, 3, 1, 0, 0],
        [3, 0, 3, 1, 1, 0],
        [3, 2, 2, 1, 0, 0],
        [3, 2, 1, 1, 0, 1],
        [1, 3, 2, 2, 0, 1],
        [0, 3, 1, 3, 0, 2],
        [0, 2, 2, 3, 1, 1],
        [0, 2, 1, 3, 1, 2],
        [1, 1, 4, 2, 1, 0],
        [0, 1, 3, 3, 1, 1],
        [3, 2, 3, 1, 0, 0],
        [2, 2, 2, 2, 0, 1],
        [2, 1, 3, 2, 1, 0],
        [2, 1, 2, 2, 1, 1],
        [2, 3, 1, 2, 0, 1],
        [2, 3, 0, 2, 0, 2],
        [0, 4, 0, 3, 0, 2],
        [2, 0, 4, 2, 1, 0],
        [0, 2, 3, 3, 1, 1],
        [2, 3, 2, 2, 0, 1],
        [1, 3, 1, 3, 0, 2],
        [1, 2, 2, 3, 1, 1],
        [1, 2, 1, 3, 1, 2],
        [0, 4, 1, 3, 0, 2],
        [2, 1, 4, 2, 1, 0],
        [1, 1, 3, 3, 1, 1],
        [1, 4, 0, 3, 0, 2],
        [4, 2, 3, 1, 0, 0],
        [1, 2, 3, 3, 1, 1],
        [0, 2, 2, 4, 1, 2],
        [1, 4, 1, 3, 0, 2],
        [0, 3, 1, 4, 1, 2],
        [3, 1, 4, 2, 1, 0],
        [3, 3, 2, 2, 0, 1],
        [0, 3, 2, 4, 1, 2],
        [2, 2, 3, 3, 1, 1],
        [2, 4, 1, 3, 0, 2],
        [1, 3, 2, 4, 1, 2],
      ],
    )

    basis_nz = basis_lie_highest_weight_nz(:A, 3, [2, 1, 1], [2, 3, 1, 2, 3, 1])

    @test issetequal(
      [only(exponents(m)) for m in monomials(basis_nz)],
      [
        [0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 1],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 1, 0, 1],
        [0, 0, 1, 1, 0, 0],
        [0, 0, 0, 1, 1, 0],
        [0, 1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0, 2],
        [0, 0, 0, 0, 1, 1],
        [0, 0, 0, 1, 1, 1],
        [0, 1, 0, 1, 0, 1],
        [0, 0, 1, 1, 1, 0],
        [0, 1, 1, 1, 0, 0],
        [0, 0, 0, 1, 0, 2],
        [0, 0, 1, 1, 0, 1],
        [0, 0, 0, 2, 0, 1],
        [0, 0, 0, 2, 1, 0],
        [0, 1, 0, 1, 1, 0],
        [0, 0, 0, 0, 1, 2],
        [0, 0, 0, 1, 1, 2],
        [0, 1, 0, 1, 0, 2],
        [0, 0, 1, 1, 1, 1],
        [0, 1, 1, 1, 0, 1],
        [0, 0, 0, 2, 1, 1],
        [0, 1, 0, 2, 0, 1],
        [0, 0, 1, 2, 1, 0],
        [1, 1, 1, 1, 0, 0],
        [0, 1, 0, 1, 1, 1],
        [0, 1, 1, 1, 1, 0],
        [0, 0, 0, 2, 0, 2],
        [0, 0, 1, 2, 0, 1],
        [0, 1, 0, 2, 1, 0],
        [0, 0, 1, 1, 0, 2],
        [0, 0, 0, 2, 1, 2],
        [0, 1, 0, 2, 0, 2],
        [0, 0, 1, 2, 1, 1],
        [0, 1, 1, 2, 0, 1],
        [1, 1, 1, 1, 0, 1],
        [0, 0, 2, 2, 1, 0],
        [0, 1, 0, 2, 1, 1],
        [0, 2, 0, 2, 0, 1],
        [0, 1, 1, 2, 1, 0],
        [1, 1, 1, 1, 1, 0],
        [0, 0, 1, 1, 1, 2],
        [0, 1, 1, 1, 0, 2],
        [0, 1, 0, 1, 1, 2],
        [0, 1, 1, 1, 1, 1],
        [0, 0, 1, 2, 0, 2],
        [0, 0, 0, 3, 0, 2],
        [0, 0, 0, 3, 1, 1],
        [0, 1, 0, 2, 1, 2],
        [0, 2, 0, 2, 0, 2],
        [0, 1, 1, 2, 1, 1],
        [1, 1, 1, 1, 1, 1],
        [0, 2, 1, 2, 0, 1],
        [0, 1, 2, 2, 1, 0],
        [0, 0, 1, 2, 1, 2],
        [0, 1, 1, 2, 0, 2],
        [1, 1, 1, 1, 0, 2],
        [0, 0, 2, 2, 1, 1],
        [0, 0, 0, 3, 1, 2],
        [0, 1, 0, 3, 0, 2],
        [0, 0, 1, 3, 1, 1],
        [1, 1, 1, 2, 0, 1],
        [0, 1, 0, 3, 1, 1],
        [1, 1, 1, 2, 1, 0],
        [0, 2, 0, 2, 1, 1],
        [0, 1, 1, 1, 1, 2],
        [0, 0, 1, 3, 0, 2],
        [0, 1, 1, 2, 1, 2],
        [1, 1, 1, 1, 1, 2],
        [0, 2, 1, 2, 0, 2],
        [0, 1, 2, 2, 1, 1],
        [0, 1, 0, 3, 1, 2],
        [0, 2, 0, 3, 0, 2],
        [0, 1, 1, 3, 1, 1],
        [1, 1, 1, 2, 1, 1],
        [1, 2, 1, 2, 0, 1],
        [1, 1, 2, 2, 1, 0],
        [0, 2, 0, 2, 1, 2],
        [0, 2, 1, 2, 1, 1],
        [0, 0, 1, 3, 1, 2],
        [0, 1, 1, 3, 0, 2],
        [1, 1, 1, 2, 0, 2],
        [0, 0, 2, 3, 1, 1],
        [0, 2, 0, 3, 1, 1],
        [0, 0, 2, 2, 1, 2],
        [0, 0, 0, 4, 1, 2],
        [0, 1, 1, 3, 1, 2],
        [1, 1, 1, 2, 1, 2],
        [0, 2, 1, 3, 0, 2],
        [1, 2, 1, 2, 0, 2],
        [0, 1, 2, 3, 1, 1],
        [1, 1, 2, 2, 1, 1],
        [0, 2, 0, 3, 1, 2],
        [0, 3, 0, 3, 0, 2],
        [0, 2, 1, 3, 1, 1],
        [1, 2, 1, 2, 1, 1],
        [0, 1, 2, 2, 1, 2],
        [0, 2, 1, 2, 1, 2],
        [0, 0, 2, 3, 1, 2],
        [0, 0, 1, 4, 1, 2],
        [1, 1, 1, 3, 0, 2],
        [0, 1, 0, 4, 1, 2],
        [1, 1, 1, 3, 1, 1],
        [0, 2, 1, 3, 1, 2],
        [1, 2, 1, 2, 1, 2],
        [0, 3, 1, 3, 0, 2],
        [0, 2, 2, 3, 1, 1],
        [0, 1, 2, 3, 1, 2],
        [1, 1, 2, 2, 1, 2],
        [0, 1, 1, 4, 1, 2],
        [1, 1, 1, 3, 1, 2],
        [1, 2, 1, 3, 0, 2],
        [1, 1, 2, 3, 1, 1],
        [0, 2, 0, 4, 1, 2],
        [1, 2, 1, 3, 1, 1],
        [0, 3, 0, 3, 1, 2],
        [0, 0, 2, 4, 1, 2],
        [0, 2, 2, 3, 1, 2],
        [0, 2, 1, 4, 1, 2],
        [1, 2, 1, 3, 1, 2],
        [1, 3, 1, 3, 0, 2],
        [1, 2, 2, 3, 1, 1],
        [0, 3, 1, 3, 1, 2],
        [0, 1, 2, 4, 1, 2],
        [1, 1, 2, 3, 1, 2],
        [0, 3, 0, 4, 1, 2],
        [1, 1, 1, 4, 1, 2],
        [0, 2, 2, 4, 1, 2],
        [1, 2, 2, 3, 1, 2],
        [0, 3, 1, 4, 1, 2],
        [1, 3, 1, 3, 1, 2],
        [1, 1, 2, 4, 1, 2],
        [1, 2, 1, 4, 1, 2],
        [0, 3, 2, 4, 1, 2],
        [1, 2, 2, 4, 1, 2],
        [1, 3, 1, 4, 1, 2],
        [1, 3, 2, 4, 1, 2],
      ],
    )
  end

  @testset "Check dimension" begin
    @testset "Monomial order $monomial_ordering" for monomial_ordering in
                                                     (:lex, :invlex, :degrevlex)
      check_dimension(:A, 3, [1, 1, 1], monomial_ordering)
      #check_dimension(:B, 3, [2,1,0], monomial_ordering)
      #check_dimension(:C, 3, [1,1,1], monomial_ordering)
      #check_dimension(:D, 4, [3,0,1,1], monomial_ordering)
      #check_dimension(:F, 4, [2,0,1,0], monomial_ordering)
      #check_dimension(:G, 2, [1,0], monomial_ordering)
      #check_dimension(:G, 2, [2,2], monomial_ordering)
    end
  end
end

@testset "Coordinate ring of Kodaira embedding" begin
  @testset "general case" begin
    mbs = @inferred basis_coordinate_ring_kodaira(
      :B, 3, [0, 0, 1], 4, [3, 2, 3, 2, 1, 2, 3, 2, 1]; monomial_ordering=:neglex
    )
    @test length(mbs) == 4

    @test dim(mbs[1][1]) == length(mbs[1][2])
    @test issetequal(monomials(mbs[1][1]), mbs[1][2])

    @test length(mbs[2][2]) > 0
    @test issubset(monomials(mbs[1][1]), monomials(mbs[2][1]))
    @test issubset(mbs[2][2], monomials(mbs[2][1]))

    @test isempty(mbs[3][2])
    @test issubset(monomials(mbs[2][1]), monomials(mbs[3][1]))

    @test isempty(mbs[4][2])
    @test issubset(monomials(mbs[3][1]), monomials(mbs[4][1]))
  end

  @testset "FFL" begin
    for case in [
      (:A, 3, [1, 0, 1], 4),
      (:B, 2, [1, 1], 4),
      (:D, 4, [1, 0, 1, 0], 4),
      (:G, 2, [1, 0], 6),
    ]
      dynkin, n, lambda, degree = case
      mbs = @inferred basis_coordinate_ring_kodaira_ffl(dynkin, n, lambda, degree)
      L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(dynkin), n, GAP.Globals.Rationals)
      gap_dim = GAP.Globals.DimensionOfHighestWeightModule(L, GAP.Obj(lambda))
      @test length(mbs) == degree
      @test dim(mbs[1][1]) == gap_dim
      @test dim(mbs[1][1]) == length(mbs[1][2])
      @test issetequal(monomials(mbs[1][1]), mbs[1][2])
      for i in 2:degree
        @test isempty(mbs[i][2])
      end
    end
  end
end

@testset "Demazure" begin
  @testset "Trivial Cases" begin
    for (type, rank, highest_weight) in [
      (:A, 3, [1, 0, 1]),
      (:A, 3, [2, 2, 2]),
      (:B, 2, [1, 1]),
      (:D, 4, [1, 0, 1, 0]),
      (:G, 2, [1, 0]),
    ], monomial_ordering in [:degrevlex, :lex, :invlex, :neglex, :wdegrevlex]

      # longest weyl group element leads to the same result as a simple module
      mb1 = basis_lie_highest_weight(type, rank, highest_weight; monomial_ordering)
      mb2 = basis_lie_demazure(type, rank, highest_weight, letters(longest_element(weyl_group(type, rank))); monomial_ordering)
      @test monomials(mb1) == monomials(mb2)

      # empty weyl group element leads to a monomialbasis with only the one element
      mb = basis_lie_demazure(type, rank, highest_weight, Int[]; monomial_ordering)
      @test monomials(mb) == Set{ZZMPolyRingElem}([one(mb.monomials_parent)])
    end
  end

  @testset "Nontrivial Cases" begin
    for (type, rank, highest_weight, weyl_group_elem) in [
      (:A, 3, [1, 0, 1], [1, 2, 1]),
      (:A, 3, [2, 2, 2], [1, 2, 1]),
      (:A, 3, [1, 2, 0], [1, 2, 1]), #example below
      #(:B, 2, [1, 0], [1, 2]), #this fails currently because operators acting upwards are missing
      #(:C, 2, [1, 1], [1, 2, 1]), #this fails currently because operators acting upwards are missing
      (:D, 4, [1, 0, 1, 0], [1, 2, 1]),
      #(:G, 2, [1, 0], [1, 2, 1]), #this fails currently because operators acting upwards are missing
    ], monomial_ordering in [:degrevlex, :lex, :invlex, :neglex, :wdegrevlex]

      mb = basis_lie_demazure(type, rank, highest_weight, weyl_group_elem; monomial_ordering)
      R = root_system(birational_sequence(mb))

      char = demazure_character(R, highest_weight, weyl_group_elem)
      
      dict = Dict{WeightLatticeElem,Int64}()
      for mon in monomials(mb)
        w = WeightLatticeElem(R, highest_weight) - Oscar.BasisLieHighestWeight.weight(mon, birational_sequence(mb))
        val = get(dict, w, 0) + 1
        dict[w] = val
      end
      @test dict == char
    end
  end

  @testset "Hand Computed Case" begin
    # Test for a specific case in which dependent on the order there are two out of three monomials in a weight space
    type = :A
    rank = 3
    highest_weight = [1, 2, 0]
    weyl_group_elem = [1, 2, 1]   

    for (monomial_ordering, expected_result) in [
      (:degrevlex, Set{Vector{Int64}}([[0, 0, 2], [1, 1, 1],])),
      (:lex, Set{Vector{Int64}}([[0, 0, 2], [1, 1, 1],])),
      (:invlex, Set{Vector{Int64}}([[2, 2, 0], [1, 1, 1],])),
      (:neglex, Set{Vector{Int64}}([[2, 2, 0], [1, 1, 1],])),
      (:wdegrevlex, Set{Vector{Int64}}([[0, 0, 2], [1, 1, 1],])),
    ]
      mb = basis_lie_demazure(type, rank, highest_weight, weyl_group_elem; monomial_ordering)
    
      R = root_system(birational_sequence(mb))
      w = WeightLatticeElem(R, highest_weight) - WeightLatticeElem(2 * simple_root(R, 1) - 2 * simple_root(R, 2))

      monomials_for_weight_w = filter(mon -> WeightLatticeElem(R, highest_weight) - Oscar.BasisLieHighestWeight.weight(mon, birational_sequence(mb)) == w, monomials(mb))

      exponent_vectors = Set{Vector{Int64}}([only(exponents(m)) for m in monomials_for_weight_w])
    
      @test exponent_vectors == expected_result
    end
  end
end
