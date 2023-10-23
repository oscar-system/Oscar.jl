
function basis_lie_highest_weight_compute(
  type::Symbol,
  rank::Int,
  highest_weight::Vector{Int},
  get_operators::Function,
  monomial_order::Union{String,Function},
)
  """
  Pseudocode:

  basis_lie_highest_weight(highest_weight)
      return compute_monomials(highest_weight)

  compute_monomials(highest_weight)
      if highest_weight was already computed 
          return old results
      if highest_weight = [0, ..., 0] or [0, ..., 1, ..., 0]
          return add_by_hand(highest_weight, {})
      else
          set_mon = {}
          go through all partitions lambda_1 + lambda_2 = highest_weight
              add compute_monomials(lambda_1) (+) compute_monomials(lambda_1) to set_mon 
          if set_mon too small
              add_by_hand(highest_weight, set_mon)
          return set_mon

  add_by_hand(highest_weight, set_mon)
      add_known_monomials(set_mon)
      go through all weightspaces that are not full
          add_new_monomials(weightspace, set_mon)
      return set_mon
    
  add_known_monomials(set_mon)
      add all monomials from set_mon to basis

  add_new_monomials(weightspace, set_mon)
      calculate monomials with weight in weightspace
      go through them one by one in monomial_order until basis is full
      return set_mon
  """
  highest_weight = convert(Vector{ZZRingElem}, highest_weight)
  # The function precomputes objects that are independent of the highest weight and that can be used in all recursion 
  # steps. Then it starts the recursion and returns the result.

  # initialization of objects that can be precomputed
  # lie_algebra of type, rank and its chevalley_basis
  lie_algebra = LieAlgebraStructure(type, rank)
  chevalley_basis = NTuple{3,Vector{GAP.Obj}}(
    GAP.Globals.ChevalleyBasis(lie_algebra.lie_algebra_gap)
  )

  # operators that are represented by our monomials. x_i is connected to operators[i]
  operators = get_operators(lie_algebra, chevalley_basis)

  weights_w = weights_for_operators(
    lie_algebra.lie_algebra_gap, chevalley_basis[3], operators
  ) # weights of the operators
  # weights_w = (weight_w->Int.(weight_w)).(weights_w)  
  weights_eps = [
    w_to_eps(type, rank, convert(Vector{QQFieldElem}, weight_w)) for weight_w in weights_w
  ] # other root system
  weights_alpha = [
    w_to_alpha(type, rank, convert(Vector{QQFieldElem}, weight_w)) for weight_w in weights_w
  ] # other root system

  asVec(v) = Oscar.GAP.gap_to_julia(GAPWrap.ExtRepOfObj(v)) # TODO
  birational_sequence = BirationalSequence(
    operators, [asVec(v) for v in operators], weights_w, weights_eps, weights_alpha
  )

  ZZx, _ = PolynomialRing(ZZ, length(operators)) # for our monomials
  monomial_order_lt = get_monomial_order_lt(monomial_order, ZZx) # less than function to sort monomials by order

  # save computations from recursions
  calc_highest_weight = Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}}(
    [ZZ(0) for i in 1:rank] => Set([ZZx(1)])
  )
  # save all highest weights, for which the Minkowski-sum did not suffice to gain all monomials
  no_minkowski = Set{Vector{ZZRingElem}}()

  # start recursion over highest_weight
  set_mon = compute_monomials(
    lie_algebra,
    birational_sequence,
    ZZx,
    highest_weight,
    monomial_order_lt,
    calc_highest_weight,
    no_minkowski,
  )

  monomials = sort(collect(set_mon); lt=monomial_order_lt)

  minkowski_gens = map(
    gen -> Int.(gen), sort(collect(no_minkowski); by=(gen -> (sum(gen), reverse(gen))))
  )

  # output
  return MonomialBasis(
    lie_algebra,
    highest_weight,
    monomial_order,
    monomials,
    minkowski_gens,
    birational_sequence,
  )
end

@doc """
basis_lie_highest_weight(
    type::Symbol,
    rank::Int, 
    highest_weight::Vector{Int};
    reduced_expression::Union{String, Vector{Union{Int, Vector{Int}}}} = "regular", 
    monomial_order::Union{String, Function} = "GRevLex", 
)

Computes a monomial basis for the highest weight module with highest weight
``highest_weight`` (in terms of the fundamental weights), for a simple Lie algebra of type
``type`` and rank ``rank``.

# Parameters
- `type`: type of liealgebra we want to investigate, one of :A, :B, :C, :D, :E, :F, :G
- `rank`: rank of liealgebra
- `highest_weight`: highest-weight
- `reduced_expression`: list of operators, either "regular" or integer array. The functionality of choosing a random longest word
                is currently not implemented, because we used https://github.com/jmichel7/Gapjm.jl to work with coxeter 
                groups need a method to obtain all non left descending elements to extend a word
- `monomial_order`: monomial order in which our basis gets defined with regards to our operators 

# Examples
```jldoctest
julia> base = BasisLieHighestWeight.basis_lie_highest_weight(:A, 2, [1, 1])
Monomial basis of a highest weight module
  of highest weight [1, 1]
  of dimension 8
  with monomial ordering degrevlex
over lie-Algebra of type A and rank 2
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [1, 0]
    [0, 1]
    [1, 1]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0]
    [0, 1]

julia> base = BasisLieHighestWeight.basis_lie_highest_weight(:A, 3, [2, 2, 3]; monomial_order = "lex")
Monomial basis of a highest weight module
  of highest weight [2, 2, 3]
  of dimension 1260
  with monomial ordering lex
over lie-Algebra of type A and rank 3
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
    [1, 1, 0]
    [0, 1, 1]
    [1, 1, 1]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]

julia> base = BasisLieHighestWeight.basis_lie_highest_weight(:A, 2, [1, 0]; reduced_expression = [1,2,1])
Monomial basis of a highest weight module
  of highest weight [1, 0]
  of dimension 3
  with monomial ordering degrevlex
over lie-Algebra of type A and rank 2
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [1, 0]
    [0, 1]
    [1, 0]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0]

julia> base = BasisLieHighestWeight.basis_lie_highest_weight(:C, 3, [1, 1, 1]; monomial_order = "lex")
Monomial basis of a highest weight module
  of highest weight [1, 1, 1]
  of dimension 512
  with monomial ordering lex
over lie-Algebra of type C and rank 3
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
    [1, 1, 0]
    [0, 1, 1]
    [1, 1, 1]
    [0, 2, 1]
    [1, 2, 1]
    [2, 2, 1]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
    [0, 1, 1]
    [1, 1, 1]
```
"""
function basis_lie_highest_weight(
  type::Symbol,
  rank::Int,
  highest_weight::Vector{Int};
  reduced_expression::Union{String,Vector{Int},Vector{GAP.GapObj},Any}="regular",
  monomial_order::Union{String,Function}="degrevlex",
)
  """
  Standard function with all options
  """
  get_operators =
    (lie_algebra, chevalley_basis) ->
      get_operators_normal(lie_algebra, chevalley_basis, reduced_expression)
  return basis_lie_highest_weight_compute(
    type, rank, highest_weight, get_operators, monomial_order
  )
end

@doc """
# Examples
```jldoctest
julia> base = BasisLieHighestWeight.basis_lie_highest_weight_lustzig(:D, 4, [1,1,1,1], [4,3,2,4,3,2,1,2,4,3,2,1])
Monomial basis of a highest weight module
  of highest weight [1, 1, 1, 1]
  of dimension 4096
  with monomial ordering oplex
over lie-Algebra of type D and rank 4
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [0, 0, 0, 1]
    [0, 0, 1, 0]
    [0, 1, 1, 1]
    [0, 1, 1, 0]
    [0, 1, 0, 1]
    [0, 1, 0, 0]
    [1, 2, 1, 1]
    [1, 1, 1, 1]
    [1, 1, 0, 1]
    [1, 1, 1, 0]
    [1, 1, 0, 0]
    [1, 0, 0, 0]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0, 0, 0]
    [0, 1, 0, 0]
    [0, 0, 1, 0]
    [0, 0, 0, 1]
    [0, 0, 1, 1]
```
"""
function basis_lie_highest_weight_lustzig(
  type::Symbol,
  rank::Int,
  highest_weight::Vector{Int},
  reduced_expression::Vector{Int};
  monomial_order::Union{String,Function}="oplex",
)
  """
  Lustzig polytope
  BasisLieHighestWeight.basis_lie_highest_weight_lustzig(:D, 4, [1,1,1,1], [4,3,2,4,3,2,1,2,4,3,2,1])
  """
  # operators = some sequence of the String / Littelmann-Berenstein-Zelevinsky polytope
  get_operators =
    (lie_algebra, chevalley_basis) ->
      get_operators_lustzig(lie_algebra, chevalley_basis, reduced_expression)
  return basis_lie_highest_weight_compute(
    type, rank, highest_weight, get_operators, monomial_order
  )
end

@doc """
# Examples
```jldoctest
julia> BasisLieHighestWeight.basis_lie_highest_weight_string(:B, 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
Monomial basis of a highest weight module
  of highest weight [1, 1, 1]
  of dimension 512
  with monomial ordering oplex
over lie-Algebra of type B and rank 3
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [0, 0, 1]
    [0, 1, 0]
    [0, 0, 1]
    [0, 1, 0]
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
    [0, 1, 0]
    [1, 0, 0]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
```
"""
function basis_lie_highest_weight_string(
  type::Symbol, rank::Int, highest_weight::Vector{Int}, reduced_expression::Vector{Int};
)
  """
  String / Littelmann-Berenstein-Zelevinsky polytope
  BasisLieHighestWeight.basis_lie_highest_weight_string(:B, 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
  BasisLieHighestWeight.basis_lie_highest_weight_string(:B, 4, [1,1,1,1], [4,3,4,3,2,3,4,3,2,1,2,3,4,3,2,1])
  BasisLieHighestWeight.basis_lie_highest_weight_string(:A, 4, [1,1,1,1], [4,3,2,1,2,3,4,3,2,3])
  """
  # reduced_expression = some sequence of the String / Littelmann-Berenstein-Zelevinsky polytope
  monomial_order = "oplex"
  get_operators =
    (lie_algebra, chevalley_basis) ->
      get_operators_normal(lie_algebra, chevalley_basis, reduced_expression)
  return basis_lie_highest_weight_compute(
    type, rank, highest_weight, get_operators, monomial_order
  )
end

@doc """
# Examples
```jldoctest
julia> BasisLieHighestWeight.basis_lie_highest_weight_fflv(:A, 3, [1,1,1])
Monomial basis of a highest weight module
  of highest weight [1, 1, 1]
  of dimension 64
  with monomial ordering oplex
over lie-Algebra of type A and rank 3
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [1, 1, 1]
    [0, 1, 1]
    [1, 1, 0]
    [0, 0, 1]
    [0, 1, 0]
    [1, 0, 0]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
```
"""
function basis_lie_highest_weight_fflv(type::Symbol, rank::Int, highest_weight::Vector{Int})
  """
  Feigin-Fourier-Littelmann-Vinberg polytope
  BasisLieHighestWeight.basis_lie_highest_weight_fflv(:A, 3, [1,1,1])
  """
  monomial_order = "oplex"
  # operators = all positive roots, reverse ordering as GAP uses
  get_operators =
    (lie_algebra, chevalley_basis) ->
      reverse(get_operators_normal(lie_algebra, chevalley_basis, "regular"))
  return basis_lie_highest_weight_compute(
    type, rank, highest_weight, get_operators, monomial_order
  )
end

@doc """
# Examples
```jldoctest
julia> BasisLieHighestWeight.basis_lie_highest_weight_nz(:C, 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
Monomial basis of a highest weight module
  of highest weight [1, 1, 1]
  of dimension 512
  with monomial ordering lex
over lie-Algebra of type C and rank 3
  where the birational sequence used consists of operators to the following weights (given as coefficients w.r.t. alpha_i):
    [0, 0, 1]
    [0, 1, 0]
    [0, 0, 1]
    [0, 1, 0]
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
    [0, 1, 0]
    [1, 0, 0]
  and the basis was generated by Minkowski sums of the bases of the following highest weight modules:
    [1, 0, 0]
    [0, 1, 0]
    [0, 0, 1]
```
"""
function basis_lie_highest_weight_nz(
  type::Symbol, rank::Int, highest_weight::Vector{Int}, reduced_expression::Vector{Int};
)
  """
  Nakashima-Zelevinsky polytope
  BasisLieHighestWeight.basis_lie_highest_weight_nz(:C, 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
  BasisLieHighestWeight.basis_lie_highest_weight_nz(:A, 4, [1,1,1,1], [4,3,2,1,2,3,4,3,2,3])
  """
  monomial_order = "lex"
  get_operators =
    (lie_algebra, chevalley_basis) ->
      get_operators_normal(lie_algebra, chevalley_basis, reduced_expression)
  return basis_lie_highest_weight_compute(
    type, rank, highest_weight, get_operators, monomial_order
  )
end

function compute_monomials(
  lie_algebra::LieAlgebraStructure,
  birational_sequence::BirationalSequence,
  ZZx::ZZMPolyRing,
  highest_weight::Vector{ZZRingElem},
  monomial_order_lt::Function,
  calc_highest_weight::Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}},
  no_minkowski::Set{Vector{ZZRingElem}},
)::Set{ZZMPolyRingElem}
  """
  This function calculates the monomial basis M_{highest_weight} recursively. The recursion saves all computed 
  results in calc_highest_weight and we first check, if we already encountered this highest weight in a prior step. 
  If this is not the case, we need to perform computations. The recursion works by using the Minkowski-sum. 
  If M_{highest_weight} is the desired set of monomials (identified by the exponents as lattice points), it is know 
  that for lambda_1 + lambda_2 = highest_weight we have M_{lambda_1} + M_{lambda_2} subseteq M_{highest_weight}. 
  The complexity grows exponentially in the size of highest_weight. Therefore, it is very helpful to obtain a part of
  M_{highest_weight} by going through all partitions of highest_weight and using the Minkowski-property. The base 
  cases of the recursion are the fundamental weights highest_weight = [0, ..., 1, ..., 0]. In this case, or if the 
  Minkowski-property did not find enough monomials, we need to perform the computations "by hand".
  """
  # simple cases
  # we already computed the highest_weight result in a prior recursion step
  if haskey(calc_highest_weight, highest_weight)
    return calc_highest_weight[highest_weight]
  elseif highest_weight == [ZZ(0) for i in 1:(lie_algebra.rank)] # we mathematically know the solution
    return Set(ZZx(1))
  end

  # calculation required
  # gap_dim is number of monomials that we need to find, i.e. |M_{highest_weight}|.
  # if highest_weight is a fundamental weight, partition into smaller summands is possible. This is the basecase of 
  # the recursion.
  highest_weight_int = convert(Vector{Int}, highest_weight)
  gap_dim = GAP.Globals.DimensionOfHighestWeightModule(
    lie_algebra.lie_algebra_gap, GAP.Obj(highest_weight_int)
  ) # fundamental weights
  if is_fundamental(highest_weight) || sum(abs.(highest_weight)) == 0
    push!(no_minkowski, highest_weight)
    set_mon = add_by_hand(
      lie_algebra,
      birational_sequence,
      ZZx,
      highest_weight,
      monomial_order_lt,
      Set{ZZMPolyRingElem}(),
    )
    push!(calc_highest_weight, highest_weight => set_mon)
    return set_mon
  else
    # use Minkowski-Sum for recursion
    set_mon = Set{ZZMPolyRingElem}()
    i = 0
    sub_weights_w = compute_sub_weights(highest_weight)
    l = length(sub_weights_w)
    # go through all partitions lambda_1 + lambda_2 = highest_weight until we have enough monomials or used all 
    # partitions
    while length(set_mon) < gap_dim && i < l
      i += 1
      lambda_1 = sub_weights_w[i]
      lambda_2 = highest_weight .- lambda_1
      mon_lambda_1 = compute_monomials(
        lie_algebra,
        birational_sequence,
        ZZx,
        lambda_1,
        monomial_order_lt,
        calc_highest_weight,
        no_minkowski,
      )
      mon_lambda_2 = compute_monomials(
        lie_algebra,
        birational_sequence,
        ZZx,
        lambda_2,
        monomial_order_lt,
        calc_highest_weight,
        no_minkowski,
      )
      # Minkowski-sum: M_{lambda_1} + M_{lambda_2} \subseteq M_{highest_weight}, if monomials get identified with 
      # points in ZZ^n
      mon_sum = Set([p * q for p in mon_lambda_1 for q in mon_lambda_2])
      union!(set_mon, mon_sum)
    end
    # check if we found enough monomials

    #println("")
    #println("highest_weight: ", highest_weight)
    #println("required monomials: ", gap_dim)
    #println("monomials from Minkowski-sum: ", length(set_mon))
    #println(set_mon)

    if length(set_mon) < gap_dim
      push!(no_minkowski, highest_weight)
      set_mon = add_by_hand(
        lie_algebra, birational_sequence, ZZx, highest_weight, monomial_order_lt, set_mon
      )
    end
    push!(calc_highest_weight, highest_weight => set_mon)
    return set_mon
  end
end

function add_known_monomials!(
  weight_w::Vector{ZZRingElem},
  set_mon_in_weightspace::Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}},
  matrices_of_operators::Vector{SMat{ZZRingElem}},
  space::Dict{Vector{ZZRingElem},Oscar.BasisLieHighestWeight.SparseVectorSpaceBasis},
  v0::SRow{ZZRingElem},
)
  """
  By using the Minkowski-sum, we know that all monomials in set_mon_in_weightspace are in our basis. Since we want to
  extend the weightspace with missing monomials, we need to calculate and add the vector of each monomial to our 
  basis.
  """
  # println("add_known_monomials")

  for mon in set_mon_in_weightspace[weight_w]
    # calculate the vector vec associated with mon
    vec = calc_vec(v0, mon, matrices_of_operators)

    # check if vec extends the basis
    if !haskey(space, weight_w)
      space[weight_w] = SparseVectorSpaceBasis([], [])
    end
    add_and_reduce!(space[weight_w], vec)
  end
end

function add_new_monomials!(
  lie_algebra::LieAlgebraStructure,
  birational_sequence::BirationalSequence,
  ZZx::ZZMPolyRing,
  matrices_of_operators::Vector{SMat{ZZRingElem}},
  monomial_order_lt::Function,
  dim_weightspace::Int,
  weight_w::Vector{ZZRingElem},
  set_mon_in_weightspace::Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}},
  space::Dict{Vector{ZZRingElem},Oscar.BasisLieHighestWeight.SparseVectorSpaceBasis},
  v0::SRow{ZZRingElem},
  set_mon::Set{ZZMPolyRingElem},
)
  """
  If a weightspace is missing monomials, we need to calculate them by trial and error. We would like to go through all
  monomials in the order monomial_order_lt and calculate the corresponding vector. If it extends the basis, we add it 
  to the result and else we try the next one. We know, that all monomials that work lay in the weyl-polytope. 
  Therefore, we only inspect the monomials that lie both in the weyl-polytope and the weightspace. Since the weyl-
  polytope is bounded these are finitely many and we can sort them and then go trough them, until we found enough. 
  """
  #println("add_new_monomials")

  # get monomials that are in the weightspace, sorted by monomial_order_lt
  poss_mon_in_weightspace = convert_lattice_points_to_monomials(
    ZZx,
    get_lattice_points_of_weightspace(
      birational_sequence.weights_alpha,
      w_to_alpha(
        lie_algebra.lie_type, lie_algebra.rank, convert(Vector{QQFieldElem}, weight_w)
      ),
    ),
  )
  #println("before sort")
  #flush(stdout)
  poss_mon_in_weightspace = sort(poss_mon_in_weightspace; lt=monomial_order_lt)
  #println("after sort")

  # check which monomials should get added to the basis
  i = 0
  if weight_w == 0 # check if [0 0 ... 0] already in basis
    i += 1
  end
  number_mon_in_weightspace = length(set_mon_in_weightspace[weight_w])
  # go through possible monomials one by one and check if it extends the basis
  while number_mon_in_weightspace < dim_weightspace
    i += 1

    mon = poss_mon_in_weightspace[i]
    if mon in set_mon
      continue
    end

    # calculate the vector vec associated with mon
    d = sz(matrices_of_operators[1])
    vec = calc_vec(v0, mon, matrices_of_operators)

    # check if vec extends the basis
    if !haskey(space, weight_w)
      space[weight_w] = SparseVectorSpaceBasis([], [])
    end
    vec_red = add_and_reduce!(space[weight_w], vec)
    if isempty(vec_red) # v0 == 0
      continue
    end

    # save monom
    number_mon_in_weightspace += 1
    push!(set_mon, mon)
  end
end

function add_by_hand(
  lie_algebra::LieAlgebraStructure,
  birational_sequence::BirationalSequence,
  ZZx::ZZMPolyRing,
  highest_weight::Vector{ZZRingElem},
  monomial_order_lt::Function,
  set_mon::Set{ZZMPolyRingElem},
)::Set{ZZMPolyRingElem}
  #println("")
  #println("")
  #println("add_by_hand", highest_weight)

  """
  This function calculates the missing monomials by going through each non full weightspace and adding possible 
  monomials manually by computing their corresponding vectors and checking if they enlargen the basis.
  """
  # initialization
  # matrices g_i for (g_1^a_1 * ... * g_k^a_k)*v
  matrices_of_operators = tensorMatricesForOperators(
    lie_algebra.lie_algebra_gap, highest_weight, birational_sequence.operators
  )
  space = Dict(ZZ(0) * birational_sequence.weights_w[1] => SparseVectorSpaceBasis([], [])) # span of basis vectors to keep track of the basis
  v0 = sparse_row(ZZ, [(1, 1)])  # starting vector v

  push!(set_mon, ZZx(1))
  # required monomials of each weightspace
  weightspaces = get_dim_weightspace(lie_algebra, highest_weight)

  # sort the monomials from the minkowski-sum by their weightspaces
  set_mon_in_weightspace = Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}}()
  for (weight_w, _) in weightspaces
    set_mon_in_weightspace[weight_w] = Set{ZZMPolyRingElem}()
  end
  for mon in set_mon
    weight_w = calc_weight(mon, birational_sequence.weights_w)
    push!(set_mon_in_weightspace[weight_w], mon)
  end

  # only inspect weightspaces with missing monomials
  weights_with_full_weightspace = Set{Vector{ZZRingElem}}()
  for (weight_w, dim_weightspace) in weightspaces
    if (length(set_mon_in_weightspace[weight_w]) == dim_weightspace)
      push!(weights_with_full_weightspace, weight_w)
    end
  end
  delete!(weightspaces, weights_with_full_weightspace)

  # The weightspaces could be calculated completely indepent (except for
  # the caching). This is not implemented, since I used the package Distributed.jl for this, which is not in the 
  # Oscar dependencies. But I plan to reimplement this. 
  # insert known monomials into basis
  for (weight_w, _) in weightspaces
    # print(".")
    add_known_monomials!(weight_w, set_mon_in_weightspace, matrices_of_operators, space, v0)
  end

  # println("")
  # calculate new monomials
  for (weight_w, dim_weightspace) in weightspaces
    # print("*")
    add_new_monomials!(
      lie_algebra,
      birational_sequence,
      ZZx,
      matrices_of_operators,
      monomial_order_lt,
      dim_weightspace,
      weight_w,
      set_mon_in_weightspace,
      space,
      v0,
      set_mon,
    )
  end

  return set_mon
end

function sub_simple_refl(word::Vector{Int}, lie_algebra_gap::GAP.Obj)::Vector{GAP.Obj}
  """
  substitute simple reflections (i,i+1), saved in dec by i, with E_{i,i+1}  
  """
  root_system = GAP.Globals.RootSystem(lie_algebra_gap)
  canonical_generators = Vector{GAP.Obj}(
    GAP.Globals.CanonicalGenerators(root_system)[1]; recursive=false
  )
  operators = [canonical_generators[i] for i in word]
  return operators
end

function get_operators_normal(
  lie_algebra::LieAlgebraStructure,
  chevalley_basis::NTuple{3,Vector{GAP.Obj}},
  reduced_expression::Union{String,Vector{Union{Int,Vector{Int}}},Vector{GAP.GapObj},Any},
)::Vector{GAP.Obj}
  """
  handles user input for operators
  "regular" for all operators
  "longest-word" for random longest-word in Weyl-group (currently not implemented)
    reduced_expression::Vector{Int} for explicit longest-word
  """
  if typeof(reduced_expression) == GAP.Obj  # If user already submitted gap-roots as operators, keep
    return reduced_expression
  elseif reduced_expression == "regular" # create standard reduced_expression, use reduced_expression as specified by GAP
    return chevalley_basis[1]
    # The functionality longest-word required Coxetergroups from Gapjm.jl (https://github.com/jmichel7/Gapjm.jl and was 
    # temporarily deleted
    # choose a random longest word. Created by extending by random not leftdescending reflections until total length is 
    # reached
    #elseif operators == "longest-word"
    #    operators = longest_weyl_word(t,n)
    #    operators = sub_simple_refl(operators, lie_algebra, n)
    #    return operators
  end

  # use user defined operator
  # Check for incorrect input:
  for x in reduced_expression
    if isa(x, Int)
      if !(1 <= x <= lie_algebra.rank)
        error(
          "Each integer in reduced_expression should be between 1 and the rank of the lie-algebra",
        )
      end
    elseif isa(x, Vector{Int})
      if !(all(1 <= i <= lie_algebra.rank for i in x))
        error(
          "All integers in each vector of reduced_expression should be between 1 and the rank of the lie-algebra",
        )
      end
    else
      error("Each item in reduced_expression needs to be an Int or Vector{Int}")
    end
  end
  # If one of the conditions is met, the algorithms works. Otherwise a warning is printed (and can be ignored).
  #if  !(is_longest_weyl_word(type, rank, reduced_expression)) && !(Set(reduced_expression) == [i for i=1:n])
  #    println("WARNING: reduced_expression may be incorrect input.")
  #end
  sanitized_reduced_expression = Vector{Union{Int,Vector{Int}}}() # creates an empty array of the desired type
  for item in reduced_expression
    if isa(item, Int)
      push!(sanitized_reduced_expression, item)
    elseif isa(item, Vector{Int})
      push!(sanitized_reduced_expression, item)
    else
      error("Wrong type")
    end
  end
  operators = get_operators_simple_reflections(
    lie_algebra, chevalley_basis, sanitized_reduced_expression
  )
  return operators
end

@doc """
    is_fundamental(highest_weight::Vector{IntegerUnion})::Bool

    returns true if ``highest_weight`` is fundamental, i.e. [0, ..., 1, ..., 0]

# Examples
```jldoctest
julia> BasisLieHighestWeight.is_fundamental([0, 1, 0])
true

julia> BasisLieHighestWeight.is_fundamental([0, 1, 1])
false
```
"""
function is_fundamental(highest_weight::Vector{<:IntegerUnion})
  hasone = false
  for i in Int.(highest_weight)
    if i == 0
      continue
    elseif i == 1
      hasone && return false
      hasone = true
    else
      return false
    end
  end
  return hasone
end

function compute_sub_weights(highest_weight::Vector{ZZRingElem})::Vector{Vector{ZZRingElem}}
  """
  returns list of weights w != 0, highest_weight with 0 <= w <= highest_weight elementwise, ordered by l_2-norm
  """
  sub_weights_w = []
  foreach(Iterators.product((0:x for x in highest_weight)...)) do i
    push!(sub_weights_w, [i...])
  end
  if isempty(sub_weights_w) || length(sub_weights_w) == 1 # case [] or [[0, ..., 0]]
    return []
  else
    popfirst!(sub_weights_w) # [0, ..., 0]
    pop!(sub_weights_w) # highest_weight
    sort!(sub_weights_w; by=x -> sum((x) .^ 2))
    return sub_weights_w
  end
end
