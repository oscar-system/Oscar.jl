function basis_lie_highest_weight_compute(
  L::LieAlgebraStructure,
  highest_weight::Vector{Int},
  operators::Vector{<:GAP.Obj},     # operators are represented by our monomials. x_i is connected to operators[i]
  monomial_ordering_symb::Symbol,
)
  # Pseudocode:

  # basis_lie_highest_weight(highest_weight)
  #     return compute_monomials(highest_weight)

  # compute_monomials(highest_weight)
  #     if highest_weight was already computed 
  #         return old results
  #     if highest_weight = [0, ..., 0] or [0, ..., 1, ..., 0]
  #         return add_by_hand(highest_weight, {})
  #     else
  #         set_mon = {}
  #         go through all partitions lambda_1 + lambda_2 = highest_weight
  #             add compute_monomials(lambda_1) (+) compute_monomials(lambda_1) to set_mon 
  #         if set_mon too small
  #             add_by_hand(highest_weight, set_mon)
  #         return set_mon

  # add_by_hand(highest_weight, set_mon)
  #     add_known_monomials(set_mon)
  #     go through all weightspaces that are not full
  #         add_new_monomials(weightspace, set_mon)
  #     return set_mon

  # add_known_monomials(set_mon)
  #     add all monomials from set_mon to basis

  # add_new_monomials(weightspace, set_mon)
  #     calculate monomials with weight in weightspace
  #     go through them one by one in monomial_ordering until basis is full
  #     return set_mon

  highest_weight = ZZ.(highest_weight)
  # The function precomputes objects that are independent of the highest weight and that can be used in all recursion 
  # steps. Then it starts the recursion and returns the result.

  weights_w = [weight(L, op) for op in operators] # weights of the operators
  weights_alpha = [w_to_alpha(L, weight_w) for weight_w in weights_w] # other root system

  asVec(v) = GAP.gap_to_julia(GAPWrap.ExtRepOfObj(v)) # TODO
  birational_sequence = BirationalSequence(
    operators, [asVec(v) for v in operators], weights_w, weights_alpha
  )

  ZZx, _ = polynomial_ring(ZZ, length(operators)) # for our monomials
  monomial_ordering = get_monomial_ordering(monomial_ordering_symb, ZZx, weights_alpha)

  # save computations from recursions
  calc_highest_weight = Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}}(
    [ZZ(0) for i in 1:rank(L)] => Set([ZZx(1)])
  )
  # save all highest weights, for which the Minkowski-sum did not suffice to gain all monomials
  no_minkowski = Set{Vector{ZZRingElem}}()

  # start recursion over highest_weight
  set_mon = compute_monomials(
    L,
    birational_sequence,
    ZZx,
    highest_weight,
    monomial_ordering,
    calc_highest_weight,
    no_minkowski,
  )
  # monomials = sort(collect(set_mon); lt=((m1, m2) -> cmp(monomial_ordering, m1, m2) < 0))
  minkowski_gens = sort(collect(no_minkowski); by=(gen -> (sum(gen), reverse(gen))))
  # output
  mb = MonomialBasis(L, highest_weight, birational_sequence, monomial_ordering, set_mon)
  set_attribute!(
    mb, :algorithm => basis_lie_highest_weight_compute, :minkowski_gens => minkowski_gens
  )
  return mb
end

function basis_coordinate_ring_kodaira_compute(
  L::LieAlgebraStructure,
  highest_weight::Vector{Int},
  degree::Int,
  operators::Vector{<:GAP.Obj},     # operators are represented by our monomials. x_i is connected to operators[i]
  monomial_ordering_symb::Symbol,
)
  # Pseudocode:

  # basis_coordinate_ring_kodaira_compute(highest_weight, degree)

  # returns all multiples of the given highest_weight up to degree 
  # such that for this ordering the monomial basis is not a Minkowski sum of bases of smaller multiples.
  #     return monomial bases for highest_weight and for each multiple up to degree, that is not the Minkowski sum 
  #       of smaller multiples, the missing monomials

  @req degree > 0 "Degree must be positive"
  highest_weight = ZZ.(highest_weight)
  # The function precomputes objects that are independent of the highest weight and that can be used in all recursion 
  # steps. Then it starts the recursion and returns the result.

  weights_w = [weight(L, op) for op in operators] # weights of the operators
  weights_alpha = [w_to_alpha(L, weight_w) for weight_w in weights_w] # other root system

  asVec(v) = GAP.gap_to_julia(GAPWrap.ExtRepOfObj(v)) # TODO
  birational_sequence = BirationalSequence(
    operators, [asVec(v) for v in operators], weights_w, weights_alpha
  )

  ZZx, _ = polynomial_ring(ZZ, length(operators)) # for our monomials
  monomial_ordering = get_monomial_ordering(monomial_ordering_symb, ZZx, weights_alpha)

  # save computations from recursions
  calc_highest_weight = Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}}(
    [ZZ(0) for i in 1:rank(L)] => Set([ZZx(1)])
  )

  # save all highest weights, for which the Minkowski-sum did not suffice to gain all monomials
  no_minkowski = Set{Vector{ZZRingElem}}()
  set_mon_k = Set{ZZMPolyRingElem}[]        # monomial basis of the module k*highest_weight
  set_mon_k_missing = Set{ZZMPolyRingElem}[]  # store the monomials that are not products of basis monomials of smaller degree
  sizehint!(set_mon_k, degree)
  sizehint!(set_mon_k_missing, degree)
  mbs = MonomialBasis[]

  # start recursion over degree
  for i in 1:degree
    set_mon_sum = Set{ZZMPolyRingElem}()
    dim_i = dim_of_simple_module(L, i * highest_weight)
    # iterate over all minkowski sums of previous steps
    for k in 1:div(i, 2)
      set_help = Set([p * q for p in set_mon_k[i - k] for q in set_mon_k[k]])
      union!(set_mon_sum, set_help)
      if length(set_mon_sum) == dim_i
        break
      end
    end
    if length(set_mon_sum) == dim_i
      @vprintln :BasisLieHighestWeight "for $(Int.(i * highest_weight)) everything is generated by smaller weights"
      set_mon = set_mon_sum
      set_mon_missing = empty(set_mon_sum)
    else
      @vprintln :BasisLieHighestWeight "for $(Int.(i * highest_weight)) we have $(length(set_mon_sum)) and need $(dim_i) monomials"
      set_mon = compute_monomials(
        L,
        birational_sequence,
        ZZx,
        i * highest_weight,
        monomial_ordering,
        calc_highest_weight,
        no_minkowski,
      )
      set_mon_missing = setdiff(set_mon, set_mon_sum)
      @vprintln :BasisLieHighestWeight "for $(Int.(i * highest_weight)) we added $(length(set_mon_missing)) monomials "

      # monomials = sort(collect(set_mon); lt=((m1, m2) -> cmp(monomial_ordering, m1, m2) < 0))
      minkowski_gens = sort(collect(no_minkowski); by=(gen -> (sum(gen), reverse(gen))))
    end

    mb = MonomialBasis(
      L, i * highest_weight, birational_sequence, monomial_ordering, set_mon
    )
    set_attribute!(mb, :algorithm => basis_coordinate_ring_kodaira_compute)
    if isempty(set_mon_missing)
      set_attribute!(
        mb,
        :minkowski_gens =>
          [k * highest_weight for k in findall(!isempty, set_mon_k_missing)],
        :missing_monomials => nothing,
      )
    else
      set_attribute!(
        mb, :minkowski_gens => minkowski_gens, :missing_monomials => set_mon_missing
      )
    end
    push!(set_mon_k, set_mon)
    push!(set_mon_k_missing, set_mon_missing)
    push!(mbs, mb)
  end

  return collect(zip(mbs, length.(set_mon_k_missing)))
end

function compute_monomials(
  L::LieAlgebraStructure,
  birational_sequence::BirationalSequence,
  ZZx::ZZMPolyRing,
  highest_weight::Vector{ZZRingElem},
  monomial_ordering::MonomialOrdering,
  calc_highest_weight::Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}},
  no_minkowski::Set{Vector{ZZRingElem}},
)
  # This function calculates the monomial basis M_{highest_weight} recursively. The recursion saves all computed 
  # results in calc_highest_weight and we first check, if we already encountered this highest weight in a prior step. 
  # If this is not the case, we need to perform computations. The recursion works by using the Minkowski-sum. 
  # If M_{highest_weight} is the desired set of monomials (identified by the exponents as lattice points), it is know 
  # that for lambda_1 + lambda_2 = highest_weight we have M_{lambda_1} + M_{lambda_2} subseteq M_{highest_weight}. 
  # The complexity grows exponentially in the size of highest_weight. Therefore, it is very helpful to obtain a part of
  # M_{highest_weight} by going through all partitions of highest_weight and using the Minkowski-property. The base 
  # cases of the recursion are the fundamental weights highest_weight = [0, ..., 1, ..., 0]. In this case, or if the 
  # Minkowski-property did not find enough monomials, we need to perform the computations "by hand".

  # simple cases
  # we already computed the highest_weight result in a prior recursion step
  if haskey(calc_highest_weight, highest_weight)
    return calc_highest_weight[highest_weight]
  elseif highest_weight == [ZZ(0) for i in 1:(L.rank)] # we mathematically know the solution
    return Set(ZZx(1))
  end
  # calculation required
  # gap_dim is number of monomials that we need to find, i.e. |M_{highest_weight}|.
  # if highest_weight is not a fundamental weight, partition into smaller summands is possible. This is the basecase of 
  # the recursion.
  gap_dim = dim_of_simple_module(L, highest_weight)
  if is_fundamental(highest_weight) || sum(abs.(highest_weight)) == 0
    push!(no_minkowski, highest_weight)
    set_mon = add_by_hand(
      L, birational_sequence, ZZx, highest_weight, monomial_ordering, Set{ZZMPolyRingElem}()
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

      if lambda_1 > lambda_2
        continue
      end

      mon_lambda_1 = compute_monomials(
        L,
        birational_sequence,
        ZZx,
        lambda_1,
        monomial_ordering,
        calc_highest_weight,
        no_minkowski,
      )
      mon_lambda_2 = compute_monomials(
        L,
        birational_sequence,
        ZZx,
        lambda_2,
        monomial_ordering,
        calc_highest_weight,
        no_minkowski,
      )
      # Minkowski-sum: M_{lambda_1} + M_{lambda_2} \subseteq M_{highest_weight}, if monomials get identified with 
      # points in ZZ^n
      mon_sum = Set([p * q for p in mon_lambda_1 for q in mon_lambda_2])
      union!(set_mon, mon_sum)
    end
    # check if we found enough monomials

    if length(set_mon) < gap_dim
      push!(no_minkowski, highest_weight)
      set_mon = add_by_hand(
        L, birational_sequence, ZZx, highest_weight, monomial_ordering, set_mon
      )
    end

    push!(calc_highest_weight, highest_weight => set_mon)
    return set_mon
  end
end

function add_known_monomials!(
  weight_w::Vector{ZZRingElem},
  set_mon_in_weightspace::Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}},
  matrices_of_operators::Vector{<:SMat{ZZRingElem}},
  space::Dict{Vector{ZZRingElem},<:SMat{QQFieldElem}},
  v0::SRow{ZZRingElem},
)
  # By using the Minkowski-sum, we know that all monomials in set_mon_in_weightspace are in our basis. Since we want to
  # extend the weightspace with missing monomials, we need to calculate and add the vector of each monomial to our 
  # basis.

  for mon in set_mon_in_weightspace[weight_w]
    # calculate the vector vec associated with mon
    vec = calc_vec(v0, mon, matrices_of_operators)

    # check if vec extends the basis
    if !haskey(space, weight_w)
      space[weight_w] = sparse_matrix(QQ)
    end
    Hecke._add_row_to_rref!(space[weight_w], change_base_ring(QQ, vec))
  end
end

function add_new_monomials!(
  L::LieAlgebraStructure,
  birational_sequence::BirationalSequence,
  ZZx::ZZMPolyRing,
  matrices_of_operators::Vector{<:SMat{ZZRingElem}},
  monomial_ordering::MonomialOrdering,
  weightspaces::Dict{Vector{ZZRingElem},Int},
  dim_weightspace::Int,
  weight_w::Vector{ZZRingElem},
  set_mon_in_weightspace::Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}},
  space::Dict{Vector{ZZRingElem},<:SMat{QQFieldElem}},
  v0::SRow{ZZRingElem},
  set_mon::Set{ZZMPolyRingElem},
  zero_coordinates::Vector{Int},
)
  # If a weightspace is missing monomials, we need to calculate them by trial and error. We would like to go through all
  # monomials in the order monomial_ordering and calculate the corresponding vector. If it extends the basis, we add it 
  # to the result and else we try the next one. We know, that all monomials that work lay in the weyl-polytope. 
  # Therefore, we only inspect the monomials that lie both in the weyl-polytope and the weightspace. Since the weyl-
  # polytope is bounded these are finitely many and we can sort them and then go trough them, until we found enough. 

  # get monomials that are in the weightspace, sorted by monomial_ordering
  poss_mon_in_weightspace = convert_lattice_points_to_monomials(
    ZZx,
    get_lattice_points_of_weightspace(
      birational_sequence.weights_alpha, w_to_alpha(L, weight_w), zero_coordinates
    ),
  )
  isempty(poss_mon_in_weightspace) && error("The input seems to be invalid.")
  poss_mon_in_weightspace = sort(
    poss_mon_in_weightspace; lt=((m1, m2) -> cmp(monomial_ordering, m1, m2) < 0)
  )

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

    # check if the weight ob each suffix is a weight of the module
    cancel = false
    for i in 1:(nvars(ZZx) - 1)
      if !haskey(
        weightspaces,
        sum(
          exp * weight for (exp, weight) in
          Iterators.drop(zip(degrees(mon), birational_sequence.weights_w), i)
        ),
      )
        cancel = true
        break
      end
    end
    if cancel
      continue
    end

    # calculate the vector vec associated with mon
    vec = calc_vec(v0, mon, matrices_of_operators)

    # check if vec extends the basis
    if !haskey(space, weight_w)
      space[weight_w] = sparse_matrix(QQ)
    end
    fl = Hecke._add_row_to_rref!(space[weight_w], change_base_ring(QQ, vec))
    if !fl
      continue
    end

    # save monom
    number_mon_in_weightspace += 1
    push!(set_mon, mon)
  end
end

function add_by_hand(
  L::LieAlgebraStructure,
  birational_sequence::BirationalSequence,
  ZZx::ZZMPolyRing,
  highest_weight::Vector{ZZRingElem},
  monomial_ordering::MonomialOrdering,
  set_mon::Set{ZZMPolyRingElem},
)
  # This function calculates the missing monomials by going through each non full weightspace and adding possible 
  # monomials manually by computing their corresponding vectors and checking if they enlargen the basis.

  # initialization
  # matrices g_i for (g_1^a_1 * ... * g_k^a_k)*v
  matrices_of_operators = tensor_matrices_of_operators(
    L, highest_weight, birational_sequence.operators
  )
  space = Dict(ZZ(0) * birational_sequence.weights_w[1] => sparse_matrix(QQ)) # span of basis vectors to keep track of the basis
  v0 = sparse_row(ZZ, [(1, 1)])  # starting vector v

  push!(set_mon, ZZx(1))
  # required monomials of each weightspace
  weightspaces = get_dim_weightspace(L, highest_weight)
  # sort the monomials from the minkowski-sum by their weightspaces
  set_mon_in_weightspace = Dict{Vector{ZZRingElem},Set{ZZMPolyRingElem}}()
  for (weight_w, _) in weightspaces
    set_mon_in_weightspace[weight_w] = Set{ZZMPolyRingElem}()
  end
  for mon in set_mon
    weight_w = weight(mon, birational_sequence.weights_w)
    push!(set_mon_in_weightspace[weight_w], mon)
  end

  # only inspect weightspaces with missing monomials
  weights_with_non_full_weightspace = Set{Vector{ZZRingElem}}()
  for (weight_w, dim_weightspace) in weightspaces
    if length(set_mon_in_weightspace[weight_w]) != dim_weightspace
      push!(weights_with_non_full_weightspace, weight_w)
    end
  end

  # The weightspaces could be calculated completely indepent (except for
  # the caching). This is not implemented, since I used the package Distributed.jl for this, which is not in the 
  # Oscar dependendencies. But I plan to reimplement this. 
  # insert known monomials into basis

  for weight_w in weights_with_non_full_weightspace
    add_known_monomials!(weight_w, set_mon_in_weightspace, matrices_of_operators, space, v0)
  end

  # identify coordinates that are trivially zero because of the action on the generator
  zero_coordinates = compute_zero_coordinates(birational_sequence, highest_weight)

  # calculate new monomials
  for weight_w in weights_with_non_full_weightspace
    dim_weightspace = weightspaces[weight_w]
    add_new_monomials!(
      L,
      birational_sequence,
      ZZx,
      matrices_of_operators,
      monomial_ordering,
      weightspaces,
      dim_weightspace,
      weight_w,
      set_mon_in_weightspace,
      space,
      v0,
      set_mon,
      zero_coordinates,
    )
  end
  return set_mon
end

function operators_by_index(
  L::LieAlgebraStructure,
  chevalley_basis::NTuple{3,Vector{GAP.Obj}},
  birational_sequence::Vector{Int},
)
  @req all(i -> 1 <= i <= number_of_positive_roots(L), birational_sequence) "Entry of birational_sequence out of bounds"

  return [chevalley_basis[1][i] for i in birational_sequence] # TODO: change to [2]
end

function operators_by_simple_roots(
  L::LieAlgebraStructure,
  chevalley_basis::NTuple{3,Vector{GAP.Obj}},
  birational_sequence::Vector{Vector{Int}},
)
  rs = root_system_gap(L)
  simple_roots = Vector{Vector{Int}}(GAP.Globals.SimpleSystem(rs))
  positive_roots = Vector{Vector{Int}}(GAP.Globals.PositiveRoots(rs))

  root_inds = Int[]
  for whgt_alpha in birational_sequence
    @req length(whgt_alpha) == rank(L) "Length mismatch"
    @req all(>=(0), whgt_alpha) "Only positive roots are allowed as input"
    root = sum(whgt_alpha .* simple_roots)
    root_ind = findfirst(==(root), positive_roots)
    @req !isnothing(root_ind) "$whgt_alpha is not a positive root"
    push!(root_inds, root_ind)
  end

  return operators_by_index(L, chevalley_basis, root_inds)
end

function operators_lusztig(
  L::LieAlgebraStructure,
  chevalley_basis::NTuple{3,Vector{GAP.Obj}},
  reduced_expression::Vector{Int},
)
  root_inds = operators_lusztig_indices(L, reduced_expression)
  return operators_by_index(L, chevalley_basis, root_inds)
end

function operators_lusztig_indices(L::LieAlgebraStructure, word::Vector{Int})
  # Computes the operators for the lusztig polytopes for a longest weyl-word 
  # reduced_expression.

  # \beta_k := s_{i_1} … s_{i_{k-1}} (\alpha_{i_k})

  # F.e. for A, 2, [1, 2, 1], we get
  # \beta_1 = \alpha_1
  # \beta_2 = \alpha_1 + \alpha_2
  # \beta_3 = \alpha_2

  rs = root_system_gap(L)

  simple_roots = GAP.Globals.SimpleSystem(rs)
  positive_roots = Vector{Vector{Int}}(GAP.Globals.PositiveRoots(rs))
  sparse_cartan_matrix = GAP.Globals.SparseCartanMatrix(GAPWrap.WeylGroup(rs))

  root_inds = Int[]

  for k in 1:length(word)
    # Calculate betas by applying simple reflections step-by-step.
    root = copy(simple_roots[word[k]])
    for j in (k - 1):-1:1
      GAP.Globals.ApplySimpleReflection(sparse_cartan_matrix, word[j], root)
    end
    root_ind = findfirst(==(Vector{Int}(root)), positive_roots)
    @req !isnothing(root_ind) "$root is not a positive root"
    push!(root_inds, root_ind)
  end
  return root_inds
end

@doc """
    is_fundamental(highest_weight::Vector{IntegerUnion}) -> Bool

Return if ``highest_weight`` is fundamental, i.e. [0, ..., 1, ..., 0].

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
  for i in highest_weight
    if iszero(i)
      continue
    elseif isone(i)
      hasone && return false
      hasone = true
    else
      return false
    end
  end
  return hasone
end

function compute_sub_weights(highest_weight::Vector{ZZRingElem})
  # returns list of weights w != 0, highest_weight with 0 <= w <= highest_weight elementwise, ordered by l_2-norm

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
