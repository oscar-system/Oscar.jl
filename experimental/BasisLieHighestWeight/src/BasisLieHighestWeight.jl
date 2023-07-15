module BasisLieHighestWeight
export basis_lie_highest_weight
export is_fundamental

using ..Oscar
using ..Oscar: GAPWrap
using Polymake

# TODO Change operators from GAP.Obj to Oscar objects
# TODO BirationalSequence needs better names for weights and operators
# TODO Fix weightnames in general
# TODO Adapt arguments in function outside of BasisLieHighestWeight.jl
# TODO Summarize function
# TODO Groundup-structure for the special functions
# TODO Bugfix cache_size != 0

# TODO Export and docstring: 
# basis_lie_highest_weight
# get_dim_weightspace
# orbit_weylgroup
# get_lattice_points_of_weightspace
# convert_lattice_points_to_monomials
# convert_monomials_to_lattice_points

# tensorMatricesForOperators
# weights_for_operators

# w_to_eps
# eps_to_w
# alpha_to_eps
# eps_to_alpha
# w_to_eps
# eps_to_w

struct LieAlgebra
    lie_type::String
    rank::Int
    lie_algebra_gap::GAP.Obj
end

struct BirationalSequence
    operators::GAP.Obj # TODO Integer 
    operators_vectors::Vector{Vector{Any}}
    weights::Vector{Vector{Int}}
    weights_eps::Vector{Vector{Int}}
end

struct MonomialBasis
    set_mon::Set{ZZMPolyRingElem}
    dimension::Int
    no_minkowski::Set{Vector{Int}}
    # polytope::polytope
end

struct BasisLieHighestWeightStructure
    lie_algebra::LieAlgebra
    birational_sequence::BirationalSequence
    highest_weight::Vector{Int}
    monomial_order::Union{String, Function}
    monomial_basis::MonomialBasis
end

include("./VectorSpaceBases.jl")
include("./NewMonomial.jl")
include("./TensorModels.jl")
include("./LieAlgebras.jl")
include("./MonomialOrder.jl")
include("./RootConversion.jl")
include("./WeylPolytope.jl")

fromGap = Oscar.GAP.gap_to_julia

@doc """
basis_lie_highest_weight(
    type::String,
    rank::Int, 
    highest_weight::Vector{Int};
    operators::Union{String, Vector{Int}} = "regular", 
    monomial_order::Union{String, Function} = "GRevLex", 
    cache_size::Int = 0,
)::BasisLieHighestWeightStructure

Computes a monomial basis for the highest weight module with highest weight
``highest_weight`` (in terms of the fundamental weights), for a simple Lie algebra of type
``type`` and rank ``rank``.

# Parameters
- `type`: type of liealgebra we want to investigate, one of "A", "B", "C", "D", "E", "F", "G"
- `rank`: rank of liealgebra
- `highest_weight`: highest-weight
- `operators`: list of operators, either "regular" or integer array. The functionality of choosing a random longest word
                is currently not implemented, because we used https://github.com/jmichel7/Gapjm.jl to work with coxeter 
                groups need a method to obtain all non left descending elements to extend a word
- `monomial_order`: monomial order in which our basis gets defined with regards to our operators 
- `cache_size`: number of computed monomials we want to cache, default is 0 

# Examples
```jldoctest
julia> BasisLieHighestWeight.basis_lie_highest_weight("A", 2, [1, 1], return_no_minkowski = true, return_operators = true)
(Set(ZZMPolyRingElem[x1*x2, x2, 1, x1*x3, x3^2, x1, x3, x2*x3]), Set([[1, 0], [0, 1]]), GAP: [ v.1, v.2, v.3 ])


julia> BasisLieHighestWeight.basis_lie_highest_weight("A", 3, [2, 2, 3], monomial_order = "Lex")
Set{ZZMPolyRingElem} with 1260 elements:
  x3*x5^2*x6^2
  x2*x3*x5^2*x6
  x4^2*x5^2*x6
  x1^2*x3^3*x5
  x2*x3*x4*x5^3*x6^2
  x1*x3*x4*x5^3*x6^2
  x1^2*x3*x4*x6
  x1*x3*x4^3
  x4^2*x5^3*x6^4
  x1*x2*x3*x5^2
  x3^2*x4^4*x5^2*x6
  x2^2*x3*x6^2
  x1*x2^2*x3^2*x5
  x1*x3*x4*x5^2
  x1^2*x2*x6
  x1*x3^2*x4*x5*x6^3
  x1^2*x2*x4*x5^2*x6^2
  x4^4*x5
  x1^2*x2*x3^2*x6
  x1*x3^2*x5^2
  x2*x3*x4*x5^3
  ⋮

julia> BasisLieHighestWeight.basis_lie_highest_weight("A", 2, [1, 0], operators = [1,2,1])
Set{ZZMPolyRingElem} with 3 elements:
  1
  x3
  x2*x3

julia> BasisLieHighestWeight.basis_lie_highest_weight("C", 3, [1, 1, 1], monomial_order = "Lex")
Set{ZZMPolyRingElem} with 512 elements:
  x1*x5*x6*x8
  x6^4
  x3*x4^2*x8
  x3*x4*x6*x7
  x8^3
  x3*x6^2
  x2*x3
  x5*x6^2*x9
  x6*x8^2*x9
  x1*x6*x7
  x5*x6*x9^2
  x6^2*x7^2*x8
  x5*x7*x8
  x4*x6^2*x7*x8^2
  x4^2*x5*x7
  x1*x5^2*x6
  x1*x6*x8
  x3*x4*x5
  x2*x4*x6^2*x7
  x4*x6*x7
  x1*x4*x7*x8^2
  ⋮
```
"""
function basis_lie_highest_weight(
        type::String,
        rank::Int, 
        highest_weight::Vector{Int};
        operators::Union{String, Vector{Int}} = "regular", 
        monomial_order::Union{String, Function} = "GRevLex", 
        cache_size::Int = 0,
    )::BasisLieHighestWeightStructure

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
    # The function precomputes objects that are independent of the highest weight and that can be used in all recursion 
    # steps. Then it starts the recursion and returns the result.

    # initialization of objects that can be precomputed
    lie_algebra_gap, chevalley_basis = create_lie_algebra(type, rank) # lie_algebra of type, rank and its chevalley_basis
    lie_algebra = LieAlgebra(type, rank, lie_algebra_gap)
    # operators that are represented by our monomials. x_i is connected to operators[i]
    operators = get_operators(lie_algebra, operators, chevalley_basis) 
    weights = weights_for_operators(lie_algebra.lie_algebra_gap, chevalley_basis[3], operators) # weights of the operators
    weights = (weight->Int.(weight)).(weights)
    weights_eps = [w_to_eps(type, rank, w) for w in weights] # other root system
    
    asVec(v) = fromGap(GAP.Globals.ExtRepOfObj(v)) # TODO
    birational_sequence = BirationalSequence(operators, [asVec(v) for v in operators], weights, weights_eps)
    
    ZZx, _ = PolynomialRing(ZZ, length(operators)) # for our monomials
    monomial_order_lt = get_monomial_order_lt(monomial_order, ZZx) # less than function to sort monomials by order
    
    # save computations from recursions
    calc_highest_weight = Dict{Vector{Int}, Set{ZZMPolyRingElem}}([0 for i in 1:rank] => Set([ZZx(1)]))
    # we save all highest weights, for which the Minkowski-sum did not suffice to gain all monomials
    no_minkowski = Set{Vector{Int}}()

    # start recursion over highest_weight
    set_mon = compute_monomials(lie_algebra, birational_sequence, ZZx, highest_weight, monomial_order_lt, calc_highest_weight, cache_size, no_minkowski)
    
    # output
    return BasisLieHighestWeightStructure(
        lie_algebra,
        birational_sequence,
        highest_weight,
        monomial_order,
        MonomialBasis(
            set_mon,
            length(set_mon),
            no_minkowski
        )
    )
end

function sub_simple_refl(word::Vector{Int}, lie_algebra_gap::GAP.Obj)::GAP.Obj
    """
    substitute simple reflections (i,i+1), saved in dec by i, with E_{i,i+1}  
    """
    root_system = GAP.Globals.RootSystem(lie_algebra_gap)
    canonical_generators = fromGap(GAP.Globals.CanonicalGenerators(root_system)[1], recursive = false)
    operators = GAP.Obj([canonical_generators[i] for i in word], recursive = false)
    return operators
end

function get_operators(lie_algebra::LieAlgebra, operators::Union{String, Vector{Int}}, 
                        chevalley_basis::GAP.Obj)::GAP.Obj
    """
    handles user input for operators
    "regular" for all operators
    "longest-word" for random longest-word in Weyl-group (currently not implemented)
    operators::Vector{Int} for explicit longest-word
    """
    # create standard operators
    if operators == "regular" # use operators as specified by GAP
        operators = chevalley_basis[1]
        return operators
    # The functionality longest-word required Coxetergroups from Gapjm.jl (https://github.com/jmichel7/Gapjm.jl and was 
    # temporarily deleted
    # choose a random longest word. Created by extending by random not leftdescending reflections until total length is 
    # reached
    #elseif operators == "longest-word"
    #    operators = longest_weyl_word(t,n)
    #    operators = sub_simple_refl(operators, lie_algebra, n)
    #    return operators
    end

    # use user defined operators
    # wrong input
    if !(typeof(operators) == Vector{Int})
        println("operators needs to be of type Vector{Int}")
        return -1
    end
    if !(all([(1 <= i <= lie_algebra.rank) for i in operators]))
        println("all values of operators need to between 1 and the rank of the lie algebra.")
    end
    # If one of the conditions is met, the algorithms works. Otherwise a warning is printed (and can be ignored).
    #if  !(is_longest_weyl_word(type, rank, operators)) && !(Set(operators) == [i for i=1:n])
    #    println("WARNING: operators may be incorrect input.")
    #end
    operators = sub_simple_refl(operators, lie_algebra.lie_algebra_gap)
    return operators
end

function compute_monomials(
    lie_algebra::LieAlgebra,
    birational_sequence::BirationalSequence,
    ZZx::ZZMPolyRing, 
    highest_weight::Vector{Int},
    monomial_order_lt::Function, 
    calc_highest_weight::Dict{Vector{Int64}, Set{ZZMPolyRingElem}},
    cache_size::Int, 
    no_minkowski::Set{Vector{Int}})::Set{ZZMPolyRingElem}
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
    elseif highest_weight == [0 for i in 1:lie_algebra.rank] # we mathematically know the solution
        return Set(ZZx(1))
    end 

    # calculation required
    # gap_dim is number of monomials that we need to find, i.e. |M_{highest_weight}|.
    # if highest_weight is a fundamental weight, partition into smaller summands is possible. This is the basecase of 
    # the recursion.
    gap_dim = GAP.Globals.DimensionOfHighestWeightModule(lie_algebra.lie_algebra_gap, GAP.Obj(highest_weight)) # fundamental weights
    if is_fundamental(highest_weight) || sum(abs.(highest_weight)) == 0
        push!(no_minkowski, highest_weight)
        set_mon = add_by_hand(lie_algebra, birational_sequence, ZZx, highest_weight, 
                                monomial_order_lt, gap_dim, Set{ZZMPolyRingElem}(), cache_size)
        push!(calc_highest_weight, highest_weight => set_mon)
        return set_mon
    else
        # use Minkowski-Sum for recursion
        set_mon = Set{ZZMPolyRingElem}()
        i = 0
        sub_weights = compute_sub_weights(highest_weight)
        l = length(sub_weights)
        # go through all partitions lambda_1 + lambda_2 = highest_weight until we have enough monomials or used all 
        # partitions
        while length(set_mon) < gap_dim && i < l
            i += 1
            lambda_1 = sub_weights[i]
            lambda_2 = highest_weight .- lambda_1
            mon_lambda_1 = compute_monomials(lie_algebra, birational_sequence, ZZx, lambda_1,
                                                monomial_order_lt, calc_highest_weight, cache_size, 
                                                no_minkowski)
            mon_lambda_2 = compute_monomials(lie_algebra, birational_sequence, ZZx, lambda_2,
                                                monomial_order_lt, calc_highest_weight, cache_size, 
                                                no_minkowski)
            # Minkowski-sum: M_{lambda_1} + M_{lambda_2} \subseteq M_{highest_weight}, if monomials get identified with 
            # points in ZZ^n
            mon_sum = Set([p*q for p in mon_lambda_1 for q in mon_lambda_2])
            union!(set_mon, mon_sum)
        end
        # check if we found enough monomials
        if length(set_mon) < gap_dim
            push!(no_minkowski, highest_weight)
            set_mon = add_by_hand(lie_algebra, birational_sequence, ZZx, highest_weight, 
                                    monomial_order_lt, gap_dim, set_mon, cache_size)
        end
        push!(calc_highest_weight, highest_weight => set_mon)
        return set_mon
    end
end

@doc """
    is_fundamental(highest_weight::Vector{Int})::Bool

    returns true if ``highest_weight`` is fundamental, i.e. [0, ..., 1, ..., 0]

# Examples
```jldoctest
julia> BasisLieHighestWeight.is_fundamental([0, 1, 0])
true

julia> BasisLieHighestWeight.is_fundamental([0, 1, 1])
false
```
"""
function is_fundamental(highest_weight::Vector{Int})::Bool
    one = false
    for i in highest_weight
        if i > 0
            if one || i > 1
                return false
            else 
                one = true
            end
        end
    end
    return one
end

function compute_sub_weights(highest_weight::Vector{Int})::Vector{Vector{Int}}
    """
    returns list of weights w != 0, highest_weight with 0 <= w <= highest_weight elementwise, ordered by l_2-norm
    """
    sub_weights = []
    foreach(Iterators.product((0:x for x in highest_weight)...)) do i
        push!(sub_weights, [i...])
    end
    if isempty(sub_weights) || length(sub_weights) == 1 # case [] or [[0, ..., 0]]
        return []
    else
        popfirst!(sub_weights) # [0, ..., 0]
        pop!(sub_weights) # highest_weight
        sort!(sub_weights, by=x->sum((x).^2))
        return sub_weights
    end
end

function add_known_monomials!(
    birational_sequence::BirationalSequence,
    ZZx::ZZMPolyRing, 
    weight::Vector{Int},
    set_mon_in_weightspace::Dict{Vector{Int64}, 
    Set{ZZMPolyRingElem}},
    matrices_of_operators::Vector{SMat{ZZRingElem}},
    calc_monomials::Dict{ZZMPolyRingElem, Tuple{SRow{ZZRingElem}, Vector{Int}}},
    space::Dict{Vector{Int64}, Oscar.BasisLieHighestWeight.SparseVectorSpaceBasis},
    v0::SRow{ZZRingElem}, 
    cache_size::Int)
    """
    By using the Minkowski-sum, we know that all monomials in set_mon_in_weightspace are in our basis. Since we want to
    extend the weightspace with missing monomials, we need to calculate and add the vector of each monomial to our 
    basis.
    """
    for mon in set_mon_in_weightspace[weight]
        # calculate the vector vec associated with mon
        if cache_size == 0
            d = sz(matrices_of_operators[1])
            vec = calc_vec(v0, mon, matrices_of_operators) 
        else
            vec = calc_new_mon!(gens(ZZx) , mon, birational_sequence.weights, matrices_of_operators, calc_monomials, space,
                                cache_size)
        end

        # check if vec extends the basis
        if !haskey(space, weight)
            space[weight] = SparseVectorSpaceBasis([], [])
        end
        add_and_reduce!(space[weight], vec)
    end
end

function add_new_monomials!(
    lie_algebra::LieAlgebra,
    birational_sequence::BirationalSequence,
    ZZx::ZZMPolyRing, 
    matrices_of_operators::Vector{SMat{ZZRingElem}},
    monomial_order_lt::Function,
    dim_weightspace::Int,
    weight::Vector{Int},
    set_mon_in_weightspace::Dict{Vector{Int64}, Set{ZZMPolyRingElem}}, 
    calc_monomials::Dict{ZZMPolyRingElem, Tuple{SRow{ZZRingElem}, Vector{Int}}},
    space::Dict{Vector{Int64}, 
    Oscar.BasisLieHighestWeight.SparseVectorSpaceBasis},
    v0::SRow{ZZRingElem}, 
    cache_size::Int,
    set_mon::Set{ZZMPolyRingElem})
    """
    If a weightspace is missing monomials, we need to calculate them by trial and error. We would like to go through all
    monomials in the order monomial_order_lt and calculate the corresponding vector. If it extends the basis, we add it 
    to the result and else we try the next one. We know, that all monomials that work lay in the weyl-polytope. 
    Therefore, we only inspect the monomials that lie both in the weyl-polytope and the weightspace. Since the weyl-
    polytope is bounded these are finitely many and we can sort them and then go trough them, until we found enough. 
    """
    
    # get monomials that are in the weightspace, sorted by monomial_order_lt
    poss_mon_in_weightspace = convert_lattice_points_to_monomials(ZZx, get_lattice_points_of_weightspace(
        birational_sequence.weights_eps, w_to_eps(lie_algebra.lie_type, lie_algebra.rank, weight), lie_algebra.lie_type))
    poss_mon_in_weightspace = sort(poss_mon_in_weightspace, lt=monomial_order_lt)

    # check which monomials should get added to the basis
    i=0
    if weight == 0 # check if [0 0 ... 0] already in basis
        i += 1
    end
    number_mon_in_weightspace = length(set_mon_in_weightspace[weight])
    # go through possible monomials one by one and check if it extends the basis
    while number_mon_in_weightspace < dim_weightspace
        i += 1

        mon = poss_mon_in_weightspace[i]
        if mon in set_mon
            continue
        end

        # calculate the vector vec associated with mon
        if cache_size == 0
            d = sz(matrices_of_operators[1])
            vec = calc_vec(v0, mon, matrices_of_operators)
        else
            vec = calc_new_mon!(gens(ZZx), mon, birational_sequence.weights, matrices_of_operators, calc_monomials, space, 
                                cache_size)
        end

        # check if vec extends the basis
        if !haskey(space, weight)
            space[weight] = SparseVectorSpaceBasis([], [])
        end
        vec_red = add_and_reduce!(space[weight], vec)
        if isempty(vec_red) # v0 == 0
            continue
        end

        # save monom
        number_mon_in_weightspace += 1
        push!(set_mon, mon)
    end
end


function add_by_hand(
    lie_algebra::LieAlgebra,
    birational_sequence::BirationalSequence,
    ZZx::ZZMPolyRing,
    highest_weight::Vector{Int},
    monomial_order_lt::Function,
    gap_dim::Int, 
    set_mon::Set{ZZMPolyRingElem},
    cache_size::Int,
    )::Set{ZZMPolyRingElem}
    """
    This function calculates the missing monomials by going through each non full weightspace and adding possible 
    monomials manually by computing their corresponding vectors and checking if they enlargen the basis.
    """
    # initialization
    # matrices g_i for (g_1^a_1 * ... * g_k^a_k)*v
    matrices_of_operators = tensorMatricesForOperators(lie_algebra.lie_algebra_gap, highest_weight, birational_sequence.operators)
    space = Dict(0*birational_sequence.weights[1] => SparseVectorSpaceBasis([], [])) # span of basis vectors to keep track of the basis
    v0 = sparse_row(ZZ, [(1,1)])  # starting vector v
    # saves the calculated vectors to decrease necessary matrix multiplicatons
    calc_monomials = Dict{ZZMPolyRingElem, Tuple{SRow{ZZRingElem}, Vector{Int}}}(ZZx(1) => (v0, 0 * birational_sequence.weights[1])) 
    push!(set_mon, ZZx(1))
    # required monomials of each weightspace
    weightspaces = get_dim_weightspace(lie_algebra, highest_weight)

    # sort the monomials from the minkowski-sum by their weightspaces
    set_mon_in_weightspace = Dict{Vector{Int}, Set{ZZMPolyRingElem}}()
    for (weight, _) in weightspaces
        set_mon_in_weightspace[weight] = Set{ZZMPolyRingElem}()
    end
    for mon in set_mon
        weight = calc_weight(mon, birational_sequence.weights)
        push!(set_mon_in_weightspace[weight], mon)
    end

    # only inspect weightspaces with missing monomials
    weights_with_full_weightspace = Set{Vector{Int}}()
    for (weight, dim_weightspace) in weightspaces
        if (length(set_mon_in_weightspace[weight]) == dim_weightspace)
            push!(weights_with_full_weightspace, weight)
        end
    end
    delete!(weightspaces, weights_with_full_weightspace)

    # The weightspaces could be calculated completely indepent (except for
    # the caching). This is not implemented, since I used the package Distributed.jl for this, which is not in the 
    # Oscar dependencies. But I plan to reimplement this. 
    # insert known monomials into basis
    for (weight, _) in weightspaces
        add_known_monomials!(birational_sequence, ZZx, weight, set_mon_in_weightspace, matrices_of_operators, 
                                calc_monomials, space, v0, cache_size)
    end 

    # calculate new monomials
    for (weight, dim_weightspace) in weightspaces
        add_new_monomials!(lie_algebra, birational_sequence, ZZx, matrices_of_operators,
                            monomial_order_lt, 
                            dim_weightspace, weight,
                            set_mon_in_weightspace, 
                            calc_monomials, space, v0, 
                            cache_size, set_mon)         
    end
    return set_mon
end

function get_dim_weightspace(
    lie_algebra::LieAlgebra, 
    highest_weight::Vector{Int}
    )::Dict{Vector{Int}, Int}
    """
    Calculates dictionary with weights as keys and dimension of corresponding weightspace as value. GAP computes the 
    dimension for all positive weights. The dimension is constant on orbits of the weylgroup, and we can therefore 
    calculate the dimension of each weightspace.
    """
    # calculate dimension for dominant weights with GAP
    root_system = GAP.Globals.RootSystem(lie_algebra.lie_algebra_gap)
    result = GAP.Globals.DominantCharacter(root_system, GAP.Obj(highest_weight))
    dominant_weights = [map(Int, item) for item in result[1]]
    dominant_weights_dim = map(Int, result[2])                                                                      
    dominant_weights = convert(Vector{Vector{Int}}, dominant_weights)
    weightspaces = Dict{Vector{Int}, Int}() 

    # calculate dimension for the rest by checking which positive weights lies in the orbit.
    for i in 1:length(dominant_weights)
        orbit_weights = orbit_weylgroup(lie_algebra, dominant_weights[i])
        dim_weightspace = dominant_weights_dim[i]
        for weight in orbit_weights
            weightspaces[highest_weight - weight] = dim_weightspace
        end
    end
    return weightspaces
end

end
export BasisLieHighestWeight
