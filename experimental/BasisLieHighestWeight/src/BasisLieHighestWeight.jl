module BasisLieHighestWeight
export LieAlgebraStructure
export BirationalSequence
export MonomialBasis
export BasisLieHighestWeightStructure
export basis_lie_highest_weight
export is_fundamental
export get_dim_weightspace
export orbit_weylgroup

using ..Oscar
using ..Oscar: GAPWrap
using Polymake

# TODO basis_lie_highest_weight_lustzig
# TODO basis_lie_highest_weight_string
# TODO basis_lie_highest_weight_feigin_fflv
# TODO basis_lie_highest_weight_nZ

# TODO write methods small when using Polymake, f.e. polyhedron instead of Polyhedron (The lowercase names are not found)

# TODO (?) Maybe export and docstring: 
# get_dim_weightspace
# orbit_weylgroup
# get_lattice_points_of_weightspace
# convert_lattice_points_to_monomials
# convert_monomials_to_lattice_points
# tensorMatricesForOperators
# weights_for_operators

# TODO Use Oscar-Lie-Algebra type instead of LieAlgebra
# TODO Data-Type for weights of Lie-Algebras? Three types, in alpha_i, w_i and eps_i, conversion is defined in RootConversion
# w_to_eps
# eps_to_w
# alpha_to_eps
# eps_to_alpha
# w_to_eps
# eps_to_w

# TODO GAPWrap-wrappers are missing for 
# ChevalleyBasis
# DimensionOfHighestWeightModule
# SimpleLieAlgebra
# Rationals
# HighestWeightModule
# List
# MatrixOfAction
# RootSystem
# CartanMatrix
# WeylGroup
# DominantCharacter
# DimensionOfHighestWeightModule
# CanonicalGenerators

struct LieAlgebraStructure
    lie_type::String
    rank::Int
    lie_algebra_gap::GAP.Obj
end

function LieAlgebraStructure(lie_type::String, rank::Int)
    return LieAlgebraStructure(lie_type, rank, create_lie_algebra(lie_type, rank))
end

function Base.show(io::IO, lie_algebra::LieAlgebraStructure)
    print(io, "Lie-Algebra of type ", lie_algebra.lie_type, " and rank ", lie_algebra.rank)
end

struct BirationalSequence
    operators::Vector{GAP.Obj} # TODO Integer 
    operators_vectors::Vector{Vector{Any}}
    weights_w::Vector{Vector{ZZRingElem}}
    weights_eps::Vector{Vector{QQFieldElem}}
end

function Base.show(io::IO, birational_sequence::BirationalSequence)
    println(io, "BirationalSequence")
    println(io, "Operators: ", birational_sequence.operators)
    print(io, "Weights in w_i:", birational_sequence.weights_w)
end

struct MonomialBasis
    ZZx::ZZMPolyRing
    set_mon::Set{ZZMPolyRingElem}
    dimension::Int
    no_minkowski::Set{Vector{Int}}
    polytope::Oscar.Polymake.BigObjectAllocated
end

function MonomialBasis(ZZx::ZZMPolyRing, set_mon::Set{ZZMPolyRingElem}, no_minkowski::Set{Vector{ZZRingElem}})
    vertices = degrees.(collect(set_mon))
    vertices_hom = transpose(reduce(hcat, [prepend!(vec, 1) for vec in vertices])) # homogenoues coordinate system
    poly = Oscar.Polymake.polytope.Polytope(POINTS=vertices_hom)
    return MonomialBasis(ZZx, set_mon, length(set_mon), no_minkowski, poly)
end

function Base.show(io::IO, monomial_basis::MonomialBasis)
    println(io, "MonomialBasis")
    println(io, "Dimension: ", monomial_basis.dimension)
    println(io, "Generators within semi-group: ", monomial_basis.no_minkowski)
    println(io, "First 10 Monomials in degrevlex: ", sort(collect(monomial_basis.set_mon), 
            lt = get_monomial_order_lt("degrevlex", monomial_basis.ZZx))[1:min(end, 10)])
    # print(io, "Volume polytope: ", monomial_basis.polytope.VOLUME)
end

struct BasisLieHighestWeightStructure
    lie_algebra::LieAlgebraStructure
    birational_sequence::BirationalSequence
    highest_weight::Vector{Int}
    monomial_order::Union{String, Function}
    monomial_basis::MonomialBasis
end

function Base.show(io::IO, base::BasisLieHighestWeightStructure)
    println(io, base.lie_algebra)
    println("")
    println(io, base.birational_sequence)
    println("")
    println(io, "Highest-weight: ", base.highest_weight)
    println(io, "Monomial-order: ", base.monomial_order)
    println("")
    print(io, base.monomial_basis)
end


include("./VectorSpaceBases.jl")
include("./NewMonomial.jl")
include("./TensorModels.jl")
include("./LieAlgebras.jl")
include("./MonomialOrder.jl")
include("./RootConversion.jl")
include("./WeylPolytope.jl")
include("./DemazureOperators.jl")
include("./WordCalculations.jl")

fromGap = Oscar.GAP.gap_to_julia

function basis_lie_highest_weight_compute(
    type::String,
    rank::Int, 
    highest_weight::Vector{Int},
    get_operators::Function,
    monomial_order::Union{String, Function},
    cache_size::Int,
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
    highest_weight = convert(Vector{ZZRingElem}, highest_weight)
    # The function precomputes objects that are independent of the highest weight and that can be used in all recursion 
    # steps. Then it starts the recursion and returns the result.

    # initialization of objects that can be precomputed
    # lie_algebra of type, rank and its chevalley_basis
    lie_algebra = LieAlgebraStructure(type, rank)
    chevalley_basis = GAP.Globals.ChevalleyBasis(lie_algebra.lie_algebra_gap)

    # operators that are represented by our monomials. x_i is connected to operators[i]
    operators = get_operators(lie_algebra, chevalley_basis)

    weights_w = weights_for_operators(lie_algebra.lie_algebra_gap, chevalley_basis[3], operators) # weights of the operators
    # weights_w = (weight_w->Int.(weight_w)).(weights_w)  
    weights_eps = [w_to_eps(type, rank, convert(Vector{QQFieldElem}, weight_w)) for weight_w in weights_w] # other root system
    
    asVec(v) = fromGap(GAPWrap.ExtRepOfObj(v)) # TODO
    birational_sequence = BirationalSequence(operators, [asVec(v) for v in operators], weights_w, weights_eps)
    
    ZZx, _ = PolynomialRing(ZZ, length(operators)) # for our monomials
    monomial_order_lt = get_monomial_order_lt(monomial_order, ZZx) # less than function to sort monomials by order
    
    # save computations from recursions
    calc_highest_weight = Dict{Vector{ZZRingElem}, Set{ZZMPolyRingElem}}([ZZ(0) for i in 1:rank] => Set([ZZx(1)]))
    # save all highest weights, for which the Minkowski-sum did not suffice to gain all monomials
    no_minkowski = Set{Vector{ZZRingElem}}()

    # start recursion over highest_weight
    set_mon = compute_monomials(lie_algebra, birational_sequence, ZZx, highest_weight, monomial_order_lt, calc_highest_weight, cache_size, no_minkowski)
    
    # output
    return BasisLieHighestWeightStructure(
        lie_algebra,
        birational_sequence,
        highest_weight,
        monomial_order,
        MonomialBasis(
            ZZx,
            set_mon,
            no_minkowski
        )
    )
end

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
julia> base = BasisLieHighestWeight.basis_lie_highest_weight("A", 2, [1, 1])
Lie-Algebra of type A and rank 2

BirationalSequence
Operators: GAP: [ v.1, v.2, v.3 ]
Weights in w_i:[[2, -1], [-1, 2], [1, 1]]

Highest-weight: [1, 1]
Monomial-order: GRevLex

MonomialBasis
Dimension: 8
Generators within semi-group: Set([[1, 0], [0, 1]])
First 10 Monomials in GRevLex: ZZMPolyRingElem[1, x3, x2, x1, x3^2, x2*x3, x1*x3, x1*x2]

julia> base = BasisLieHighestWeight.basis_lie_highest_weight("A", 3, [2, 2, 3], monomial_order = "Lex")
Lie-Algebra of type A and rank 3

BirationalSequence
Operators: GAP: [ v.1, v.2, v.3, v.4, v.5, v.6 ]
Weights in w_i:[[2, -1, 0], [-1, 2, -1], [0, -1, 2], [1, 1, -1], [-1, 1, 1], [1, 0, 1]]

Highest-weight: [2, 2, 3]
Monomial-order: Lex

MonomialBasis
Dimension: 1260
Generators within semi-group: Set([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
First 10 Monomials in GRevLex: ZZMPolyRingElem[1, x6, x5, x4, x3, x2, x1, x6^2, x5*x6, x4*x6]

julia> base = BasisLieHighestWeight.basis_lie_highest_weight("A", 2, [1, 0], operators = [1,2,1])
Lie-Algebra of type A and rank 2

BirationalSequence
Operators: GAP: [ v.1, v.2, v.1 ]
Weights in w_i:[[2, -1], [-1, 2], [2, -1]]

Highest-weight: [1, 0]
Monomial-order: GRevLex

MonomialBasis
Dimension: 3
Generators within semi-group: Set([[1, 0]])
First 10 Monomials in GRevLex: ZZMPolyRingElem[1, x3, x2*x3]

julia> base = BasisLieHighestWeight.basis_lie_highest_weight("C", 3, [1, 1, 1], monomial_order = "Lex")
Lie-Algebra of type C and rank 3

BirationalSequence
Operators: GAP: [ v.1, v.2, v.3, v.4, v.5, v.6, v.7, v.8, v.9 ]
Weights in w_i:[[2, -1, 0], [-1, 2, -1], [0, -2, 2], [1, 1, -1], [-1, 0, 1], [1, -1, 1], [-2, 2, 0], [0, 1, 0], [2, 0, 0]]

Highest-weight: [1, 1, 1]
Monomial-order: Lex

MonomialBasis
Dimension: 512
Generators within semi-group: Set([[0, 1, 0], [1, 0, 0], [0, 0, 1], [1, 1, 1], [0, 1, 1]])
First 10 Monomials in GRevLex: ZZMPolyRingElem[1, x9, x8, x7, x6, x5, x4, x3, x2, x1]
```
"""
function basis_lie_highest_weight(
        type::String,
        rank::Int, 
        highest_weight::Vector{Int};
        reduced_expression::Union{String, Vector{Int}, Vector{GAP.GapObj}, Any} = "regular", 
        monomial_order::Union{String, Function} = "degrevlex", 
        cache_size::Int = 0,
    )::BasisLieHighestWeightStructure
    """
    Standard function with all options
    """
    get_operators = (lie_algebra, chevalley_basis) -> get_operators_normal(lie_algebra, chevalley_basis, reduced_expression)
    return basis_lie_highest_weight_compute(type, rank, highest_weight, get_operators, monomial_order, cache_size)
end

@doc """
# Examples
```jldoctest
julia> base = BasisLieHighestWeight.basis_lie_highest_weight_lustzig("D", 4, [1,1,1,1], [4,3,2,4,3,2,1,2,4,3,2,1])
Lie-Algebra of type D and rank 4

BirationalSequence
Operators: Union{Bool, Int64, GAP.FFE, GAP_jll.GapObj}[GAP: v.4, GAP: v.3, GAP: v.10, GAP: v.6, GAP: v.7, GAP: v.2, GAP: v.12, GAP: v.11, GAP: v.9, GAP: v.8, GAP: v.5, GAP: v.1]
Weights in w_i:Vector{ZZRingElem}[[0, -1, 0, 2], [0, -1, 2, 0], [-1, 0, 1, 1], [-1, 1, 1, -1], [-1, 1, -1, 1], [-1, 2, -1, -1], [0, 1, 0, 0], [1, -1, 1, 1], [1, 0, -1, 1], [1, 0, 1, -1], [1, 1, -1, -1], [2, -1, 0, 0]]

Highest-weight: [1, 1, 1, 1]
Monomial-order: oplex

MonomialBasis
Dimension: 4096
Generators within semi-group: Set([[0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 1], [0, 0, 1, 0], [1, 0, 0, 0]])
First 10 Monomials in degrevlex: ZZMPolyRingElem[1, x12, x11, x10, x9, x8, x7, x6, x5, x4]
```
"""
function basis_lie_highest_weight_lustzig(
    type::String,
    rank::Int, 
    highest_weight::Vector{Int},
    reduced_expression::Vector{Int};    
    monomial_order::Union{String, Function} = "oplex", 
    cache_size::Int = 0,
    )::BasisLieHighestWeightStructure
    """
    Lustzig polytope
    BasisLieHighestWeight.basis_lie_highest_weight_lustzig("D", 4, [1,1,1,1], [4,3,2,4,3,2,1,2,4,3,2,1])
    """
    # operators = some sequence of the String / Littelmann-Berenstein-Zelevinsky polytope
    get_operators = (lie_algebra, chevalley_basis) -> get_operators_lustzig_nz(lie_algebra, chevalley_basis, reduced_expression)
    return basis_lie_highest_weight_compute(type, rank, highest_weight, get_operators, monomial_order, cache_size)
end

@doc """
# Examples
```jldoctest
julia> BasisLieHighestWeight.basis_lie_highest_weight_string("B", 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
Lie-Algebra of type B and rank 3

BirationalSequence
Operators: Union{Bool, Int64, GAP.FFE, GAP_jll.GapObj}[GAP: v.3, GAP: v.2, GAP: v.3, GAP: v.2, GAP: v.1, GAP: v.2, GAP: v.3, GAP: v.2, GAP: v.1]
Weights in w_i:Vector{ZZRingElem}[[0, -1, 2], [-1, 2, -2], [0, -1, 2], [-1, 2, -2], [2, -1, 0], [-1, 2, -2], [0, -1, 2], [-1, 2, -2], [2, -1, 0]]

Highest-weight: [1, 1, 1]
Monomial-order: oplex

MonomialBasis
Dimension: 512
Generators within semi-group: Set([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
First 10 Monomials in degrevlex: ZZMPolyRingElem[1, x5, x2, x1, x5*x6, x2*x5, x1*x5, x2*x3, x1*x2, x5*x6*x7]
```
"""
function basis_lie_highest_weight_string(
    type::String,
    rank::Int, 
    highest_weight::Vector{Int},
    reduced_expression::Vector{Int};
    cache_size::Int = 0,
    )::BasisLieHighestWeightStructure
    """
    String / Littelmann-Berenstein-Zelevinsky polytope
    BasisLieHighestWeight.basis_lie_highest_weight_string("B", 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
    BasisLieHighestWeight.basis_lie_highest_weight_string("B", 4, [1,1,1,1], [4,3,4,3,2,3,4,3,2,1,2,3,4,3,2,1])
    """
    # reduced_expression = some sequence of the String / Littelmann-Berenstein-Zelevinsky polytope
    monomial_order = "oplex"
    get_operators = (lie_algebra, chevalley_basis) -> get_operators_normal(lie_algebra, chevalley_basis, reduced_expression)
    return basis_lie_highest_weight_compute(type, rank, highest_weight, get_operators, monomial_order, cache_size)
end

@doc """
# Examples
```jldoctest
julia> BasisLieHighestWeight.basis_lie_highest_weight_fflv("A", 3, [1,1,1])
Lie-Algebra of type A and rank 3

BirationalSequence
Operators: Union{Bool, Int64, GAP.FFE, GAP_jll.GapObj}[GAP: v.1, GAP: v.2, GAP: v.3, GAP: v.4, GAP: v.5, GAP: v.6]
Weights in w_i:Vector{ZZRingElem}[[2, -1, 0], [-1, 2, -1], [0, -1, 2], [1, 1, -1], [-1, 1, 1], [1, 0, 1]]

Highest-weight: [1, 1, 1]
Monomial-order: oplex

MonomialBasis
Dimension: 64
Generators within semi-group: Set([[0, 1, 0], [1, 0, 1], [0, 0, 1], [1, 0, 0]])
First 10 Monomials in degrevlex: ZZMPolyRingElem[1, x6, x5, x4, x3, x2, x1, x5*x6, x2*x6, x4*x5]
```
"""
function basis_lie_highest_weight_fflv(
    type::String,
    rank::Int,
    highest_weight::Vector{Int};
    cache_size::Int = 0,
    )::BasisLieHighestWeightStructure
    """
    Feigin-Fourier-Littelmann-Vinberg polytope
    BasisLieHighestWeight.basis_lie_highest_weight_fflv("A", 3, [1,1,1])
    """
    monomial_order = "oplex"
    # operators = all positive roots, same ordering as GAP uses
    get_operators = (lie_algebra, chevalley_basis) -> get_operators_normal(lie_algebra, chevalley_basis, "regular")
    return basis_lie_highest_weight_compute(type, rank, highest_weight, get_operators, monomial_order, cache_size)
end

@doc """
# Examples
```jldoctest
julia> BasisLieHighestWeight.basis_lie_highest_weight_nz("C", 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
Lie-Algebra of type C and rank 3

BirationalSequence
Operators: Union{Bool, Int64, GAP.FFE, GAP_jll.GapObj}[GAP: v.3, GAP: v.5, GAP: v.7, GAP: v.2, GAP: v.8, GAP: v.6, GAP: v.9, GAP: v.4, GAP: v.1]
Weights in w_i:Vector{ZZRingElem}[[0, -2, 2], [-1, 0, 1], [-2, 2, 0], [-1, 2, -1], [0, 1, 0], [1, -1, 1], [2, 0, 0], [1, 1, -1], [2, -1, 0]]

Highest-weight: [1, 1, 1]
Monomial-order: lex

MonomialBasis
Dimension: 512
Generators within semi-group: Set([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
First 10 Monomials in degrevlex: ZZMPolyRingElem[1, x9, x8, x7, x6, x5, x4, x3, x2, x1]
```
"""
function basis_lie_highest_weight_nz(
    type::String,
    rank::Int, 
    highest_weight::Vector{Int},
    reduced_expression::Vector{Int};
    cache_size::Int = 0,
    )::BasisLieHighestWeightStructure
    """
    Nakashima-Zelevinsky polytope
    BasisLieHighestWeight.basis_lie_highest_weight_nz("C", 3, [1,1,1], [3,2,3,2,1,2,3,2,1])
    """
    monomial_order = "lex"
    get_operators = (lie_algebra, chevalley_basis) -> get_operators_lustzig_nz(lie_algebra, chevalley_basis, reduced_expression)
    return basis_lie_highest_weight_compute(type, rank, highest_weight, get_operators, monomial_order, cache_size)
end

function sub_simple_refl(word::Vector{Int}, lie_algebra_gap::GAP.Obj)::Vector{GAP.Obj}
    """
    substitute simple reflections (i,i+1), saved in dec by i, with E_{i,i+1}  
    """
    root_system = GAP.Globals.RootSystem(lie_algebra_gap)
    canonical_generators = fromGap(GAP.Globals.CanonicalGenerators(root_system)[1], recursive = false)
    operators = [canonical_generators[i] for i in word]
    return operators
end

function get_operators_normal(
    lie_algebra::LieAlgebraStructure, 
    chevalley_basis::GAP.Obj,
    operators::Union{String, Vector{Int}, Vector{GAP.GapObj}, Any}, 
    )::Vector{GAP.Obj}
    """
    handles user input for operators
    "regular" for all operators
    "longest-word" for random longest-word in Weyl-group (currently not implemented)
    operators::Vector{Int} for explicit longest-word
    """
    #if typeof(operators) == GAP.Obj  # If user already submitted gap-roots as operators, keep
    if typeof(operators) != String && typeof(operators) != Vector{Int}
        return operators
    elseif operators == "regular" # create standard operators, use operators as specified by GAP
        return [v for v in chevalley_basis[1]]
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
    lie_algebra::LieAlgebraStructure,
    birational_sequence::BirationalSequence,
    ZZx::ZZMPolyRing, 
    highest_weight::Vector{ZZRingElem},
    monomial_order_lt::Function, 
    calc_highest_weight::Dict{Vector{ZZRingElem}, Set{ZZMPolyRingElem}},
    cache_size::Int, 
    no_minkowski::Set{Vector{ZZRingElem}})::Set{ZZMPolyRingElem}
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
    elseif highest_weight == [ZZ(0) for i in 1:lie_algebra.rank] # we mathematically know the solution
        return Set(ZZx(1))
    end 

    # calculation required
    # gap_dim is number of monomials that we need to find, i.e. |M_{highest_weight}|.
    # if highest_weight is a fundamental weight, partition into smaller summands is possible. This is the basecase of 
    # the recursion.
    highest_weight_int = convert(Vector{Int}, highest_weight)
    gap_dim = GAP.Globals.DimensionOfHighestWeightModule(lie_algebra.lie_algebra_gap, GAP.Obj(highest_weight_int)) # fundamental weights
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
        sub_weights_w = compute_sub_weights(highest_weight)
        l = length(sub_weights_w)
        # go through all partitions lambda_1 + lambda_2 = highest_weight until we have enough monomials or used all 
        # partitions
        while length(set_mon) < gap_dim && i < l
            i += 1
            lambda_1 = sub_weights_w[i]
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

        #println("")
        #println("highest_weight: ", highest_weight)
        #println("required monomials: ", gap_dim)
        #println("monomials from Minkowski-sum: ", length(set_mon))
        #println(set_mon)

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
function is_fundamental(highest_weight::Vector{ZZRingElem})::Bool
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
        sort!(sub_weights_w, by=x->sum((x).^2))
        return sub_weights_w
    end
end

function add_known_monomials!(
    birational_sequence::BirationalSequence,
    ZZx::ZZMPolyRing, 
    weight_w::Vector{ZZRingElem},
    set_mon_in_weightspace::Dict{Vector{ZZRingElem}, 
    Set{ZZMPolyRingElem}},
    matrices_of_operators::Vector{SMat{ZZRingElem}},
    calc_monomials::Dict{ZZMPolyRingElem, Tuple{SRow{ZZRingElem}, Vector{Int}}},
    space::Dict{Vector{ZZRingElem}, Oscar.BasisLieHighestWeight.SparseVectorSpaceBasis},
    v0::SRow{ZZRingElem}, 
    cache_size::Int)
    """
    By using the Minkowski-sum, we know that all monomials in set_mon_in_weightspace are in our basis. Since we want to
    extend the weightspace with missing monomials, we need to calculate and add the vector of each monomial to our 
    basis.
    """
    for mon in set_mon_in_weightspace[weight_w]
        # calculate the vector vec associated with mon
        if cache_size == 0
            d = sz(matrices_of_operators[1])
            vec = calc_vec(v0, mon, matrices_of_operators) 
        else
            vec = calc_new_mon!(gens(ZZx) , mon, birational_sequence.weights_w, matrices_of_operators, calc_monomials, space,
                                cache_size)
        end

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
    set_mon_in_weightspace::Dict{Vector{ZZRingElem}, Set{ZZMPolyRingElem}}, 
    calc_monomials::Dict{ZZMPolyRingElem, Tuple{SRow{ZZRingElem}, Vector{Int}}},
    space::Dict{Vector{ZZRingElem}, Oscar.BasisLieHighestWeight.SparseVectorSpaceBasis},
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
        birational_sequence.weights_eps, w_to_eps(lie_algebra.lie_type, lie_algebra.rank, convert(Vector{QQFieldElem}, weight_w)), lie_algebra.lie_type))
    poss_mon_in_weightspace = sort(poss_mon_in_weightspace, lt=monomial_order_lt)

    # check which monomials should get added to the basis
    i=0
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
        if cache_size == 0
            d = sz(matrices_of_operators[1])
            vec = calc_vec(v0, mon, matrices_of_operators)
        else
            vec = calc_new_mon!(gens(ZZx), mon, birational_sequence.weights_w, matrices_of_operators, calc_monomials, space, 
                                cache_size)
        end

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
    gap_dim::Int, 
    set_mon::Set{ZZMPolyRingElem},
    cache_size::Int,
    )::Set{ZZMPolyRingElem}

    # println("add_by_hand", highest_weight)
    set_mon_temp = copy(set_mon)

    """
    This function calculates the missing monomials by going through each non full weightspace and adding possible 
    monomials manually by computing their corresponding vectors and checking if they enlargen the basis.
    """
    # initialization
    # matrices g_i for (g_1^a_1 * ... * g_k^a_k)*v
    matrices_of_operators = tensorMatricesForOperators(lie_algebra.lie_algebra_gap, highest_weight, birational_sequence.operators)    
    space = Dict(ZZ(0)*birational_sequence.weights_w[1] => SparseVectorSpaceBasis([], [])) # span of basis vectors to keep track of the basis
    v0 = sparse_row(ZZ, [(1,1)])  # starting vector v
    # saves the calculated vectors to decrease necessary matrix multiplicatons
    calc_monomials = Dict{ZZMPolyRingElem, Tuple{SRow{ZZRingElem}, Vector{Int}}}(ZZx(1) => (v0, ZZ(0) * birational_sequence.weights_w[1])) 
    push!(set_mon, ZZx(1))
    # required monomials of each weightspace
    weightspaces = get_dim_weightspace(lie_algebra, highest_weight)

    # sort the monomials from the minkowski-sum by their weightspaces
    set_mon_in_weightspace = Dict{Vector{ZZRingElem}, Set{ZZMPolyRingElem}}()
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
        add_known_monomials!(birational_sequence, ZZx, weight_w, set_mon_in_weightspace, matrices_of_operators, 
                                calc_monomials, space, v0, cache_size)
    end 

    # println("")
    # calculate new monomials
    for (weight_w, dim_weightspace) in weightspaces
        # print("*")
        add_new_monomials!(lie_algebra, birational_sequence, ZZx, matrices_of_operators,
                            monomial_order_lt, 
                            dim_weightspace, weight_w,
                            set_mon_in_weightspace, 
                            calc_monomials, space, v0, 
                            cache_size, set_mon)       
    end

    #println("")
    #println("highest_weight: ", highest_weight)
    #println("added-by-hand: ", [mon for mon in set_mon if !(mon in set_mon_temp)])

    return set_mon
end

end
export BasisLieHighestWeight
