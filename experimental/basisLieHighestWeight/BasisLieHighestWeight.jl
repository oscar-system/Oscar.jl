module BasisLieHighestWeight
export basisLieHighestWeight2
export is_fundamental

using Polymake
using Distributed
using Markdown
#using CUDA
#using CuArrays

include("./NewMonomial.jl")

#include("./VectorSpaceBases.jl") #--- bekommt gerade noch ZZ, Short und TVEC aus VectorSpaceBases

fromGap = Oscar.GAP.gap_to_julia

@doc Markdown.doc"""
    basisLieHighestWeight2(t, n, hw; ops = "regular", known_monomials = [], monomial_order = "GRevLex", cache_size::Int = 1000000, parallel::Bool = false, return_no_minkowski::Bool = false, return_ops::Bool = false)

Compute a monomial basis for the highest weight module with highest weight
``hw`` (in terms of the fundamental weights), for a simple Lie algebra of type
``t`` and rank ``n``.

# Parameters
- `t`: Explain
- `n`: Explain
- `hw`: Explain

# Examples
```jldoctest
julia> 1+1
2
```
"""
function basisLieHighestWeight2(t, n, hw; ops = "regular", known_monomials = [], monomial_order = "GRevLex", cache_size::Int = 1000000, parallel::Bool = false, return_no_minkowski::Bool = false, return_ops::Bool = false)
    # The function precomputes objects that are independent of the highest weight and can be used in all recursion steps. Then it starts the recursion and returns the result.

    # initialization
    L, CH = lieAlgebra(t, n) # L Lie Algebra of Type t,n and CH Chevalleybasis of L
    ops = get_ops(t, n, ops, L, CH) # operators that are represented by our monomials. x_i is connected to ops[i]
    wts = weightsForOperators(L, CH[3], ops) # weights of the operators
    wts = (v->Int.(v)).(wts)
    wts_eps = [w_to_eps(t, n, w) for w in wts] # this is another way of representing the weights and equivalent to wts, but some functionality is easier with those.
    calc_hw = Dict{Vector{Int}, Set{Vector{Int}}}([0 for i=1:n] => Set([[0 for i=1:n]]))
    no_minkowski = Set() # we save all highest weights, for which the Minkowski-sum did not suffice to gain all monomials

    # start recursion over hw
    res = compute_monomials(t, n, L, hw, ops, wts, wts_eps, monomial_order, calc_hw, cache_size, parallel, no_minkowski)
    
    # output
    if return_no_minkowski && return_ops
        return res, no_minkowski, ops
    elseif return_no_minkowski
        return res, no_minkowski
    elseif return_ops
        return res, ops
    else
        return res
    end
end

function sub_simple_refl(word, L, n)
    """
    substitute simple reflections (i,i+1), saved in dec by i, with E_{i,i+1}  
    """
    R = GAP.Globals.RootSystem(L)
    CG = fromGap(GAP.Globals.CanonicalGenerators(R)[1], recursive = false)
    ops = GAP.Obj([CG[i] for i in word], recursive = false)
    return ops
end

function get_ops(t,n, ops, L, CH)
    """
    handles user input for ops
    "regular" for all operators
    "longest-word" for random longest-word in Weyl-group
    ops::Vector{Int} for explicit longest-word
    """
    # create standard ops
    if ops == "regular" # use ops as specified by GAP
        ops = CH[1]
        return ops
    #elseif ops == "longest-word" # choose a random longest word. Created by extending by random not leftdescending reflections until total length is reached
    #    ops = longest_weyl_word(t,n)
    #    ops = sub_simple_refl(ops, L, n)
    #    return ops
    end

    # use user defined ops
    # wrong input
    if !(typeof(ops) == Vector{Int})
        println("ops needs to  of type Vector{Int}")
        return -1
    end
    if !(all([(1 <= i <= n) for i in ops]))
        println("all values of ops need to between 1 and the rank of the lie algebra.")
    end
    # If one of the conditions is met, the algorithms works for sure. Otherwise a warning is printed (and can be ignored).
    #if  !(is_longest_weyl_word(t, n, ops)) && !(Set(ops) == [i for i=1:n])
    #    println("WARNING: ops may be incorrect input.")
    #end
    ops = sub_simple_refl(ops, L, n)
    return ops
end

function compute_monomials(t, n, L, hw, ops, wts, wts_eps, monomial_order, calc_hw, cache_size::Int, parallel::Bool, no_minkowski)
    """
    This function calculates the monomial basis M_{hw} recursively. The recursion saves all computed results in calc_hw and we first check, if we already
    encountered this highest weight in a prior step. If this is not the case, we need to perform computations. The recursion works by using the
    Minkowski-sum. If M_{hw} is the desired set of monomials (identified by the exponents as lattice points), it is know that for lambda_1 + lambda_2 = hw
    we have M_{lambda_1} + M_{lambda_2} subseteq M_{hw}. The complexity grows exponentially in the size of hw. Therefore, it is very helpful to obtain
    a part of M_{hw} by going through all partitions of hw and using the Minkowski-property. The base cases of the recursion are the fundamental weights
    hw = [0, ..., 1, ..., 0]. In this case, or if the Minkowski-property did not find enough monomials, we need to perform the computations "by hand".
    """
    #println("compute_monomials: ", t, n, L, hw, ops, wts, wts_eps, monomial_order, calc_hw)
    # simple cases
    if haskey(calc_hw, hw) # we already computed the result in a prior recursion step
        return calc_hw[hw]
    elseif hw == [0 for i=1:n] # we mathematically know the solution
        return Set([0 for i=1:n])
    end
    
    # calculation required
    gapDim = GAP.Globals.DimensionOfHighestWeightModule(L, GAP.Obj(hw)) # number of monomials that we need to find, i.e. |M_{hw}|.
    # fundamental weights
    if is_fundamental(hw) # if hw is a fundamental weight, no partition into smaller summands is possible. This is the basecase of the recursion.
        push!(no_minkowski, hw)
        set_mon = add_by_hand(t, n, L, hw, ops, wts, wts_eps, monomial_order, gapDim, Set(), cache_size, parallel)
        push!(calc_hw, hw => set_mon)
        return set_mon
    else
        # use Minkowski-Sum for recursion
        set_mon = Set()
        i = 0
        sub_weights = compute_sub_weights(hw)
        l = length(sub_weights)
        # go through partitions lambda_1 + lambda_2 = hw until we have monomials enough monomials or used all partitions
        while length(set_mon) < gapDim && i < l 
            i += 1
            lambda_1 = sub_weights[i]
            lambda_2 = hw .- lambda_1
            mon_lambda_1 = compute_monomials(t, n, L, lambda_1, ops, wts, wts_eps, monomial_order, calc_hw, cache_size, parallel, no_minkowski)
            mon_lambda_2 = compute_monomials(t, n, L, lambda_2, ops, wts, wts_eps, monomial_order, calc_hw, cache_size, parallel, no_minkowski)
            mon_sum = Set([p .+ q for p in mon_lambda_1 for q in mon_lambda_2]) # Minkowski-sum: M_{lambda_1} + M_{lambda_2} subseteq M_{hw}
            union!(set_mon, mon_sum)
        end
        # check if we found enough monomials
        if length(set_mon) < gapDim
            push!(no_minkowski, hw)
            set_mon = add_by_hand(t, n, L, hw, ops, wts, wts_eps, monomial_order, gapDim, set_mon, cache_size, parallel)
        end
        push!(calc_hw, hw => set_mon)
        return set_mon
    end
end

function is_fundamental(hw)
    """
    returns true if hw is fundamental weight, i.e. hw = [0, ..., 1, ..., 0] (or hw = [0, ..., 0])
    """
    one = false
    for i in hw
        if i > 0
            if one || i > 1
                return false
            else 
                one = true
            end
        end
    end
    return true
end

function order_sub_weights(hw, sub_weights)
    """
    sort weights incrasing by their 2-norm
    """
    sort!(sub_weights, by=x->sum((x).^2))
end

function compute_sub_weights(hw)
    """
    returns list of weights w != 0 with 0 <= w <= hw elementwise
    """
    sub_weights = []
    foreach(Iterators.product((0:x for x in hw)...)) do i
        push!(sub_weights, [i...])
    end
    popfirst!(sub_weights) # [0, ..., 0]
    pop!(sub_weights) # hw
    order_sub_weights(hw, sub_weights)
    return sub_weights
end

function add_known_monomials!(weightspace, set_mon_in_weightspace, m, wts, mats, calc_monomials, space, e, cache_size)
    #println("add_known_monomials: ", weightspace, set_mon_in_weightspace, m, wts, mats, calc_monomials, space, e)
    #println("")
    """
    By using the Minkowski-sum, we know that all monomials in set_mon_in_weightspace are in our basis. Since we want to extend the weightspacse with
    missing monomials, we go need to calculate and add the vector of each monomial to our basis.
    """
    weight = weightspace[1]
    for mon in set_mon_in_weightspace[weight]
        #vec, wt = calc_new_mon!(mon, m, wts, mats, calc_monomials, space, e, cache_size)
        d = sz(mats[1])
        v0 = sparse_row(ZZ, [(1,1)]) # starting vector v
        vec = calc_vec(v0, mon, mats)
        wt = calc_wt(mon, wts)  
        #println("vec:" , vec)
        #println("wt: ", wt)
        if !haskey(space, wt)
            space[wt] = nullSpace()
        end
        addAndReduce!(space[wt], vec)
    end
    GC.gc()
end

function add_new_monomials!(t, n, mats, wts, monomial_order, weightspace, wts_eps, set_mon_in_weightspace, calc_monomials, space, e, cache_size, set_mon, m)
    """
    If a weightspace is missing monomials, we need to calculate them by trial and error. We would like to go through all monomials in the order monomial_order
    and calculate the corresponding vector. If it extends the basis, we add it to the result and else we try the next one. We know, that all monomials that work
    lie in the weyl-polytope. Therefore, we only inspect the monomials that lie both in the weyl-polytope and the weightspace. Since the weyl-polytope is bounded
    these are finitely many, we can sort them and then go trough them, until we found enough. 
    """
    #println("add_new_monomials: ", weightspace, set_mon_in_weightspace, m, wts, mats, calc_monomials, space, e)
    #println("")
    #println("memory: ", Int(Base.Sys.free_memory()) / 2^20)
    weight = weightspace[1]
    dim_weightspace = weightspace[2]
    
    # get monomials from weyl-polytope that are in the weightspace
    poss_mon_in_weightspace = get_monomials_of_weightspace(wts_eps, w_to_eps(t, n, weight), t)
    lt_function = lt_monomial_order(monomial_order)
    poss_mon_in_weightspace = sort(poss_mon_in_weightspace, lt=lt_function)
    
    # check which monomials should get added to the basis
    i=0
    if weight == 0 # check if [0 0 ... 0] already in basis
        i += 1
    end
    number_mon_in_weightspace = length(set_mon_in_weightspace[weight])
    # go through possible monomials one by one and check if it extends the basis
    while number_mon_in_weightspace < dim_weightspace
        i += 1

        # get a new mon
        if t in ["A", "G"] # necessary because of the structure of get_monomials_of_weightspace_An
            mon = convert(Vector{Int64}, poss_mon_in_weightspace[i][2:end])
        else
            mon = convert(Vector{Int64}, poss_mon_in_weightspace[i])
        end
        if mon in set_mon
            continue
        end

        # calculate the vector and check if it extends the basis
        #vec, _ = calc_new_mon!(mon, m, wts, mats, calc_monomials, space, e, cache_size)
        d = sz(mats[1])
        v0 = sparse_row(ZZ, [(1,1)]) # starting vector v
        vec = calc_vec(v0, mon, mats)
        #println("vec:" , vec)
        if !haskey(space, weight)
            space[weight] = nullSpace()
        end
        vec_red = addAndReduce!(space[weight], vec)
        if isempty(vec_red) # v0 == 0
            continue
        end

        # save monom
        number_mon_in_weightspace += 1
        push!(set_mon, mon)
    end

    # Weirdly, removing this causes memory problems. The lines should not be necessary since the objects are not referenced
    # and the garbage collector should run automatically. If I find out what causes this, it will be removed.
    weight = 0
    dim_weightspace = 0
    lt_function = 0
    number_mon_in_weightspace = 0
    poss_mon_in_weightspace = 0
    GC.gc()
end


function add_by_hand(t, n, L, hw, ops, wts, wts_eps, monomial_order, gapDim, set_mon, cache_size::Int, parallel::Bool)
    """
    This function calculates the missing monomials by going through each non full weightspace and adding possible monomials
    manually by computing their corresponding vectors and checking if they enlargen the basis.
    """
    #println("add_by_hand: ", t, n, L, hw, ops, wts, wts_eps, monomial_order, gapDim, set_mon)
    #println("")
    #println("")
    #println("add_by_hand: ", hw)
    # initialization
    mats = tensorMatricesForOperators(L, hw, ops) # matrices g_i for (g_1^a_1 * ... * g_k^a_k)*v
    m = length(mats)
    e = [(1:m .== i) for i in 1:m] # e_i
    space = Dict(0*wts[1] => nullSpace()) # span of basis vectors to keep track of the basis
    d = sz(mats[1])
    v0 = sparse_row(ZZ, [(1,1)])  # starting vector v
    calc_monomials = Dict{Vector{Int}, Tuple{TVec, Vector{Int}}}([0 for i in 1:m] => (v0, 0 * wts[1])) # saves the calculated vectors to decrease necessary matrix multiplicatons
    push!(set_mon, [0 for i in 1:m])
    weightspaces = get_dim_weightspace(t,n, L, hw) # required monomials of each weightspace

    # sort the monomials from the minkowski-sum by their weightspaces
    set_mon_in_weightspace = Dict{Vector{Int}, Set{Vector{Int}}}()
    for (weight, _) in weightspaces
        set_mon_in_weightspace[weight] = Set()
    end
    for mon in set_mon
        weight = calc_wt(mon, wts)
        push!(set_mon_in_weightspace[weight], mon)
    end

    # only inspect weightspaces with missing monomials
    full_weightspaces = zeros(Bool, length(weightspaces))
    for i=1:length(weightspaces)
        full_weightspaces[i] = (length(set_mon_in_weightspace[weightspaces[i][1]]) == weightspaces[i][2])
    end
    deleteat!(weightspaces, full_weightspaces)
    weightsapces = full_weightspaces
    

    # use parallel computations if parallel=true. The weightspaces could be calculated completely indepent. We just save
    # the computed monomials.
    # insert known monomials into basis
    if parallel
        @distributed for weightspace in weightspaces
            add_known_monomials!(weightspace, set_mon_in_weightspace, m, wts, mats, calc_monomials, space, e, cache_size)
        end
    else
        for weightspace in weightspaces
            #println("known memory: ", Int(Base.Sys.free_memory()) / 2^20)
            #println("size space: ", Int(Base.summarysize(space)) / 2^20)
            #println("size calc_monomials: ", Int(Base.summarysize(calc_monomials)) / 2^20)
            add_known_monomials!(weightspace, set_mon_in_weightspace, m, wts, mats, calc_monomials, space, e, cache_size)
        end 
    end

    # calculate new monomials
    if parallel
        @distributed for weightspace in weightspaces
            add_new_monomials!(t, n, mats, wts, monomial_order, weightspace, wts_eps, set_mon_in_weightspace, calc_monomials, space, e, cache_size, set_mon, m)
        end 
    else
        for weightspace in weightspaces
            #println("varinfo: ", InteractiveUtils.varinfo(MB3))
            #println("new memory: ", Int(Base.Sys.free_memory()) / 2^20)
            #println("size space: ", Int(Base.summarysize(space)) / 2^20)
            #println("size calc_monomials: ", Int(Base.summarysize(calc_monomials)) / 2^20)
            add_new_monomials!(t, n, mats, wts, monomial_order, weightspace, wts_eps, set_mon_in_weightspace, calc_monomials, space, e, cache_size, set_mon, m)
        end
    end
    #println("calc_monomials: ", length(calc_monomials))
    return set_mon
end

function get_dim_weightspace(t, n, L, hw)
    """
    calculates matrix, first row weights, second row dimension of corresponding weightspace
    GAP computes the dimension for all positive weights. The dimension is constant on orbits of the weylgroup,
    and we can therefore calculate the dimension of each weightspace
    """
    # calculate dimension for dominant weights with GAP
    R = GAP.Globals.RootSystem(L)
    W = fromGap(GAP.Globals.DominantCharacter(R, GAP.Obj(hw)))
    dominant_weights = W[1]
    dominant_weights_dim = W[2]
    dim_weightspace = []

    # calculate dimension for the rest by checking which positive weights lies in the orbit.
    for i=1:length(dominant_weights)
        orbit_weights = orbit_weylgroup(t,n, dominant_weights[i])
        dim = dominant_weights_dim[i]
        for weight in eachrow(orbit_weights)
            append!(dim_weightspace, [[hw-eps_to_w(t,n,weight), dim]])
        end
    end
    return dim_weightspace
end



end
