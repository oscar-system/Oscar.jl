# Add your new types, functions, and methods here.


module OrigamiHelper

using ..GAP

function __init__()
    mod_p = "https://ag-weitze-schmithusen.github.io/ModularGroup/PackageInfo.g"
    ori_p = "https://ag-weitze-schmithusen.github.io/Origami/PackageInfo.g"
    GAP.Packages.install(mod_p)
    GAP.Packages.install(ori_p)
    GAP.Packages.load("Origami")
end

end

include("types.jl")
include("canonical.jl")
include("deck_group.jl")
include("action.jl")
include("special_origami.jl")
include("generate.jl")
include("homology_action.jl")
include("cyclic_torus_covers.jl")
include("homology.jl")
include("normal_origami.jl")
include("systoles.jl")

@doc raw"""
    origami(h::PermGroupElem, v::PermGroupElem)

blablabla
# Examples
```jldoctest
julia> h = @perm (1,2)
(1,2)

julia> v = @perm (1,2)
(1,2)

julia> o = origami(h,v)
Origami ((1,2),(1,2))
```

"""
function origami(h::PermGroupElem, v::PermGroupElem)
    deg = max(degree(h), degree(v))
    return Origami(GAP.Globals.Origami(GapObj(h), GapObj(v)), h, v, deg)
    # TODO check for transitivity? GAP already does this and this was so slow
    # for huge origamis
    #perm_group = permutation_group(deg, [h, v])
    #if transitivity(perm_group) > 0
        
    #end
end

# ugly workaround
function from_GAP_origami(o::GapObj)
    d = GAP.Globals.DegreeOrigami(o)::Int
    G = symmetric_group(d)
    h = PermGroupElem(G, GAP.Globals.HorizontalPerm(o))
    v = PermGroupElem(G, GAP.Globals.VerticalPerm(o))
    return origami(h, v)
end

function origami_disconnected(h::PermGroupElem, v::PermGroupElem, d::Integer)
    return Origami(GAP.Globals.OrigamiNC(GapObj(h), GapObj(v), d), h, v, d)
end

function Base.:(==)(a::Origami, b::Origami)
    # TODO rewrite this? for now use Gap equality
    return (a.h == b.h) && (a.v == b.v)
end

function Base.hash(o::Origami, h::UInt=0x000000000)
    return hash(o.h, hash(o.v, hash(o.d, h)))
end

function horizontal_perm(o::Origami)
    return o.h
end

function vertical_perm(o::Origami)
    return o.v
end

function degree(o::Origami)
    return o.d
end

function Base.show(io::IO, o::Origami)
    h = horizontal_perm(o)
    v = vertical_perm(o)
    d = degree(o)
    print(io, "Origami ($(h),$(v), $(d))")
end

GapObj(O::Origami) = O.o

function stratum(o::Origami)
    # cannot use comm from Oscar because of different definitions
    h = horizontal_perm(o)
    v = vertical_perm(o)
    commutator = h * v * h^-1 * v^-1
    cycle_struc = cycle_structure(commutator)
    different_cycle_count = length(cycle_struc)
    stratum::Vector{Integer} = []
    for i in 1:different_cycle_count
        cycle_type_entry = cycle_struc[i]
        cycle_type_length = cycle_struc[i][1]
        if(cycle_type_entry[1] == 1)
            # ignore cycles of length 1
            continue
        end
        for j in 1:cycle_type_entry[2]
            push!(stratum, ZZ(cycle_type_length - 1))
        end
    end
    return stratum
end

function genus(o::Origami)
    return ZZ((sum(stratum(o)) + 2) * 0.5)
end

function veech_group(O::Origami)
    # TODO use hashtables or not? Implement modular subgroup?
    GAP.Globals.ComputeVeechGroupWithHashTables(GapObj(O))
end

function index_monodromy_group(o::Origami)
    GAP.Globals.IndexOfMonodromyGroup(GapObj(o))
end

function sum_of_lyapunov_exponents(o::Origami)
    gap_obj = GAP.Globals.SumOfLyapunovExponents(GapObj(o))
    return GAP.gap_to_julia(gap_obj)
end

function translations(o::Origami)
    h = horizontal_perm(o)
    v = vertical_perm(o)
    G = normalform_conjugators(o)

    O = (i -> origami(i^-1 * h * i, i^-1 * v * i)).(G)
    list = [cperm([degree(o)])]

    l = length(O)
    for i in 1:l
        positions = findall(item -> item == O[i], O)
        if length(positions) != 1
            for j in positions
                push!(list, G[i] * G[j]^-1)
            end
        end
    end

    return collect(Set(list))
end

function is_hyperelliptic(o::Origami)
    # check whether -1 is in the veech group
	if !veech_group_is_even(o) 
        return false
    end

    n = 2
    n_inv = 1 / 2
	x = horizontal_perm(o)
	y = vertical_perm(o)
	g = genus(o)

    L = point_reflections(o)
	L = filter(i -> order(i) == 2, L)

	if L == []
        return false
    end

    # helper to find all elements in x that are not in y
    function list_diff(x, y)
        set_1 = Set(x)
        set_2 = Set(y)
        diff = setdiff(set_1, set_2)
        return collect(diff)
    end

    degree_list = collect(1:degree(o))
    for sigma in L
        # fixpoints
        b = 0
        b = b + length(list_diff(degree_list, moved_points(sigma)))
        b = b + length(list_diff(degree_list, moved_points(sigma*x)))
        b = b + length(list_diff(degree_list, moved_points(sigma*y)))

        for i in degree_list
            if i^(sigma * x^-1 * y^-1) == i^(y*x*(x*y)^-1)
                b = b + 1
            end
        end

        # TODO can n_inv yields floating point errors?
        if n_inv * (g - 1 - b * 0.5) + 1 == 0
            return true
        end
    end

    return false
end

function cylinder_structure(o::Origami)
    # TODO returns Vector{Any}, maybe cast to Integer?
    gap_obj = GAP.Globals.CylinderStructure(GapObj(o))
    return GAP.gap_to_julia(gap_obj)
end

function veech_group_and_orbit(o::Origami)
    return GAP.Globals.VeechGroupAndOrbit(GapObj(o))
end

function veech_group_is_even(o::Origami)
    return GAP.Globals.VeechGroupIsEven(GapObj(o))
end

function are_equivalent(o1::Origami, o2::Origami)
    return GAP.Globals.OrigamisEquivalent(GapObj(o1), GapObj(o2))
end

# TODO why is this not in canonical.jl?
function normalform_conjugators(o::Origami)
    x = horizontal_perm(o)
    y = vertical_perm(o)
    n = degree(o)
    G = []

    # Starting from each of the vertices found above, do a breadth-first search
	# and list the vertices in the order they appear.
	# This defines a permutation l with which we conjugate x and y.
	# From the resulting list of pairs of permutations (all of which are by
	# definition simultaneously conjugated to (x,y)) we choose the
	# lexicographically smallest one as the canonical form.
    for i in 1:n
        L = fill(0, n)
        seen = fill(false, n)
        Q = [i]
        seen[i] = true
        numSeen = 1
        L[i] = 1
        while numSeen < n
            v = popfirst!(Q)
            wx = v^x
            wy = v^y
            if !seen[wx]
                push!(Q, wx)
                seen[wx] = true
                numSeen = numSeen + 1
                L[wx] = numSeen
            end
            if !seen[wy]
                push!(Q, wy)
                seen[wy] = true
                numSeen = numSeen + 1
                L[wy] = numSeen
            end
        end
        push!(G, L)
    end

    return perm.(G)
end

function point_reflections(o::Origami)
    if !veech_group_is_even(o)
        throw("VeechGroup must contain -1")
    end

    h = horizontal_perm(o)
    v = vertical_perm(o)
    o1 = origami(h^-1, v^-1)
    h1 = horizontal_perm(o1)
    v1 = vertical_perm(o1)
    G = normalform_conjugators(o)
    G1 = normalform_conjugators(o1)

    f(i) = origami(i^-1 * h * i, i^-1 * v * i)
    f1(i) = origami(i^-1 * h1 * i, i^-1 * v1 * i)

    # origamis derived from the permutations above
    O = f.(G)
    # we need to calculate these to test find k s.t.
    # sigma_i *origami *sigma_i^-1=delta_k(i)*origami_1*delta_k(i)^-1
    O1 = f1.(G1)

    # fitting the permuations together
    result::Vector{PermGroupElem} = []

    l = length(O)
    for i in 1:l
        index = findfirst(item -> item == O[i], O1)
        push!(result, G[i] * (G1[index])^-1)
    end

    # Remove duplicates - TODO better solution?
    return collect(Set(result))
end

function automorphisms(o::Origami)
    return [[translations(o),1], [point_reflections(o),-1]]
end

export origami, veech_group, GapObj, vertical_perm, horizontal_perm, stratum,
        index_monodromy_group, sum_of_lyapunov_exponents, translations,
        is_hyperelliptic, cylinder_structure, veech_group_and_orbit,
        veech_group_is_even, are_equivalent, normalform_conjugators,
        point_reflections, automorphisms, in_deck_group, deck_group, is_normal,
        origami_disconnected, action_s, action_t, action_s_inv, action_t_inv,
        action_sl2, random_origami, staircase_origami, CylinderDiagram, compute_rays,
        all_combinations, realizable_lengths, origami_from_cylinder_coordinates, product_gray_code,
        cylinders, system_of_equations, origamis_in_h11, cylinder_diagrams_h11,
        possible_lengths_and_heights, partition_degree, realizable_lengths_of_cylinder_diagram,
        x_origami, elevator_origami, homology, non_taut_part_of_homology,
        action_of_matrix_on_non_taut, shadow_veech_group, homology_to_string,
        action_of_matrix_on_homology, generalized_cyclic_torus_cover, comb_origami, cyclic_torus_cover_origamiS,
        cyclic_torus_cover_origamiL, base_change_l_to_s, translation_group_on_homology_of_tn, action_of_t_on_homology_of_tn,
        action_of_s_on_homology_of_tn, action_of_matrix_on_homology_of_tn, symplectic_basis_of_homology, has_spin_structure,
        spin_parity, normal_stored_origami, as_permutation_pepresentation, all_normal_origamis_by_degree,
        all_normal_origamis_from_group, read_cylinder_diagrams, split_diagrams, parse_cycle, origamis, extract_permutation,
        from_GAP_origami, systolic_ratio, systolic_ratio_bigger_one_over_pi_in_h11
