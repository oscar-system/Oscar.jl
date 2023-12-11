# Add your new types, functions, and methods here.

export origami, veech_group, GapObj, vertical_perm, horizontal_perm, stratum, index_monodromy_group,
        sum_of_lyapunov_exponents, translations, is_hyperelliptic, cylinder_structure, veech_group_and_orbit,
        veech_group_is_even, are_equivalent

module OrigamiHelper

using ..GAP

function __init__()
    GAP.Packages.install("https://ag-weitze-schmithusen.github.io/ModularGroup/PackageInfo.g")
    GAP.Packages.install("https://ag-weitze-schmithusen.github.io/Origami/PackageInfo.g")
    GAP.Packages.load("Origami")
end

end

include("types.jl")

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
    print(io, "Origami ($(horizontal_perm(o)),$(vertical_perm(o)), $(degree(o)))")
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
    return GAP.Globals.TranslationsOfOrigami(GapObj(o))
end

function is_hyperelliptic(o::Origami)
    return GAP.Globals.IsHyperelliptic(GapObj(o))
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