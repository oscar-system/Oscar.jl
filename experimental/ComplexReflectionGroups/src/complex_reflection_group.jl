# This file implements explicit models of complex reflection groups.
#
# References:
#
# * Lehrer, G. I., & Taylor, D. E. (2009). Unitary reflection groups (Vol. 20, p. viii). Cambridge University Press, Cambridge.
# 
# * Thiel, U. (2014). On restricted rational Cherednik algebras. TU Kaiserslautern.
#
# [Ulrich Thiel](https://ulthiel.com/math), 2023 

export complex_reflection_group
export is_complex_reflection_group
export complex_reflection_group_type
export complex_reflection_group_model

function complex_reflection_group(G::ComplexReflectionGroupType; model=nothing)
    
    # this will be the list of matrix groups corresponding to the components of G
    component_groups = MatrixGroup[]

    # list of models
    modellist = []

    for C in components(G)

        t = C.type[1]
        gens = nothing
        Cmodel = nothing

        # we now create the list "gens" of matrix generators for C in the selected model
        if length(t) == 1
            return true
        elseif length(t) == 3
            (m,p,n) = t

            # Symmetric group case needs special care
            if m == 1 && p == 1
                if model === nothing
                    Cmodel = :CHEVIE
                else
                    Cmodel = model
                end
                if Cmodel === :CHEVIE
                    # the generators correspond to the transpositions r_i=(i,i+1)
                    # but note that we consider the irreducible representation of 
                    # dimension n-1, so it's not simply the transposition matrix
                    matspace = matrix_space(QQ, n-1, n-1)
                    gens = elem_type(matspace)[]
                    for i = 1:n-1
                        ri = matspace(1)
                        ri[i,i] = -1
                        if i > 1
                            ri[i,i-1] = 1
                        end
                        if i < n-1
                            ri[i,i+1] = 1
                        end
                        push!(gens, ri)
                    end
                end
            else
                if model === nothing
                    Cmodel = :Magma
                else
                    Cmodel = model
                end
                if Cmodel === :Magma
                    # See Lehrer & Taylor (2009), p 35-36
                    # What's implemented in Magma deviates slightly from Lehrer & Taylor, 
                    # see the remarks below.
                    K, z = cyclotomic_field(m)
                    matspace = matrix_space(K, n, n)

                    #the transpositions
                    transp = elem_type(matspace)[]
                    for i = 1:n-1
                        r = matspace(1)
                        r[i,i] = 0
                        r[i+1,i+1] = 0
                        r[i,i+1] = 1
                        r[i+1,i] = 1
                        push!(transp, r)
                    end

                    #the matrix t
                    t = matspace(1)
                    t[n,n] = z #LT would take t[1,1]=z here

                    #the matrix s
                    if n > 1
                        s = t^-1*transp[n-1]*t #LT would take transp[1] here
                    end

                    #now, create the list of generators
                    #LT would order the transpositions to the end
                    if m == p
                        gens = vcat(transp, [s])
                    elseif p == 1
                        gens = vcat(transp, [t])
                    else
                        gens = vcat(transp, [s], [t^p])
                    end
                end
            end
        end

        # create matrix group from the generators
        if gens === nothing
            throw(ArgumentError("Specified model not found"))
        end
        matgrp = matrix_group(gens)

        # set attributes that are already known from type
        set_attribute!(matgrp, :order, order(C))
        set_attribute!(matgrp, :is_complex_reflection_group, true)
        set_attribute!(matgrp, :complex_reflection_group_type, C)
        set_attribute!(matgrp, :complex_reflection_group_model, Cmodel)

        # add to list
        push!(component_groups, matgrp)
        push!(modellist, Cmodel)
    end

    if length(component_groups) == 1
        return component_groups[1]
    else
        matgrp = direct_product(component_groups)

        set_attribute!(matgrp, :order, order(G))
        set_attribute!(matgrp, :is_complex_reflection_group, true)
        set_attribute!(matgrp, :complex_reflection_group_type, G)
        set_attribute!(matgrp, :complex_reflection_group_model, modellist)
    end

end

# Convenience constructors
function complex_reflection_group(i::Int; model=nothing)
    return complex_reflection_group(ComplexReflectionGroupType(i); model=model)
end

function complex_reflection_group(t::Tuple{Int}; model=nothing)
    return complex_reflection_group(ComplexReflectionGroupType(t); model=model)
end

function complex_reflection_group(m::Int, p::Int, n::Int; model=nothing)
    return complex_reflection_group(ComplexReflectionGroupType(m,p,n); model=model)
end

function complex_reflection_group(t::Tuple{Int,Int,Int}; model=nothing)
    return complex_reflection_group(ComplexReflectionGroupType(t); model=model)
end

function is_complex_reflection_group(G::MatrixGroup)
    if has_attribute(G, :is_complex_reflection_group)
        return get_attribute(G, :is_complex_reflection_group)
    end
    return false 
    #this should be upgraded later to work with a general matrix group
end

function complex_reflection_group_type(G::MatrixGroup)
    if has_attribute(G, :complex_reflection_group_type)
        return get_attribute(G, :complex_reflection_group_type)
    end
    return nothing
    # this should be upgraded later to work with a general matrix group
end

function complex_reflection_group_model(G::MatrixGroup)
    if has_attribute(G, :complex_reflection_group_model)
        return get_attribute(G, :complex_reflection_group_model)
    end
    return nothing
end