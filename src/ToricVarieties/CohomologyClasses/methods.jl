@doc Markdown.doc"""
    integrate(c::CohomologyClass)

Integrate the cohomolgy class `c` over the normal
toric variety `toric_variety(c)`.

# Examples
```jldoctest
julia> dP3 = del_pezzo(3)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> (x1,e1,x2,e3,x3,e2) = gens(cohomology_ring(dP3))
6-element Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}:
 x1
 e1
 x2
 e3
 x3
 e2

julia> c = CohomologyClass(dP3, e3*e3 + e3)
A cohomology class on a normal toric variety given by e3 + e2^2

julia> integrate(c)
-1

julia> F3 = hirzebruch_surface(3)
A normal, non-affine, smooth, projective, gorenstein, non-fano, 2-dimensional toric variety without torusfactor

julia> (x1,x2,x3,x4) = gens(cohomology_ring(F3))
4-element Vector{MPolyQuoElem{MPolyElem_dec{fmpq, fmpq_mpoly}}}:
 t1
 x1
 t2
 x2

julia> c = CohomologyClass(F3, x1*x2 + x3*x4)
A cohomology class on a normal toric variety given by 2//3*x2^2

julia> integrate(c)
2
```
"""
function integrate(c::CohomologyClass)::fmpq
    # can only integrate if the variety is simplicial, complete
    if !issimplicial(toric_variety(c)) || !iscomplete(toric_variety(c))
        throw(ArgumentError("Integration only supported over complete and simplicial toric varieties."))
    end

    # if the intersection form is known, we can use it
    if has_attribute(toric_variety(c),:_intersection_form_via_exponents)
        intersection_dict = _intersection_form_via_exponents(toric_variety(c))
        coeffs = coefficients(c)
        expos = exponents(c)
        integral = 0
        for i in 1:nrows(expos)
            if expos[i,:] in keys(intersection_dict)
                integral += coeffs[i] * intersection_dict[expos[i,:]]
            end
        end
        return integral
    end
    
    # otherwise, proceed "by hand"
    if istrivial(c)
        return 0
    end
    poly = polynomial(c)
    dict = homogeneous_components(poly)
    elem = parent(poly).R.D([dim(toric_variety(c))])
    if !(elem in keys(dict))
        return 0
    end
    top_form = dict[elem]
    if iszero(top_form)
        return 0
    end
    n = leading_coefficient(top_form.f)
    m = leading_coefficient(polynomial(volume_form(toric_variety(c))).f)
    return fmpq(n//m)
end
export integrate
