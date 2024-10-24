@doc raw"""
    integrate(c::CohomologyClass; check::Bool = true)

Integrate the cohomolgy class `c` over the normal
toric variety `toric_variety(c)`.

The theory underlying this method requires that the toric variety
in question is simplicial and complete. The check of completeness
may take a long time to complete. If desired, this can be switched
off by setting the optional argument `check` to the value `false`.

# Examples
```jldoctest
julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> (x1, x2, x3, e1, e2, e3) = gens(cohomology_ring(dP3))
6-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 x1
 x2
 x3
 e1
 e2
 e3

julia> c = cohomology_class(dP3, e3*e3 + e3)
Cohomology class on a normal toric variety given by e3^2 + e3

julia> integrate(c)
-1

julia> F3 = hirzebruch_surface(NormalToricVariety, 3)
Normal toric variety

julia> (x1, x2, x3, x4) = gens(cohomology_ring(F3))
4-element Vector{MPolyQuoRingElem{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 t1
 x1
 t2
 x2

julia> c = cohomology_class(F3, x1*x2 + x3*x4)
Cohomology class on a normal toric variety given by 2//3*x2^2

julia> integrate(c)
2
```

The following example constructs the Fano variety 2-36
(cf. https://www.fanography.info/2-36) and verifies
that the triple self-intersection number of its anticanonical
bundle is 62.

# Examples
```jldoctest
julia> e1 = [1,0,0];

julia> e2 = [0,1,0];

julia> e3 = [0,0,1];

julia> m = 2;

julia> ray_generators = [e1, -e1, e2, e3, - e2 - e3 - m * e1];

julia> max_cones = incidence_matrix([[1,3,4], [1,3,5], [1,4,5], [2,3,4], [2,3,5], [2,4,5]]);

julia> X = normal_toric_variety(max_cones, ray_generators; non_redundant = true)
Normal toric variety

julia> cox_ring(X)
Multivariate polynomial ring in 5 variables over QQ graded by
  x1 -> [1 0]
  x2 -> [1 2]
  x3 -> [0 -1]
  x4 -> [0 -1]
  x5 -> [0 -1]

julia> cohomology_ring(X)
Quotient
  of multivariate polynomial ring in 5 variables over QQ graded by
    x1 -> [1]
    x2 -> [1]
    x3 -> [1]
    x4 -> [1]
    x5 -> [1]
  by ideal (x1 - x2 - 2*x5, x3 - x5, x4 - x5, x1*x2, x3*x4*x5)

julia> integrate(cohomology_class(anticanonical_divisor(X))^3)
62

julia> integrate(cohomology_class(anticanonical_divisor_class(X))^3)
62
```
"""
function integrate(c::CohomologyClass; check::Bool = true)
    # can only integrate if the variety is simplicial, complete
    if check
      @req is_simplicial(toric_variety(c)) && is_complete(toric_variety(c)) "Integration only supported over complete and simplicial toric varieties"
    end

    # if the intersection form is known, we can use it
    if has_attribute(toric_variety(c), :_intersection_form_via_exponents)
        intersection_dict = _intersection_form_via_exponents(toric_variety(c))
        coeffs = coefficients(c)
        expos = exponents(c)
        integral = zero(QQ)
        for i in 1:nrows(expos)
            if expos[i, :] in keys(intersection_dict)
                integral += coeffs[i] * intersection_dict[expos[i, :]]
            end
        end
        return integral::QQFieldElem
    end
    
    # otherwise, proceed "by hand"
    if is_trivial(c)
        return zero(QQ)
    end
    poly = polynomial(c)
    dict = homogeneous_components(poly)
    elem = base_ring(parent(poly)).D([dim(toric_variety(c))])
    if !(elem in keys(dict))
        return zero(QQ)
    end
    top_form = dict[elem]
    if iszero(top_form)
        return zero(QQ)
    end
    n = AbstractAlgebra.leading_coefficient(top_form.f)
    m = AbstractAlgebra.leading_coefficient(polynomial(volume_form(toric_variety(c))).f)
    return n//m
end
