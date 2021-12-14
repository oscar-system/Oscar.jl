@doc Markdown.doc"""
    intersect(TV1, TV2)

Intersect two tropical varieties.

# Examples
```jldoctest
julia> T = tropical_numbers(min)
Tropical ring (min)

julia> Txy,(x,y) = T["x","y"]
(Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalNumbersElem{typeof(min)}}[x, y])

julia> f1 = x+y+1
x + y + (1)

julia> f2 = x^2+y^2+T(-6)
x^2 + y^2 + (-6)

julia> hyp1 = TropicalHypersurface(f1)
A min tropical hypersurface embedded in 2-dimensional Euclidian space

julia> hyp2 = TropicalHypersurface(f2)
A min tropical hypersurface embedded in 2-dimensional Euclidian space

julia> tv12 = intersect(hyp1, hyp2)
A min tropical variety of dimension 1 embedded in 2-dimensional Euclidian space
```
"""
function intersect(TV1::TropicalVarietySupertype{M, EMB}, TV2::TropicalVarietySupertype{M, EMB}) where {M, EMB}
    pm_tv1 = pm_object(TV1)
    pm_tv2 = pm_object(TV2)
    result = Polymake.fan.PolyhedralComplex(Polymake.fan.common_refinement(pm_tv1, pm_tv2))
    result = polyhedral_complex_workaround(result)
    return TropicalVariety{M, EMB}(result)
end


@doc Markdown.doc"""
    stably_intersect(TV1, TV2)

# Examples
```jldoctest
```
"""
function stably_intersect(TV1::TropicalVarietySupertype{M, EMB}, TV2::TropicalVarietySupertype{M, EMB}) where {M, EMB}
    pm_tv1 = pm_object(TV1)
    pm_tv2 = pm_object(TV2)
    result = Polymake.fan.common_refinement(pm_tv1, pm_tv2)
    k = dim(TV1) + dim(TV2) - ambient_dim(TV1)
    result = Polymake.fan.PolyhedralComplex(Polymake.fan.k_skeleton(result, k+1))
    result = polyhedral_complex_workaround(result)
    return TropicalVariety{M, EMB}(result)
end
export stably_intersect
