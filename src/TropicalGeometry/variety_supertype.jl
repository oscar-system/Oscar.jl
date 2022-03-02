###
# Tropical variety supertype in Oscar (not for public use)
# ===================================
###


###
# 0. Definition
# -------------
# M = typeof(min) or typeof(max):
#   min or max convention
# EMB = true or false:
#   embedded or abstract
# see also: variety.jl, hypersurface.jl, curve.jl, linear_space.jl
###

abstract type TropicalVarietySupertype{M,EMB} end
function pm_object(v::TropicalVarietySupertype)
  return v.polymakeTV
end


###
# 1. Basic constructions
# ----------------------
###

@doc Markdown.doc"""
    intersect(T1, T2)

Intersect two tropical varieties.

# Examples
```jldoctest
julia> RR = tropical_semiring(min)
Tropical ring (min)

julia> S,(x,y) = RR["x","y"]
(Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}[x, y])

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
function intersect(T1::TropicalVarietySupertype{M, EMB}, T2::TropicalVarietySupertype{M, EMB}) where {M, EMB}
    return TropicalVariety{M, EMB}(intersect(T1.polyhedralComplex,T2.polyhedralComplex))
end


@doc Markdown.doc"""
    intersect_stably(T1, T2)

# Examples
```jldoctest
```
"""
function stably_intersect(T1::TropicalVarietySupertype{M, EMB}, T2::TropicalVarietySupertype{M, EMB}) where {M, EMB}
    return TropicalVariety{M, EMB}(intersect_stably(T1.polyhedralComplex,T2.polyhedralComplex))
end
export intersect_stably



###
# 2. Basic properties
# -------------------
###

@doc Markdown.doc"""
    ambient_dim(T::TropicalVariety{M, EMB})
    ambient_dim(T::TropicalCurve{M, EMB})
    ambient_dim(T::TropicalHypersurface{M, EMB})
    ambient_dim(T::TropicalLinearSpace{M, EMB})

Returns the ambient dimension of `T` if it is embedded. Returns an error otherwise.

# Examples
A tropical hypersurface in RR^n is of ambient dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> ambient_dim(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function ambient_dim(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    if !EMB
        error("ambient_dim: tropical variety not embedded")
    end

    return ambient_dim(T.polyhedralComplex)
end



@doc Markdown.doc"""
    codim(T::TropicalVariety{M, EMB})
    codim(T::TropicalCurve{M, EMB})
    codim(T::TropicalHypersurface{M, EMB})
    codim(T::TropicalLinearSpace{M, EMB})

Returns the codimension of `T`.

# Examples
A tropical hypersurface in RR^n is always of dimension n-1
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> codim(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function codim(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return codim(T.polyhedralComplex)
end



@doc Markdown.doc"""
    dim(T::TropicalVariety{M, EMB})
    dim(T::TropicalCurve{M, EMB})
    dim(T::TropicalHypersurface{M, EMB})
    dim(T::TropicalLinearSpace{M, EMB})

Returns the dimension of `T`.

# Examples
A tropical hypersurface in RR^n is always of dimension n-1
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> dim(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function dim(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return dim(T.polyhedralComplex)
end



@doc Markdown.doc"""
    f_vector(T::TropicalVariety{M, EMB})
    f_vector(T::TropicalCurve{M, EMB})
    f_vector(T::TropicalHypersurface{M, EMB})
    f_vector(T::TropicalLinearSpace{M, EMB})

Returns the f-Vector of `T`.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> f_vector(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function f_vector(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return f_vector(T.polyhedralComplex)
end



@doc Markdown.doc"""
    lineality_dim(T::TropicalVariety{M, EMB})
    lineality_dim(T::TropicalCurve{M, EMB})
    lineality_dim(T::TropicalHypersurface{M, EMB})
    lineality_dim(T::TropicalLinearSpace{M, EMB})

Returns the dimension of the lineality space of `T` if it is embedded. Returns an error otherwise.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y;

julia> tropicalAndAffineLine = TropicalHypersurface(f);

julia> lineality_dim(tropicalAndAffineLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function lineality_dim(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    if !EMB
        error("lineality_dim: tropical variety not embedded")
    end

    return lineality_dim(T.polyhedralComplex)
end



@doc Markdown.doc"""
    lineality_space(T::TropicalVariety{M, EMB})
    lineality_space(T::TropicalCurve{M, EMB})
    lineality_space(T::TropicalHypersurface{M, EMB})
    lineality_space(T::TropicalLinearSpace{M, EMB})

Returns the lineality space of `T` if it is embedded. Returns an error otherwise.

# Examples
A tropical hypersurface in RR^n is of lineality spaceension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y;

julia> tropicalAndAffineLine = TropicalHypersurface(f);

julia> lineality_space(tropicalAndAffineLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function lineality_space(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    if !EMB
        error("lineality_space: tropical variety not embedded")
    end

    return lineality_space(T.polyhedralComplex)
end



@doc Markdown.doc"""
    maximal_polyhedra(T::TropicalVariety{M, EMB})
    maximal_polyhedra(T::TropicalCurve{M, EMB})
    maximal_polyhedra(T::TropicalHypersurface{M, EMB})
    maximal_polyhedra(T::TropicalLinearSpace{M, EMB})

Returns the maximal polyhedra of `T`.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> maximal_polyhedra(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function maximal_polyhedra(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return maximal_polyhedra(T.polyhedralComplex)
end



@doc Markdown.doc"""
    nmaximal_polyhedra(T::TropicalVariety{M, EMB})
    nmaximal_polyhedra(T::TropicalCurve{M, EMB})
    nmaximal_polyhedra(T::TropicalHypersurface{M, EMB})
    nmaximal_polyhedra(T::TropicalLinearSpace{M, EMB})

Returns the number of maximal polyhedra of `T`.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> nmaximal_polyhedra(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
# todo: do maximal polyhedra at infinity count?
```
"""
function nmaximal_polyhedra(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return nmaximal_polyhedra(T.polyhedralComplex)
end



@doc Markdown.doc"""
    npolyhedra(T::TropicalVariety{M, EMB})
    npolyhedra(T::TropicalCurve{M, EMB})
    npolyhedra(T::TropicalHypersurface{M, EMB})
    npolyhedra(T::TropicalLinearSpace{M, EMB})

Returns the number of polyhedra of `T`.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> npolyhedra(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
# todo: do polyhedra at infinity count?
```
"""
function npolyhedra(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return npolyhedra(T.polyhedralComplex)
end



@doc Markdown.doc"""
    nvertices(T::TropicalVariety{M, EMB})
    nvertices(T::TropicalCurve{M, EMB})
    nvertices(T::TropicalHypersurface{M, EMB})
    nvertices(T::TropicalLinearSpace{M, EMB})

Returns the number of vertices of `T`.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> nvertices(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
# todo: do vertices at infinity count?
```
"""
function nvertices(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return nvertices(T.polyhedralComplex)
end



@doc Markdown.doc"""
    polyhedra(T::TropicalVariety{M, EMB})
    polyhedra(T::TropicalCurve{M, EMB})
    polyhedra(T::TropicalHypersurface{M, EMB})
    polyhedra(T::TropicalLinearSpace{M, EMB})

Returns the polyhedra of `T`.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> polyhedra(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
# todo: do vertices at infinity count?
```
"""
function polyhedra(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return polyhedra(T.polyhedralComplex)
end



@doc Markdown.doc"""
    ispure(T::TropicalVariety{M, EMB})
    ispure(T::TropicalCurve{M, EMB})
    ispure(T::TropicalHypersurface{M, EMB})
    ispure(T::TropicalLinearSpace{M, EMB})

Return true if `T` is a pure polyhedral complex, false otherwise.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> ispure(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function ispure(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return ispure(T.polyhedralComplex)
end



@doc Markdown.doc"""
    issimplicial(T::TropicalVariety{M, EMB})
    issimplicial(T::TropicalCurve{M, EMB})
    issimplicial(T::TropicalHypersurface{M, EMB})
    issimplicial(T::TropicalLinearSpace{M, EMB})

Returns true if `T` is a simplicial polyhedral complex, false otherwise.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> issimplicial(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function issimplicial(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    return issimplicial(T.polyhedralComplex)
end



@doc Markdown.doc"""
    vertices(T::TropicalVariety{M, EMB})
    vertices(T::TropicalCurve{M, EMB})
    vertices(T::TropicalHypersurface{M, EMB})
    vertices(T::TropicalLinearSpace{M, EMB})

Returns the vertices of `T`, which are points in euclidean space if T is embedded or elements in an ordered set otherwise.

# Examples
The vertices of a plane tropical line, plane tropical honeycomb quadric, and plane tropical honeycomb cubic
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f1 = x+y+1;

julia> tropicalLine = TropicalHypersurface(f1);

julia> vertices(tropicalLine)

julia> f2 = 1*x^2+x*y+1*y^2+x+y+1;

julia> tropicalQuadric = TropicalHypersurface(f1);

julia> vertices(tropicalQuadric)

julia> f3 = x^3+x*y^2+x^2*y+y^3+x^2+x*y+y^2+x+y+1;

julia> tropicalCubic = TropicalHypersurface(f3);

julia> vertices(tropicalCubic)
```
"""
function vertices(as::Type{PointVector{S}}, T::TropicalVarietySupertype{M,EMB}) where {S,M,EMB}
    return vertices(as,T.polyhedralComplex)
end

function vertices(T::TropicalVarietySupertype{M, EMB}) where {M,EMB}
    return vertices(T.polyhedralComplex)
end



@doc Markdown.doc"""
    weights(T::TropicalVariety{M, EMB})
    weights(T::TropicalCurve{M, EMB})
    weights(T::TropicalHypersurface{M, EMB})
    weights(T::TropicalLinearSpace{M, EMB})

Returns the weights of `T`.

# Examples
A tropical hypersurface in RR^n is of lineality dimension n
```jldoctest
julia> RR = tropical_semiring(min);

julia> S,(x,y) = RR["x","y"];

julia> f = x+y+1;

julia> tropicalLine = TropicalHypersurface(f);

julia> weights(tropicalLine)
# todo: add examples for varieties, curves and linear spaces
```
"""
function weights(T::TropicalVarietySupertype{M,EMB}) where {M,EMB}
    if !has_attribute(T,:weights)
        error("weights: no weights attributed")
    end
    return get_attribute(T,:weights)
end



# @doc Markdown.doc"""
#     PolyhedralComplex(TV::TropicalVarietySupertype)

# Return the underlying polyhedral complex.

# # Examples
# ```jldoctest
# julia> RR = tropical_semiring(min)
# Tropical ring (min)

# julia> S,(x,y) = RR["x","y"]
# (Multivariate Polynomial Ring in x, y over Tropical ring (min), AbstractAlgebra.Generic.MPoly{Oscar.TropicalSemiringElem{typeof(min)}}[x, y])

# julia> f = x+y+1
# x + y + (1)

# julia> hyp = TropicalHypersurface(f)
# A min tropical hypersurface embedded in 2-dimensional Euclidian space

# julia> pc = PolyhedralComplex(hyp)
# A polyhedral complex in ambient dimension 3

# julia> vertices(pc)
# 4-element SubObjectIterator{Union{PointVector{Polymake.Rational}, RayVector{Polymake.Rational}}}:
#  [0, -1, -1]
#  [0, 1, 0]
#  [0, 0, 1]
#  [0, 1, 1]

# julia> for v in vertices(pc)
#        println(typeof(v))
#        end
# RayVector{Polymake.Rational}
# RayVector{Polymake.Rational}
# RayVector{Polymake.Rational}
# PointVector{Polymake.Rational}
# ```
# """
# function PolyhedralComplex(TV::TropicalVarietySupertype{M, EMB}) where {M,EMB}
#     return PolyhedralComplex(pm_object(TV))
# end
