################################################################################
################################################################################
## Cyclic Quotient Singularit struct
################################################################################
################################################################################

struct CyclicQuotientSingularity <: AbstractNormalToricVariety
    polymakeNTV::Polymake.BigObject
end
export CyclicQuotientSingularity


################################################################################
################################################################################
## Constructor
################################################################################
################################################################################

@doc Markdown.doc"""
    CyclicQuotientSingularity(n::Int64, q::Int64)

Return the cyclic quotient singularity for the parameters $n$ and $q$, with
$0<q<n$ and $q,n$ coprime.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7,5)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> isaffine(cqs)
true

julia> issmooth(cqs)
false
```
"""
function CyclicQuotientSingularity(n::Int64, q::Int64)
    n > 0 || error("n (=$(n)) must be positive")
    q > 0 || error("q (=$(q)) must be positive")
    q < n || error("q must be smaller than n (q=$(q) >= n=$(n))")
    gcd(n,q)==1 || error("n and q must be coprime (gcd=$(gcd(n,q)))")
    pmntv = Polymake.fulton.CyclicQuotient(N=n, Q=q)
    return CyclicQuotientSingularity(pmntv)
end


@doc Markdown.doc"""
    continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)

Return the Hirzebruch-Jung continued fraction associated with the cyclic
quotient singularity, i.e. the Hirzebruch-Jung continued fraction corresponding
to $n/q$.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1,c_2,\ldots,c_n]$ is $r([c_1,c_2,\ldots,c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots,c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7,5)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> cf = continued_fraction_hirzebruch_jung(cqs)
3-element Vector{Int64}:
 2
 2
 3

julia> ecf = cf[1]-1/(cf[2]-Rational(1,cf[3]))
7//5
```
"""
continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity) = Vector{Int64}(pm_ntv(cqs).CONTINUED_FRACTION)
export continued_fraction_hirzebruch_jung


@doc Markdown.doc"""
    dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)

Return the dual Hirzebruch-Jung continued fraction associated with the cyclic
quotient singularity, i.e. the Hirzebruch-Jung continued fraction corresponding
to $q/(n-q)$.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1,c_2,\ldots,c_n]$ is $r([c_1,c_2,\ldots,c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots,c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7,5)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> dcf = dual_continued_fraction_hirzebruch_jung(cqs)
2-element Vector{Int64}:
 4
 2

julia> edcf = dcf[1] - Rational(1,dcf[2])
7//2
```
"""
dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity) = Vector{Int64}(pm_ntv(cqs).DUAL_CONTINUED_FRACTION)
export dual_continued_fraction_hirzebruch_jung


@doc Markdown.doc"""
    continued_fraction_hirzebruch_jung_to_rational(v::Vector{Int64})

Return the rational number corresponding to a Hirzebruch-Jung continued
fraction given as a vector of (positive) integers.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1,c_2,\ldots,c_n]$ is $r([c_1,c_2,\ldots,c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots,c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7,5)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> v = continued_fraction_hirzebruch_jung(cqs)
3-element Vector{Int64}:
 2
 2
 3

julia> continued_fraction_hirzebruch_jung_to_rational(v)
7//5
```
"""
function continued_fraction_hirzebruch_jung_to_rational(v::Vector{Int64})
    return Rational(Polymake.fulton.cf2rational(v))
end
export continued_fraction_hirzebruch_jung_to_rational


@doc Markdown.doc"""
    rational_to_continued_fraction_hirzebruch_jung(r::Rational)

Encode a (positive) rational number as a Hirzebruch-Jung continued fraction,
i.e. find the Hirzebruch-Jung continued fraction corresponding to the given
rational number.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1,c_2,\ldots,c_n]$ is $r([c_1,c_2,\ldots,c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots,c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> r = Rational(2464144958, 145732115)
2464144958//145732115

julia> cf = rational_to_continued_fraction_hirzebruch_jung(r)
7-element Vector{Int64}:
 17
 11
 23
 46
 18
 19
 37

julia> continued_fraction_hirzebruch_jung_to_rational(cf)
2464144958//145732115

julia> r == continued_fraction_hirzebruch_jung_to_rational(cf)
true
```
"""
function rational_to_continued_fraction_hirzebruch_jung(r::Rational)
    return Vector{Int64}(Polymake.fulton.rational2cf(r))
end
export rational_to_continued_fraction_hirzebruch_jung
