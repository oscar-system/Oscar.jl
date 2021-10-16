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
    continued_fraction(cqs::CyclicQuotientSingularity)

Return the continued fraction associated with the cyclic quotient singularity,
i.e. the continued fraction whose evaluation is $n/q$.

The evaluation of a continued fraction $[c_1,c_2,/ldots,c_n]$ is
$e([c_1,c_2,/ldots])\ =\ c_1-\frac{1}{e([c_2,\ldots,c_n])}$
where $e([c_n]) = c_n$.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7,5)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> cf = continued_fraction(cqs)
3-element Vector{Int64}:
 2
 2
 3

julia> ecf = cf[1]-1/(cf[2]-Rational(1,cf[3]))
7//5
```
"""
continued_fraction(cqs::CyclicQuotientSingularity) = Vector{Int64}(pm_ntv(cqs).CONTINUED_FRACTION)
export continued_fraction


@doc Markdown.doc"""
    dual_continued_fraction(cqs::CyclicQuotientSingularity)

Return the dual continued fraction associated with the cyclic quotient
singularity, i.e. the continued fraction whose evaluation is $q/(n-q)$.

The evaluation of a continued fraction $[c_1,c_2,/ldots,c_n]$ is
$e([c_1,c_2,/ldots])\ =\ c_1-\frac{1}{e([c_2,\ldots,c_n])}$
where $e([c_n]) = c_n$.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7,5)
A normal toric variety corresponding to a polyhedral fan in ambient dimension 2

julia> dcf = dual_continued_fraction(cqs)
2-element Vector{Int64}:
 4
 2

julia> edcf = dcf[1] - Rational(1,dcf[2])
7//2
```
"""
dual_continued_fraction(cqs::CyclicQuotientSingularity) = Vector{Int64}(pm_ntv(cqs).DUAL_CONTINUED_FRACTION)
export dual_continued_fraction
