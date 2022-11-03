################################################################################
################################################################################
## Cyclic Quotient Singularit struct
################################################################################
################################################################################

@attributes mutable struct CyclicQuotientSingularity <: AbstractNormalToricVariety
    polymakeNTV::Polymake.BigObject
end
export CyclicQuotientSingularity

function Base.show(io::IO, cqs::CyclicQuotientSingularity)
    n = pm_object(cqs).N
    q = pm_object(cqs).Q
    print(io, "The cyclic quotient singularity Y($(n), $(q))")
end



################################################################################
################################################################################
## Constructor
################################################################################
################################################################################

@doc Markdown.doc"""
    CyclicQuotientSingularity(n::fmpz, q::fmpz)

Return the cyclic quotient singularity for the parameters $n$ and $q$, with
$0<q<n$ and $q, n$ coprime.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7, 5)
The cyclic quotient singularity Y(7, 5)

julia> is_affine(cqs)
true

julia> is_smooth(cqs)
false
```
"""
function CyclicQuotientSingularity(n::fmpz, q::fmpz)
    n > 0 || error("n (=$(n)) must be positive")
    q > 0 || error("q (=$(q)) must be positive")
    q < n || error("q must be smaller than n (q=$(q) >= n=$(n))")
    gcd(n, q) == 1 || error("n and q must be coprime (gcd=$(gcd(n, q)))")
    pmntv = Polymake.fulton.CyclicQuotient(N=convert(Polymake.Integer, n), Q=convert(Polymake.Integer, q))
    return CyclicQuotientSingularity(pmntv, Dict())
end
CyclicQuotientSingularity(n::Int64, q::Int64) = CyclicQuotientSingularity(fmpz(n), fmpz(q))


@doc Markdown.doc"""
    continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)

Return the Hirzebruch-Jung continued fraction associated with the cyclic
quotient singularity, i.e. the Hirzebruch-Jung continued fraction corresponding
to $n/q$.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1, c_2,\ldots, c_n]$ is $r([c_1, c_2,\ldots, c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots, c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7, 5)
The cyclic quotient singularity Y(7, 5)

julia> cf = continued_fraction_hirzebruch_jung(cqs)
3-element Vector{fmpz}:
 2
 2
 3

julia> ecf = cf[1]-1//(cf[2]-fmpq(1, cf[3]))
7//5
```
"""
@attr Vector{fmpz} function continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
    return Vector{fmpz}(pm_object(cqs).CONTINUED_FRACTION)
end
export continued_fraction_hirzebruch_jung


@doc Markdown.doc"""
    dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)

Return the dual Hirzebruch-Jung continued fraction associated with the cyclic
quotient singularity, i.e. the Hirzebruch-Jung continued fraction corresponding
to $q/(n-q)$.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1, c_2,\ldots, c_n]$ is $r([c_1, c_2,\ldots, c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots, c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7, 5)
The cyclic quotient singularity Y(7, 5)

julia> dcf = dual_continued_fraction_hirzebruch_jung(cqs)
2-element Vector{fmpz}:
 4
 2

julia> edcf = dcf[1] - fmpq(1, dcf[2])
7//2
```
"""
@attr Vector{fmpz} function dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
    return Vector{fmpz}(pm_object(cqs).DUAL_CONTINUED_FRACTION)
end
export dual_continued_fraction_hirzebruch_jung


@doc Markdown.doc"""
    continued_fraction_hirzebruch_jung_to_rational(v::Vector{fmpz})

Return the rational number corresponding to a Hirzebruch-Jung continued
fraction given as a vector of (positive) integers.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1, c_2,\ldots, c_n]$ is $r([c_1, c_2,\ldots, c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots, c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> cqs = CyclicQuotientSingularity(7, 5)
The cyclic quotient singularity Y(7, 5)

julia> v = continued_fraction_hirzebruch_jung(cqs)
3-element Vector{fmpz}:
 2
 2
 3

julia> continued_fraction_hirzebruch_jung_to_rational(v)
7//5
```
"""
function continued_fraction_hirzebruch_jung_to_rational(v::Vector{fmpz})
    return convert(fmpq, Polymake.fulton.cf2rational(convert(Vector{Polymake.Integer}, v)))
end
export continued_fraction_hirzebruch_jung_to_rational


@doc Markdown.doc"""
    rational_to_continued_fraction_hirzebruch_jung(r::fmpq)

Encode a (positive) rational number as a Hirzebruch-Jung continued fraction,
i.e. find the Hirzebruch-Jung continued fraction corresponding to the given
rational number.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1, c_2,\ldots, c_n]$ is $r([c_1, c_2,\ldots, c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots, c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> r = fmpq(2464144958, 145732115)
2464144958//145732115

julia> cf = rational_to_continued_fraction_hirzebruch_jung(r)
7-element Vector{fmpz}:
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
function rational_to_continued_fraction_hirzebruch_jung(r::fmpq)
    cf = continued_fraction(r)
    z = fmpz[]
    n = length(cf)
    for i in 1:n
       cfi = cf[i]
       if iseven(i)
          cfi < 2^30 || @warn "blowing up your memory"
          while (cfi -= 1) > 0
             push!(z, fmpz(2))
          end
       else
          push!(z, cfi + (1 < i < n) + (1 < n))
       end
    end
    z
end
export rational_to_continued_fraction_hirzebruch_jung
