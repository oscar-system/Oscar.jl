################################################################################
################################################################################
## Cyclic Quotient Singularit struct
################################################################################
################################################################################

function Base.show(io::IO, cqs::CyclicQuotientSingularity)
    n = pm_object(cqs).N
    q = pm_object(cqs).Q
    print(io, "Cyclic quotient singularity Y($(n), $(q))")
end



################################################################################
################################################################################
## Constructor
################################################################################
################################################################################

@doc raw"""
    cyclic_quotient_singularity(n::ZZRingElem, q::ZZRingElem)

Return the cyclic quotient singularity for the parameters $n$ and $q$, with
$0<q<n$ and $q, n$ coprime.

# Examples
```jldoctest
julia> cqs = cyclic_quotient_singularity(7, 5)
Cyclic quotient singularity Y(7, 5)

julia> is_affine(cqs)
true

julia> is_smooth(cqs)
false
```
"""
function cyclic_quotient_singularity(n::T, q::T) where {T <: IntegerUnion}
    n > 0 || error("n (=$(n)) must be positive")
    q > 0 || error("q (=$(q)) must be positive")
    q < n || error("q must be smaller than n (q=$(q) >= n=$(n))")
    gcd(n, q) == 1 || error("n and q must be coprime (gcd=$(gcd(n, q)))")
    pmntv = Polymake.fulton.CyclicQuotient(N=convert(Polymake.Integer, n), Q=convert(Polymake.Integer, q))
    return CyclicQuotientSingularity(pmntv, Dict())
end


@doc raw"""
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
julia> cqs = cyclic_quotient_singularity(7, 5)
Cyclic quotient singularity Y(7, 5)

julia> cf = continued_fraction_hirzebruch_jung(cqs)
3-element Vector{ZZRingElem}:
 2
 2
 3

julia> ecf = cf[1]-1//(cf[2]-QQFieldElem(1, cf[3]))
7//5
```
"""
@attr Vector{ZZRingElem} function continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
    return Vector{ZZRingElem}(pm_object(cqs).CONTINUED_FRACTION)
end


@doc raw"""
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
julia> cqs = cyclic_quotient_singularity(7, 5)
Cyclic quotient singularity Y(7, 5)

julia> dcf = dual_continued_fraction_hirzebruch_jung(cqs)
2-element Vector{ZZRingElem}:
 4
 2

julia> edcf = dcf[1] - QQFieldElem(1, dcf[2])
7//2
```
"""
@attr Vector{ZZRingElem} function dual_continued_fraction_hirzebruch_jung(cqs::CyclicQuotientSingularity)
    return Vector{ZZRingElem}(pm_object(cqs).DUAL_CONTINUED_FRACTION)
end


@doc raw"""
    continued_fraction_hirzebruch_jung_to_rational(v::Vector{ZZRingElem})

Return the rational number corresponding to a Hirzebruch-Jung continued
fraction given as a vector of (positive) integers.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1, c_2,\ldots, c_n]$ is $r([c_1, c_2,\ldots, c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots, c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> cqs = cyclic_quotient_singularity(7, 5)
Cyclic quotient singularity Y(7, 5)

julia> v = continued_fraction_hirzebruch_jung(cqs)
3-element Vector{ZZRingElem}:
 2
 2
 3

julia> continued_fraction_hirzebruch_jung_to_rational(v)
7//5
```
"""
function continued_fraction_hirzebruch_jung_to_rational(v::Vector{ZZRingElem})
    return convert(QQFieldElem, Polymake.fulton.cf2rational(convert(Vector{Polymake.Integer}, v)))
end


@doc raw"""
    rational_to_continued_fraction_hirzebruch_jung(r::QQFieldElem)

Encode a (positive) rational number as a Hirzebruch-Jung continued fraction,
i.e. find the Hirzebruch-Jung continued fraction corresponding to the given
rational number.

The rational number corresponding to a Hirzebruch-Jung continued fraction
$[c_1, c_2,\ldots, c_n]$ is $r([c_1, c_2,\ldots, c_n])\ =\
c_1-\frac{1}{r([c_2,\ldots, c_n])}$ where $r([c_n]) = c_n$.  Note that this is
differs in sign from what is commonly known as continued fraction.

# Examples
```jldoctest
julia> r = QQFieldElem(2464144958, 145732115)
2464144958//145732115

julia> cf = rational_to_continued_fraction_hirzebruch_jung(r)
7-element Vector{ZZRingElem}:
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
function rational_to_continued_fraction_hirzebruch_jung(r::QQFieldElem)
    cf = continued_fraction(r)
    z = ZZRingElem[]
    n = length(cf)
    for i in 1:n
       cfi = cf[i]
       if iseven(i)
          cfi < 2^30 || @warn "blowing up your memory"
          while (cfi -= 1) > 0
             push!(z, ZZRingElem(2))
          end
       else
          push!(z, cfi + (1 < i < n) + (1 < n))
       end
    end
    z
end
