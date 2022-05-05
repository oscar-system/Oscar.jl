```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["galois.md"]
```

# Galois Theory

Let `K` be a finite (separable) field extension of `k`. Then, in contrast to most of 
the literature we distinguish two concepts

 - the *automorphism group*
 - the *Galois group*
 
The automorphism group deals with the actual automorphism of `K` fixing `k`
and thus is, in general trivial. Access is via two constructions:

 - a list of all automorphisms (usually only the identity)
 - the group of automorphisms, returned as an abstract group and
   a map linking group elements to actual automorphisms

On the other hand, the Galois group is isomorphic to the automorphism
group of the normal closure and is explicitly given as a group of
permutations of the roots  of the defining polynomial. Thus even
in the case of `K` over `k` being normal, elements of the
Galois group do not immediately give automorphisms at all.

Currently, the computation of Galois groups is possible for
 
  - `K` a simple extension of the rationals (`AnticNumberField`)
  - `K` a simple extension of an `AnticNumberField` 
  - `K` a finite extension of the rational function field over the
     rationals. In this case the
     monodromy group can be computed as well, ie. the automorphism group over the
     complex numbers.
  - `f` a polynomial over the rationals, or an `AnticNumberField`

Independently of the Galois group, subfields, that is intermediate fields
between `K` and `k` can be computed as well.

## Automorphism Group

The automorphisms are computed using various specialised factoring
algorithms: lifting the roots of the defining polynomial in the
given field modulo suitable prime ideal powers and
recovering the true roots from this information.

The main information is included in the number field chaper, see for example

  - [`automorphisms(::Hecke.NumFieldMor)`](@ref)
  - [`automorphism_group(::NumField)`](@ref)
  - [`automorphism_group(::NumField, ::NumField)`](@ref)

## Subfields

The main information is included in the number field chaper, see

  - [`subfields(K::SimpleNumField; degree::Int = -1)`](@ref)
  - [`Hecke.principal_subfields(K::SimpleNumField)`](@ref)
  - [`subfields(FF::Generic.FunctionField{fmpq})`](@ref)

By setting `set_verbose_level(:Subfields, n::Int)` to 1 or 2
information about the progress can be obtained.

## Galois Group

The computation of Galois groups follows Stauduhars algorithm with many
improvements, see ... for an overview.

The entrire computation can also be thought of finding a description of the
splitting field of the polynomial. In fact, the information returned
can be used to verify any algebraic identity between the roots, and
find explicit subfields of the splitting field as well.

Information about the progress is available via
 
 - `set_verbose_level(:GaloisGroup, n::Int)`
 - `set_verbose_level(:GaloisInvariants, n::Int)`

```@docs
galois_group(K::AnticNumberField, extra::Int = 5; useSubfields::Bool = true, pStart::Int = 2*degree(K), prime::Int = 0)
galois_group(f::PolyElem{<:FieldElem})
```

Over the rational function field, we can also compute the monodromy group:
```jldoctest galqt; setup = :(using Oscar, Random ; Random.seed!(1))
julia> Qt, t = RationalFunctionField(QQ, "t");

julia> Qtx, x = Qt["x"];

julia> F, a = function_field(x^6 + 108*t^2 + 108*t + 27);

julia> subfields(F)
4-element Vector{Any}:
 (Function Field over Rational Field with defining polynomial a^2 + 108*t^2 + 108*t + 27, _a^3)
 (Function Field over Rational Field with defining polynomial a^3 - 54*t - 27, (-1//12*_a^4 + (3//2*t + 3//4)*_a)//(t + 1//2))
 (Function Field over Rational Field with defining polynomial a^3 + 54*t + 27, (1//12*_a^4 + (3//2*t + 3//4)*_a)//(t + 1//2))
 (Function Field over Rational Field with defining polynomial a^3 - 108*t^2 - 108*t - 27, -_a^2)

julia> galois_group(F)
(Group([ (), (1,5)(2,3)(4,6), (1,3,4)(2,5,6) ]), Galois Context for s^6 + 108*t^2 + 540*t + 675)

julia> G, C, k = galois_group(F, overC = true)
(Group([ (1,3,4)(2,5,6) ]), Galois Context for s^6 + 108*t^2 + 540*t + 675, Number field over Rational Field with defining polynomial x^2 + 12*x + 24336)

```
So, while the splitting field over `Q(t)` has degree `6`, the galois group there
is isomorphic to the `S(3)` or `D(3)` (on 6 points), the splitting field
over `C(t)` is only of degree `3`. Here the group collapses to a cyclic group
of degree `3`, the algebraic closure of `Q` in the splitting field is the
quadratic field returned last. It can be seen to be isomorphic to a cyclotomic field:

```jldoctest galqt
julia> isisomorphic(k, cyclotomic_field(3)[1])
(true, Map with following data
Domain:
=======
Number field over Rational Field with defining polynomial x^2 + 12*x + 24336
Codomain:
=========
Cyclotomic field of order 3)

```

The information returned consists always at least of a group `G` and a `GaloisCtx`: `C`.
Jointly, they can be used to further work with the information:

```@docs
roots(C::Oscar.GaloisGrp.GaloisCtx{Hecke.qAdicRootCtx}, pr::Int)
Oscar.GaloisGrp.upper_bound
Oscar.GaloisGrp.isinteger
Oscar.GaloisGrp.resolvent(C::Oscar.GaloisGrp.GaloisCtx, G::PermGroup, U::PermGroup, extra::Int = 5)
```

To illustrate:
```jldoctest galois1; setup = :(using Oscar, Random ; Random.seed!(1))
julia> Qx, x = QQ["x"];

julia> f = (x^2-2)*(x^2-3);

julia> G, C = galois_group(f)
(Group([ (1,2), (3,4) ]), Galois Context for x^4 - 5*x^2 + 6 and prime 11)

julia> r = roots(C, 5)
4-element Vector{qadic}:
 5*11^0 + 2*11^1 + 6*11^2 + 8*11^3 + 11^4 + O(11^5)
 6*11^0 + 8*11^1 + 4*11^2 + 2*11^3 + 9*11^4 + O(11^5)
 (10*11^0 + 4*11^1 + 4*11^2 + 10*11^3 + 8*11^4 + O(11^5))*a + 2*11^0 + 6*11^1 + 4*11^2 + 3*11^3 + 9*11^4 + O(11^5)
 (11^0 + 6*11^1 + 6*11^2 + 2*11^4 + O(11^5))*a + 9*11^0 + 4*11^1 + 6*11^2 + 7*11^3 + 11^4 + O(11^5)

julia> r[1]^2
3*11^0 + O(11^5)

julia> r[3]^2
2*11^0 + O(11^5)
```
To illustrate the use as a splitting field, we will prove that `r[1]^2` is actually
an integer - and that `r[1]+r[3]` is not.

Any multivariate polynomial in four variables and with integer coefficients defines via evaluation at the
roots an element in the splitting field. In case the evaluation is 
actually an integer, this can be proven with the tools provided.

```jldoctest galois1
julia> I, s = PolynomialRing(ZZ, 4);

julia> s[1]^2
x1^2

```
Next, we need a bound for the evaluation as a complex number, and compute the precision
neccessary:

```jldoctest galois1
julia> B = Oscar.GaloisGrp.upper_bound(C, s[1]^2)
(x <= 36)

julia> pr = Oscar.GaloisGrp.bound_to_precision(C, B)
7

```
Finally, we evalute the polynomial at the roots and verify that the exact value is `3`:
```jldoctest galois1
julia> evaluate(s[1]^2, roots(C, 7))
3*11^0 + O(11^7)

julia> Oscar.GaloisGrp.isinteger(C, B, ans)
(true, 3)
```
Now, to show that `r[1] + r[3]` is not an integer:
```jldoctest galois1
julia> B = Oscar.GaloisGrp.upper_bound(C, s[1] + s[3])
(x <= 12)

julia> Oscar.GaloisGrp.isinteger(C, B, evaluate(s[1] + s[3], roots(C, 7)))
(false, 0)
```
More interestingly, we can use this to find the minimal polynomial of `r[1] + r[3]`.
Generically, the Galois-conjugates of `r[1]+r[3]` should be the `G`-orbit
of `s[1]+s[3]` evaluated at the roots.

Once the orbit is known, the coefficients of the minimal polynomial are just the elementary
symmetric functions evaluated at the roots: 

```jldoctest galois1
julia> o = collect(orbit(G, s[1]+s[3]))
4-element Vector{fmpz_mpoly}:
 x1 + x3
 x1 + x4
 x2 + x4
 x2 + x3

julia> for i=1:4
         B = Oscar.GaloisGrp.upper_bound(C, elementary_symmetric, o, i)
         pr = Oscar.GaloisGrp.bound_to_precision(C, B)
         co = [evaluate(x, roots(C, pr)) for x = o]
         println(i, ": ", Oscar.GaloisGrp.isinteger(C, B, elementary_symmetric(co, i)))
       end
1: (true, 0)
2: (true, -10)
3: (true, 0)
4: (true, 1)
```

So, `x^4-10x^2+1` should be the minimal polynomial to $\sqrt 3 + \sqrt 2$ - which it is.

In the case of computations over the rational function field, both the
precision and the bound are more complicated - but can be used in the same way:
Here, the roots are power series with `q`-adic coefficients, thus the precision
has to cover both the precision of the coefficient as well as the number
of terms in the series. Similarly, in this context, an `isinteger`
is now a polynomial with integer coefficients. Thus the bound needs
to bound the degree as well as the coefficient size.

```jldoctest
julia> Qt,t = RationalFunctionField(QQ, "t");

julia> Qtx, x = Qt["x"];

julia> F, a = function_field(x^3+t+2);

julia> G, C = galois_group(F);

julia> describe(G)
"S3"

julia> _, s = slpoly_ring(ZZ, 3);

julia> B = Oscar.GaloisGrp.upper_bound(C, prod(s))
(x <= (9261, 2, 1))

julia> pr = Oscar.GaloisGrp.bound_to_precision(C, B)
(2, 2)

julia> Oscar.GaloisGrp.isinteger(C, B, evaluate(prod(s), roots(C, pr)))
(true, -t - 2)
```

```@docs
galois_quotient(C::Oscar.GaloisGrp.GaloisCtx, Q::PermGroup)
galois_quotient(C::Oscar.GaloisGrp.GaloisCtx, d::Int)
galois_quotient(C::Oscar.GaloisGrp.GaloisCtx, d::Int, n::Int)
galois_quotient(f::PolyElem, p::Vector{Int})
fixed_field(GC::Oscar.GaloisGrp.GaloisCtx, U::PermGroup, extra::Int = 5)
minpoly(C::Oscar.GaloisGrp.GaloisCtx, I, extra::Int = 5)
```

```@docs
Oscar.GaloisGrp.cauchy_ideal(f::PolyElem{<:FieldElem})
Oscar.GaloisGrp.galois_ideal(C::Oscar.GaloisGrp.GaloisCtx, extra::Int = 5)
```
