```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Action polynomial rings](@id actionpolyring)

In Oscar we provide the action polynomial interface via the abstract types `ActionPolyRing{T} <: Ring`
and `ActionPolyRingElem{T} <: RingElem`. The type parameter `T` is the element type of the coefficient
ring. All concrete subtypes use the functionality of [universal polynomials](@ref "Universal polynomial ring")
from the AbstractAlgebra package for polynomial arithmetic, as well as maintaining variables and adding
new ones on demand. Currently, there are two concrete subtypes available, namely `DifferencePolyRing{T}` and
`DifferentialPolyRing{T}` with element types `DifferencePolyRingElem{T}` and `DifferentialPolyRingElem{T}`.
See [difference polynomial rings](@ref differencepolyring)
and [differential polynomial rings](@ref differentialpolyring) for their unique functionality.

Each action polynomial ring maintains a sorted list of currently tracked jet variables,
that can be accessed and extended by a number of methods, see, e.g. the section [Element Constructors](@ref element_constructors_apr).
The jet variables are sorted with respect to a user-defined [ranking](@ref actionpolyranking).

!!! note "Tracked jet variables"
    The set of valid jet variables of an action polynomial ring depend only on the integers ``m`` and
    ``n`` and are thus known at the time of construction. For reasons of efficiency, we keep the list of
    tracked jet variables as short as possible and track jet variables only, if necessary. The list of
    currently tracked jet variables is obtained using
    [`gens`](@ref gens(apr::ActionPolyRing)).

Additionally, each action polynomial ring maintains a vector of pairwise commuting action maps associated with it. 
These action maps are implemented via the new abstract type `ActionMap{D} <: Map{D,D,Any,Any}`, with type parameter
`D <: Ring`.

Depending on the type of the action polynomial ring, `ActionMap` branches into two core abstract subtypes:
- `ActionShift`: Represents shift operators used in difference polynomial rings.
- `ActionDerivation`: Represents derivations used in differential polynomial rings.

---

For detailed information on constructing specific action polynomial rings, see:
- [Difference polynomial rings](@ref differencepolyring_construction)
- [Differential polynomial rings](@ref differentialpolyring_construction)

## [Specifying jet variables](@id specifying_jet_variables)

Recall that a jet variable is of the form ``(u_i)_J`` with ``i \in \lbrace 1, \ldots, m \rbrace`` and
``J \in \mathbb{Z}_{\geq 0}^n``. There are four ways in which jet variables can be specified as an input for methods,
which can be found below. The first two do not require the jet variable in question to be tracked, the last two do.
- By a tuple consisting of the index `i` and the jet `J`.
- By passing the tuple as separate arguments, starting with the index.
- By passing the index of the jet variable in the list of the currently tracked jet variables.
- By immediately passing the jet variable as an element of an action polynomial ring.

!!! note 
    For many methods, e.g. [`degree`](@ref degree(p::ActionPolyRingElem, i::Int, jet::Vector{Int})) or
    [`derivative`](@ref derivative(p::ActionPolyRingElem, i::Int, jet::Vector{Int})) we provide all the above
    versions, but only record one in this documentation for readability. Usually, we choose the second version
    from the above list.

## [Element Constructors](@id element_constructors_apr)

`(A::ActionPolyRing)()` returns the zero polynomial of the action polynomial ring `A`.
`(A::ActionPolyRing)(a::T) where {T<:RingElement}` returns `a` as an element of `A`, if possible.
This can be used for creating constant polynomials.

!!! warning
    The next three methods take input arguments specifying a jet variable or a vector of jet variables. After calling one of them,
    the provided action polynomial ring `A` will track all jet variables specified.

```@docs
gen(A::ActionPolyRing, i::Int, jet::Vector{Int})
getindex(A::ActionPolyRing, i::Int, jet::Vector{Int})
gens(A::ActionPolyRing, jet_idxs::Vector{Tuple{Int, Vector{Int}}})
```

!!! note "Creating polynomials"
    Polynomials can be created by applying the usual arithmetic operations, such as `+`, `-`, `*`, and `^`,
    to jet variables.

## [Generators and variables](@id generators_and_variables_apr)

```@docs
gen(A::ActionPolyRing, i::Int)
gens(A::ActionPolyRing)
is_gen(p::ActionPolyRingElem)
var_index(p::ActionPolyRingElem)
vars(p::ActionPolyRingElem)
leader(p::ActionPolyRingElem)
```

We also provide the usual `ngens` and `nvars` methods that respectively return the number of currently tracked jet variables. Finally,
we provide the `getindex` method below to allow for fast access to jet variables from existing ones.

```@docs
getindex(var::ActionPolyRingElem, index_shift::Int...)
```

## [Basic methods for action polynomial rings](@id basic_functionality_apr)

```@docs
zero(A::ActionPolyRing)
one(A::ActionPolyRing)
n_action_indeterminates(A::ActionPolyRing)
action_indeterminates(A::DifferencePolyRing)
action_maps(A::DifferencePolyRing)
action_map(A::DifferencePolyRing, i::Int)
n_action_maps(A::DifferencePolyRing)
```

## [Iterators](@id iterators_apr)

The following iterators are available for elements of action polynomial rings. The entries across the different iterators are
guaranteed to match. Moreover, the order of the entries of the iterators depends only on the ranking of the action polynomial
ring, leading with the most significant entry.

```@docs
coefficients(p::ActionPolyRingElem)
exponents(p::ActionPolyRingElem)
monomials(p::ActionPolyRingElem)
terms(p::ActionPolyRingElem)
```

### Iterator based methods

The following methods are based on the [iterators](@ref iterators_apr) for elements of action polynomial rings.

---

Basic access to entries of the iterators:

```@docs
coeff(p::ActionPolyRingElem, i::Int)
exponent_vector(p::ActionPolyRingElem, i::Int)
monomial(p::ActionPolyRingElem, i::Int)
term(p::ActionPolyRingElem, i::Int)
```

---

Access to the first and last entries:

```@docs
leading_coefficient(p::ActionPolyRingElem)
leading_monomial(p::ActionPolyRingElem)
leading_term(p::ActionPolyRingElem)

trailing_coefficient(p::ActionPolyRingElem)
trailing_monomial(p::ActionPolyRingElem)
trailing_term(p::ActionPolyRingElem)
```

---

Other useful methods:

```@docs
is_monomial(p::ActionPolyRingElem)
is_term(p::ActionPolyRingElem)

initial(p::ActionPolyRingElem)
length(p::ActionPolyRingElem)
tail(p::ActionPolyRingElem)
```

## [Miscellaneous](@id miscellaneous_apr)
In this subsection, we enumerate methods that might be useful but primarily exists, because they already do
for other polynomial types.

### Constant polynomials

```@docs
is_constant(p::ActionPolyRingElem)
constant_coefficient(p::ActionPolyRingElem)
```

### [Degree](@id degree_apr)

```@docs
degree(p::ActionPolyRingElem, i::Int, jet::Vector{Int})
degrees(p::ActionPolyRingElem)
total_degree(p::ActionPolyRingElem)
```

### [Derivative](@id derivative_apr)

```@docs
derivative(p::ActionPolyRingElem, i::Int, jet::Vector{Int})
```

### [Discriminant and resultant](@id discriminant_resultant_apr)

```@docs
discriminant(p::ActionPolyRingElem)
resultant(f::ActionPolyRingElem, g::ActionPolyRingElem, i::Int, jet::Vector{Int})
```

### [Evaluation](@id evaluation_apr)

The following function allows evaluation of a polynomial at all its variables. The result is always
in the ring that a product of a coefficient and one of the values belongs to, i.e. if all the values
are in the coefficient ring, the result of the evaluation will be too.

```@docs
evaluate(a::ActionPolyRingElem{T}, vals::Vector{V}) where {T <: RingElement, V <: RingElement}
```

The following functions allow evaluation of a polynomial at some of its variables. Note that the
result will be a product of values and an element of the polynomial ring, i.e. even if all the
values are in the coefficient ring and all variables are given values, the result will be a
constant polynomial, not a coefficient.

```@docs
evaluate(a::ActionPolyRingElem{T}, vars::Vector{Int}, vals::Vector{V}) where {T <: RingElement, V <: RingElement}
evaluate(a::PolyT, vars::Vector{PolyT}, vals::Vector{V}) where {PolyT <: ActionPolyRingElem, V <: RingElement}
```

### Univariate polynomials

```@docs
is_univariate(p::ActionPolyRingElem)
to_univariate(R::PolyRing{T}, p::ActionPolyRingElem{T}) where {T <: RingElement}
to_univariate(p::ActionPolyRingElem)
univariate_coefficients(p::ActionPolyRingElem, i::Int, jet::Vector{Int})
univariate_leading_coefficient(p::ActionPolyRingElem, i::Int, jet::Vector{Int})
```

## [Polynomial reduction methods](@id polynomial_reduction_methods_apr)

The following two methods `pseudorem` and `pseudodivrem` for the pseudo-division of an action polynomial ``p`` by another action polynomial
``q`` form the backbone of most reduction methods. Recall that if ``s`` is the pseudo-quotient and the pseudo-remainder ``r`` of ``p`` by ``q``,
we have the identity
```math
\operatorname{init}(q)^a p = s \cdot q + r,
```
with ``\operatorname{deg}_{v}(r) < \operatorname{deg}_{v}(q)`` or ``r = 0``, ``v = \operatorname{ld}(q)`` and some ``a \in \mathbb{N}_0``.
In order to avoid coefficient swell, both methods specifically return the values for ``r`` and ``s``, where the exponent ``a`` is minimal. For
additional flexibility, both methods also allow the user to specify with respect to which jet variable the pseudo-division should be performed, so they
are not just restricted to pseudo-division by the leader of the second input. However, if no such jet variable is specified, the leader is used by default.

```@docs
pseudorem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem}
pseudodivrem(p::PolyT, q::PolyT, i::Int, jet::Vector{Int}) where {PolyT <: ActionPolyRingElem}
```

---

We provide the following methods for reducing the action polynomial ``p`` with respect to the action polynomial ``q`` and to verify reducedness:

```@docs
is_partially_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
is_reduced(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}

partially_reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
reduce(p::PolyT, q::PolyT) where {PolyT <: ActionPolyRingElem}
```

We also provide similar methods for the set-related notions of reducedness:

```@docs
is_partially_reduced(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
is_reduced(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
is_autoreduced(S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}

partially_reduce(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
reduce(p::PolyT, S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
autoreduce(S::Vector{PolyT}) where {PolyT <: ActionPolyRingElem}
```
