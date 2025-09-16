```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Action polynomial rings

An *action polynomial ring* over the commutative ring ``R`` is a polynomial ring
```math
S = R[\, (u_i)_J \mid i \in \lbrace 1, \ldots, m \rbrace, J \in \mathbb{Z}_{\geq 0}^n ]
```
in the countably infinitely many *jet variables* ``(u_i)_J``, equipped with ``n`` commuting
``R``-linear *action maps* ``\sigma_1, \ldots, \sigma_n``, where ``m`` and ``n`` are positive
integers. The symbols ``u_1, \ldots, u_m`` are called elementary symbols, the multiindices
``J \in \mathbb{Z}_{\geq 0}^n`` are called jets.

The ``j``-th action map has the property that, when applied to a jet variable, it increments
the ``j``-th entry of its jet by one. Depending on the setting it also has further properties,
e.g. it is multiplicative for difference polynomial rings and it is a derivation for differential
polynomial rings.

---
In Oscar we provide the action polynomial interface via the abstract types `ActionPolyRing{T} <: Ring`
and `ActionPolyRingElem{T} <: RingElem`. The type parameter `T` is the element type of the coefficient
ring. All concrete subtypes use the functionality of [universal polynomials](@ref "Universal polynomial")
from the AbstractAlgebra package for polynomial arithmetic, as well as maintaining variables and adding
new ones on demand. Any action polynomial ring maintains a sorted list of currently tracked jet variables,
that can be accessed and extended by a number of methods, see, e.g., [Element Constructors](@ref). The jet
variables are sorted with respect to a user-defined [ranking](@ref actionpolyranking).

!!! note "Tracked jet variables"
    The set of valid jet variables of an action polynomial ring depend only on the integers ``m`` and
    ``n`` and are thus known at the time of construction. For reasons of efficiency, we keep the list of
    tracked jet variables as short as possible and track jet variables only, if necessary. The list of
    currently tracked jet variables is obtained, using
    [`gens`](@ref gens(apr::ActionPolyRing)).

Currently, there are two concrete subtypes available, namely `DifferencePolyRing{T}` and
`DifferentialPolyRing{T}` with element types `DifferencePolyRingElem{T}` and `DifferentialPolyRingElem{T}`.
See [difference polynomial rings](@ref differencepolyring)
and [differential polynomial rings](@ref differentialpolyring) for their unique functionality.

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

`(A::ActionPolyRing)()` returns the zero polynomial of the action polynomial ring `A`.\
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
    Polynomials can be created using by applying the usual arithmetic operations, such as `+`, `-`, `*`, `^`
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

We also provide the usual `ngens` and `nvars` methods that respectively return the number of currently tracked jet variables.

## [Basic methods for action polynomial rings](@id basic_functionality_apr)

```@docs
zero(A::ActionPolyRing)
one(A::ActionPolyRing)
n_elementary_symbols(A::ActionPolyRing)
elementary_symbols(A::DifferencePolyRing)
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
In this subsection we enumerate methods that might be useful but primarily exists, because they already do
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

### Univariate polynomials

```@docs
is_univariate(A::ActionPolyRing)
is_univariate(p::ActionPolyRingElem)
to_univariate(R::PolyRing{T}, p::ActionPolyRingElem{T}) where {T <: RingElement}
to_univariate(p::ActionPolyRingElem)
```

