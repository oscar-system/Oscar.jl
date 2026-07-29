# Introduction

This project aims to provide functionality for what we call action polynomial rings. This phrase refers to a generalized
framework that allows for an algorithmic treatment of both difference polynomial rings and differential polynomial rings.
In the future, further similar algebraic structure might be covered as well.

## [Definitions and notations](@id definitionsandnotations_apr)

For a positive integer ``k`` we denote the set ``\{1, \ldots, k\}`` by ``\underline{k}``.
An *action polynomial ring* over the commutative ring ``R`` is a polynomial ring
```math
S = R[\, (u_i)_J \mid i \in \underline{m}, \, J \in \mathbb{Z}_{\geq 0}^n ]
```
in the countably infinitely many *jet variables* ``(u_i)_J``, equipped with ``n`` commuting
``R``-linear *action maps* ``\sigma_1, \ldots, \sigma_n``, where ``m`` and ``n`` are positive
integers. The symbols ``u_1, \ldots, u_m`` are called *action indeterminates*, the multiindices
``J \in \mathbb{Z}_{\geq 0}^n`` are called *jets*. Note that depending on the specific setting, the 
action maps and action indeterminates get a more specific name, e.g. if ``S`` is a difference polynomial ring then
the ``\sigma_j`` are called *shift operators* and the ``u_i`` are called *difference indeterminates*. If instead,
``S`` is a differential polynomial ring, then the ``\sigma_j`` are called *derivatives* and the ``u_i``
are called *differential indeterminates*.

We put ``\Delta \coloneqq \{ \sigma_1, \ldots, \sigma_n\}`` and denote the set of monomials in ``\Delta`` by
```math
\operatorname{Mon}(\Delta) = \{\sigma_1^{a_1} \ldots \sigma_n^{a_n} \mid a_j \in \mathbb{N}_0, \, j \in \underline{n}\}.
```
The ``j``-th action map ``\sigma_j`` has the property that, when applied to a jet variable, it increments the 
``j``-th entry of its jet by one. Depending on the specific setting it also has further properties,
e.g. it is multiplicative for difference polynomial rings and it is a derivation for differential
polynomial rings.\
Given an action polynomial ``p \in S`` we call an expression of the form ``\Theta(p)``, with some
``\Theta \in \operatorname{Mon}(\Delta)`` an *action transform* of ``p``. It is called proper, if
``\Theta \neq 1``. Note that an action transform of a jet variable is itself a jet variable. Specifically,
if ``\Theta = \sigma_1^{a_1} \ldots \sigma_n^{a_n}``, with ``a = (a_1, \ldots, a_n) \in \mathbb{N}_0^n``, then
```math
\Theta\left( (u_i)_J \right) = (u_i)_{J + a}. 
```

### [Rankings](@id actionpolyranking)

Most algorithms require the jet variables of the action polynomial ring ``S`` to be totally ordered. This is achieved using
so-called *rankings*. A total ordering ``>`` on the set of jet variables of ``S`` is called a ranking, if the following two conditions are satisfied:
- ``\Theta(v) > v`` for each ``1 \neq \Theta \in \operatorname{Mon}(\Delta)`` and each jet variable ``v``
- If ``v`` and ``w`` are jet variables in ``S`` with ``v > w`` then ``\Theta(v) > \Theta(w)`` for all ``\Theta \in \operatorname{Mon}(\Delta)``

Analogous to how monomial orderings in a multivariate polynomial ring in ``n`` variables correspond to total
orderings of ``\mathbb{N}_0^n`` that respect multiplication by monomials, rankings of jet variables correspond to total orderings
of ``X := \underline{m} \times \mathbb{N}_0^n`` that respect the action maps. Specifically, one can identify the jet variable
``(u_i)_J \in S`` with the pair ``(i, J) \in X``. Consequently, by Dickson's Lemma, each properly descending chain of jet
variables that are ordered with respect to some ranking is necessarily finite.

---

Let ``p \in S \setminus R`` be an action polynomial equipped with a ranking ``>``. The *leader* of ``p``, denoted by ``\operatorname{ld}(p)``,
is the ``>``-largest (or highest ranked) jet variable occurring in ``p``. If we regard ``p`` as a univariate polynomial in
``\operatorname{ld}(p)`` then its coefficients are itself action polynomials in lower ranked jet variables, which we once again can regard as
univariate polynomials in their respective leaders and so on. The *initial* of ``p``, denoted by ``\operatorname{init}(p)``, is the
leading coefficient of ``p`` regarded as a univariate polynomial in ``\operatorname{ld}(p)``. The initial of ``p`` is very important
in the context of [polynomial reduction](@ref polynomial_reduction_apr).

If ``c \in R \setminus \{0\}`` is a nonzero constant polynomial, we define its leader as the multiplicative identity of ``S``. The initial
of ``c`` is just ``c`` itself. The zero polynomial has no leader and consequently also no initial.

### [Riquier rankings](@id actionpolyriquierranking)

The rankings we use are called *Riquier rankings*. By definition, these are the rankings of ``X`` that extend to a ranking of
``\{1\} \times \mathbb{N}^{m+n}``. Equivalently, we have that ``(u_i)_I > (u_i)_J`` implies ``(u_j)_I > (u_j)_J`` for all
``i,j \in \underline{m}`` and all ``I,J \in \mathbb{N}_0^{n}``. Another equivalent condition and the one used for implementation
is the following:

There exists a positive integer ``s`` and an ``s \times (m+n)`` real matrix ``M`` such that the total ordering
of the jet variables defined by ``X`` coincides with the ordering obtained from the [matrix ordering](@ref "Matrix Orderings")
on ``\mathbb{N}_0^{m+n}`` defined by ``M``.

For this construction, we identify a jet variable ``(u_i)_J \in A`` with ``(e_i, J) \in \mathbb{N}_0^{m+n}``,
where ``e_i`` is the ``i``-th unit row and restrict ourselves to integer matrices ``M``. In this context, we call ``M``
a *Riquier matrix*.

!!! note "Integer Riquier matrix"
    Not all Riquier rankings are obtained from integral Riquier matrices. However, this is the case if we only require a total
    ordering of a finite subset of ``X``. Thus, only considering integer matrices is sufficient for practical use.

#### Ritt ordering

While it is possible to define a total ordering of the set of jet variables of the action polynomial ring ``S``, there is no obvious way
to do the same for the set of all action polynomials. However, one can still compare two action polynomials using the so-called
*Ritt ordering* (associated to a given ranking ``<`` on ``S``). Given ``p, q \in S``, we say that *``p`` is smaller than  ``q``* with
respect to Ritt ordering, if one of the following conditions holds:
- ``p \in R`` but ``q \notin R``
- ``p, q \notin R`` and ``\operatorname{ld}(p) < \operatorname{ld}(q)``
- ``p, q \notin R``, ``v \coloneqq \operatorname{ld}(p) = \operatorname{ld}(q)`` and ``\operatorname{deg}_v(p) < \operatorname{deg}_v(q)``

Additionally, we make the convention that the zero polynomial is smaller with respect to Ritt ordering than any other non-zero polynomial.

!!! note "Incomparable action polynomials"
    In any action polynomial ring there exist polynomials ``p,q \in S`` such that ``p \neq q`` but neither ``p`` is less than ``q`` with
    respect to Ritt ordering nor the other way around. In this case, ``p`` and ``q`` are said to be *incomparable* (with respect to Ritt
    ordering). Clearly, non-constant ``p`` and ``q`` are incomparable if and only if they have the same leader as well as
    the same degree in that leader.

### [Polynomial reduction](@id polynomial_reduction_apr)

Let ``p, q \in S = R[\, (u_i)_J \mid i \in \lbrace 1, \ldots, m \rbrace, J \in \mathbb{Z}_{\geq 0}^n]`` with ``q \notin R`` be
action polynomials and let ``v`` denote the leader of ``q``, so that both ``p`` and ``q`` can be regarded as univariate polynomials
in the jet variable ``v``. For the implementations of the concepts discussed below, see [polynomial reduction methods](@ref polynomial_reduction_methods_apr).

#### [Pseudo-division for action polynomials](@id pseudodivision_for_action_polynomials)

Classical long-division of ``p`` by ``q`` will usually fail since ``\operatorname{init}(q)`` need not be invertible in ``S``. This can be fixed by
pre-multiplying the dividend ``p`` by the initial of ``q`` to some high enough power. In fact, it is easy to see that one can always
choose ``a = \max(0, \operatorname{deg}_{v}(p) - \operatorname{deg}_{v}(q) + 1)`` as this power. However, this choice is usually not minimal.

The process described above is called *pseudo-division* of ``p`` by ``q``. The *pseudo-remainder* of ``p`` by ``q`` is an
action polynomial ``r`` with ``\operatorname{deg}_{v}(r) < \operatorname{deg}_{v}(q)`` or ``r = 0`` such that there exists an
integer ``a \in \mathbb{N}_0`` and some ``s \in S``, the so-called *pseudo-quotient* of ``p`` by ``q`` such that
```math
\operatorname{init}(q)^a p = s \cdot q + r.
```
Note that for each fixed exponent ``a``, both ``s`` and ``r`` are unique (if they exist). In our implementations and in further discussions,
if we refer to the pseudo-remainder or pseudo-quotient of ``p`` by ``q``, we mean the values of ``s`` and ``r`` that correspond to the above
identity where the exponent ``a`` is minimal.

Slightly generalising the above, we also allow the divisor ``q`` to be a nonzero constant polynomial ``c`` by regarding ``1`` as the leader
of ``c`` and by putting ``\operatorname{init}(c) = c``. Thus, using the same formula as above, the pseudo-remainder of ``p`` by ``c`` is
always equal to zero. The pseudo-quotient is ``p/c`` (with ``a = 0``) if ``p`` is divisible by ``c``; otherwise it is ``p`` (with ``a = 1``).

#### [Notions of reducedness for action polynomials](@id notions_of_reducedness_for_action_polynomials)

Let ``p, q \in S``, ``q \notin R`` and ``v = \operatorname{ld}(q)``. We have the following notions of reducedness for action polynomials:
- ``p`` is *partially reduced* with respect to ``q``, if ``\operatorname{deg}_{w}(p) < \operatorname{deg}_{w}(\Theta(q))`` 
    for all proper action transforms ``w = \Theta(v)``, ``1 \neq \Theta \in \operatorname{Mon}(\Delta)``
- ``p`` is *reduced* with respect to ``q``, if ``\operatorname{deg}_{w}(p) < \operatorname{deg}_{w}(\Theta(q))``
    for all action transforms ``w = \Theta(v)``, ``\Theta \in \operatorname{Mon}(\Delta)``

!!! note "partially reduced differential polynomials"
    In the case where ``S`` is a differential polynomial, the degree of ``\Theta(q)`` in each proper transform ``w`` of
    ``q`` is one. Thus, being partially reduced amounts to the condition that ``p`` contains no proper derivative
    of ``v`` as a variable.

Also allowing ``q`` to be a nonzero constant polynomial ``c``, the only polynomial that is partially reduced or reduced with respect to ``c``
is the zero polynomial.

---

There also exist the following set based notions of reducedness for action polynomials. For this, let ``p \in S`` and ``M \subset S \setminus \{0\}``.
- ``p`` is called *partially reduced* with respect to ``M``, if ``p`` is partially reduced with respect to all elements of ``M``
- ``p`` is called *reduced* with respect to ``M``, if ``p`` is reduced with respect to all elements of ``M``
- The set ``M`` is *autoreduced*, if for any two distinct ``p_1, p_2 \in M``, ``p_1`` is reduced with respect to ``p_2``

Note that for our implementation we use vectors of action polynomials instead of sets of action polynomials. In view of this, we add
the requirement that an autoreduced set must be ordered with respect to Ritt ordering.

## Content

This project is structured as follows:

- The section [Action polynomial rings](@ref actionpolyring) contains generic methods that work for all subtypes of the action polynomial ring interface.
- The sections [Difference polynomial rings](@ref differencepolyring) and [Differential polynomial rings](@ref differentialpolyring) contain functionality that is unique to difference polynomial rings and differential polynomial rings respectively.
- The section [Rankings](@ref actionpolyranking) contains the functionality for rankings of jet variables of action polynomial rings.

## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what this means.

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Tobias Braun](https://www.math.rwth-aachen.de/homes/Tobias.Braun/)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

