```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# [Quotients](@id quotient)

Quotient groups in OSCAR can be defined using the instruction `quo` in two ways.

* Quotients by normal subgroups.
```@docs
quo(G::GAPGroup, N::GAPGroup)
```

* Quotients by elements.
```@docs
quo(G::T, elements::Vector{S}) where T <: GAPGroup where S <: GAPGroupElem
```
This is the typical way to build finitely presented groups.

  **Example:**
```jldoctest
julia> F = @free_group(2);

julia> G,_=quo(F,[f1^2,f2^3,(f1*f2)^2]);

julia> is_finite(G)
true

julia> is_isomorphic(G,symmetric_group(3))
true
```
Similarly to the subgroups, the output consists of a pair (`Q`,`p`), where `Q` is the quotient group and `p` is the projection homomorphism of `G` into `Q`.

```@docs
maximal_abelian_quotient
```
