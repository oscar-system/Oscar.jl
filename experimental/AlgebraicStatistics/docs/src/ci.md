```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Conditional independence statements

Conditional independence (CI) statements over a ground set ``N`` are
triples of pairwise disjoint subsets ``I, J, K \subseteq N`` denoted
as ``[I \mathrel{⫫} J \mid K]``. The ground set indexes objects under
consideration and the CI statement asserts that once the objects in
``K`` are "controlled" (conditioned on, in statistical language), the
objects in ``I`` reveal no information about (are independent of) the
objects in ``J``.

The functionality documented here deals with CI statements are combinatorial
objects. Collections of CI statements are often used to state *Markov
properties* of graphical models in statistics and are ultimately used
to define ideals. Their interpretations as polynomial equations depend
on the ambient ring ([`markov_ring`](@ref) or [`gaussian_ring`](@ref)).

```@docs
ci_stmt(I::Vector{<:VarName}, J::Vector{<:VarName}, K::Vector{<:VarName})
@CI_str(str)
Base.:(==)(lhs::CIStmt, rhs::CIStmt)
Base.hash(stmt::CIStmt, h::UInt)
ci_statements(random_variables::Vector{<:VarName})
make_elementary(stmt::CIStmt; semigaussoid=false)
```
