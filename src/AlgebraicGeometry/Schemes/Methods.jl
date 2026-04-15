@doc raw"""
    base_change(phi::Any, X::Scheme)

For a `Scheme` ``X`` over a `base_ring` ``𝕜`` and a map ``φ : 𝕜 → R`` 
we compute ``X' = X ×ₖ Spec(R)`` and return a pair `(X', f)` where 
``f : X' → X`` is the canonical morphism.

!!! note
    We do not restrict `phi` to be of type `Map` so that one can also use coercion, anonymous functions, etc. 
"""
function base_change(phi::Any, X::Scheme)
  error("base_change not implemented for input of type $(typeof(X))")
end

@doc raw"""
    base_change(phi::Any, f::SchemeMor;
        domain_map::SchemeMor, codomain_map::SchemeMor
      )

For a morphism ``f : X → Y`` with both ``X`` and ``Y`` defined over a 
`base_ring` ``𝕜`` and a map ``φ : 𝕜 → R`` return a triple `(a, F, b)` 
where ``a : X' → X`` is the morphism from `base_change(phi, X)`, 
``b : Y' → Y`` the one for ``Y``, and ``F : X' → Y'`` the induced 
morphism on those fiber products.

!!! note
    We do not restrict `phi` to be of type `Map` so that one can also use coercion, anonymous functions, etc. 

!!! note 
    The morphisms ``a`` and ``b`` can be passed as the optional arguments `domain_map` and `codomain_map`, respectively. 
"""
function base_change(phi::Any, f::SchemeMor;
    domain_map::SchemeMor, codomain_map::SchemeMor
  )
  error("base_change not implemented for input of type $(typeof(f))")
end

@doc raw"""
    irreducible_components(X::Scheme)
    
Return the irreducible components of ``X`` with their reduced structure. 

Note that the irreducible components are defined over the base ring of ``X``.
An irreducible component may split into several irreducible components after a base change, e.g. a field extension. 
"""
irreducible_components(X::Scheme) = error("Not implemented")
