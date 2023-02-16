export affine_normal_toric_variety

export normal_toric_variety

########################################################################
# ToricSpec                                                            #
########################################################################

underlying_scheme(X::ToricSpec) = X.X
affine_normal_toric_variety(X::ToricSpec) = X.antv
#antv(X::ToricSpec) = affine_normal_toric_variety(X)
antv = affine_normal_toric_variety
cone(X::ToricSpec) = cone(antv(X))
dual_cone(X::ToricSpec) = X.dual_cone
hilbert_basis(X::ToricSpec) = X.hb

@Markdown.doc """
    torus_inclusions(X::ToricSpec)

For an affine toric variety ``X`` this returns a list `l` 
containing the inclusions ``Tʳ⁽ⁱ⁾ ↪ X`` of the different 
tori. 
"""
function torus_inclusions(X::ToricSpec)::Vector{<:AbsSpecMor}
  #TODO: Fill in 
end

@Markdown.doc """
    torus_action(X::ToricSpec)

For an affine toric variety ``X`` with a dense open torus ``T``
this returns a quintuple of morphisms `(pT, pX, incX, mult)` 
consisting of 

 * the projection ``T × X → T`` of the product with the torus ``T`` to ``T``
 * the projection ``T × X → X``
 * the inclusion ``X ↪ T × X`` taking ``x`` to ``(1, x)``
 * the group action ``T × X → X``.
"""
function torus_action(X::ToricSpec)::AbsSpecMor
  #TODO: Fill in 
end

function Base.show(io::IO, X::ToricSpec) 
  print(io, "Spec of a toric variety with cone spanned by $(rays(cone(X)))")
end

########################################################################
# ToricCoveredScheme                                                   #
########################################################################

underlying_scheme(X::ToricCoveredScheme) = X.X
normal_toric_variety(X::ToricCoveredScheme) = X.ntv
ntv = normal_toric_variety
fan(X::ToricCoveredScheme) = fan(normal_toric_variety(X))
