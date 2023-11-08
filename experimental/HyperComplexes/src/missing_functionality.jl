# The commented lines are also contained in the DoubleComplexes 
# and hence removed here to suppress warnings.
#morphism_type(::Type{T}) where {T<:ModuleFP} = ModuleFPHom{T, T}
morphism_type(::Type{T}) where {T<:FreeMod} = FreeModuleHom{T, T}

#typ(C::ComplexOfMorphisms) = C.typ
#tensor_product(v::Vector{<:ModuleFP}) = tensor_product(v...)
