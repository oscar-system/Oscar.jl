
# A union type for all algebraic extension of the field of rational numbers
QQAlgField = Union{NumField, QQField, QQBarField, QQAbField}
QQAlgFieldElem = Union{NumFieldElem, QQFieldElem, QQBarFieldElem, QQAbElem}

include("ComplexReflectionGroupType.jl")
include("complex_reflection_group.jl")
include("unitary_matrices.jl")
include("symplectic_reflection_group.jl")
include("ComplexReflection.jl")