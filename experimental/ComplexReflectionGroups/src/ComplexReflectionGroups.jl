
# A union type for all algebraic extension of the field of rational numbers
QQAlgField = Union{NumField, QQField, QQBarField, QQAbField}
QQAlgFieldElem = Union{NumFieldElem, QQFieldElem, QQBarFieldElem, QQAbElem}

include("dual_vector_space.jl")
include("is_root_of_unity.jl")
include("hermitian_things.jl")

include("ComplexReflection.jl")

include("ComplexReflectionGroupType.jl")

include("complex_reflection_group.jl")
include("complex_reflection_group_LT.jl")
include("complex_reflection_group_Magma.jl")
include("complex_reflection_group_CHEVIE.jl")

include("symplectic_reflection_group.jl")