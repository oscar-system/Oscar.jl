
# A union type for all algebraic extensions of the field of rational numbers
const QQAlgField = Union{NumField, QQField, QQBarField, QQAbField}
const QQAlgFieldElem = Union{NumFieldElem, QQFieldElem, QQBarFieldElem, QQAbElem}

# Imports (for stuff from experimental)
import Oscar.LieAlgebras: coroot #no conflict, just same function name

# Project files
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

# Exports
export canonical_pairing
export codegrees
export coexponents
export ComplexReflection
export ComplexReflectionGroupType
export complex_reflection
export complex_reflections
export complex_reflection_group
export complex_reflection_group_cartan_matrix
export complex_reflection_group_dual
export complex_reflection_group_model
export complex_reflection_group_type
export components
export coroot_form
export coxeter_number
export degrees
export eigenvalue
export hyperplane
export hyperplane_basis
export hyperplane_inclusion
export is_complex_reflection
export is_complex_reflection_group
export is_complex_reflection_with_data
export is_coxeter_group
export is_imprimitive
export is_irreducible
export is_primitive
export is_pseudo_real
export is_orthogonal
export is_rational
export is_real
export is_root_of_unity_with_data
export is_root_of_unity
export is_spetsial
export is_symplectic_reflection_group
export is_weyl_group
export is_well_generated
export is_unitary
export linear_form
export number_of_components
export number_of_hyperplanes
export number_of_reflections
export number_of_reflection_classes
export root_line
export root_line_inclusion
export symplectic_doubling
export symplectic_reflection_group
export unitary_reflection

# Aliases
@alias n_components number_of_components
@alias n_hyperplanes number_of_hyperplanes
@alias n_reflections number_of_reflections
@alias n_reflection_classes number_of_reflection_classes
