
# A union type for all algebraic extensions of the field of rational numbers
const QQAlgField = Union{NumField, QQField, QQBarField, QQAbField}
const QQAlgFieldElem = Union{NumFieldElem, QQFieldElem, QQBarFieldElem, QQAbFieldElem}

# Imports (for stuff from experimental)
import Oscar.LieAlgebras: coroot #no conflict, just same function name

# Project files
include("dual_vector_space.jl")
include("is_root_of_unity.jl")
include("hermitian_things.jl")

include("ComplexReflection.jl")

include("ComplexReflectionGroupType.jl")
include("reflection_data.jl")

include("complex_reflection_group.jl")
include("complex_reflection_group_LT.jl")
include("complex_reflection_group_Magma.jl")
include("complex_reflection_group_CHEVIE.jl")
include("reflection_orbits.jl")

include("symplectic_reflection_group.jl")

# Exports
export ComplexReflection
export ComplexReflectionClass
export ComplexReflectionGroupType
export ComplexReflectionHyperplaneOrbit
export ReflectionClassType
export ReflectionHyperplaneOrbitType
export canonical_pairing
export codegrees
export coexponents
export complex_reflection
export complex_reflection_group
export complex_reflection_group_cartan_matrix
export complex_reflection_group_component_embeddings
export complex_reflection_group_dual
export complex_reflection_group_dual_source
export complex_reflection_group_dual_source_model
export complex_reflection_group_model
export complex_reflection_group_type
export complex_reflections
export component_index
export components
export coroot_form
export coxeter_number
export degrees
export distinguished_reflection
export distinguished_reflections
export eigenvalue
export hyperplane
export hyperplane_basis
export hyperplane_inclusion
export hyperplane_orbit
export invariant_hermitian_form
export is_complex_reflection
export is_complex_reflection_group
export is_complex_reflection_with_data
export is_coxeter_group
export is_equivalent
export is_imprimitive
export is_irreducible
export is_orthogonal
export is_primitive
export is_pseudo_real
export is_rational
export is_real
export is_root_of_unity
export is_root_of_unity_with_data
export is_spetsial
export is_symplectic_reflection_group
export is_unitary
export is_well_generated
export is_weyl_group
export linear_form
export local_orbit_index
export number_of_components
export number_of_hyperplanes
export number_of_reflection_classes
export number_of_reflections
export orbit_size
export pointwise_stabilizer_order
export reflection_class_type
export reflection_classes
export reflection_exponent
export reflection_hyperplane_orbit_type
export reflection_hyperplane_orbits
export reflection_hyperplanes
export reflection_library
export representative_word
export root_line
export root_line_inclusion
export set_reflection_hyperplane_orbit_marking!
export symplectic_doubling
export symplectic_doubling_block_embeddings
export symplectic_doubling_source
export symplectic_doubling_source_model
export symplectic_doubling_source_type
export symplectic_form
export symplectic_reflection_group
export unitary_reflection

# Aliases
@alias n_components number_of_components
@alias n_hyperplanes number_of_hyperplanes
@alias n_reflections number_of_reflections
@alias n_reflection_classes number_of_reflection_classes
