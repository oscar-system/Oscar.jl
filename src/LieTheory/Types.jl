###############################################################################
#
#   Root systems, root space elements, weight lattice elements
#
###############################################################################

@doc raw"""
    RootSystem

Type for abstract root systems.

See [`root_system(::ZZMatrix)`](@ref) for the constructor.
"""
@attributes mutable struct RootSystem
  cartan_matrix::ZZMatrix # (generalized) Cartan matrix
  positive_roots::Vector #::Vector{RootSpaceElem} (cyclic reference)
  positive_roots_map::Dict{QQMatrix,Int}
  positive_coroots::Vector #::Vector{DualRootSpaceElem} (cyclic reference)
  positive_coroots_map::Dict{QQMatrix,Int}
  weyl_group::Any #::WeylGroup (cyclic reference)
  weight_lattice::Any #::WeightLattice (cyclic reference)

  # optional:
  type::Vector{Tuple{Symbol,Int}}
  type_ordering::Vector{Int}

  function RootSystem(mat::ZZMatrix; check::Bool=true, detect_type::Bool=true)
    check && @req is_cartan_matrix(mat) "Requires a generalized Cartan matrix"

    pos_roots, pos_coroots, refl = _positive_roots_and_reflections(mat)
    finite = count(refl .== 0) == nrows(mat)

    R = new(mat)
    R.positive_roots = map(r -> RootSpaceElem(R, r), pos_roots)
    R.positive_roots_map = Dict(
      (coefficients(root), ind) for
      (ind, root) in enumerate(R.positive_roots::Vector{RootSpaceElem})
    )
    R.positive_coroots = map(r -> DualRootSpaceElem(R, r), pos_coroots)
    R.positive_coroots_map = Dict(
      (coefficients(root), ind) for
      (ind, root) in enumerate(R.positive_coroots::Vector{DualRootSpaceElem})
    )
    R.weyl_group = WeylGroup(finite, refl, R)
    R.weight_lattice = WeightLattice(R)

    detect_type && is_finite(weyl_group(R)) && assure_root_system_type(R)
    return R
  end
end

@doc raw"""
    RootSpaceElem

Type for roots and linear combinations thereof.
"""
struct RootSpaceElem
  root_system::RootSystem
  vec::QQMatrix # the coordinate (row) vector with respect to the simple roots

  @doc raw"""
      RootSpaceElem(R::RootSystem, vec::QQMatrix) -> RootSpaceElem

  Construct a root space element in the root system `R` with the given coefficien vector w.r.t. the simple roots of `R`.

  `vec` must be a row vector of the same length as the rank of `R`.
  """
  function RootSpaceElem(R::RootSystem, vec::QQMatrix)
    @req size(vec) == (1, rank(R)) "Invalid dimension"
    return new(R, vec)
  end
end

@doc raw"""
    DualRootSpaceElem

Type for coroots and linear combinations thereof.
"""
mutable struct DualRootSpaceElem
  root_system::RootSystem
  vec::QQMatrix # the coordinate (row) vector with respect to the simple coroots

  @doc raw"""
      DualRootSpaceElem(R::RootSystem, vec::QQMatrix) -> DualRootSpaceElem

  Construct a dual root space element in the root system `R` with the given coefficien vector w.r.t. the simple coroots of `R`.

  `vec` must be a row vector of the same length as the rank of `R`.
  """
  function DualRootSpaceElem(root_system::RootSystem, vec::QQMatrix)
    @req size(vec) == (1, rank(root_system)) "Invalid dimension"
    return new(root_system, vec)
  end
end

###############################################################################
#
#   Weight lattices
#
###############################################################################

@doc raw"""
    WeightLattice <: AbstractAlgebra.AdditiveGroup

Type for weight lattices, parent type of `WeightLatticeElem`.
"""
@attributes mutable struct WeightLattice <: AbstractAlgebra.AdditiveGroup
  root_system::RootSystem

  function WeightLattice(root_system::RootSystem)
    return new(root_system)
  end
end

@doc raw"""
    WeightLatticeElem <: AbstractAlgebra.AdditiveGroupElem

Type for weights and linear combinations thereof, elem type of `WeightLattice`.
"""
mutable struct WeightLatticeElem <: AbstractAlgebra.AdditiveGroupElem
  parent_lat::WeightLattice
  vec::ZZMatrix # the coordinate (row) vector with respect to the fundamental weights

  @doc raw"""
      WeightLatticeElem(P::WeightLattice, vec::ZZMatrix) -> WeightLatticeElem

  Construct a weight lattice element in `P` with the given coefficients w.r.t. the fundamental weights of corresponding root system.

  `vec` must be a row vector of the same length as the rank of `P`.
  """
  function WeightLatticeElem(P::WeightLattice, vec::ZZMatrix)
    @req size(vec) == (1, rank(P)) "Invalid dimension"
    return new(P, vec)
  end
end

###############################################################################
#
#   Weyl groups
#
###############################################################################

@doc raw"""
    WeylGroup <: Group

Type for Weyl groups of root systems.

See [`weyl_group(::RootSystem)`](@ref) for the constructor.
"""
@attributes mutable struct WeylGroup <: Group
  finite::Bool              # finite indicates whether the Weyl group is finite
  refl::Matrix{UInt}        # see _positive_roots_and_reflections
  root_system::RootSystem   # root_system is the RootSystem from which the Weyl group was constructed

  function WeylGroup(finite::Bool, refl::Matrix{UInt}, root_system::RootSystem)
    return new(finite, refl, root_system)
  end
end

@doc raw"""
    WeylGroupElem <: GroupElem

Type for elements of Weyl groups.
"""
struct WeylGroupElem <: GroupElem
  parent::WeylGroup     # parent group
  word::Vector{UInt8}   # short revlex normal form of the word

  function WeylGroupElem(W::WeylGroup, word::Vector{<:Integer}; normalize::Bool=true)
    if !normalize
      if word isa Vector{UInt8}
        return new(W, word)
      else
        return new(W, UInt8.(word))
      end
    end

    @req all(1 <= i <= ngens(W) for i in word) "word contains invalid generators"
    x = new(W, sizehint!(UInt8[], length(word)))
    for s in word
      rmul!(x, s)
    end

    return x
  end
end

const WeylIteratorNoCopyState = Tuple{WeightLatticeElem,WeylGroupElem}

@doc raw"""
    ReducedExpressionIterator

Iterator for reduced expressions of a Weyl group element.

See [`reduced_expressions(::WeylGroupElem)`](@ref) for the constructor.
"""
struct ReducedExpressionIterator
  el::WeylGroupElem         # the Weyl group element for which we a searching reduced expressions
  up_to_commutation::Bool   # if true and say s1 and s3 commute, we only list s3*s1 and not s1*s3
end

struct WeylIteratorNoCopy
  weight::WeightLatticeElem # dominant weight
  weyl_group::WeylGroup

  function WeylIteratorNoCopy(wt::WeightLatticeElem)
    return new(conjugate_dominant_weight(wt), weyl_group(root_system(wt)))
  end
end

struct WeylOrbitIterator
  nocopy::WeylIteratorNoCopy

  function WeylOrbitIterator(wt::WeightLatticeElem)
    return new(WeylIteratorNoCopy(wt))
  end
end
