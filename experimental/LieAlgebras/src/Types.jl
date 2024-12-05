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

    pos_roots, pos_coroots, refl = positive_roots_and_reflections(mat)
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
mutable struct RootSpaceElem
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
@attributes mutable struct WeylGroup <: AbstractAlgebra.Group
  finite::Bool              # finite indicates whether the Weyl group is finite
  refl::Matrix{UInt}        # see positive_roots_and_reflections
  root_system::RootSystem   # root_system is the RootSystem from which the Weyl group was constructed

  function WeylGroup(finite::Bool, refl::Matrix{UInt}, root_system::RootSystem)
    return new(finite, refl, root_system)
  end
end

@doc raw"""
    WeylGroupElem <: GroupElem

Type for elements of Weyl groups.
"""
struct WeylGroupElem <: AbstractAlgebra.GroupElem
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

###############################################################################
#
#   Lie algebras
#
###############################################################################

abstract type LieAlgebra{C<:FieldElem} <: AbstractAlgebra.Set end

abstract type LieAlgebraElem{C<:FieldElem} <: AbstractAlgebra.SetElem end

###############################################################################
# AbstractLieAlgebra

@attributes mutable struct AbstractLieAlgebra{C<:FieldElem} <: LieAlgebra{C}
  R::Field
  dim::Int
  struct_consts::Matrix{<:SRow{C}}
  s::Vector{Symbol}

  # only set if known
  root_system::RootSystem
  chevalley_basis::NTuple{3,Vector{<:LieAlgebraElem{C}}}

  function AbstractLieAlgebra{C}(
    R::Field,
    struct_consts::Matrix{<:SRow{C}},
    s::Vector{Symbol};
    check::Bool=true,
  ) where {C<:FieldElem}
    @assert struct_consts isa Matrix{sparse_row_type(R)} "Invalid structure constants type."
    (n1, n2) = size(struct_consts)
    @req n1 == n2 "Invalid structure constants dimensions."
    dimL = n1
    @req length(s) == dimL "Invalid number of basis element names."
    if check
      @req all(
        r -> all(e -> parent(last(e)) === R, r), struct_consts
      ) "Invalid structure constants."
      @req all(iszero, struct_consts[i, i] for i in 1:dimL) "Not anti-symmetric."
      @req all(
        iszero, struct_consts[i, j] + struct_consts[j, i] for i in 1:dimL, j in 1:dimL
      ) "Not anti-symmetric."
      @req all(
        iszero,
        begin
          row = sparse_row(R)
          for (k, k_val) in struct_consts[i, j]
            row = addmul!(row, k_val, struct_consts[k, l])
          end
          for (k, k_val) in struct_consts[j, l]
            row = addmul!(row, k_val, struct_consts[k, i])
          end
          for (k, k_val) in struct_consts[l, i]
            row = addmul!(row, k_val, struct_consts[k, j])
          end
          row
        end for i in 1:dimL, j in 1:dimL, l in 1:dimL
      ) "Jacobi identity does not hold."
    end
    return new{C}(R, dimL, struct_consts, s)
  end
end

struct AbstractLieAlgebraElem{C<:FieldElem} <: LieAlgebraElem{C}
  parent::AbstractLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
# LinearLieAlgebra

@attributes mutable struct LinearLieAlgebra{C<:FieldElem} <: LieAlgebra{C}
  R::Field
  n::Int  # the n of the gl_n this embeds into
  dim::Int
  basis::Vector{<:MatElem{C}}
  s::Vector{Symbol}

  # only set if known
  root_system::RootSystem
  chevalley_basis::NTuple{3,Vector{<:LieAlgebraElem{C}}}

  function LinearLieAlgebra{C}(
    R::Field,
    n::Int,
    basis::Vector{<:MatElem{C}},
    s::Vector{Symbol};
    check::Bool=true,
  ) where {C<:FieldElem}
    @req all(b -> size(b) == (n, n), basis) "Invalid basis element dimensions."
    @req length(s) == length(basis) "Invalid number of basis element names."
    @req eltype(basis) == dense_matrix_type(R) "Invalid basis element type."
    L = new{C}(R, n, length(basis), basis, s)
    if check
      @req all(b -> all(e -> parent(e) === R, b), basis) "Invalid matrices."
      # TODO: make work
      # for xi in basis(L), xj in basis(L)
      #   @req (xi * xj) in L
      # end
    end
    return L
  end
end

struct LinearLieAlgebraElem{C<:FieldElem} <: LieAlgebraElem{C}
  parent::LinearLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
# DirectSumLieAlgebra

@attributes mutable struct DirectSumLieAlgebra{C<:FieldElem} <: LieAlgebra{C}
  R::Field
  dim::Int
  summands::Vector{<:LieAlgebra{C}}
  s::Vector{Symbol}

  # only set if known
  root_system::RootSystem
  chevalley_basis::NTuple{3,Vector{<:LieAlgebraElem{C}}}

  function DirectSumLieAlgebra{C}(
    R::Field,
    summands::Vector{<:LieAlgebra{C}},
  ) where {C<:FieldElem}
    @req all(x -> coefficient_ring(x) == R, summands) "Summands must have the same coefficient ring."
    totaldim = sum(dim, summands; init=0)
    s = [Symbol("$(x)^($(i))") for (i, S) in enumerate(summands) for x in symbols(S)]
    L = new{C}(R, totaldim, summands, s)
    # TODO: glue root systems if all summands have one
    return L
  end
end

struct DirectSumLieAlgebraElem{C<:FieldElem} <: LieAlgebraElem{C}
  parent::DirectSumLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Lie algebra constructions
#
###############################################################################

###############################################################################
# LieSubalgebra

@attributes mutable struct LieSubalgebra{C<:FieldElem,LieT<:LieAlgebraElem{C}}
  base_lie_algebra::LieAlgebra{C}
  gens::Vector{LieT}
  basis_elems::Vector{LieT}
  basis_matrix::MatElem{C}

  function LieSubalgebra{C,LieT}(
    L::LieAlgebra{C}, gens::Vector{LieT}; is_basis::Bool=false
  ) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
    @req all(g -> parent(g) === L, gens) "Parent mismatch."
    @req L isa parent_type(LieT) "Parent type mismatch."
    if is_basis
      basis_elems = gens
      basis_matrix = if length(gens) == 0
        matrix(coefficient_ring(L), 0, dim(L), C[])
      else
        matrix(coefficient_ring(L), [coefficients(g) for g in gens])
      end
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
    else
      basis_matrix = matrix(coefficient_ring(L), 0, dim(L), C[])
      gens = unique(g for g in gens if !iszero(g))
      left = copy(gens)
      while !isempty(left)
        g = pop!(left)
        can_solve(basis_matrix, _matrix(g); side=:left) && continue
        for row in eachrow(basis_matrix)
          push!(left, g * L(row))
        end
        basis_matrix = vcat(basis_matrix, _matrix(g))
        rank = rref!(basis_matrix)
        @assert rank == nrows(basis_matrix) # otherwise the continue above would've triggered
      end
      basis_elems = L.(eachrow(basis_matrix))
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
    end
  end

  function LieSubalgebra{C,LieT}(
    L::LieAlgebra{C}, gens::Vector; kwargs...
  ) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
    return LieSubalgebra{C,LieT}(L, Vector{LieT}(map(L, gens)); kwargs...)
  end
end

###############################################################################
# LieAlgebraIdeal

@attributes mutable struct LieAlgebraIdeal{C<:FieldElem,LieT<:LieAlgebraElem{C}}
  base_lie_algebra::LieAlgebra{C}
  gens::Vector{LieT}
  basis_elems::Vector{LieT}
  basis_matrix::MatElem{C}

  function LieAlgebraIdeal{C,LieT}(
    L::LieAlgebra{C}, gens::Vector{LieT}; is_basis::Bool=false
  ) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
    @req all(g -> parent(g) === L, gens) "Parent mismatch."
    L::parent_type(LieT)
    if is_basis
      basis_elems = gens
      basis_matrix = if length(gens) == 0
        matrix(coefficient_ring(L), 0, dim(L), C[])
      else
        matrix(coefficient_ring(L), [coefficients(g) for g in gens])
      end
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
    else
      basis_matrix = matrix(coefficient_ring(L), 0, dim(L), C[])
      gens = unique(g for g in gens if !iszero(g))
      left = copy(gens)
      while !isempty(left)
        g = pop!(left)
        can_solve(basis_matrix, _matrix(g); side=:left) && continue
        for b in basis(L)
          push!(left, b * g)
        end
        basis_matrix = vcat(basis_matrix, _matrix(g))
        rank = rref!(basis_matrix)
        @assert rank == nrows(basis_matrix) # otherwise the continue above would've triggered
      end
      basis_elems = L.(eachrow(basis_matrix))
      return new{C,LieT}(L, gens, basis_elems, basis_matrix)
    end
  end

  function LieAlgebraIdeal{C,LieT}(
    L::LieAlgebra{C}, gens::Vector; kwargs...
  ) where {C<:FieldElem,LieT<:LieAlgebraElem{C}}
    return LieAlgebraIdeal{C,LieT}(L, Vector{LieT}(map(L, gens)); kwargs...)
  end
end

###############################################################################
# LieAlgebraHom

@attributes mutable struct LieAlgebraHom{T1<:LieAlgebra,T2<:LieAlgebra,MatT<:MatElem} <:
                           Map{T1,T2,Hecke.HeckeMap,LieAlgebraHom}
  header::MapHeader{T1,T2}
  matrix::MatT

  inverse_isomorphism::LieAlgebraHom{T2,T1}

  function LieAlgebraHom(
    L1::LieAlgebra, L2::LieAlgebra, imgs::Vector{<:LieAlgebraElem}; check::Bool=true
  )
    @req coefficient_ring(L1) === coefficient_ring(L2) "Coefficient rings must be the same" # for now at least
    @req all(x -> parent(x) === L2, imgs) "Images must lie in the codomain"
    @req length(imgs) == dim(L1) "Number of images must match dimension of domain"

    mat = zero_matrix(coefficient_ring(L2), dim(L1), dim(L2))
    for (i, img) in enumerate(imgs)
      mat[i, :] = _matrix(img)
    end
    return LieAlgebraHom(L1, L2, mat; check)
  end

  function LieAlgebraHom(L1::LieAlgebra, L2::LieAlgebra, mat::MatElem; check::Bool=true)
    @req coefficient_ring(L1) === coefficient_ring(L2) "Coefficient rings must be the same" # for now at least
    @req size(mat) == (dim(L1), dim(L2)) "Matrix size must match dimensions of domain and codomain"
    @req mat isa MatElem{elem_type(coefficient_ring(L2))} "Matrix must be over coefficient ring of codomain"
    h = new{typeof(L1),typeof(L2),typeof(mat)}()
    h.matrix = mat
    h.header = MapHeader(L1, L2)
    if check
      @req is_welldefined(h) "Not a homomorphism"
    end
    return h
  end
end

###############################################################################
#
#   Lie algebra modules
#
###############################################################################

@attributes mutable struct LieAlgebraModule{C<:FieldElem,LieT<:LieAlgebraElem{C}} <:
                           AbstractAlgebra.Set
  L::LieAlgebra{C} # parent_type(LieT)
  dim::Int
  transformation_matrices::Vector{<:MatElem{C}} # Vector{dense_matrix_type(C)}
  s::Vector{Symbol}

  function LieAlgebraModule{C}(
    L::LieAlgebra{C},
    dimV::Int,
    transformation_matrices::Vector{<:MatElem{C}},
    s::Vector{Symbol};
    check::Bool=true,
  ) where {C<:FieldElem}
    @req dimV == length(s) "Invalid number of basis element names."
    @req dim(L) == length(transformation_matrices) "Invalid number of transformation matrices."
    @req all(m -> size(m) == (dimV, dimV), transformation_matrices) "Invalid transformation matrix dimensions."

    V = new{C,elem_type(L)}(L, dimV, transformation_matrices, s)
    if check
      @req all(m -> all(e -> parent(e) === coefficient_ring(V), m), transformation_matrices) "Invalid transformation matrix entries."
      for xi in basis(L), xj in basis(L), v in basis(V)
        @req (xi * xj) * v == xi * (xj * v) - xj * (xi * v) "Transformation matrices do not define a module."
      end
    end
    return V
  end
end

struct LieAlgebraModuleElem{C<:FieldElem,LieT<:LieAlgebraElem{C}} <: AbstractAlgebra.SetElem
  parent::LieAlgebraModule{C,LieT}
  mat::MatElem{C}
end

###############################################################################
# LieAlgebraModuleHom

@attributes mutable struct LieAlgebraModuleHom{T1<:LieAlgebraModule,T2<:LieAlgebraModule} <:
                           Map{T1,T2,Hecke.HeckeMap,LieAlgebraModuleHom}
  header::MapHeader{T1,T2}
  matrix::MatElem

  inverse_isomorphism::LieAlgebraModuleHom{T2,T1}

  function LieAlgebraModuleHom(
    V1::LieAlgebraModule,
    V2::LieAlgebraModule,
    imgs::Vector{<:LieAlgebraModuleElem};
    check::Bool=true,
  )
    @req base_lie_algebra(V1) === base_lie_algebra(V2) "Lie algebras must be the same" # for now at least
    @req all(x -> parent(x) === V2, imgs) "Images must lie in the codomain"
    @req length(imgs) == dim(V1) "Number of images must match dimension of domain"

    mat = zero_matrix(coefficient_ring(V2), dim(V1), dim(V2))
    for (i, img) in enumerate(imgs)
      mat[i, :] = _matrix(img)
    end
    return LieAlgebraModuleHom(V1, V2, mat; check)
  end

  function LieAlgebraModuleHom(
    V1::LieAlgebraModule, V2::LieAlgebraModule, mat::MatElem; check::Bool=true
  )
    @req base_lie_algebra(V1) === base_lie_algebra(V2) "Lie algebras must be the same" # for now at least
    @req size(mat) == (dim(V1), dim(V2)) "Matrix size must match dimensions of domain and codomain"
    h = new{typeof(V1),typeof(V2)}()
    h.matrix = mat::dense_matrix_type(coefficient_ring(V2))
    h.header = MapHeader(V1, V2)
    if check
      @req is_welldefined(h) "Not a homomorphism"
    end
    return h
  end
end
