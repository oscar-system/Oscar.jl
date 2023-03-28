@attributes mutable struct LieAlgebraAbstractModule{C<:RingElement} <: LieAlgebraModule{C}
  L::LieAlgebra{C}
  dim::Int
  transformation_matrices::Vector{MatElem{C}}
  s::Vector{Symbol}

  function LieAlgebraAbstractModule{C}(
    L::LieAlgebra{C},
    dimV::Int,
    transformation_matrices::Vector{<:MatElem{C}},
    s::Vector{Symbol};
    cached::Bool=true,
    check::Bool=true,
  ) where {C<:RingElement}
    return get_cached!(
      LieAlgebraAbstractModuleDict, (L, dimV, transformation_matrices, s), cached
    ) do
      @req dimV == length(s) "Invalid number of basis element names."
      @req dim(L) == length(transformation_matrices) "Invalid number of transformation matrices."
      @req all(
        m -> size(m) == (dimV, dimV), transformation_matrices
      ) "Invalid transformation matrix dimensions."

      V = new{C}(L, dimV, transformation_matrices, s)
      if check
        for xi in basis(L), xj in basis(L), v in basis(V)
          @req bracket(xi, xj) * v ==
            xi * (xj * v) - xj * (xi * v) "Transformation matrices do not define a module."
        end
      end
      V
    end::LieAlgebraAbstractModule{C}
  end
end

const LieAlgebraAbstractModuleDict = CacheDictType{
  Tuple{LieAlgebra,Int,Vector{MatElem},Vector{Symbol}},LieAlgebraAbstractModule
}()

struct LieAlgebraAbstractModuleElem{C<:RingElement} <: LieAlgebraModuleElem{C}
  parent::LieAlgebraAbstractModule{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraAbstractModuleElem{C}}) where {C<:RingElement} =
  LieAlgebraAbstractModule{C}

elem_type(::Type{LieAlgebraAbstractModule{C}}) where {C<:RingElement} =
  LieAlgebraAbstractModuleElem{C}

parent(v::LieAlgebraAbstractModuleElem{C}) where {C<:RingElement} = v.parent

base_ring(V::LieAlgebraAbstractModule{C}) where {C<:RingElement} =
  base_ring(base_liealgebra(V))

base_liealgebra(V::LieAlgebraAbstractModule{C}) where {C<:RingElement} = V.L

dim(V::LieAlgebraAbstractModule{C}) where {C<:RingElement} = V.dim

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  if has_attribute(V, :show) && get_attribute(V, :show) isa Function
    get_attribute(V, :show)(io, V)
  else
    print(io, "AbstractModule of ")
    print(IOContext(io, :compact => true), base_liealgebra(V))
  end
end

function show_standard_module(io::IO, V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  print(io, "StdModule of ")
  print(IOContext(io, :compact => true), base_liealgebra(V))
end

function show_exterior_power(io::IO, V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th exterior power of ")
  print(IOContext(io, :compact => true), get_attribute(V, :inner_module))
end

function show_symmetric_power(io::IO, V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th symmetric power of ")
  print(IOContext(io, :compact => true), get_attribute(V, :inner_module))
end

function show_tensor_power(io::IO, V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  print(io, "$(get_attribute(V, :power))-th tensor power of ")
  print(IOContext(io, :compact => true), get_attribute(V, :inner_module))
end

function symbols(V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  return V.s
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraAbstractModule{C})(
  a::Vector{T}
) where {T<:LieAlgebraAbstractModuleElem{C}} where {C<:RingElement}
  @req is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V) "Only implemented for power modules."
  @req length(a) == get_attribute(V, :power) "Length of vector does not match power."
  @req all(x -> parent(x) == get_attribute(V, :inner_module), a) "Incompatible modules."
  mat = zero_matrix(base_ring(V), 1, dim(V))
  if is_exterior_power(V)
    for (i, _inds) in enumerate(get_attribute(V, :ind_map)),
      inds in Combinatorics.permutations(_inds)

      sgn = Combinatorics.levicivita(sortperm(inds))
      mat[1, i] += sgn * prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
  elseif is_symmetric_power(V)
    for (i, _inds) in enumerate(get_attribute(V, :ind_map)),
      inds in unique(Combinatorics.permutations(_inds))

      mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
  elseif is_tensor_power(V)
    for (i, inds) in enumerate(get_attribute(V, :ind_map))
      mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
  end
  return LieAlgebraAbstractModuleElem{C}(V, mat)
end

###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_by_basisindex(
  V::LieAlgebraAbstractModule{C}, i::Int
) where {C<:RingElement}
  return (V.transformation_matrices[i])::dense_matrix_type(C)
end

###############################################################################
#
#   Attribute getters
#
###############################################################################

function is_standard_module(V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :standard_module
end

function is_exterior_power(V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :exterior_power
end

function is_symmetric_power(V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :symmetric_power
end

function is_tensor_power(V::LieAlgebraAbstractModule{C}) where {C<:RingElement}
  return get_attribute(V, :type, :fallback) == :tensor_power
end

###############################################################################
#
#   Constructor
#
###############################################################################

function abstract_module(
  L::LieAlgebra{C},
  dimV::Int,
  transformation_matrices::Vector{<:MatElem{C}},
  s::Vector{<:Union{AbstractString,Char,Symbol}}=[Symbol("v_$i") for i in 1:dimV];
  cached::Bool=true,
  check::Bool=true,
) where {C<:RingElement}
  return LieAlgebraAbstractModule{C}(
    L, dimV, transformation_matrices, Symbol.(s); cached, check=check
  )
end

function abstract_module(
  L::LieAlgebra{C},
  dimV::Int,
  struct_consts::Matrix{SRow{C}},
  s::Vector{<:Union{AbstractString,Char,Symbol}}=[Symbol("v_$i") for i in 1:dimV];
  cached::Bool=true,
  check::Bool=true,
) where {C<:RingElement}
  @req dim(L) == size(struct_consts, 1) "Invalid structure constants dimensions."
  @req dimV == size(struct_consts, 2) "Invalid structure constants dimensions."
  @req dimV == length(s) "Invalid number of basis element names."

  transformation_matrices = [zero_matrix(base_ring(L), dimV, dimV) for _ in 1:dim(L)]
  for i in 1:dim(L), j in 1:dimV
    transformation_matrices[i][:, j] = transpose(dense_row(struct_consts[i, j], dimV))
  end

  return LieAlgebraAbstractModule{C}(
    L, dimV, transformation_matrices, Symbol.(s); cached, check=check
  )
end

function highest_weight_module(
  L::LieAlgebra{C}, weight::Vector{Int}; cached::Bool=true
) where {C<:RingElement}
  struct_consts = liealgebra_highest_weight_module_struct_consts_gap(L, weight)
  dimV = size(struct_consts, 2)
  V = abstract_module(L, dimV, struct_consts; cached, check=false)
  set_attribute!(V, :highest_weight => weight)
  return V
end

function standard_module(L::LinearLieAlgebra{C}; cached::Bool=true) where {C<:RingElement}
  dim_std_V = L.n
  transformation_matrices = matrix_repr_basis(L)
  s = [Symbol("v_$(i)") for i in 1:dim_std_V]
  std_V = LieAlgebraAbstractModule{elem_type(base_ring(L))}(
    L, dim_std_V, transformation_matrices, s; cached, check=false
  )
  set_attribute!(std_V, :type => :standard_module, :show => show_standard_module)
  return std_V
end

function exterior_power(
  V::LieAlgebraAbstractModule{C}, k::Int; cached::Bool=true
) where {C<:RingElement}
  L = base_liealgebra(V)
  dim_pow_V = binomial(dim(V), k)
  ind_map = collect(Combinatorics.combinations(1:dim(V), k))

  T = tensor_power(V, k)
  basis_change_E2T = zero_matrix(base_ring(V), dim(T), dim_pow_V)
  basis_change_T2E = zero_matrix(base_ring(V), dim_pow_V, dim(T))
  T_ind_map = get_attribute(T, :ind_map)
  for (i, _inds) in enumerate(ind_map), inds in Combinatorics.permutations(_inds)
    sgn = Combinatorics.levicivita(sortperm(inds))
    j = findfirst(==(inds), T_ind_map)
    basis_change_E2T[j, i] = sgn//factorial(k)
    basis_change_T2E[i, j] = sgn
  end
  transformation_matrices = map(1:dim(L)) do i
    basis_change_T2E * transformation_matrix_by_basisindex(T, i) * basis_change_E2T
  end

  if k == 1
    s = symbols(V)
  else
    if is_standard_module(V)
      parentheses = identity
    else
      parentheses = x -> "($x)"
    end
    s = [
      Symbol(join(parentheses.(s), " ∧ ")) for
      s in Combinatorics.combinations(symbols(V), k)
    ]
  end

  pow_V = LieAlgebraAbstractModule{C}(
    L, dim_pow_V, transformation_matrices, s; cached, check=false
  )
  set_attribute!(
    pow_V,
    :type => :exterior_power,
    :power => k,
    :inner_module => V,
    :ind_map => ind_map,
    :show => show_exterior_power,
  )
  return pow_V
end

function symmetric_power(
  V::LieAlgebraAbstractModule{C}, k::Int; cached::Bool=true
) where {C<:RingElement}
  L = base_liealgebra(V)
  dim_pow_V = binomial(dim(V) + k - 1, k)
  ind_map = collect(Combinatorics.with_replacement_combinations(1:dim(V), k))

  T = tensor_power(V, k)
  basis_change_S2T = zero_matrix(base_ring(V), dim(T), dim_pow_V)
  basis_change_T2S = zero_matrix(base_ring(V), dim_pow_V, dim(T))
  T_ind_map = get_attribute(T, :ind_map)
  for (i, _inds) in enumerate(ind_map), inds in Combinatorics.permutations(_inds)
    sgn = Combinatorics.levicivita(sortperm(inds))
    j = findfirst(==(inds), T_ind_map)
    basis_change_S2T[j, i] += 1//factorial(k)
    basis_change_T2S[i, j] = 1
  end
  transformation_matrices = map(1:dim(L)) do i
    basis_change_T2S * transformation_matrix_by_basisindex(T, i) * basis_change_S2T
  end

  if k == 1
    s = symbols(V)
  else
    if is_standard_module(V)
      parentheses = identity
    else
      parentheses = x -> "($x)"
    end
    s = [
      Symbol(
        join(
          (
            begin
              e = count(==(i), inds)
              if e == 1
                parentheses(s)
              else
                "$(parentheses(s))^$e"
              end
            end for (i, s) in enumerate(symbols(V)) if in(i, inds)
          ),
          "*",
        ),
      ) for inds in ind_map
    ]
  end

  pow_V = LieAlgebraAbstractModule{C}(
    L, dim_pow_V, transformation_matrices, s; cached, check=false
  )
  set_attribute!(
    pow_V,
    :type => :symmetric_power,
    :power => k,
    :inner_module => V,
    :ind_map => ind_map,
    :show => show_symmetric_power,
  )
  return pow_V
end

function tensor_power(
  V::LieAlgebraAbstractModule{C}, k::Int; cached::Bool=true
) where {C<:RingElement}
  L = base_liealgebra(V)
  dim_pow_V = dim(V)^k
  ind_map = reverse.(collect(ProductIterator(1:dim(V), k)))

  transformation_matrices = map(1:dim(L)) do i
    y = transformation_matrix_by_basisindex(V, i)
    sum(reduce(kronecker_product, (j == i ? y : one(y) for j in 1:k)) for i in 1:k)
  end

  if k == 1
    s = symbols(V)
  else
    if is_standard_module(V)
      parentheses = identity
    else
      parentheses = x -> "($x)"
    end
    s = [
      Symbol(join(parentheses.(s), " ⊗ ")) for s in reverse.(ProductIterator(symbols(V), k))
    ]
  end

  pow_V = LieAlgebraAbstractModule{C}(
    L, dim_pow_V, transformation_matrices, s; cached, check=false
  )
  set_attribute!(
    pow_V,
    :type => :tensor_power,
    :power => k,
    :inner_module => V,
    :ind_map => ind_map,
    :show => show_tensor_power,
  )
  return pow_V
end
