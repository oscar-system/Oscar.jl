## Generic functionality
@doc raw"""
    total_complex(D::AbsDoubleComplexOfMorphisms)

Construct the total complex of the double complex `D`. 

Note that `D` needs to be reasonably bounded for this to work so that the strands
``⨁ ᵢ₊ⱼ₌ₖ Dᵢⱼ`` are finite for every `k`. Moreover, the generic code uses the internal 
function `_direct_sum`. See the docstring of that function to learn more.
"""
function total_complex(D::AbsDoubleComplexOfMorphisms)
  is_bounded(D) || error("computation of total complexes is only implemented for bounded double complexes")
  vertical_direction(D) == horizontal_direction(D) || error("horizontal and vertical typs must be the same")
  if vertical_direction(D) == :chain
    return _total_chain_complex(D)
  else
    return _total_cochain_complex(D)
  end
end

function _total_chain_complex(dc::AbsDoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  r1 = horizontal_range(dc)
  r2 = vertical_range(dc)
  r_tot = (first(r1) + first(r2)):-1:(last(r1) + last(r2))
  new_chains = ChainType[]
  new_maps = MorphismType[]

  # First round for initialization
  index_pairs = [I for I in Iterators.product(r1, r2) if sum(I) == first(r_tot)]
  summands = [dc[I] for I in index_pairs]
  #new_chain, inc, pr = direct_sum(summands...; task=:both)
  new_chain, inc, pr = _direct_sum(summands)
  last_inc = inc
  last_pr = pr
  last_chain = new_chain
  last_index_pairs = index_pairs
  push!(new_chains, new_chain)
  for k in r_tot
    k == first(r_tot) && continue
    index_pairs = [I for I in Iterators.product(r1, r2) if sum(I) == k]
    summands = [dc[I] for I in index_pairs]
    #new_chain, inc, pr = direct_sum(summands...; task=:both)
    new_chain, inc, pr = _direct_sum(summands)
    @assert all(f->domain(f) === new_chain, pr)
    @assert all(f->codomain(f) === new_chain, inc)
    push!(new_chains, new_chain)
    # assemble the map
    boundary_map = hom(last_chain, new_chain, elem_type(new_chain)[zero(new_chain) for i in 1:ngens(last_chain)])
    for (l, I) in enumerate(last_index_pairs)
      p = last_pr[l]
      @assert domain(p) === last_chain

      if I[2] > lower_bound(dc)
        vert = compose(p, vertical_map(dc, I))
        m = findfirst(x->x==(I[1], I[2]-1), index_pairs)
        @assert m !== nothing
        @assert codomain(vert) === dc[I[1], I[2]-1] === domain(inc[m])
        boundary_map = boundary_map + compose(vert, inc[m])
      end

      if I[1] > left_bound(dc)
        horz = compose(p, horizontal_map(dc, I))
        n = findfirst(x->x==(I[1]-1, I[2]), index_pairs)
        @assert n !== nothing
        @assert codomain(horz) === dc[I[1]-1, I[2]] === domain(inc[n])
        boundary_map = boundary_map + (-1)^I[2]*compose(horz, inc[n])
      end
    end
    push!(new_maps, boundary_map)
    last_index_pairs = index_pairs
    last_chain = new_chain
    last_inc = inc
    last_pr = pr
  end
  return ComplexOfMorphisms(ChainType, new_maps, seed=last(r_tot))
end

function _total_cochain_complex(dc::AbsDoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  error("total complex of double cochain complexes currently not implemented")
end

### Missing functionality for complexes
typ(C::ComplexOfMorphisms) = C.typ
is_complete(C::ComplexOfMorphisms) = C.complete


### cached homology

function kernel(c::HyperComplex{ChainType}, p::Int, i::Tuple) where {ChainType <: SparseFPModule}
  if !isdefined(c, :kernel_cache) 
    c.kernel_cache = Dict{Tuple{Tuple, Int}, Map}()
  end
  if haskey(c.kernel_cache, (i, p))
    inc = c.kernel_cache[(i, p)]
    return domain(inc), inc
  end

  if !can_compute_map(c, p, i)
    K, inc = sub(c[i], gens(c[i]))
    c.kernel_cache[(i, p)] = inc
    return K, inc
  end

  @assert domain(map(c, p, i)) === c[i]
  K, inc = kernel(map(c, p, i))
  c.kernel_cache[(i, p)] = inc
  return K, inc
end

function boundary(c::HyperComplex{ChainType}, p::Int, i::Tuple) where {ChainType <: SparseFPModule}
  I = collect(i)
  prev = I + (direction(c, p) == :chain ? 1 : -1)*[k==p ? 1 : 0 for k in 1:dim(c)]
  Prev = Tuple(prev)

  if !isdefined(c, :boundary_cache) 
    c.boundary_cache = Dict{Tuple{Tuple, Int}, Map}()
  end
  if haskey(c.boundary_cache, (i, p))
    inc = c.boundary_cache[(i, p)]
    return domain(inc), inc
  end

  if !can_compute_map(c, p, Prev) 
    !can_compute_index(c, Prev) || error("map can not be computed")
    Im, inc = sub(c[i], elem_type(c[i])[])
    @assert codomain(inc) === c[i]
    c.boundary_cache[(i, p)] = inc
    return Im, inc
  end

  Im, inc = image(map(c, p, Prev))
  @assert codomain(inc) === c[i]
  c.boundary_cache[(i, p)] = inc
  return Im, inc
end

function homology(c::HyperComplex{ChainType}, p::Int, i::Tuple) where {ChainType <: SparseFPModule}
  if !isdefined(c, :homology_cache) 
    c.homology_cache = Dict{Tuple{Tuple, Int}, Map}()
  end
  if haskey(c.homology_cache, (i, p))
    pr = c.homology_cache[(i, p)]
    return codomain(pr), pr
  end

  H, pr = quo(kernel(c, p, i)[1], boundary(c, p, i)[1])
  c.homology_cache[(i, p)] = pr
  return H, pr
end

#function Base.show(io::IO, ::MIME"text/plain", c::AbsHyperComplex)
function Base.show(io::IO, c::AbsHyperComplex)
  if dim(c) == 1 
    return _print_standard_complex(io, c)
  end
  print(io, "Hyper complex of dimension $(dim(c))")
end

function Base.show(io::IO, c::HomComplex)
  io = pretty(io)
  println(io, "Hom complex from")
  print(io, Indent())
  print(io, "$(domain(c))")
  println(io, Dedent())
  println(io, "to")
  print(io, Indent())
  print(io, "$(codomain(c))")
end

function Base.show(io::IO, c::HCTensorProductComplex)
  io = pretty(io)
  println(io, "Tensor product of the complexes")
  print(io, Indent())
  for a in factors(c)
    println(io, "$(a)")
  end
  print(io, Dedent())
end

function Base.show(io::IO, c::LinearStrandComplex)
  io = pretty(io)
  io = IOContext(io, :compact => true)
  println(io, "Linear strand of")
  println(io, Indent(), "$(original_complex(c))")
  print(io, Dedent(), "of degree $(degree(c))")
end

function Base.show(io::IO, c::LinearStrandComplementComplex)
  io = pretty(io)
  io = IOContext(io, :compact => true)
  println(io, "Quotient of")
  println(io, Indent(), "$(original_complex(c))")
  print(io, Dedent(), "by its linear strand of degree $(degree(original_strand(c)))")
end

function _print_standard_complex(io::IO, c::AbsHyperComplex)
  io = IOContext(io, :compact => true)
  if has_lower_bound(c, 1) && has_upper_bound(c, 1)
    lb = lower_bound(c, 1)
    ub = upper_bound(c, 1)
    while !can_compute_index(c, lb)
      lb = lb+1
    end
    str = "$(c[lb])" 
    for i in lb+1:ub
      !can_compute_index(c, i) && break
      if direction(c, 1) == :chain
        str = str * " <-- " 
      else
        str = str * " --> " 
      end
      str = str * "$(c[i])"
    end
    println(io, str)
    return
  end

  if has_lower_bound(c, 1)
    lb = lower_bound(c, 1)
    while !can_compute_index(c, lb)
      lb = lb+1
    end
    str = "$(c[lb])"
    for i in lb+1:lb+3
      if !can_compute_index(c, (i,))
        println(io, str)
        return
      end

      if direction(c, 1) == :chain
        str = str * " <-- "
      else
        str = str * " --> "
      end
      str = str * "$(c[i])"
    end
    if direction(c, 1) == :chain
      str = str * " <-- ..."
    else
      str = str * " --> ..."
    end
    println(io, str)
    return
  end
  
  if has_upper_bound(c, 1)
    ub = upper_bound(c, 1)
    while !can_compute_index(c, ub)
      ub = ub-1
    end
    str = "$(c[ub])"
    for i in ub-1:-1:ub-3
      if !can_compute_index(c, (i,))
        println(io, str)
        return
      end

      if direction(c, 1) == :chain
        str = " <-- " * str
      else
        str = " --> " * str
      end
      str = "$(c[i])" * str
    end
    if direction(c, 1) == :chain
      str = "... <-- " * str
    else
      str = "... --> " * str
    end
    println(io, str)
    return
  end

  # If no bounds are known, we do not know where to start printing, so we give up.
  print(io, "Hyper complex of dimension $(dim(c)) with no known bounds")
  return 
end
 
function Base.show(io::IO, c::ZeroDimensionalComplex)
  print(io, "Zero-dimensional complex given by $(c[()])")
end

function Base.show(io::IO, c::SimpleFreeResolution)
  has_upper_bound(c) && return _free_show(io, c)
  return _print_standard_complex(io, c)
end

function _free_show(io::IO, C::AbsHyperComplex)
  # copied and adapted from src/Modules/UngradedModules/FreeResolutions.jl
  name_mod = String[]
  rank_mod = Int[]

  rng = upper_bound(C, 1):-1:lower_bound(C, 1)
  arr = ("<--", "--")

  R = Nemo.base_ring(C[first(rng)])
  R_name = AbstractAlgebra.get_name(R)
  if isnothing(R_name)
    R_name = "($R)"
  end
 
  for i=reverse(rng)
    M = C[i]
    M_name = AbstractAlgebra.get_name(M)
    if isnothing(M_name)
      M_name = "$R_name^$(rank(M))"
    end
    push!(name_mod, M_name)
    push!(rank_mod, rank(M))
  end

  io = IOContext(io, :compact => true)
  if C isa SimpleFreeResolution
    print(io, "Free resolution")
    print(io, " of ", C.M)
  end
  print(io, "\n")

  pos = 0
  pos_mod = Int[]
  
  for i=1:length(name_mod)
    print(io, name_mod[i])
    push!(pos_mod, pos)
    pos += length(name_mod[i])
    if i < length(name_mod)
      print(io, " ", arr[1], arr[2], " ")
      pos += length(arr[1]) + length(arr[2]) + 2
    end
  end

  print(io, "\n")
  len = 0
  for i=1:length(name_mod)
    if i>1
      print(io, " "^(pos_mod[i] - pos_mod[i-1]-len))
    end
    print(io, reverse(rng)[i])
    len = length("$(reverse(rng)[i])")
  end
end

### Koszul contraction
# Given a free `R`-module `F`, a morphism φ: F → R, and an element `v` in ⋀ ᵖ F, 
# compute the contraction φ(v) ∈ ⋀ ᵖ⁻¹ F.
# Note: For this internal method φ is only represented as a dense (column) vector.
# Warning: If the user provides their own parent, they need to make sure that things 
#          are compatible. As this is an internal function, no sanity checks are done.
function _contract(
    v::FreeModElem{T}, phi::Vector{T}; 
    parent::FreeMod{T}=begin
      success, F0, p = _is_exterior_power(Oscar.parent(v))
      @req success "parent is not an exterior power"
      exterior_power(F0, p-1)
    end
  ) where {T}
  success, F0, p = _is_exterior_power(Oscar.parent(v))
  @req success "parent is not an exterior power"
  @assert length(phi) == ngens(F0) "lengths are incompatible"
  result = zero(parent)
  n = ngens(F0)
  for (i, ind) in enumerate(OrderedMultiIndexSet(p, n))
    is_zero(v[i]) && continue
    for j in 1:p
      I = deleteat!(copy(indices(ind)), j)
      new_ind = OrderedMultiIndex(I, n)
      result = result + (-1)^j * v[i] * phi[ind[j]] * parent[linear_index(new_ind)]
    end
  end
  return result
end


