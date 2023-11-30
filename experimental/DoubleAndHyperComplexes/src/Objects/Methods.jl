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

function kernel(c::HyperComplex{ChainType}, p::Int, i::Tuple) where {ChainType <: ModuleFP}
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

function boundary(c::HyperComplex{ChainType}, p::Int, i::Tuple) where {ChainType <: ModuleFP}
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

function homology(c::HyperComplex{ChainType}, p::Int, i::Tuple) where {ChainType <: ModuleFP}
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
