### Generic functionality
function total_complex(dc::DoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  is_bounded(dc) || error("computation of total complexes is only implemented for bounded double complexes")
  vertical_typ(dc) == horizontal_typ(dc) || error("horizontal and vertical typs must be the same")
  if vertical_typ(dc) == :chain
    return _total_chain_complex(dc)
  else
    return _total_cochain_complex(dc)
  end
end

function _total_chain_complex(dc::DoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  r1 = horizontal_range(dc)
  r2 = vertical_range(dc)
  r_tot = (first(r1) + first(r2)):-1:(last(r1) + last(r2))
  new_chains = ChainType[]
  new_maps = MorphismType[]

  # First round for initialization
  index_pairs = [I for I in Iterators.product(r1, r2) if sum(I) == first(r_tot)]
  summands = [dc[I] for I in index_pairs]
  new_chain, inc, pr = direct_sum(summands...; task=:both)
  last_inc = inc
  last_pr = pr
  last_chain = new_chain
  last_index_pairs = index_pairs
  push!(new_chains, new_chain)
  for k in r_tot
    k == first(r_tot) && continue
    index_pairs = [I for I in Iterators.product(r1, r2) if sum(I) == k]
    summands = [dc[I] for I in index_pairs]
    new_chain, inc, pr = direct_sum(summands...; task=:both)
    @assert all(f->domain(f) === new_chain, pr)
    @assert all(f->codomain(f) === new_chain, inc)
    push!(new_chains, new_chain)
    # assemble the map
    boundary_map = hom(last_chain, new_chain, elem_type(new_chain)[zero(new_chain) for i in 1:ngens(last_chain)])
    for (l, I) in enumerate(last_index_pairs)
      p = last_pr[l]
      @assert domain(p) === last_chain

      if I[2] > vertical_lower_bound(dc)
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

function _total_cochain_complex(dc::DoubleComplexOfMorphisms{ChainType, MorphismType}) where {ChainType, MorphismType}
  error("total complex of double cochain complexes currently not implemented")
end

### Missing functionality for complexes
typ(C::ComplexOfMorphisms) = C.typ
is_complete(C::ComplexOfMorphisms) = C.complete

