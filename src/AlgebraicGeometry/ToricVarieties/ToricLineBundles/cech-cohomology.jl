export _toric_cech_complex, all_cohomologies_via_cech

@doc raw"""
    all_cohomologies_via_cech(tl::ToricLineBundle)

Compute line bundle cohomology via Cech cohomology.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l1 = toric_line_bundle(X, [1])
Toric line bundle on a normal toric variety

julia> all_cohomologies_via_cech(l1)
5-element Vector{ZZRingElem}:
 3
 0
 0

julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> l2 = toric_line_bundle(dP3, [-3,-2,-2,-2])
Toric line bundle on a normal toric variety

julia> all_cohomologies_via_cech(l2)
```
"""
function all_cohomologies_via_cech(tl::ToricLineBundle)
  our_free_modules, our_mapping_matrices = _toric_cech_complex(tl)
  rks_modules = [rank(m) for m in our_free_modules]
  X = toric_variety(tl)
  @req length(rks_modules) >= dim(X) + 1 "Inconsistency encountered"
  alternating_sum = 0
  for k in dim(X) + 2:length(rks_modules)
    alternating_sum += rks_modules[k] * (-1)^k
  end
  if alternating_sum < 0
    alternating_sum = alternating_sum * (-1)
  end
  dim_kernel_last_map = rks_modules[dim(X) + 1] - alternating_sum

  @req dim_kernel_last_map >= 0 "Inconsistency encountered"
  #return [dim_kernel_last_map, rks]

  # Drop all mapping matrices but the ones at position 1 to dim(X)
  rks_maps = [rank(our_mapping_matrices[k]) for k in 1:dim(X)]
  push!(rks_maps, alternating_sum)
  #return rks_modules[1:dim(X) + 1], rks_maps

  # Compute cohomologies
  coho = Int[]
  push!(coho, rks_modules[1] - rks_maps[1])
  for k in 2:dim(X)+1
    push!(coho, rks_modules[k] - rks_maps[k] - rks_maps[k-1])
  end
  #return rks_modules[1:dim(X) + 1], rks_maps, coho
  return coho
  #return rks
end

@doc raw"""
    _toric_cech_complex(tl::ToricLineBundle)

Construct the toric Cech complex.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l1 = toric_line_bundle(X, [1])
Toric line bundle on a normal toric variety

julia> cech_cohomologies(l1);

julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> l2 = toric_line_bundle(dP3, [-3,-2,-2,-2])
Toric line bundle on a normal toric variety

julia> cech_cohomologies(l2);

julia> v = dP3 * dP3
Normal toric variety

julia> l3 = toric_line_bundle(v, [-i for i in 1:8])
Toric line bundle on a normal toric variety

julia> all_cohomologies(l3)
5-element Vector{ZZRingElem}:
 0
 0
 196
 119
 0
```
"""
function _toric_cech_complex(tl::ToricLineBundle)

  # Extract essential information
  X = toric_variety(tl)
  @req is_complete(X) "Toric variety must be complete!"

  # Compute the hyperplane arrangement
  a_plane = matrix(ZZ, n_rays(X), 1, coefficients(toric_divisor(tl)))

  sc = cone_from_inequalities(matrix(ZZ, [[-1; zeros(Int, dim(X))]]))
  H = Polymake.fan.HyperplaneArrangement( HYPERPLANES = [a_plane matrix(ZZ, rays(X))], SUPPORT=sc.pm_cone)
  chamber_signs = matrix(ZZ,H.CHAMBER_SIGNATURES)

  # Compute the maximal chambers
  pff = Polymake.fan.PolyhedralComplex(POINTS=H.CHAMBER_DECOMPOSITION.RAYS, INPUT_CONES=H.CHAMBER_DECOMPOSITION.MAXIMAL_CONES)
  max_polys = maximal_polyhedra(polyhedral_complex(pff))

  # Keep only bounded chambers, but remember their sign
  bounded_max_polys = Dict()
  for i in 1:length(max_polys)
    if is_bounded(max_polys[i])
      bounded_max_polys[chamber_signs[i,:]] = max_polys[i]
    end
  end

  # Prepare information, which we use to iterate over the Cech complex and identify the relevant lattice points
  RI = ray_indices(maximal_cones(X))
  ray_index_list = map(row -> findall(!iszero, collect(row)), eachrow(RI))

  # Length
  cech_length = n_maximal_cones(X)
  #cech_length = dim(X)

  # Now iterate over the Cech complex
  cech_complex_points = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}[]
  #cech_complex_maps = Vector{Any}(undef, cech_length)
  cech_complex_maps = Vector{Any}(undef, dim(X))
  cech_complexes = Vector{FreeMod}(undef, cech_length+1)
  comb_dict = Dict(); d_k = 0

  for k in 0:cech_length
    polyhedron_dict = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}()
    combs = collect(combinations(n_maximal_cones(X), k+1))
    for i in 1:length(combs)
      if k == dim(X)
        list_of_lattice_points = unique(vcat(values(cech_complex_points[k])...))
      else
        list_of_lattice_points = PointVector{ZZRingElem}[]
        generating_ray_indices = reduce(intersect, ray_index_list[combs[i]])
        remaining_points = Dict{PointVector{ZZRingElem}, Any}()

        for (sign_list, p) in bounded_max_polys

          # Verify if p's lattice points contribute
          any(x -> x < 0, sign_list[generating_ray_indices]) && continue
          
          # Interior points always contribute, but boundary points will have to be processed later
          append!(list_of_lattice_points, interior_lattice_points(p))
          for pt in boundary_lattice_points(p)
            haskey(remaining_points, pt) && continue
            remaining_points[pt] = sign_list
          end
        end

        #Process remaining points
        for (pt, sign_list) in remaining_points
          pt_sign = matrix(QQ, rays(X)) * pt + a_plane
          if all(j -> (sign_list[j] < 0 && pt_sign[j] < 0) || (sign_list[j] >= 0 && pt_sign[j] >= 0), 1:length(sign_list))
            push!(list_of_lattice_points, pt)
          elseif count(==(0), pt_sign) > dim(X) && all(x -> x >= 0, pt_sign[generating_ray_indices, :])
            push!(list_of_lattice_points, pt)
          end
        end
      end

      polyhedron_dict[combs[i]] = list_of_lattice_points
    end
    cech_complexes[k+1] = FreeMod(QQ, sum(length.(values(polyhedron_dict))))

    # Initialize comb_dict before using it
    offset = 0
    for i in 1:length(combs)
      this_comb = combs[i]
      pts = polyhedron_dict[this_comb]
      inner_dict = Dict(pt => j + offset for (j, pt) in enumerate(pts))
      comb_dict[this_comb] = inner_dict
      offset += length(pts)
    end
    
    # Compute Cech differential maps
    if k > 0 && k <= dim(X)
      n_rows = rank(cech_complexes[k])
      n_cols = rank(cech_complexes[k+1])
      d_k = zero_matrix(QQ, n_rows, n_cols)    
      col_idx = 1
      for comb in combs
        pts = polyhedron_dict[comb]
        sub_combs = collect(combinations(comb, k))  # returns lex-sorted k-subsets
        for pt in pts
          for j in 1:length(sub_combs)
            sub_comb = sub_combs[j]
            if haskey(comb_dict, sub_comb) && haskey(comb_dict[sub_comb], pt)
              row_idx = comb_dict[sub_comb][pt]
              sign = (-1)^(j + 1)  # Čech sign convention (1-based indexing)
              d_k[row_idx, col_idx] = sign
            end
          end
          col_idx += 1
        end
      end
      cech_complex_maps[k] = d_k
    end

    push!(cech_complex_points, polyhedron_dict)
  end
  cech_complexes[cech_length+1] = FreeMod(QQ, 0)
  #cech_complex_maps[cech_length] = matrix(QQ, zeros(QQ, rank(cech_complexes[cech_length+1]), 0))

  #return cech_complexes, cech_complex_maps, cech_complex_points
  return cech_complexes, cech_complex_maps

end
