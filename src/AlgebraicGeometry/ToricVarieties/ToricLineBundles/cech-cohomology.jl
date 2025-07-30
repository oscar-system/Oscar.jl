export _toric_cech_complex, all_cohomologies_via_cech, tester


@doc raw"""
    tester(tl::ToricLineBundle)

Compare line bundle cohomology from Cech cohomology with cohomCalg result.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> [tester(toric_line_bundle(X, [1])) for k in -10:10]

julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> l2 = toric_line_bundle(dP3, [3,3,3,3])
Toric line bundle on a normal toric variety


```
"""
function tester(tl::ToricLineBundle)
  return all_cohomologies_via_cech(tl) == all_cohomologies(tl)
end


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

julia> l2 = toric_line_bundle(dP3, [3,3,3,3])
Toric line bundle on a normal toric variety

julia> all_cohomologies_via_cech(l2)
```
"""
function all_cohomologies_via_cech(tl::ToricLineBundle)
  our_maps = _toric_cech_complex(tl)
  X = toric_variety(tl)
  alternating_sum = abs(sum(binomial(n_maximal_cones(X), k) * (-1)^k for k in dim(X)+2:n_maximal_cones(X)))
  alternating_sum *= Int(ncols(our_maps[end]) // binomial(n_maximal_cones(X), dim(X)+1))
  @req ncols(our_maps[end]) >= alternating_sum "Inconsistency encountered"
  rks_maps = push!(rank.(our_maps), alternating_sum)
  return [nrows(our_maps[1]) - rks_maps[1]; [nrows(our_maps[k]) - rks_maps[k] - rks_maps[k-1] for k in 2:dim(X)]; ncols(our_maps[dim(X)]) - rks_maps[dim(X) + 1] - rks_maps[dim(X)]]
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

  # Compute the maximal chambers
  pff = Polymake.fan.PolyhedralComplex(POINTS=H.CHAMBER_DECOMPOSITION.RAYS, INPUT_CONES=H.CHAMBER_DECOMPOSITION.MAXIMAL_CONES)
  bounded_max_polys = filter(is_bounded, maximal_polyhedra(polyhedral_complex(pff)))

  # Find all lattice points in bounded_max_polys, excluding duplicates
  list_of_contributing_lattice_points = unique(vcat([lattice_points(p) for p in bounded_max_polys]...))
  
  # Prepare information, which we use to iterate over the Cech complex and identify the relevant lattice points
  RI = ray_indices(maximal_cones(X))
  ray_index_list = map(row -> findall(!iszero, collect(row)), eachrow(RI))

  # Now iterate over the Cech complex
  cech_complex_points = Vector{Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}}(undef, dim(X) + 1)
  cech_complex_maps = Vector{QQMatrix}(undef, dim(X))
  comb_dict = Dict{Combination{Int64}, Dict{PointVector{ZZRingElem}, Int64}}()
  d_k = 0
  previous_number_of_generators = 0
  for k in 0:dim(X)

    # Find the contributing lattice points
    polyhedron_dict = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}()
    combs = collect(combinations(n_maximal_cones(X), k+1))
    for i in 1:length(combs)
      list_of_lattice_points = PointVector{ZZRingElem}[]
      generating_ray_indices = reduce(intersect, ray_index_list[combs[i]])
      for pt in list_of_contributing_lattice_points
        signs = matrix(QQ, rays(X)) * pt + a_plane
        all(x -> x >= 0, signs[generating_ray_indices, :]) && push!(list_of_lattice_points, pt)
      end
      polyhedron_dict[combs[i]] = list_of_lattice_points
    end

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
    if k > 0
      n_rows = previous_number_of_generators
      n_cols = sum(length.(values(polyhedron_dict)))
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
    
    cech_complex_points[k+1] = polyhedron_dict
    previous_number_of_generators = sum(length.(values(polyhedron_dict)))

  end

  # Return the result
  return cech_complex_maps

end
