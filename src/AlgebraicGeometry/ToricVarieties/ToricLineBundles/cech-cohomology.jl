#The following function is called, when algorithm="chambers" is specified for all_cohomologies
function _all_cohomologies_via_cech(tl::ToricLineBundle)
  our_maps = transpose.(matrix.(_toric_cech_complex(tl)))
  X = toric_variety(tl)
  alternating_sum = abs(sum(binomial(n_maximal_cones(X), k) * (-1)^k for k in dim(X)+2:n_maximal_cones(X); init = 0))
  alternating_sum *= Int(ncols(our_maps[end]) // binomial(n_maximal_cones(X), dim(X)+1))
  @req ncols(our_maps[end]) >= alternating_sum "Inconsistency encountered"
  rks_maps = push!(rank.(our_maps), alternating_sum)
  return ZZRingElem.([nrows(our_maps[1]) - rks_maps[1]; [nrows(our_maps[k]) - rks_maps[k] - rks_maps[k-1] for k in 2:dim(X)]; ncols(our_maps[dim(X)]) - rks_maps[dim(X) + 1] - rks_maps[dim(X)]])
end

function _toric_cech_complex(tl::ToricLineBundle)
  # Extract essential information
  X = toric_variety(tl)
  @req is_complete(X) "Toric variety must be complete!"

  # Compute the hyperplane arrangement
  a_plane = matrix(ZZ, n_rays(X), 1, coefficients(toric_divisor(tl)))
  sc = cone_from_inequalities(matrix(ZZ, [[-1; zeros(Int, dim(X))]]))
  H = Polymake.fan.HyperplaneArrangement( HYPERPLANES = [a_plane matrix(ZZ, rays(X))], SUPPORT=sc.pm_cone)

  # Classify the regions obtained from the hyperlane arrangement
  pff = Polymake.fan.PolyhedralComplex(POINTS=H.CHAMBER_DECOMPOSITION.RAYS, INPUT_CONES=H.CHAMBER_DECOMPOSITION.MAXIMAL_CONES)
  max_polys = maximal_polyhedra(polyhedral_complex(pff))
  other_polys = filter(is_bounded, vcat([polyhedra_of_dim(polyhedral_complex(pff), i) for i in 0:dim(polyhedral_complex(pff))-1]...))
  chamber_signs = 2*matrix(ZZ, H.CHAMBER_SIGNATURES).-1
  sign_of_chamber = Dict(chamber_signs[i,:] => interior_lattice_points(p) for (i, p) in enumerate(max_polys) if is_bounded(p))

  # For the other regions, we don't get their signs for free and have to calculate them
  for p in other_polys
    if dim(p) == 0
      sign_list = ZZ.(sign.(matrix(QQ, rays(X)) * vertices(p)[1] + a_plane)[:,1])
    else
      sign_list = ZZ.(sign.(matrix(QQ, rays(X)) * 1//n_vertices(p) * sum(vertices(p)) + a_plane)[:,1])
    end
    sign_of_chamber[sign_list] = interior_lattice_points(p)
  end

  # Prepare information, which we use to iterate over the Cech complex and identify the relevant lattice points
  RI = ray_indices(maximal_cones(X))
  ray_index_list = map(row -> findall(!iszero, collect(row)), eachrow(RI))

  # Now iterate over the Cech complex
  cech_complex_points = Vector{Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}}(undef, dim(X) + 1)
  cech_complex_maps = Vector{Any}(undef, dim(X))
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
      for (signs, pts) in sign_of_chamber
        all(x -> x >= 0, signs[generating_ray_indices]) && append!(list_of_lattice_points, pts)
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
      our_sparse_matrix = sparse_matrix(ZZ, 0, previous_number_of_generators)
      col_idx = 1
      for comb in combs
        pts = polyhedron_dict[comb]
        sub_combs = collect(combinations(comb, k))  # returns lex-sorted k-subsets
        for pt in pts
          position_vector = Vector{Int}()
          value_vector = Vector{ZZRingElem}()
          for j in 1:length(sub_combs)
            sub_comb = sub_combs[j]
            if haskey(comb_dict, sub_comb) && haskey(comb_dict[sub_comb], pt)
              row_idx = comb_dict[sub_comb][pt]
              sign = (-1)^(j + 1)  # Čech sign convention (1-based indexing)
              push!(position_vector, row_idx)
              push!(value_vector, sign)
            end
          end
          push!(our_sparse_matrix, sparse_row(ZZ, position_vector, value_vector))
          col_idx += 1
        end
      end
      cech_complex_maps[k] = our_sparse_matrix
    end
    cech_complex_points[k+1] = polyhedron_dict
    previous_number_of_generators = sum(length.(values(polyhedron_dict)))

  end
  return cech_complex_maps
end
