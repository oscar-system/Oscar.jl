export cech_cohomologies, cech_cohomologies2

@doc raw"""
    structure_sheaf(v::NormalToricVarietyType)

Construct the structure sheaf of a normal toric variety.

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
function cech_cohomologies(tl::ToricLineBundle)

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

  # Length
  #cech_length = n_maximal_cones(X)
  cech_length = dim(X)

  # Now iterate over the Cech complex
  cech_complex_points = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}[]
  cech_complex_maps = Vector{Any}(undef, cech_length+1)
  cech_complexes = Vector{FreeMod}(undef, cech_length+1)
  comb_dict = Dict(); d_k = 0
  for k in 0:cech_length

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

    # Create the modules in the Cech complex
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
    if k > 0
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

  # Append one final map to zero. Why is this needed?
  cech_complexes[cech_length+1] = FreeMod(QQ, 0)
  cech_complex_maps[cech_length+1] = matrix(QQ, zeros(QQ, rank(cech_complexes[cech_length+1]), 0))

  # Return the result
  return cech_complexes, cech_complex_maps, cech_complex_points
end

function cech_cohomologies2(tl::ToricLineBundle)

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
  
  # Prepare information, which we use to iterate over the Cech complex and identify the relevant lattice points
  RI = ray_indices(maximal_cones(X))
  ray_index_list = map(row -> findall(!iszero, collect(row)), eachrow(RI))

  # Length
  #cech_length = n_maximal_cones(X)
  cech_length = dim(X)

  # Now iterate over the Cech complex
  cech_complex_points = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}[]
  cech_complex_maps = Vector{Any}(undef, cech_length)
  cech_complexes = Vector{FreeMod}(undef, cech_length+1)
  comb_dict = Dict(); d_k = 0

  #for k in 0:cech_length
  for k in 0:1
    polyhedron_dict = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}()
    combs = collect(combinations(n_maximal_cones(X), k+1))
    for i in 1:length(combs)
      list_of_lattice_points = PointVector{ZZRingElem}[]
      generating_ray_indices = reduce(intersect, ray_index_list[combs[i]])
      for p in bounded_max_polys
        
        # Verify if p's lattice points contribute
        sign_list = matrix(QQ, rays(X)) * 1//n_vertices(p) * sum(vertices(p)) + a_plane
        any(x -> x < 0, sign_list[generating_ray_indices, :]) && continue
        
        # Find out which lattice points contribute
        append!(list_of_lattice_points, interior_lattice_points(p))
        for pt in boundary_lattice_points(p)
          pt_sign = matrix(QQ, rays(X)) * pt + a_plane
          if all(i -> (sign_list[i] < 0 && pt_sign[i] < 0) || (sign_list[i] >= 0 && pt_sign[i] >= 0), 1:length(sign_list))
            push!(list_of_lattice_points, pt)
          end
        end
      end

      #We are also interested in vertices that are 'almost regions', i.e., if a hyperplane were moved slightly then a new region could show up here 
      for pt in vertices(polyhedral_complex(pff))
        pt =  ZZRingElem.(Int.(pt))
        pt_sign = matrix(QQ, rays(X)) * pt + a_plane
        count(==(0), pt_sign) <= dim(X) && continue
        any(x -> x < 0, pt_sign[generating_ray_indices, :]) && continue
        push!(list_of_lattice_points, pt)
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
    if k > 0
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
  cech_complex_maps[cech_length] = matrix(QQ, zeros(QQ, rank(cech_complexes[cech_length+1]), 0))

  return cech_complexes, cech_complex_maps, cech_complex_points
end
