export cech_cohomologies

@doc raw"""
    structure_sheaf(v::NormalToricVarietyType)

Construct the structure sheaf of a normal toric variety.

# Examples
```jldoctest
julia> X = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> l1 = toric_line_bundle(X, [1])
Toric line bundle on a normal toric variety

julia> cech_cohomologies(l1)

julia> dP3 = del_pezzo_surface(NormalToricVariety, 3)
Normal toric variety

julia> l2 = toric_line_bundle(dP3, [-3,-2,-2,-2])
Toric line bundle on a normal toric variety

julia> cech_cohomologies(l2)
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

  # Prepare information, which we use to iterate over the Cech complex and identify the relevant lattice points
  RI = ray_indices(maximal_cones(X))
  ray_index_list = map(row -> findall(!iszero, collect(row)), eachrow(RI))

  # Now iterate over the Cech complex
  cech_complex_points = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}[]
  cech_complex_maps = Vector{Any}(undef, dim(X)+1)
  cech_complexes = Vector{FreeMod}(undef, dim(X)+2)
  comb_dict = Dict(); d_k = 0

  for k in 0:dim(X)
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
          if !any(x -> x < 0, sign_list[generating_ray_indices, :])
            push!(list_of_lattice_points, pt)
          end
        end
      end
      polyhedron_dict[combs[i]] = unique(list_of_lattice_points)
    end
    cech_complexes[k+1] = FreeMod(QQ, sum(length.(values(polyhedron_dict))))

    if k > 0
      d_k = matrix(QQ, zeros(QQ, rank(cech_complexes[k]), rank(cech_complexes[k+1])))
      counter = 1
      for i in 1:length(combs)
        sub_combs = collect(combinations(combs[i], k))
        for pt in polyhedron_dict[combs[i]]
          my_sign = 1
          for j in length(sub_combs):-1:1
            if haskey(comb_dict[sub_combs[j]], pt)
              d_k[comb_dict[sub_combs[j]][pt], counter] = my_sign
              my_sign *= -1
            end
          end
          counter += 1
        end
      end
      cech_complex_maps[k] = d_k
    end

    #Perhaps the most incriminating line
    comb_dict = Dict(combs[i] => Dict(polyhedron_dict[combs[i]][j] => j + sum(length.(values(Dict(key => polyhedron_dict[key] for key in combs[1:i-1]))))  for j in 1:length(polyhedron_dict[combs[i]])) for i in 1:length(combs))
    # return comb_dict
    push!(cech_complex_points, polyhedron_dict)
  end
  cech_complexes[dim(X)+2] = FreeMod(QQ, 0)
  cech_complex_maps[dim(X)+1] = matrix(QQ, zeros(QQ, rank(cech_complexes[dim(X)+1]), 0))

  return cech_complexes, cech_complex_maps, cech_complex_points
end
