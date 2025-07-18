export cech_cohomologies

#=
=#
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
  cech_complex = Dict{Combination{Int64}, Vector{PointVector{ZZRingElem}}}[]
  for k in 0:dim(X)
    polyhedron_dict = Dict{Combination{Int64}, Vector{PointVector{ZZRingElem}}}()
    for l in combinations(n_maximal_cones(X), k+1)
      list_of_lattice_points = PointVector{ZZRingElem}[]
      generating_ray_indices = reduce(intersect, ray_index_list[l])
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
      polyhedron_dict[l] = unique(list_of_lattice_points)
    end
    push!(cech_complex, polyhedron_dict)
  end

  return cech_complex
end
