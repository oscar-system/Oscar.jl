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
  comb_dict = Dict(); d_k = 0; d_k_extended = 0; old_length = 0

  for k in 0:dim(X)
    polyhedron_dict = Dict{Vector{Int64}, Vector{PointVector{ZZRingElem}}}()
    combs = combinations(n_maximal_cones(X), k+1)
    if k > 0
      d_k = matrix(QQ, zeros(QQ, old_length, length(combs)))
      d_k_extended = matrix(QQ, zeros(QQ, 0, length(combs)))
      d_k_final = matrix(QQ, zeros(QQ, rank(cech_complexes[k]), 0))
    end
    for i in 1:length(combs)
      if k > 0 
        sub_combs = combinations(combs[i], k)
        d_k[[comb_dict[sub_combs[j]] for j in length(sub_combs):-1:1], i] = [(-1).^(t+mod(length(sub_combs), 2)) for t in length(sub_combs):-1:1]#This perhaps should have multiple copies of each +-1 for each lattice point corresponding to sub_combs(j)?
      end
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
    if k > 0
      for j in 1:old_length
        sub_comb = findfirst(x[2]==j for x in comb_dict)
        repeat_count = length(cech_complex_points[k][sub_comb])
        d_k_extended = vcat(d_k_extended, matrix(QQ,vcat([d_k[j,:] for _ in 1:repeat_count])))
      end
      for i in 1:length(combs)
        repeat_count = length(polyhedron_dict[combs[i]])
        d_k_final = hcat(d_k_final, matrix(QQ, hcat([d_k_extended[:,i] for _ in 1:repeat_count]...)))
      end
      cech_complex_maps[k] = d_k_final
    end
    comb_dict = Dict(combs[i] => i for i in 1:length(combs))
    cech_complexes[k+1] = FreeMod(QQ, sum(length.(values(polyhedron_dict))))
    push!(cech_complex_points, polyhedron_dict)
    old_length = length(combs)
  end
  cech_complexes[dim(X)+2] = FreeMod(QQ, 0)
  cech_complex_maps[dim(X)+1] = matrix(QQ, zeros(QQ, rank(cech_complexes[dim(X)+1]), 0))

  return cech_complexes, cech_complex_maps, cech_complex_points
end
