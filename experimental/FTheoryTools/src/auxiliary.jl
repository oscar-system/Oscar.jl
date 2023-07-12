################################################################
# 1: Construct auxiliary base space
################################################################

function _auxiliary_base_space(auxiliary_base_variable_names::Vector{String}, auxiliary_base_grading::Matrix{Int64}, d::Int)
  
  # Find candidate base spaces
  candidates = normal_toric_varieties_from_glsm(matrix(ZZ, auxiliary_base_grading))
  @req length(candidates) > 0 "Could not find a full star triangulation"
  f = fan(candidates[1])
  @req dim(f) >= d "Cannot construct an auxiliary base space of the desired dimension"
  
  # Construct one base space of desired dimension
  fan_rays = matrix(ZZ, rays(f))
  fan_rays = [vec([Int(a) for a in fan_rays[k,:]]) for k in 1:nrows(fan_rays)]
  fan_max_cone_matrix = matrix(ZZ, cones(f))
  new_max_cones = Vector{Int}[]
  for k in 1:nrows(fan_max_cone_matrix)
    indices = findall(x -> x == 1, vec(fan_max_cone_matrix[k,:]))
    if length(indices) == d
      push!(new_max_cones, indices)
    end
  end
  auxiliary_base_space = normal_toric_variety(fan_rays, new_max_cones; non_redundant = true)
  
  # Set attributes for this base space
  set_coordinate_names(auxiliary_base_space, auxiliary_base_variable_names)
  G1 = free_abelian_group(ncols(auxiliary_base_grading))
  G2 = free_abelian_group(nrows(auxiliary_base_grading))
  grading_of_cox_ring = hom(G1, G2, transpose(matrix(ZZ, auxiliary_base_grading)))
  set_attribute!(auxiliary_base_space, :map_from_torusinvariant_weil_divisor_group_to_class_group, grading_of_cox_ring)
  set_attribute!(auxiliary_base_space, :class_group, G2)
  set_attribute!(auxiliary_base_space, :torusinvariant_weil_divisor_group, G1)
  
  return auxiliary_base_space
end

################################################################
# 2: Construct ambient space from given base
################################################################

_ambient_space_from_base(base::ToricCoveredScheme) = _ambient_space_from_base(underlying_toric_variety(base))

_ambient_space_from_base(base::ToricCoveredScheme, fiber_ambient_space::ToricCoveredScheme, D1::ToricDivisorClass, D2::ToricDivisorClass) = _ambient_space_from_base(underlying_toric_variety(base), underlying_toric_variety(fiber_ambient_space), D1, D2)

function _ambient_space_from_base(base::AbstractNormalToricVariety)
  fiber_ambient_space = weighted_projective_space(NormalToricVariety, [2,3,1])
  D1 = 2 * anticanonical_divisor_class(base)
  D2 = 3 * anticanonical_divisor_class(base)
  set_coordinate_names(fiber_ambient_space, ["x", "y", "z"])
  return _ambient_space(base, fiber_ambient_space, D1, D2)
end

function _ambient_space(base::AbstractNormalToricVariety, fiber_ambient_space::AbstractNormalToricVariety, D1::ToricDivisorClass, D2::ToricDivisorClass)
  
  # Consistency checks
  @req ((toric_variety(D1) === base) && (toric_variety(D2) === base)) "The divisors must belong to the base space"
  
  # Extract information about the toric base
  base_rays = matrix(ZZ, rays(base))
  base_cones = matrix(ZZ, ray_indices(maximal_cones(base)))
  
  # Extract information about the fiber ambient space
  fiber_rays = matrix(ZZ, rays(fiber_ambient_space))
  fiber_cones = matrix(ZZ, ray_indices(maximal_cones(fiber_ambient_space)))
  
  # Compute the u-matrix
  base_weights = transpose(vcat([elem.coeff for elem in cox_ring(base).d]))
  m1 = transpose(vcat([divisor_class(D1).coeff, divisor_class(D2).coeff]))
  m2 = fiber_rays[1:2,:]
  u_matrix = solve(base_weights,(-1)*m1*m2)
  
  # Form the rays of the toric ambient space
  new_base_rays = hcat(base_rays, u_matrix)
  new_fiber_rays = hcat(zero_matrix(ZZ, nrows(fiber_rays), ncols(base_rays)), fiber_rays)
  ambient_space_rays = vcat(new_base_rays, new_fiber_rays)
  ambient_space_rays = vcat([[k for k in ambient_space_rays[i,:]] for i in 1:nrows(ambient_space_rays)]...)
  
  # Construct the incidence matrix for the maximal cones of the ambient space
  ambient_space_max_cones = []
  for i in 1:nrows(base_cones)
    for j in 1:nrows(fiber_cones)
      push!(ambient_space_max_cones, [k for k in hcat([b for b in base_cones[i,:]], [c for c in fiber_cones[j,:]])])
    end
  end
  ambient_space_max_cones = IncidenceMatrix(vcat(ambient_space_max_cones...))
  
  # Construct the ambient space
  ambient_space = normal_toric_variety(polyhedral_fan(ambient_space_rays, ambient_space_max_cones; non_redundant = true))
  
  # Compute torusinvariant weil divisor group and the class group
  ambient_space_torusinvariant_weil_divisor_group = free_abelian_group(nrows(ambient_space_rays))
  ambient_space_class_group = free_abelian_group(nrows(base_weights) + rank(class_group(fiber_ambient_space)))
  
  # Construct grading matrix of ambient space
  ambient_space_grading = zero_matrix(ZZ,rank(ambient_space_torusinvariant_weil_divisor_group),rank(ambient_space_class_group))
  for i in 1:ncols(base_weights)
    for j in 1:nrows(base_weights)
      ambient_space_grading[i,j] = base_weights[j,i]
    end
  end
  fiber_weights = transpose(vcat([elem.coeff for elem in cox_ring(fiber_ambient_space).d]))
  for i in 1:ncols(fiber_weights)
    for j in 1:nrows(fiber_weights)
      ambient_space_grading[i + nrows(base_rays),j + nrows(base_weights)] = fiber_weights[j,i]
    end
  end
  for i in 1:ncols(divisor_class(D1).coeff)
    ambient_space_grading[1 + nrows(base_rays),i] = divisor_class(D1).coeff[i]
  end
  for i in 1:ncols(divisor_class(D2).coeff)
    ambient_space_grading[2 + nrows(base_rays),i] = divisor_class(D2).coeff[i]
  end
  
  # Construct the grading map for the ambient space
  ambient_space_grading = hom(ambient_space_torusinvariant_weil_divisor_group, ambient_space_class_group, ambient_space_grading)
  
  set_coordinate_names(ambient_space, vcat([string(k) for k in gens(cox_ring(base))], [string(k) for k in gens(cox_ring(fiber_ambient_space))]))
  set_attribute!(ambient_space, :map_from_torusinvariant_weil_divisor_group_to_class_group, ambient_space_grading)
  set_attribute!(ambient_space, :class_group, ambient_space_class_group)
  set_attribute!(ambient_space, :torusinvariant_weil_divisor_group, ambient_space_torusinvariant_weil_divisor_group)
  
  # Return the constructed space
  return ambient_space
  
end


################################################################
# 3: Construct the Weierstrass/Tate ambient space
################################################################

_weierstrass_ambient_space_from_base(base::ToricCoveredScheme) = _weierstrass_ambient_space_from_base(underlying_toric_variety(base))

function _weierstrass_ambient_space_from_base(base::AbstractNormalToricVariety)
  
  # Extract information about the toric base
  base_rays = matrix(ZZ, rays(base))
  base_cones = matrix(ZZ, ray_indices(maximal_cones(base)))
  
  # Construct the rays for the fibre ambient space
  xray = [0 for i in 1:ncols(base_rays)+2]
  yray = [0 for i in 1:ncols(base_rays)+2]
  zray = [0 for i in 1:ncols(base_rays)+2]
  xray[ncols(base_rays)+1] = 1
  yray[ncols(base_rays)+2] = 1
  zray[ncols(base_rays)+1] = -2
  zray[ncols(base_rays)+2] = -3
  
  # Construct the rays of the toric ambient space
  ambient_space_rays = hcat([r for r in base_rays], [-2 for i in 1:nrows(base_rays)], [-3 for i in 1:nrows(base_rays)])
  ambient_space_rays = vcat(ambient_space_rays, transpose(xray), transpose(yray), transpose(zray))
  
  # Construct the incidence matrix for the maximal cones of the ambient space
  ambient_space_max_cones = []
  for i in 1:nrows(base_cones)
    push!(ambient_space_max_cones, [k for k in hcat([b for b in base_cones[i,:]], [1 1 0])])
    push!(ambient_space_max_cones, [k for k in hcat([b for b in base_cones[i,:]], [1 0 1])])
    push!(ambient_space_max_cones, [k for k in hcat([b for b in base_cones[i,:]], [0 1 1])])
  end
  ambient_space_max_cones = IncidenceMatrix(vcat(ambient_space_max_cones...))
  
  # Construct and return the ambient space
  ambient_space = normal_toric_variety(polyhedral_fan(ambient_space_rays, ambient_space_max_cones; non_redundant = true))
  set_coordinate_names(ambient_space, vcat([string(k) for k in gens(cox_ring(base))], ["x", "y", "z"]))
  return ambient_space
  
end


################################################################
# 4: Construct the Weierstrass polynomial
################################################################

function _weierstrass_sections(base::AbstractNormalToricVariety)
  f = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
  g = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])
  return [f, g]
end

function _weierstrass_polynomial(base::AbstractNormalToricVariety, S::MPolyRing)
  (f, g) = _weierstrass_sections(base)
  return _weierstrass_polynomial(f, g, S)
end

function _weierstrass_polynomial(f::MPolyRingElem, g::MPolyRingElem, S::MPolyRing)
  x, y, z = gens(S)[ngens(S)-2:ngens(S)]
  ring_map = hom(parent(f), S, gens(S)[1:ngens(S)-3])
  return x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6
end


################################################################
# 4: Construct the Tate polynomial
################################################################

function _tate_sections(base::AbstractNormalToricVariety)
  a1 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base))])
  a2 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^2)])
  a3 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^3)])
  a4 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^4)])
  a6 = sum([rand(Int) * b for b in basis_of_global_sections(anticanonical_bundle(base)^6)])
  return [a1, a2, a3, a4, a6]
end

function _tate_polynomial(base::AbstractNormalToricVariety, S::MPolyRing)
  (a1, a2, a3, a4, a6) = _tate_sections(base)
  return _tate_polynomial([a1, a2, a3, a4, a6], S)
end

function _tate_polynomial(ais::Vector{<:MPolyRingElem}, S::MPolyRing)
  x, y, z = gens(S)[ngens(S)-2:ngens(S)]
  ring_map = hom(parent(ais[1]), S, gens(S)[1:ngens(S)-3])
  (a1, a2, a3, a4, a6) = [ring_map(k) for k in ais]
  return x^3 - y^2 - x*y*z*a1 + x^2*z^2*a2 - y*z^3*a3 + x*z^4*a4 + z^6*a6
end


################################################################
# 5: A base space for efficient testing
################################################################

@doc raw"""
    sample_toric_variety()

This method constructs a 3-dimensional toric variety, which we
use for efficient testing of the provided functionality.
"""
function sample_toric_variety()
  rays = [-1 -1 -1; -1 -1 0; -1 -1 1; -1 -1 2; -1 -1 3; -1 -1 4;
          -1 -1 5; -1 0 -1; -1 0 0; -1 0 1; -1 0 2; -1 0 3; -1 0 4;
          -1 1 -1; -1 1 0; -1 1 1; -1 1 2; -1 1 3; -1 2 -1; -1 2 0;
          -1 2 1; -1 2 2; -1 3 -1; -1 3 0; -1 3 1; -1 4 -1; -1 4 0;
          -1 5 -1; 0 -1 -1; 0 -1 0; 0 -1 1; 0 -1 2; 0 0 -1; 0 0 1;
          0 1 -1; 0 1 0; 0 2 -1; 1 -1 -1]
  cones = IncidenceMatrix([[36, 37, 38], [35, 37, 38], [34, 36, 38],
          [33, 35, 38], [32, 34, 38], [31, 32, 38], [30, 31, 38],[29, 33, 38],
          [29, 30, 38], [27, 28, 37], [26, 28, 37], [26, 27, 28], [25, 36, 37],
          [25, 27, 37], [24, 26, 27], [24, 25, 27], [23, 26, 37], [23, 24, 26],
          [22, 34, 36], [22, 25, 36], [21, 24, 25], [21, 22, 25], [20, 23, 24],
          [20, 21, 24], [19, 35, 37], [19, 23, 37], [19, 20, 23], [18, 32, 34],
          [18, 22, 34], [17, 21, 22], [17, 18, 22], [16, 20, 21], [16, 17, 21],
          [15, 19, 20], [15, 16, 20], [14, 33, 35], [14, 19, 35], [14, 15, 19],
          [13, 18, 32], [12, 17, 18], [12, 13, 18], [11, 16, 17], [11, 12, 17],
          [10, 15, 16], [10, 11, 16], [9, 14, 15], [9, 10, 15], [8, 29, 33],
          [8, 14, 33], [8, 9, 14], [7, 13, 32], [6, 12, 13], [6, 7, 32], [6, 7, 13],
          [5, 11, 12], [5, 6, 32], [5, 6, 12], [4, 31, 32], [4, 10, 11], [4, 5, 32],
          [4, 5, 11], [3, 30, 31], [3, 9, 10], [3, 4, 31], [3, 4, 10], [2, 29, 30],
          [2, 8, 9], [2, 3, 30], [2, 3, 9], [1, 8, 29], [1, 2, 29], [1, 2, 8]])
  return normal_toric_variety(polyhedral_fan(rays, cones))
end

@doc raw"""
    sample_toric_scheme()

This method constructs a 3-dimensional toric variety, which we
use for efficient testing of the provided functionality.
"""
sample_toric_scheme() = toric_covered_scheme(sample_toric_variety())


################################################################
# 6: Check if an ideal/subvariety is nontrivial
################################################################

_is_nontrivial(id::MPolyIdeal{T}, irr::MPolyIdeal{T}) where {T<:MPolyRingElem} = !is_one(id) && !is_one(saturation(id, irr))


################################################################
# 7: Compute singularity Kodaira type and refined Tate type
################################################################

_count_factors(poly::QQMPolyRingElem) = mapreduce(p -> p[end], +, absolute_primary_decomposition(ideal([poly])))

_string_from_factor_count(poly::QQMPolyRingElem, string_list::Vector{String}) = string_list[_count_factors(poly)]

function _kodaira_type(id::MPolyIdeal{T}, f::T, g::T, d::T, ords::Tuple{Int64, Int64, Int64}) where {T<:MPolyRingElem}
  f_ord = ords[1]
  g_ord = ords[2]
  d_ord = ords[3]
    
  if d_ord == 0
    kod_type = "I_0"
  elseif d_ord == 1 && f_ord == 0 && g_ord == 0
    kod_type = "I_1"
  elseif d_ord == 2 && g_ord == 1 && f_ord >= 1
    kod_type = "II"
  elseif d_ord == 3 && f_ord == 1 && g_ord >= 2
    kod_type = "III"
  elseif d_ord == 9 && f_ord == 3 && g_ord >= 5
    kod_type = "III^*"
  elseif d_ord == 10 && g_ord == 5 && f_ord >= 4
    kod_type = "II^*"
  elseif d_ord >= 12 && f_ord >= 4 && g_ord >= 6
    kod_type = "Non-minimal"
  else
    R = parent(f)
    S, (_psi, ) = polynomial_ring(QQ, ["_psi"; [string(v) for v in gens(R)]], cached = false)
    ring_map = hom(R, S, gens(S)[2:end])
    poly_f = ring_map(f)
    poly_g = ring_map(g)
    poly_d = ring_map(d)
    locus = ring_map(gens(id)[1])
    
    if f_ord == 0 && g_ord == 0
      monodromy_poly = _psi^2 + divexact(evaluate(9 * poly_g, [locus], [0]), evaluate(2 * poly_f, [locus], [0]))
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I_$d_ord", "Split I_$d_ord"])
    elseif d_ord == 4 && g_ord == 2 && f_ord >= 2
      monodromy_poly = _psi^2 - evaluate(divexact(poly_g, locus^2), [locus], [0])
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split IV", "Split IV"])
    elseif d_ord == 6 && f_ord >= 2 && g_ord >= 3
      monodromy_poly =  _psi^3 + _psi * evaluate(divexact(poly_f, locus^2), [locus], [0]) + evaluate(divexact(poly_g, locus^3), [locus], [0])
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I^*_0", "Semi-split I^*_0", "Split I^*_0"])
    elseif f_ord == 2 && g_ord == 3 && d_ord >= 7 && d_ord % 2 == 1
      monodromy_poly = _psi^2 + divexact(evaluate(divexact(poly_d, locus^d_ord) * divexact(2 * poly_f, locus^2)^3, [locus], [0]), 4 * evaluate(divexact(9 * poly_g, locus^3), [locus], [0])^3)
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I^*_$(d_ord - 6)", "Split I^*_$(d_ord - 6)"])
    elseif f_ord == 2 && g_ord == 3 && d_ord >= 8 && d_ord % 2 == 0
      monodromy_poly = _psi^2 + divexact(evaluate(divexact(poly_d, locus^d_ord) * divexact(2 * poly_f, locus^2)^2, [locus], [0]), evaluate(divexact(9 * poly_g, locus^3), [locus], [0])^2)
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I^*_$(d_ord - 6)", "Split I^*_$(d_ord - 6)"])
    elseif d_ord == 8 && g_ord == 4 && f_ord >= 3
      monodromy_poly = _psi^2 - evaluate(divexact(poly_g, locus^4), [locus], [0])
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split IV^*", "Split IV^*"])
    else
      kod_type = "Unrecognized"
    end
  end
  
  return kod_type
end


################################################################
# 8: Blowups
################################################################

function _blowup_global(id::MPolyIdeal{QQMPolyRingElem}, center::MPolyIdeal{QQMPolyRingElem}, irr::MPolyIdeal{QQMPolyRingElem}, sri::MPolyIdeal{QQMPolyRingElem}, lin::MPolyIdeal{QQMPolyRingElem}; index::Integer = 1)
  # @warn "The function _blowup_global is experimental; absence of bugs and proper results are not guaranteed"
  
  R = base_ring(id)
  center_size = ngens(center)
  
  # Various sanity checks
  @req (!is_zero(center)) "The blowup center must be non-empty"
  
  # @req is_subset(id, center) "The ideal of the blowup center must contain the ideal to be blown up"
  @req base_ring(irr) == R "The given irrelevant ideal must share the base ring of the ideal to be blown up"
  @req base_ring(sri) == R "The given Stanley–Reisner ideal must share the base ring of the ideal to be blown up"
  @req ngens(base_ring(lin)) == ngens(R) "The base ring of ideal of linear relations must have the same number of generators as the base ring of the ideal to be blown up"
  
  # Make sure the ideal of linear relations has the same base ring as the others
  lin = ideal(map(hom(base_ring(lin), R, collect(1:ngens(R))), gens(lin)))
  
  # Create new base ring for the blown up ideal and a map between the rings
  S, S_gens = polynomial_ring(QQ, [string("e_", index); [string("b_", index, "_", i) for i in 1:center_size]; [string(v) for v in gens(R)]], cached = false)
  (_e, new_coords...) = S_gens[1:center_size + 1]
  ring_map = hom(R, S, S_gens[center_size + 2:end])
  
  # Compute the total transform
  center_gens_S = map(ring_map, gens(center))
  total_transform = ideal(map(ring_map, gens(id))) + ideal([new_coords[i] * _e - center_gens_S[i] for i in 1:center_size])
  
  # Compute the exceptional locus and strict transform, checking for crepancy
  # Could alternatively replace _e with center_gens_S in the exceptional locus here, then take the
  # primary decomposition and remove parts whose saturation by the irrelevant ideal is the whole ring
  exceptional_ideal = total_transform + ideal([_e])
  strict_transform, exceptional_factor = saturation_with_index(total_transform, exceptional_ideal)
  crepant = (exceptional_factor == center_size - 1)
  
  # Compute the new irrelevant ideal, SRI, and ideal of linear relations
  # These may need to be changed after reintroducing e
  new_irr = ideal(map(ring_map, gens(irr))) * ideal(new_coords)
  new_sri = ideal(map(ring_map, gens(sri))) + ideal([prod(new_coords)])
  new_lin = ideal(map(ring_map, gens(lin))) + ideal([g - new_coords[end] for g in new_coords[1:end - 1]])
  
  return total_transform, strict_transform, exceptional_ideal, crepant, new_irr, new_sri, new_lin, S, S_gens, ring_map
end
_blowup_global(id::T, center::T, irr::T, sri::T, lin::MPolyIdeal{QQMPolyRingElem}; index::Integer = 1) where {T<:MPolyIdeal{<:MPolyRingElem}} = _blowup_global(ideal(map(g -> g.f, gens(id))), ideal(map(g -> g.f, gens(center))), ideal(map(g -> g.f, gens(irr))), ideal(map(g -> g.f, gens(sri))), lin, index = index)


function _blowup_global_sequence(id::MPolyIdeal{QQMPolyRingElem}, centers::Vector{<:Vector{<:Integer}}, irr::MPolyIdeal{QQMPolyRingElem}, sri::MPolyIdeal{QQMPolyRingElem}, lin::MPolyIdeal{QQMPolyRingElem}; index::Integer = 1)
  # @warn "The function _blowup_global_sequence is experimental; absence of bugs and proper results are not guaranteed"
  
  (cur_strict_transform, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, cur_index) = (id, irr, sri, lin, base_ring(id), gens(base_ring((id))), index)
  crepant = true
  ring_map = hom(cur_S, cur_S, cur_S_gens) # Identity map
  
  exceptionals = MPolyIdeal{<:MPolyRingElem}[]
  for center in centers
    @req all(ind -> 1 <= ind <= length(cur_S_gens), center) "The given indices for the center generators are out of bounds"
    
    (_, cur_strict_transform, cur_ex, cur_crep, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, cur_ring_map) = _blowup_global(cur_strict_transform, ideal(map(ind -> cur_S_gens[ind], center)), cur_irr, cur_sri, cur_lin, index = cur_index)
    
    map!(cur_ring_map, exceptionals, exceptionals)
    push!(exceptionals, cur_ex)
    
    crepant = crepant && cur_crep
    
    ring_map = compose(ring_map, cur_ring_map)
    
    cur_index += 1
  end
  
  return cur_strict_transform, exceptionals, crepant, cur_irr, cur_sri, cur_lin, cur_S, cur_S_gens, ring_map
end
_blowup_global_sequence(id::T, centers::Vector{<:Vector{<:Integer}}, irr::T, sri::T, lin::MPolyIdeal{QQMPolyRingElem}; index::Integer = 1) where {T<:MPolyIdeal{<:MPolyRingElem}} = _blowup_global_sequence(ideal(map(g -> g.f, gens(id))), centers, ideal(map(g -> g.f, gens(irr))), ideal(map(g -> g.f, gens(sri))), lin, index = index)
