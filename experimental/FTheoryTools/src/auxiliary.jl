################################################################
# 1: Construct ambient space from given base
################################################################

function _ambient_space(base::NormalToricVariety, fiber_amb_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass})
  @req all(D -> toric_variety(D) === base, fiber_twist_divisor_classes) "The divisors must belong to the (same) base space"
  
  # Extract information about the toric base
  b_rays = matrix(ZZ, rays(base))
  b_cones = matrix(ZZ, ray_indices(maximal_cones(base)))
  b_grades = reduce(vcat, [elem.coeff for elem in cox_ring(base).d])
  b_var_names = [string(k) for k in gens(cox_ring(base))]
  
  # Extract information about the fiber ambient space
  f_rays = matrix(ZZ, rays(fiber_amb_space))
  f_cones = matrix(ZZ, ray_indices(maximal_cones(fiber_amb_space)))
  f_grades = reduce(vcat, [elem.coeff for elem in cox_ring(fiber_amb_space).d])
  f_var_names = [string(k) for k in gens(cox_ring(fiber_amb_space))]
  
  # Extract coefficients of divisors D1, D2 and compute u_matrix
  fiber_twist_divisor_classes_coeffs = [divisor_class(D).coeff for D in fiber_twist_divisor_classes]
  m1 = reduce(vcat, fiber_twist_divisor_classes_coeffs)
  m2 = transpose(f_rays)
  u_matrix = solve(b_grades, (-1)*m2*m1; side = :left)
  
  # Form toric ambient space
  a_rays = zero_matrix(ZZ, nrows(b_rays) + nrows(f_rays), ncols(b_rays) + ncols(f_rays))
  a_rays[1:nrows(b_rays), 1:ncols(b_rays)] = b_rays
  a_rays[1:nrows(b_rays), 1+ncols(b_rays):ncols(a_rays)] = transpose(u_matrix)
  a_rays[1+nrows(b_rays):nrows(a_rays), 1+ncols(b_rays):ncols(a_rays)] = f_rays
  a_cones = [hcat([b for b in b_cones[i:i,:]], [c for c in f_cones[j:j,:]]) for i in 1:nrows(b_cones), j in 1:nrows(f_cones)]
  a_space = normal_toric_variety(IncidenceMatrix(vcat(a_cones...)), a_rays; non_redundant = true)
  set_coordinate_names(a_space, vcat(b_var_names, f_var_names))
  
  # Compute divisor group and the class group of a_space
  a_space_divisor_group = free_abelian_group(nrows(a_rays))
  a_space_class_group = free_abelian_group(ncols(b_grades) + torsion_free_rank(class_group(fiber_amb_space)))
  
  # Compute grading of Cox ring of a_space
  a_space_grading = zero_matrix(ZZ, torsion_free_rank(a_space_divisor_group), torsion_free_rank(a_space_class_group))
  a_space_grading[1:nrows(b_grades), 1:ncols(b_grades)] = b_grades
  a_space_grading[1+nrows(b_rays):nrows(b_rays) + nrows(f_grades), 1+ncols(b_grades):ncols(b_grades) + ncols(f_grades)] = f_grades
  for k in 1:length(fiber_twist_divisor_classes_coeffs)
    a_space_grading[k+nrows(b_rays), 1:ncols(fiber_twist_divisor_classes_coeffs[k])] = fiber_twist_divisor_classes_coeffs[k]
  end
  
  # Set important attributes of a_space and return it
  a_space_grading = hom(a_space_divisor_group, a_space_class_group, a_space_grading)
  set_attribute!(a_space, :map_from_torusinvariant_weil_divisor_group_to_class_group, a_space_grading)
  set_attribute!(a_space, :class_group, a_space_class_group)
  set_attribute!(a_space, :torusinvariant_weil_divisor_group, a_space_divisor_group)
  return a_space
end


################################################################
# 2: Construct the Weierstrass polynomial
################################################################

function _weierstrass_sections(base::NormalToricVariety)
  return [generic_section(anticanonical_bundle(base)^4), generic_section(anticanonical_bundle(base)^6)]
end

function _weierstrass_polynomial(base::NormalToricVariety, S::MPolyRing)
  (f, g) = _weierstrass_sections(base)
  return _weierstrass_polynomial(f, g, S)
end

function _weierstrass_polynomial(f::MPolyRingElem, g::MPolyRingElem, S::MPolyRing)
  x, y, z = gens(S)[ngens(S)-2:ngens(S)]
  ring_map = hom(parent(f), S, gens(S)[1:ngens(S)-3])
  return x^3 - y^2 + ring_map(f)*x*z^4 + ring_map(g)*z^6
end


################################################################
# 3: Construct the Tate polynomial
################################################################

function _tate_sections(base::NormalToricVariety)
  a1 = generic_section(anticanonical_bundle(base))
  a2 = generic_section(anticanonical_bundle(base)^2)
  a3 = generic_section(anticanonical_bundle(base)^3)
  a4 = generic_section(anticanonical_bundle(base)^4)
  a6 = generic_section(anticanonical_bundle(base)^6)
  return [a1, a2, a3, a4, a6]
end

function _tate_polynomial(base::NormalToricVariety, S::MPolyRing)
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
# 4: A base space for efficient testing
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
  return normal_toric_variety(cones, rays)
end


################################################################
# 5: Check if an ideal/subvariety is nontrivial
################################################################

_is_nontrivial(id::MPolyIdeal{T}, irr::MPolyIdeal{T}) where {T<:MPolyRingElem} = !is_one(id) && !is_one(saturation(id, irr))


################################################################
# 6: Compute singularity Kodaira type and refined Tate type
################################################################

_count_factors(poly::QQMPolyRingElem) = mapreduce(p -> p[end], +, absolute_primary_decomposition(ideal([poly])))

_string_from_factor_count(poly::QQMPolyRingElem, string_list::Vector{String}) = string_list[_count_factors(poly)]

function _kodaira_type(id::MPolyIdeal{T}, f::T, g::T, d::T, ords::Tuple{Int64, Int64, Int64}) where {T<:MPolyRingElem}
  f_ord = ords[1]
  g_ord = ords[2]
  d_ord = ords[3]
    
  # Check for cases where there cannot be Tate monodromy
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
    # Create new ring with auxiliary variable to construct the monodromy polynomial
    R = parent(f)
    S, (_psi, ) = polynomial_ring(QQ, ["_psi"; [string(v) for v in gens(R)]], cached = false)
    ring_map = hom(R, S, gens(S)[2:end])
    poly_f = ring_map(f)
    poly_g = ring_map(g)
    poly_d = ring_map(d)
    locus = ring_map(gens(id)[1])
    
    # Compute monodromy polynomial and check factorization for remaining cases
    if f_ord == 0 && g_ord == 0
      g_quotient = divrem(9 * poly_g, locus)[2]
      f_quotient = divrem(2 * poly_f, locus)[2]
      quotient_val = div(g_quotient, f_quotient)

      monodromy_poly = _psi^2 + quotient_val
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I_$d_ord", "Split I_$d_ord"])
    elseif d_ord == 4 && g_ord == 2 && f_ord >= 2
      g_quotient = divrem(div(poly_g, locus^2), locus)[2]

      monodromy_poly = _psi^2 - g_quotient
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split IV", "Split IV"])
    elseif d_ord == 6 && f_ord >= 2 && g_ord >= 3
      f_quotient = divrem(div(poly_f, locus^2), locus)[2]
      g_quotient = divrem(div(poly_g, locus^3), locus)[2]
      
      monodromy_poly = _psi^3 + _psi * f_quotient + g_quotient
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I^*_0", "Semi-split I^*_0", "Split I^*_0"])
    elseif f_ord == 2 && g_ord == 3 && d_ord >= 7
      d_quotient = div(poly_d, locus^d_ord)
      f_quotient = div(2 * poly_f, locus^2)
      g_quotient = div(9 * poly_g, locus^3)
      num_quotient = divrem(d_quotient * f_quotient^(2 + d_ord % 2), locus)[2]
      den_quotient = divrem(4 * g_quotient^(2 + d_ord % 2), locus)[2]
      quotient_val = div(num_quotient, den_quotient)

      monodromy_poly = _psi^2 + quotient_val
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I^*_$(d_ord - 6)", "Split I^*_$(d_ord - 6)"])
    elseif d_ord == 8 && g_ord == 4 && f_ord >= 3
      g_quotient = divrem(div(poly_g, locus^4), locus)[2]

      monodromy_poly = _psi^2 - g_quotient
      kod_type = _string_from_factor_count(monodromy_poly, ["Non-split IV^*", "Split IV^*"])
    else
      kod_type = "Unrecognized"
    end
  end
  
  return kod_type
end


##################################################################
# 7: Blowups (old helper function, to be used for family of bases)
##################################################################

function _blowup_global(id::MPolyIdeal{QQMPolyRingElem}, center::MPolyIdeal{QQMPolyRingElem}, irr::MPolyIdeal{QQMPolyRingElem}, sri::MPolyIdeal{QQMPolyRingElem}, lin::MPolyIdeal{<:MPolyRingElem}; index::Integer = 1)
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
_blowup_global(id::T, center::T, irr::T, sri::T, lin::T; index::Integer = 1) where {T<:MPolyIdeal{<:MPolyRingElem}} = _blowup_global(ideal(map(g -> g.f, gens(id))), ideal(map(g -> g.f, gens(center))), ideal(map(g -> g.f, gens(irr))), ideal(map(g -> g.f, gens(sri))), lin, index = index)


function _blowup_global_sequence(id::MPolyIdeal{QQMPolyRingElem}, centers::Vector{<:Vector{<:Integer}}, irr::MPolyIdeal{QQMPolyRingElem}, sri::MPolyIdeal{QQMPolyRingElem}, lin::MPolyIdeal{<:MPolyRingElem}; index::Integer = 1)
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
_blowup_global_sequence(id::T, centers::Vector{<:Vector{<:Integer}}, irr::T, sri::T, lin::T; index::Integer = 1) where {T<:MPolyIdeal{<:MPolyRingElem}} = _blowup_global_sequence(ideal(map(g -> g.f, gens(id))), centers, ideal(map(g -> g.f, gens(irr))), ideal(map(g -> g.f, gens(sri))), lin, index = index)



###########################################################################
# 8: Constructing a generic sample for models over not-fully specified spaces
###########################################################################

function _construct_generic_sample(base_grading::Matrix{Int64}, base_vars::Vector{String}, d::Int)
  base_space = family_of_spaces(polynomial_ring(QQ, base_vars, cached = false)[1], base_grading, d)
  ambient_space_vars = vcat(base_vars, ["x", "y", "z"])
  coordinate_ring_ambient_space = polynomial_ring(QQ, ambient_space_vars, cached = false)[1]
  ambient_space_grading = zero_matrix(Int, nrows(base_grading)+1,ncols(base_grading)+3)
  ambient_space_grading[1:nrows(base_grading),1:ncols(base_grading)] = base_grading
  ambient_space_grading[1,ncols(base_grading)+1:ncols(base_grading)+2] = [2; 3]
  ambient_space_grading[nrows(base_grading) + 1,ncols(base_grading) + 1:ncols(base_grading) + 3] = [2; 3; 1]
  ambient_space = family_of_spaces(coordinate_ring_ambient_space, ambient_space_grading, d+2)
  return [coordinate_ring(base_space), base_space, ambient_space]
end


function _construct_generic_sample(base_grading::Matrix{Int64}, base_vars::Vector{String}, d::Int, fiber_ambient_space::NormalToricVariety, fiber_twist_divisor_classes::ZZMatrix)
  base_space = family_of_spaces(polynomial_ring(QQ, base_vars, cached = false)[1], base_grading, d)
  ambient_space_vars = vcat(base_vars, coordinate_names(fiber_ambient_space))
  coordinate_ring_ambient_space = polynomial_ring(QQ, ambient_space_vars, cached = false)[1]
  w = Matrix{Int64}(reduce(vcat, [k.coeff for k in cox_ring(fiber_ambient_space).d]))
  z_block = zeros(Int64, ncols(w), ncols(base_grading))
  D_block = hcat([[Int(fiber_twist_divisor_classes[k,l]) for k in 1:nrows(fiber_twist_divisor_classes)] for l in 1:ncols(fiber_twist_divisor_classes)]...)
  ambient_space_grading = [base_grading D_block; z_block w']
  ambient_space = family_of_spaces(coordinate_ring_ambient_space, ambient_space_grading, d+dim(fiber_ambient_space))
  return [coordinate_ring(ambient_space), base_space, ambient_space]
end



###########################################################################
# 9: Evaluating a string to a polynomial
###########################################################################

function _eval_poly(E::Expr, vars)
  @assert E.head == :call
  if E.args[1] == :+
    return reduce(+, (_eval_poly(E.args[i], vars) for i in 2:length(E.args)))
  elseif E.args[1] == :*
    return reduce(*, (_eval_poly(E.args[i], vars) for i in 2:length(E.args)))
  elseif E.args[1] == :-
    if length(E.args) == 2
      return -_eval_poly(E.args[2], vars)
    else
      @assert length(E.args) == 3
      return _eval_poly(E.args[2], vars) - _eval_poly(E.args[3], vars)
    end
  elseif E.args[1] == :^
    return _eval_poly(E.args[2], vars)^_eval_poly(E.args[3], vars)
  elseif E.args[1] == ://
    @assert E.args[2] isa Number && E.args[3] isa Number
    return E.args[2]//E.args[3]
  end
end

function _eval_poly(E::Symbol, vars)
  return vars[E]
end

function _eval_poly(E::Number, vars)
  return E
end

function eval_poly(s::String, R)
  if (R isa PolyRing || R isa MPolyRing)
    symR = symbols(R) # Symbol[]
    genR = gens(R)
  else
    symR = []
    genR = []
  end

  return R(_eval_poly(Meta.parse(s), Dict(symR[i] => genR[i] for i in 1:length(symR))))
end

eval_poly(n::Number, R) = R(n)

# Example
# julia> Qx, (x1, x2) = QQ["x1", "x2"];
#
# julia> eval_poly("-x1 - 3//5*x2^3 + 5 - 3", Qx)
# -x1 - 3//5*x2^3 + 2



##########################################
### 10 strict_transform helpers
##########################################

_strict_transform(bd::AbsCoveredSchemeMorphism, II::AbsIdealSheaf; coordinate_name = "e") = strict_transform(bd, II)

function _strict_transform(bd::ToricBlowdownMorphism, II::ToricIdealSheafFromCoxRingIdeal; coordinate_name = "e")
  center_ideal = ideal_in_cox_ring(center(bd))
  if (ngens(ideal_in_cox_ring(II)) != 1) || (all(x -> x in gens(base_ring(center_ideal)), gens(center_ideal)) == false)
    return strict_transform(bd, II)
  end
  S = cox_ring(domain(bd))
  _e = eval_poly(coordinate_name, S)
  images = MPolyRingElem[]
  for v in gens(S)
    v == _e && continue
    if string(v) in [string(k) for k in gens(ideal_in_cox_ring(center(bd)))]
      push!(images, v * _e)
    else
      push!(images, v)
    end
  end
  ring_map = hom(cox_ring(codomain(bd)), S, images)
  total_transform = ring_map(ideal_in_cox_ring(II))
  exceptional_ideal = total_transform + ideal([_e])
  strict_transform, exceptional_factor = saturation_with_index(total_transform, exceptional_ideal)
  return ideal_sheaf(domain(bd), strict_transform)
end
