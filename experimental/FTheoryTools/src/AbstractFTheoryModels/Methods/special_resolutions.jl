function _resolve_ambient_space_1511_03209(resolved_model::GlobalTateModel)
  # Additional blowup 1:
  # Additional blowup 1:
  # Ambient space has the following rays
  #x: [0,0,0,-3,1]
  #y: [0,0,0,2,-1]
  #z: [0,0,0,0,1]
  # We add the ray m1: (0,0,0,1,0). This looks like y + z = 2 * m1.
  # So naively, I think of this as blowing up y^2 = z^2 = 0 and introducing the variable m1. For the strict transform, we thus do
  # y^2 -> y^2 * m1
  # z^2 -> z^2 * m1
  # y * z -> y * z * m1
  as = ambient_space(resolved_model)
  bl = domain(blow_up(as, [0, 0, 0, 1, 0]; coordinate_name="m1"))
  f = hypersurface_equation(resolved_model)
  my_mons = collect(monomials(f))
  pos_1 = findfirst(k -> k == "y", [string(a) for a in gens(coordinate_ring(as))])
  pos_2 = findfirst(k -> k == "z", [string(a) for a in gens(coordinate_ring(as))])
  exp_list = [collect(exponents(m))[1] for m in my_mons]
  my_exps = [[k[pos_1], k[pos_2]] for k in exp_list]
  @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
  m_power = [Int(1//2 * sum(a)) for a in my_exps]
  overall_factor = minimum(m_power)
  new_coeffs = collect(coefficients(f))
  new_exps = [
    vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)
  ]
  my_builder = MPolyBuildCtx(coordinate_ring(bl))
  for a in 1:length(new_exps)
    push_term!(my_builder, new_coeffs[a], new_exps[a])
  end
  new_tate_polynomial = finish(my_builder)
  model_bl = GlobalTateModel(
    explicit_model_sections(resolved_model),
    model_section_parametrization(resolved_model),
    new_tate_polynomial,
    base_space(resolved_model),
    bl,
  )
  set_attribute!(model_bl, :partially_resolved, true)
  # Additional blowup 2:
  # Additional blowup 2:
  # Ambient space has the following rays:
  # x: [0,0,0,-3,1]
  # y: [0,0,0,2,-1]
  # z: [0,0,0,0,1]
  # m1: [0,0,0,1,0]
  # We add the ray m2: (0,0,0,-2,1). This looks like 2 * x + z = 3 * m2. For the strict transform, we thus do
  # x^3 -> x^3 * m2^2
  # z^3 -> z^3 * m2
  as = ambient_space(model_bl)
  bl = domain(blow_up(as, [0, 0, 0, -2, 1]; coordinate_name="m2"))
  f = hypersurface_equation(model_bl)
  my_mons = collect(monomials(f))
  pos_1 = findfirst(k -> k == "x", [string(a) for a in gens(coordinate_ring(as))])
  pos_2 = findfirst(k -> k == "z", [string(a) for a in gens(coordinate_ring(as))])
  exp_list = [collect(exponents(m))[1] for m in my_mons]
  my_exps = [[k[pos_1], k[pos_2]] for k in exp_list]
  @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
  m_power = [Int(2//3 * a[1] + 1//3 * a[2]) for a in my_exps]
  overall_factor = minimum(m_power)
  new_coeffs = collect(coefficients(f))
  new_exps = [
    vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)
  ]
  my_builder = MPolyBuildCtx(coordinate_ring(bl))
  for a in 1:length(new_exps)
    push_term!(my_builder, new_coeffs[a], new_exps[a])
  end
  new_tate_polynomial = finish(my_builder)
  model_bl2 = GlobalTateModel(
    explicit_model_sections(model_bl),
    model_section_parametrization(model_bl),
    new_tate_polynomial,
    base_space(model_bl),
    bl,
  )
  set_attribute!(model_bl2, :partially_resolved, true)
  # Additional blowup 3:
  # Additional blowup 3:
  # Ambient space has the following rays:
  # x: [0,0,0,-3,1]
  # y: [0,0,0,2,-1]
  # z: [0,0,0,0,1]
  # m1: [0,0,0,1,0]
  # m2: [0,0,0,-2,1]
  # We add the ray m3: (0,0,0,-1,1). This looks like m2 + z = 2 * m3. For the strict transform, we thus do
  # m2^2 -> m2^2 * m3
  # z^2 -> z^2 * m3
  as = ambient_space(model_bl2)
  bl = domain(blow_up(as, [0, 0, 0, -1, 1]; coordinate_name="m3"))
  f = hypersurface_equation(model_bl2)
  my_mons = collect(monomials(f))
  pos_1 = findfirst(k -> k == "m2", [string(a) for a in gens(coordinate_ring(as))])
  pos_2 = findfirst(k -> k == "z", [string(a) for a in gens(coordinate_ring(as))])
  exp_list = [collect(exponents(m))[1] for m in my_mons]
  my_exps = [[k[pos_1], k[pos_2]] for k in exp_list]
  @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
  m_power = [Int(1//2 * sum(a)) for a in my_exps]
  overall_factor = minimum(m_power)
  new_coeffs = collect(coefficients(f))
  new_exps = [
    vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)
  ]
  my_builder = MPolyBuildCtx(coordinate_ring(bl))
  for a in 1:length(new_exps)
    push_term!(my_builder, new_coeffs[a], new_exps[a])
  end
  new_tate_polynomial = finish(my_builder)
  model_bl3 = GlobalTateModel(
    explicit_model_sections(model_bl2),
    model_section_parametrization(model_bl2),
    new_tate_polynomial,
    base_space(model_bl2),
    bl,
  )
  set_attribute!(model_bl3, :partially_resolved, true)

  # We confirm that after these steps, we achieve what we desire.
  @req is_smooth(ambient_space(model_bl3)) "Ambient space not yet smooth. Please inform the authors!"
  @req is_homogeneous(hypersurface_equation(model_bl3)) "Strict transform is not homogeneous. Please inform the authors!"
  return model_bl3
end
