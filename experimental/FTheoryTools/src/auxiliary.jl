################################################################
# 1: Construct ambient space from given base
################################################################

function _ambient_space(base::NormalToricVariety, fiber_amb_space::NormalToricVariety, fiber_twist_divisor_classes::Vector{ToricDivisorClass})
  @req all(D -> toric_variety(D) === base, fiber_twist_divisor_classes) "The divisors must belong to the (same) base space"
  
  # Extract information about the toric base
  b_rays = matrix(ZZ, rays(base))
  b_cones = matrix(ZZ, ray_indices(maximal_cones(base)))
  b_grades = reduce(vcat, [elem.coeff for elem in coordinate_ring(base).d])
  b_var_names = symbols(coordinate_ring(base))
  
  # Extract information about the fiber ambient space
  f_rays = matrix(ZZ, rays(fiber_amb_space))
  f_cones = matrix(ZZ, ray_indices(maximal_cones(fiber_amb_space)))
  f_grades = reduce(vcat, [elem.coeff for elem in coordinate_ring(fiber_amb_space).d])
  f_var_names = symbols(coordinate_ring(fiber_amb_space))
  
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
  a_space_class_group = free_abelian_group(ncols(b_grades) + torsion_free_rank(class_group_with_map(fiber_amb_space)[1]))
  
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
  set_attribute!(a_space, :class_group, codomain(a_space_grading))
  set_attribute!(a_space, :torusinvariant_weil_divisor_group, a_space_divisor_group)
  return a_space
end


################################################################
# 2: Construct the Weierstrass polynomial
################################################################

function _weierstrass_sections(base::NormalToricVariety; rng::AbstractRNG = Random.default_rng())
  return [generic_section(anticanonical_bundle(base)^4; rng), generic_section(anticanonical_bundle(base)^6; rng)]
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

function _tate_sections(base::NormalToricVariety; rng::AbstractRNG = Random.default_rng())
  a1 = generic_section(anticanonical_bundle(base); rng)
  a2 = generic_section(anticanonical_bundle(base)^2; rng)
  a3 = generic_section(anticanonical_bundle(base)^3; rng)
  a4 = generic_section(anticanonical_bundle(base)^4; rng)
  a6 = generic_section(anticanonical_bundle(base)^6; rng)
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
# 4: Compute singularity Kodaira type and refined Tate type
################################################################

_count_factors(poly::QQMPolyRingElem) = mapreduce(p -> p[end], +, absolute_primary_decomposition(ideal([poly])))

_string_from_factor_count(poly::QQMPolyRingElem, string_list::Vector{String}) = string_list[_count_factors(poly)]

function _kodaira_type(id::MPolyIdeal{<:MPolyRingElem}, ords::Tuple{Int64, Int64, Int64}, w::WeierstrassModel; rng::AbstractRNG = Random.default_rng())
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
  elseif d_ord == 6 && f_ord >= 2 && g_ord >= 3
    # For type I_0^* singularities, we have to rely on the old method for now,
    # which is not always dependable

    f = weierstrass_section_f(w)
    g = weierstrass_section_g(w)
    d = discriminant(w)

    # Create new ring with auxiliary variable to construct the monodromy polynomial
    R = parent(f)
    S, (_psi,), _old_gens = polynomial_ring(QQ, [:_psi], symbols(R); cached = false)
    ring_map = hom(R, S, _old_gens)
    poly_f = ring_map(f)
    poly_g = ring_map(g)
    locus = ring_map(gens(id)[1])

    f_quotient = divrem(div(poly_f, locus^2), locus)[2]
    g_quotient = divrem(div(poly_g, locus^3), locus)[2]
    
    monodromy_poly = _psi^3 + _psi * f_quotient + g_quotient
    kod_type = _string_from_factor_count(monodromy_poly, ["Non-split I^*_0", "Semi-split I^*_0", "Split I^*_0"])
  else
    # If the base is arbitrary, we tune the model over projective space of the
    # appropriate dimension. This allows us to use the same algorithm for all
    # cases. The choice of projective space here is an attempt to minimize the
    # chances of accidental gauge enhancement
    if !is_base_space_fully_specified(w)
      # Build the new concrete base, and get the anticanonical and hyperplane
      # bundles. We choose the hyperplane bundle for all gauge loci over the
      # concrete base as an additional measure to avoid accidental gauge
      # enhancement
      concrete_base = projective_space(NormalToricVariety, dim(base_space(w)))
      KBar = anticanonical_bundle(concrete_base)
      hyperplane_bundle = toric_line_bundle(torusinvariant_prime_divisors(concrete_base)[1])

      # Get the grading matrix and the coordinates of the arbitrary base
      grading = weights(base_space(w))
      base_coords_symbols = symbols(coordinate_ring(base_space(w)))
      @req (length(base_coords_symbols) == length(grading[1, :])) "The number of columns in the weight matrix does not match the number of base coordinates"

      # Choose explicit sections for all parameters of the model,
      # and then put the model over the concrete base using these data
      concrete_data = merge(Dict(string(base_coords_symbols[i]) => generic_section(KBar^grading[1, i] * prod(hyperplane_bundle^grading[j, i] for j in 2:length(grading[:, 1])); rng) for i in eachindex(base_coords_symbols)), Dict("base" => concrete_base))
      w = put_over_concrete_base(w, concrete_data; rng)

      # We also need to determine the gauge locus over the new base
      # by using the explicit forms of all of the sections chosen above
      list_of_sections = [concrete_data[string(base_coords_symbols[i])] for i in eachindex(base_coords_symbols)]
      id = ideal([evaluate(p, list_of_sections) for p in gens(id)])
    end

    f = weierstrass_section_f(w)
    g = weierstrass_section_g(w)
    d = discriminant(w)

    # For now, we explicitly require that the gauge ideal is principal
    @req (ngens(id) == 1) "Gauge ideal is not principal"

    # Over concrete bases, we randomly reduce the polynomials defining the gauge
    # divisor to only two variables so that the is_radical check is faster. This
    # could give an incorrect result (radical or not), so we actually try this
    # five times and see if we get agreement among all of the results
    num_gens = ngens(parent(f))
    gauge2s, f2s, g2s, d2s = [], [], [], []
    for _ in 1:5
      coord_inds = randperm(rng, num_gens)[1:end-2]
      rand_ints = rand(rng, -100:100, num_gens - 2)

      push!(gauge2s, evaluate(forget_decoration(gens(id)[1]), coord_inds, rand_ints))
      push!(f2s, evaluate(forget_decoration(f), coord_inds, rand_ints))
      push!(g2s, evaluate(forget_decoration(g), coord_inds, rand_ints))
      push!(d2s, evaluate(forget_decoration(d), coord_inds, rand_ints))
    end

    # Check monodromy conditions for remaining cases.
    # Default to split when there is disagreement among the five attempts,
    # because this approach seems to skew toward accidentally identifying
    # a singularity as non-split
    if f_ord == 0 && g_ord == 0
      quotients = []
      for i in eachindex(gauge2s)
        push!(quotients, quotient(ideal([9 * g2s[i], gauge2s[i]]), ideal([2 * f2s[i], gauge2s[i]])))
      end

      kod_type = if all(is_radical, quotients) "Non-split I_$d_ord" else "Split I_$d_ord" end
    elseif d_ord == 4 && g_ord == 2 && f_ord >= 2
      quotients = []
      for i in eachindex(gauge2s)
        push!(quotients, quotient(ideal([g2s[i]]), ideal([gauge2s[i]^2])) + ideal([gauge2s[i]]))
      end

      kod_type = if all(is_radical, quotients) "Non-split IV" else "Split IV" end
    elseif f_ord == 2 && g_ord == 3 && d_ord >= 7
      quotients = []
      if d_ord % 2 == 0
        for i in eachindex(gauge2s)
          push!(quotients, quotient(ideal([4 // 81 * (d2s[i] * f2s[i]^2) / gauge2s[i]^(d_ord + 4), gauge2s[i]]), ideal([g2s[i]^2 / gauge2s[i]^6, gauge2s[i]])))
        end
      else
        for i in eachindex(gauge2s)
          push!(quotients, quotient(ideal([2 // 729 * (d2s[i] * f2s[i]^3) / gauge2s[i]^(d_ord + 6), gauge2s[i]]), ideal([g2s[i]^3 / gauge2s[i]^9, gauge2s[i]])))
        end
      end

      kod_type = if all(is_radical, quotients) "Non-split I^*_$(d_ord - 6)" else "Split I^*_$(d_ord - 6)" end
    elseif d_ord == 8 && g_ord == 4 && f_ord >= 3
      quotients = []
      for i in eachindex(gauge2s)
        push!(quotients, quotient(ideal([g2s[i]]), ideal([gauge2s[i]^4])) + ideal([gauge2s[i]]))
      end

      kod_type = if all(is_radical, quotients) "Non-split IV^*" else "Split IV^*" end
    else
      kod_type = "Unrecognized"
    end
  end
  
  return kod_type
end


###########################################################################
# 5: Constructing a generic sample for models over not-fully specified spaces
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
  coordinate_ring_ambient_space, _ = polynomial_ring(QQ, base_vars, coordinate_names(fiber_ambient_space); cached=false)
  w = Matrix{Int64}(reduce(vcat, [k.coeff for k in coordinate_ring(fiber_ambient_space).d]))
  z_block = zeros(Int64, ncols(w), ncols(base_grading))
  D_block = hcat([[Int(fiber_twist_divisor_classes[k,l]) for k in 1:nrows(fiber_twist_divisor_classes)] for l in 1:ncols(fiber_twist_divisor_classes)]...)
  ambient_space_grading = [base_grading D_block; z_block w']
  ambient_space = family_of_spaces(coordinate_ring_ambient_space, ambient_space_grading, d+dim(fiber_ambient_space))
  return [coordinate_ring(ambient_space), base_space, ambient_space]
end



###########################################################################
# 6: Evaluating a string to a polynomial
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
# julia> Qx, (x1, x2) = QQ[:x1, :x2];
#
# julia> eval_poly("-x1 - 3//5*x2^3 + 5 - 3", Qx)
# -x1 - 3//5*x2^3 + 2



###########################################################################
# 7: Apply a function to the innermost nested structure
###########################################################################

function deepmap(f, x)
    if x isa AbstractArray
        return map(e -> deepmap(f, e), x)
    else
        return f(x)
    end
end



###########################################################################
# 8: Macro for function generation
###########################################################################

macro define_model_attribute_getter(arg_expr, doc_example="", doc_link="", attr_name=nothing)
  if !(arg_expr isa Expr && arg_expr.head == :tuple && length(arg_expr.args) == 2)
    error("Expected input like: (function_name, ReturnType)")
  end

  fname_expr = arg_expr.args[1]
  rettype_expr = arg_expr.args[2]
  fname = fname_expr isa Symbol ? fname_expr : error("function_name is not a symbol")

  # Determine attribute name symbol: use attr_name if provided, else function name
  attr_sym = attr_name === nothing ? fname : (attr_name isa Symbol ? attr_name : error("attr_name must be a Symbol if provided"))
  sym = QuoteNode(attr_sym)  
  
  msg = "No $(replace(string(fname), '_' => ' ')) known for this model"

  # Conditionally build doc section
  doc_link_section = isempty(doc_link) ? "" : "\n\n$doc_link"
  examples_section = isempty(doc_example) ? "" : "\n\n# Examples\n$doc_example"  
  default_doc = """
      $(fname)(m::AbstractFTheoryModel)
  
  Return `$(fname)` of the F-theory model if known, otherwise throw an error.$doc_link_section$examples_section
  """

  return quote
    @doc $default_doc
    function $(esc(fname))(m::AbstractFTheoryModel)
      @req has_attribute(m, $sym) $msg
      return get_attribute(m, $sym)::$(esc(rettype_expr))
    end
  end
end
