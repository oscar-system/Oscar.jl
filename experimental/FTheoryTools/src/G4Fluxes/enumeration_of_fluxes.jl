#################################################
# (1) Read in the model (~ 8 minutes)
#################################################

using Oscar
using ProgressMeter

h_res = load("/home/i/1511.03209/ResolvedModelFiles/1511-03209-res.mrdi")





#################################################
# (2) Prepare the toric ambient space
#################################################

const amb = ambient_space(h_res);
const pt = tate_polynomial(h_res);
# The following line is important but takes a long time to execute...
@req is_homogeneous(pt) == true "Tate polynomial is not homogeneous!"
const S = parent(pt);
const gS = gens(S);
set_attribute!(amb, :cox_ring, S)
set_attribute!(amb, :coordinate_names, string.(gS))
const mnf = Oscar._minimal_nonfaces(amb);
const indices = [Set(Vector{Int}(Polymake.row(mnf,i))) for i in 1:Polymake.nrows(mnf)];
function can_be_ignored(my_set)
  for k in 1:length(indices)
    if is_subset(indices[k], my_set)
      return true
    end
  end
  return false
end
stanley_reisner_ideal(amb);
ideal_of_linear_relations(amb);
cohomology_ring(amb, check = false);
const c_ds = [k.f for k in gens(cohomology_ring(amb))];
const sr_ideal_gens = gens(stanley_reisner_ideal(amb));
const sr_ideal_exponents = [collect(exponents(sr_ideal_gens[k]))[1] for k in 1:length(sr_ideal_gens)];
const scaling_relations = transpose(vcat([c.coeff for c in S.d]...));
@req rank(scaling_relations) == 305 "The rank of the scaling relations is reduced!"





#################################################
# (3) Work out a basis of H^(2,2)(amb) = C^{2104}
#################################################

# (3.1) Rule out variable by use of the linear relations
# gs = gens(ideal_of_linear_relations(amb));
# Via gs[1] we get w1: Position 2
# Via gs[2] we get w2: Position 3
# Via gs[3] we get w0: Position 1
# Via gs[4] we get e99_12: Position 308
# Via gs[5] we get e99_4: Position 300
# So should ignore those!
good_positions = [k for k in 1:length(gS)];
deleteat!(good_positions, 308);
deleteat!(good_positions, 300);
deleteat!(good_positions, 3);
deleteat!(good_positions, 2);
deleteat!(good_positions, 1);


# (3.2) Work out the basis
g4_basis = [];
g4_basis_indices = [];
for i in 1:length(good_positions)
  for j in i:length(good_positions)
    il = [good_positions[i], good_positions[j]]
    can_be_ignored(Set(il)) && continue
    push!(g4_basis, c_ds[il[1]] * c_ds[il[2]])
    push!(g4_basis_indices, il)
  end
end
@req length(g4_basis) == 1990 "Error in computation of H^(2,2)(ambient_space)!!!"





#################################################
# (4) Introduce a sophisticated intersection
# (4) product method, so that we can just loop
# (4) with this function and get the intersection
# (4) numbers that we seek
#################################################


# (4.1) A function to reduce the tate polynomial
function reduce_tate_polynomial(mon::MPolyRingElem)

  # Set variables to zero in Tate polynomial and check for degenerate case
  exp_list = collect(exponents(mon))[1];
  vanishing_vars_pos = sort(findall(k -> exp_list[k] != 0, 1:length(gS)))
  new_pt = pt;
  for m in vanishing_vars_pos
    new_pt = divrem(new_pt, gS[m])[2];
  end
  if is_constant(new_pt)
    return [new_pt, [], [], zero(matrix_algebra(ZZ, 0))]
  end
  
  # Identify all variables that can be set to one
  set_to_one_list = Int[]
  for pos_seek in vanishing_vars_pos
    for my_exps in sr_ideal_exponents
      if my_exps[pos_seek] != 0
        non_zero_positions = findall(k -> k != 0, my_exps)
        if length(non_zero_positions) == 2
          other_variable_index = non_zero_positions[1]
          if other_variable_index == pos_seek
            other_variable_index = non_zero_positions[2]
          end
          push!(set_to_one_list, other_variable_index)
        end
      end
    end
  end
  for my_exps in sr_ideal_exponents
    non_zero_positions = findall(k -> k != 0, my_exps)
    if length(non_zero_positions) == 3
      for k1 in 1:length(vanishing_vars_pos)
        for k2 in k1+1:length(vanishing_vars_pos)
          if my_exps[vanishing_vars_pos[k1]] != 0 && my_exps[vanishing_vars_pos[k2]] != 0
            other_variable_index = non_zero_positions[1]
            if Set([vanishing_vars_pos[k1], vanishing_vars_pos[k2]]) == Set([non_zero_positions[1], non_zero_positions[2]])
              other_variable_index = non_zero_positions[3]
            end
            if Set([vanishing_vars_pos[k1], vanishing_vars_pos[k2]]) == Set([non_zero_positions[1], non_zero_positions[3]])
              other_variable_index = non_zero_positions[2]
            end
            push!(set_to_one_list, other_variable_index)
          end
        end
      end
    end
  end
  set_to_one_list = sort(unique(set_to_one_list));

  # Compute remaining variables and scaling relations, thereby verifying that indeed we can set the above variables all to one
  remaining_vars_pos = [k for k in 1:length(gS) if k ∉ vcat(set_to_one_list, vanishing_vars_pos)]
  remaining_vars = [gS[k] for k in remaining_vars_pos]

  # Extract remaining Stanley-Reisner ideal relations
  sr_reduced = Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
  my_set = Set(set_to_one_list)
  for k in 1:length(indices)
    if isdisjoint(my_set, indices[k])
      push!(sr_reduced, sr_ideal_gens[k])
    end
  end
  images = [m in vanishing_vars_pos ? one(S) : gS[m] for m in 1:length(gS)]
  my_ring_map = hom(S, S, images)
  sr_reduced = unique([my_ring_map(g) for g in sr_reduced])

  # Extract the remaining scaling relations
  @req 305 >= length(set_to_one_list) "Not enough scaling relations to set all of the desired variables to one"
  prepared_scaling_relations = hcat(vcat([scaling_relations[:, k] for k in set_to_one_list], [scaling_relations[:, k] for k in remaining_vars_pos])...);
  prepared_scaling_relations = hnf(matrix(ZZ, prepared_scaling_relations));
  reduced_scaling_relations = prepared_scaling_relations[length(set_to_one_list) + 1: nrows(prepared_scaling_relations), length(set_to_one_list) + 1 : ncols(prepared_scaling_relations)]
  @req rank(reduced_scaling_relations) == rank(scaling_relations) - length(set_to_one_list) "Found inconsistency in extracting remaining scaling relations"

  # Set all variables to one that we can set to one
  images = [k in remaining_vars_pos ? gS[k] : one(S) for k in 1:length(gS)];
  pt_reduced = hom(S, S, images)(new_pt);

  # Return the result
  return [pt_reduced, sr_reduced, remaining_vars, reduced_scaling_relations]
end


# (4.2) Function to find linearly equivalent expression for the case that one variable appears trice.
function rationally_equivalent_cycle(mon::MPolyRingElem)

  # Extract the exponents
  exps = collect(exponents(mon))[1]

  # Identify positions of the single and triple variable
  power_variable = findfirst(k -> k > 1, exps)
  other_variables = findall(k -> k == 1, exps)
  if power_variable === nothing
    index = rand(1:length(other_variables))
    power_variable = other_variables[index]
  end
  filter!(x -> x != power_variable, other_variables)
  @req length(other_variables) + 1 <= 5 "Found too many variables -- will likely not find a suitable relation!"

  # Let us simplify the problem by extracting the entries in the columns of single_variables and double_variables of the linear relation matrix
  linear_relations_matrix = matrix(ZZ, rays(amb));
  simpler_matrix = linear_relations_matrix[vcat(other_variables, power_variable), :]
  simpler_matrix = matrix(QQ, simpler_matrix)
  b = [[0] for k in 1:length(other_variables) + 1]
  b[end] = [1]
  b = matrix(QQ, b)
  A = solve(simpler_matrix, b; side = :right)
  
  # Now form the relation in case...
  employed_relation = (-1) * sum([QQ(l) for l in linear_relations_matrix[:,k]] * A[k] for k in 1:5)
  employed_relation[power_variable] = 0
  employed_relation = sum(employed_relation[k] * c_ds[k] for k in 1:length(employed_relation))
  
  # Form an equivalent linear combination of transverse intersections
  return mon / c_ds[power_variable] * employed_relation

end


# (4.3) Sophisticated intersection product function
#intersection_dict = Dict{String, Int}()
intersection_dict = load("/home/i/1511.03209/intersection_dict.mrdi")
#special_intersection_dict = Dict{String, Int}()
special_intersection_dict = load("/home/i/1511.03209/special_intersection_dict.mrdi")

function intersection_from_equivalent_cycle(mon::MPolyRingElem)
  intersect_numb = 0
  equivalent_linear_combination = rationally_equivalent_cycle(mon)
  coeffs_list = collect(coefficients(equivalent_linear_combination))
  mons_list = collect(monomials(equivalent_linear_combination))
  for k in 1:length(mons_list)
    intersect_numb += coeffs_list[k] * sophisticated_intersection_product(mons_list[k])
  end
  @req is_integer(intersect_numb) "Should have expected to find only integer intersection numbers..."
  intersection_dict[string(mon)] = ZZ(intersect_numb)
  return ZZ(intersect_numb)
end

function sophisticated_intersection_product(mon::MPolyRingElem)

  # (A) Have we computed this intersection number in the past? If so, just use that result...
  # (A) Have we computed this intersection number in the past? If so, just use that result...
  if haskey(intersection_dict, string(mon))
    return intersection_dict[string(mon)]
  end


  # (B) Get the indices of the variables that we try to intersect
  # (B) Get the indices of the variables that we try to intersect
  exps = collect(exponents(mon))[1]
  variable_pos = findall(k -> k != 0, exps)
  if can_be_ignored(Set(variable_pos))
    intersection_dict[string(mon)] = 0
    return 0
  end


  # (C) Let us deal with transverse intersection
  # (C) Let us deal with transverse intersection
  if length(variable_pos) == 4

    # C.1 Work out the intersection locus in detail.
    pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations = reduce_tate_polynomial(mon)

    # C.2 If pt == 0, then we are likely not looking at a transverse intersection. So take an equivalent cycle and try again...
    if is_zero(pt_reduced)
      return intersection_from_equivalent_cycle(mon)
    end

    # C.3 If pt is constant and non-zero, then the intersection is trivial.
    if is_constant(pt_reduced) && is_zero(pt_reduced) == false
      intersection_dict[string(mon)] = 0
      return 0
    end

    # C.4 Helper function for the cases below
    function has_one_and_rest_zero(vec)
      return count(==(1), vec) == 1 && all(x -> x == 0 || x == 1, vec)
    end

    # C.5 Cover a case that seems to appear frequently for our investigation:
    # pt_reduced of the form a * x + b * y for non-zero number a,b and remaining variables x, y subject to a reduced SR generator x * y and scaling relation [1,1].
    # This will thus always give exactly one solution (x = 1, y = -a/b), and so the intersection number is one.
    if length(gs_reduced) == 1 && length(remaining_vars) == 2
      mons_list = collect(monomials(pt_reduced))
      if length(mons_list) == 2
        if all(x -> x != 0, collect(coefficients(pt_reduced)))
          exps_list = [collect(exponents(k))[1] for k in mons_list]
          if has_one_and_rest_zero(exps_list[1]) && has_one_and_rest_zero(exps_list[2])
            if gs_reduced[1] == remaining_vars[1] * remaining_vars[2]
              if reduced_scaling_relations == matrix(ZZ, [[1,1]])
                intersection_dict[string(mon)] = 1
                return 1
              end
            end
          end
        end
      end
    end

    # C.6 Cover a case that seems to appear frequently for our investigation:
    # pt_reduced of the form a * x^n for non-zero number a and remaining variables x, y subject to a reduced SR generator x * y and scaling relation [*, != 0].
    # While this only gives the solution [0:1], I think that we must count it with multiplicity n. So n solutions.

    # Acutally, to avoid this subtle case. Let me try by allowing only n = 1... Maybe we can avoid those weird cases with n higher than one...

    # But this is a somewhat more tricky case it seems...
    if length(gs_reduced) == 1 && length(remaining_vars) == 2
      mons_list = collect(monomials(pt_reduced))
      if length(mons_list) == 1 && collect(coefficients(pt_reduced))[1] != 0
        list_of_exps = collect(exponents(mons_list[1]))[1]
        number_of_zeros = count(==(0), list_of_exps)
        if number_of_zeros == length(list_of_exps) - 1
          highest_power = list_of_exps[findfirst(x -> x > 0, list_of_exps)]
          if gs_reduced[1] == remaining_vars[1] * remaining_vars[2]
            if reduced_scaling_relations[1,2] != 0
              if highest_power == 1
                intersection_dict[string(mon)] = highest_power
                return highest_power
              end
            end
          end
        end
      end
    end
    
    # C.7 Cover a case that seems to appear frequently for our investigation. It looks as follows:
    # pt_reduced = -5700*w8*w10
    # remaining_vars = MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w8, w10]
    # gs_reduced = MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[w8*w10]
    # reduced_scaling_relations = [1 1]
    # This gives exactly two solutions, namely [0:1] and [1:0].
    if length(gs_reduced) == 1 && length(remaining_vars) == 2
      mons_list = collect(monomials(pt_reduced))
      if length(mons_list) == 1 && mons_list[1] == remaining_vars[1] * remaining_vars[2]
        if gs_reduced[1] == remaining_vars[1] * remaining_vars[2]
          if reduced_scaling_relations == matrix(ZZ, [[1,1]])
            intersection_dict[string(mon)] = 2
            return 2
          end
        end
      end
    end

    # C.8 Check if this was covered in our special cases
    if haskey(special_intersection_dict, string([pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations]))
      numb = special_intersection_dict[string([pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations])]
      intersection_dict[string(mon)] = numb
      return numb
    end

    # C.9 In all other cases, proceed via a rationally equivalent cycle
    println("")
    println("FOUND CASE THAT CANNOT YET BE DECIDED!")
    println("$pt_reduced")
    println("$remaining_vars")
    println("$gs_reduced")
    println("$mon")
    println("$reduced_scaling_relations")
    println("TRYING WITH EQUIVALENT CYCLE")
    println("")
    numb = intersection_from_equivalent_cycle(mon)
    special_intersection_dict[string([pt_reduced, gs_reduced, remaining_vars, reduced_scaling_relations])] = numb
    return numb
  end


  # (D) Deal with self-intersection cases.
  # (D) Deal with self-intersection cases.
  if length(variable_pos) < 4 && length(variable_pos) >= 1
    return intersection_from_equivalent_cycle(mon)
  end


  # (E) Ignore the worst self-intersections (for now).
  # (E) Ignore the worst self-intersections (for now).
  if length(variable_pos) == 0
    println("WEIRD! THIS SHOULD NEVER HAPPEN!")
    println(mon)
    println("")
  end

end





#################################################
# (5) Identify the constraints from verticality
# (5) constraints with the zero section.
# (5) Should take about 15 minutes.
#################################################

#=
function verticality_constraints_from_zero_section()
  constraint_matrix = Vector{Vector{Int64}}()
  pos_zero_section = findfirst(x -> x == "z", string.(gS));
  println("")
  @showprogress 1 "Processing..." for k in 1:101
    push!(constraint_matrix, [sophisticated_intersection_product(g4_basis[i] * c_ds[pos_zero_section] * c_ds[k]) for i in 1:length(g4_basis)])
  end
  println("")
  return matrix(ZZ, hcat(constraint_matrix...));
end
@time verticality_constraints_from_zero_section_matrix = verticality_constraints_from_zero_section();
@req rank(verticality_constraints_from_zero_section_matrix) == 98 "SOME INTERSECTION NUMBERS MUST HAVE CHANGED. RANK OF CONSTRAINT MATRIX SHOULD BE 98!"





#################################################
# (6) Identify the constraints from all other
# (6) verticality constraints.
# (6) Takes roughly: 7679 seconds (~ 2h)
#################################################

function other_verticality_constraints()
  constraint_matrix = Vector{Vector{Int64}}()
  for i in 1:101
    println("-------------------------------")
    println("Working with base divisor $i (of 101)\n")
    @showprogress 1 "Processing..." for j in 1:101
      push!(constraint_matrix, [sophisticated_intersection_product(g4_basis[k] * c_ds[i] * c_ds[j]) for k in 1:length(g4_basis)])
    end
    println("\n")
  end
  return matrix(ZZ, hcat(constraint_matrix...));
end
@time other_verticality_constraints_matrix = other_verticality_constraints();
@req rank(other_verticality_constraints_matrix) == 229 "SOME INTERSECTION NUMBERS MUST HAVE CHANGED. RANK OF SECOND CONSTRAINT MATRIX SHOULD BE 229!"





#################################################
# (7) Work out the vertical G4-fluxes
#################################################

# Combine all constraints into one matrix: 1990 columns and 102 * 101 = 10.302 rows:
total_verticality_constraint_matrix = transpose(hcat(verticality_constraints_from_zero_section_matrix, other_verticality_constraints_matrix));
n, N = nullspace(total_verticality_constraint_matrix);
@req n == 1712 "Different number of basis elements of vertical G4 fluxes found!"
@req rank(N) == 1712 "Internal inconsistency!"
basis_of_vertical_g4_fluxes = [sum(N[i, k] * g4_basis[i] for i in 1:1990) for k in 1:1712];
=#
basis_of_vertical_g4_fluxes = load("/home/i/1511.03209/basis_of_vertical_g4_fluxes =.mrdi")





#################################################
# (8) Towards the quantization condition
# (8) First time (from scratch): 9400 seconds.
# (8) Second time only a little bit faster...
#################################################

#=
# Work out all products Di * Dj of two toric divisors.
# Maybe I could use g4_basis from above, but I am not sure if the linear relations are ok to be invoked for the quantization checks.
h22_basis = [];
for i in 1:length(gS)
  for j in i:length(gS)
    il = [i, j]
    can_be_ignored(Set(il)) && continue
    push!(h22_basis, c_ds[il[1]] * c_ds[il[2]])
  end
end


# Compute intersection product of a polynomial
function sophisticated_intersection_product_from_poly(poly::MPolyRingElem)
  intersect_numb = 0
  coeffs_list = collect(coefficients(poly))
  mons_list = collect(monomials(poly))
  for k in 1:length(mons_list)
    intersect_numb += coeffs_list[k] * sophisticated_intersection_product(mons_list[k])
  end
  @req is_integer(intersect_numb) "Should have expected to find only integer intersection numbers..."
  return ZZ(intersect_numb)
end


# Iteration over all necessary intersection numbers to work on the quantization conditions
function quantization_constraints()
  constraint_matrix = Vector{Vector{Int64}}()
  for i in 1:length(basis_of_vertical_g4_fluxes)
    println("----------------------------")
    println("Working with basis element $i (of $(length(basis_of_vertical_g4_fluxes)))")
    condition = Int[]
    @showprogress 1 "Processing..." for j in 1:length(h22_basis)
      push!(condition, sophisticated_intersection_product_from_poly(basis_of_vertical_g4_fluxes[i] * h22_basis[j]))
    end
    push!(constraint_matrix, condition)
  end
  return matrix(ZZ, hcat(constraint_matrix...));
end
@time towards_quantization_conditions_matrix = quantization_constraints();
@req rank(towards_quantization_conditions_matrix) == 224 "Some intersection numbers must have changed!"
=#
towards_quantization_conditions_matrix = load("/home/i/1511.03209/towards_quantization_conditions_matrix.mrdi");





#################################################
# (9) Identify the quantization conditions as
# (9) polynomials, and prepare to solve them.
#################################################

# Figure out the non-trivial quantization conditions
non_zero_quantization_conditions_matrix = Vector{Vector{ZZRingElem}}()
for k in 1:nrows(towards_quantization_conditions_matrix)
  row = towards_quantization_conditions_matrix[k,:]
  zero_count = count(==(0), row)
  if zero_count != length(row)
    push!(non_zero_quantization_conditions_matrix, row)
  end
end
@req length(non_zero_quantization_conditions_matrix) == 1173 "Some intersection numbers must have changed!"


# A general vertical flux is given as \sum_i{a_i * G_i} with G_i the basis of vertical fluxes computed above and ai a rational number.
# Now express the quantiziation constraints as linear polynomials in the coefficients ai.
# We also introduce integral quantities Ni that we will use to solve the quantization constraints in a moment.
my_variable_names = vcat( [join(["a", string(k)]) for k in 1:length(non_zero_quantization_conditions_matrix[1])],
                          [join(["N", string(k)]) for k in 1:length(non_zero_quantization_conditions_matrix) + 50]);

# To simplify the code below, I introduce 50 additional variables Ni here.
# After the first round of solving the quantization constraints, there should be 16 constraints remaining.
# But of course, if I made a mistake somewhere, then this number should change. To solve those constraints,
# I want 16 "fresh" variables in the following polynomial rings. With room to spare, I thus introduce 50
# additional variables.
const R1, R1_gens = polynomial_ring(QQ, my_variable_names);
quantization_constraints_as_poly = [];
for k in 1:length(non_zero_quantization_conditions_matrix)
  push!(quantization_constraints_as_poly, sum(non_zero_quantization_conditions_matrix[k][m] * R1_gens[m] for m in 1:length(non_zero_quantization_conditions_matrix[1])))
end
const fixed_quantization_constraints_as_poly = quantization_constraints_as_poly;


# Introduce a small helper function to reduce a rational coefficient in the quantization conditions.
function reduce_numerator(r::Rational{Int64})
  n = numerator(r)
  d = denominator(r)
  if n % d == 0
    return QQ(0)
  else
    return QQ((n % d)//d)
  end
end


function reduce_numerator(r::QQFieldElem)
  n = numerator(r)
  d = denominator(r)
  if n % d == 0
    return QQ(0)
  else
    return QQ((n % d)//d)
  end
end





#################################################
# (10) First round of solving the quantization
# (10) constraints. This take a LOOONG time,
# (10) while it should be instantly. There seems
# (10) a lot of room for optimization!
# (10) Start at 15:50PM. End around 18:30PM
# (10) So takes almost 3 hours!
#################################################

new_quantization_constraints = deepcopy(fixed_quantization_constraints_as_poly);
my_replacement_dict = Dict{typeof(R1_gens[1]), typeof(R1_gens[1])}();
my_replacement_dict_integral_quantities = Dict{typeof(R1_gens[1]), typeof(R1_gens[1])}();
offset = length(non_zero_quantization_conditions_matrix[1]);

for k in 1:length(new_quantization_constraints)

  println("---------------------------------------")
  println("Working on constraint $k (of $(length(new_quantization_constraints)))")
  println("---------------------------------------\n")
  println("$(new_quantization_constraints[k])")
  println("")

  # Analyze the k-th quantization constraint, to find the rational coefficient that we want to express in terms of other rational coefficients and integral quantities
  list_of_exps = collect(exponents(new_quantization_constraints[k]));
  list_of_exps = [[my_exp[l] for l in 1:length(non_zero_quantization_conditions_matrix[1])] for my_exp in list_of_exps];
  list_of_exps = [findfirst(l -> l == 1, my_exp) for my_exp in list_of_exps];
  list_of_exps = filter!(x -> x !== nothing, list_of_exps);

  # Initialize a value for good measure, then identify the true value...
  new_linear_combination = zero(R1)
  if length(list_of_exps) == 0

    # This might become tricky...
    println("Constraint among inserted integer variables encountered\n")
    my_coeffs = collect(coefficients(new_quantization_constraints[k]))
    if all(is_integer, my_coeffs)
      println("Can be ignored\n")
      new_quantization_constraints[k] = zero(R1)
      continue
    end

    # First, remove any integer multiple of an Ni, as those do not matter.
    new_coeffs = collect(coefficients(new_quantization_constraints[k]));
    new_mons = collect(monomials(new_quantization_constraints[k]));
    new_constraint = zero(R1)
    for l in 1:length(new_coeffs)
      if is_integer(new_coeffs[l]) == false
        new_constraint += reduce_numerator(new_coeffs[l]) * new_mons[l]
      end
    end
    println("Simplified constraint: $(new_constraint)\n")
    new_quantization_constraints[k] = new_constraint;

    # Extract key data about the constraint
    new_coeffs = collect(coefficients(new_constraint));
    new_mons = collect(monomials(new_constraint));
    denominator_list = [denominator(l) for l in new_coeffs]
    numerator_list = [numerator(l) for l in new_coeffs]
    
    # Bring all coefficients to the same denominator
    denominator_list_updated = [lcm(denominator_list) for l in 1:length(denominator_list)]
    numerator_list_updated = [numerator_list[l] * ZZ(denominator_list_updated[l] // denominator_list[l]) for l in 1:length(numerator_list)]
    println("New denominators: $denominator_list_updated")
    println("New numerators: $numerator_list_updated\n")

    # Try to solve the constraint. Alternatively, at least simplify the constraint...
    if findfirst(x -> (x == 1 || x == -1), numerator_list_updated) !== nothing

      # E.g. of the form 1/2 N1 + 1/2 N2 + 1/2 N3...
      println("We can solve this constraint\n")

      my_pos = findfirst(x -> (x == 1 || x == -1), numerator_list_updated)
      var_to_be_replaced = new_mons[my_pos];
      new_linear_combination = (-1) * sum(numerator_list_updated[l] * new_mons[l] for l in 1:length(numerator_list_updated))
      new_linear_combination += numerator_list_updated[my_pos] * new_mons[my_pos] + denominator_list_updated[1] * R1_gens[offset + k]
      new_linear_combination = numerator_list_updated[my_pos] * new_linear_combination
      my_replacement_dict_integral_quantities[var_to_be_replaced] = new_linear_combination;
      var_pos_to_be_replaced = findfirst(x -> x == var_to_be_replaced, R1_gens);

    elseif length(numerator_list_updated) == 1 && length(denominator_list_updated) == 1

      # E.g. of the form 3/17 * N1 or 13/25 * N1
      println("We can solve this constraint\n")

      var_to_be_replaced = new_mons[1]
      multies = [denominator_list_updated[1]]
      for l in 1:denominator_list_updated[1] - 1
        if numerator_list_updated[1] * l % denominator_list_updated[1] == 0
          push!(multies, l)
        end
      end
      multies = minimum(multies)
      new_linear_combination = multies * R1_gens[offset + k]
      my_replacement_dict_integral_quantities[var_to_be_replaced] = new_linear_combination;
      var_pos_to_be_replaced = findfirst(x -> x == var_to_be_replaced, R1_gens);

    else

      println("We cannot (yet) solve this constraint\n")
      continue

    end
  
  else

    var_pos_to_be_replaced = minimum(list_of_exps);
    var_to_be_replaced = R1_gens[var_pos_to_be_replaced];
    coeff_of_var_to_be_replaced = collect(coefficients(new_quantization_constraints[k]))[findfirst(l -> l == var_pos_to_be_replaced, list_of_exps)];
    Q_linear_replacement_combination = new_quantization_constraints[k] - coeff_of_var_to_be_replaced * var_to_be_replaced - R1_gens[offset + k];
    new_coeffs = [(-1) * l // coeff_of_var_to_be_replaced for l in collect(coefficients(Q_linear_replacement_combination))];
    new_mons = collect(monomials(Q_linear_replacement_combination));
    new_linear_combination = sum(new_coeffs[l] * new_mons[l] for l in 1:length(new_mons));
    my_replacement_dict[var_to_be_replaced] = new_linear_combination;

  end

  # Inform what we are doing
  println("Variable to be replaced: $(var_to_be_replaced)\n")
  println("Image: $(new_linear_combination)\n")

  # Create ring map to update the constraints with the relation in question.
  images = deepcopy(R1_gens);
  @req var_pos_to_be_replaced !== nothing "Inconsistency in determining position of variable to be replaced"
  if haskey(my_replacement_dict, var_to_be_replaced)
    images[var_pos_to_be_replaced] = my_replacement_dict[var_to_be_replaced];
  else
    images[var_pos_to_be_replaced] = my_replacement_dict_integral_quantities[var_to_be_replaced];
  end
  my_ring_map = hom(R1, R1, images);

  # Update the constraints
  @showprogress 1 "Updating list of constraints..." for l in 1:length(new_quantization_constraints)
    list_of_exps = collect(exponents(new_quantization_constraints[l]));
    list_of_exps = [list_of_exps[m][var_pos_to_be_replaced] for m in 1:length(list_of_exps)];
    monomial_position = findfirst(x -> x == 1, list_of_exps)
    if monomial_position !== nothing
      new_quantization_constraints[l] = my_ring_map(new_quantization_constraints[l])
    end
  end

  # Update the replacements
  println("Updating replacement dicts...\n")
  for (key, value) in my_replacement_dict
    my_replacement_dict[key] = my_ring_map(value)
  end
  for (key, value) in my_replacement_dict_integral_quantities
    my_replacement_dict[key] = my_ring_map(value)
  end

  # Display the updated constraint
  println("New constraint $k:\n$(new_quantization_constraints[k])\n")

end


# Prepare and filter out all the remaining non-trivial quantization constraints
remaining_quantization_constraints = []
for k in 1:length(new_quantization_constraints)
  new_coeffs = collect(coefficients(new_quantization_constraints[k]))
  if all(is_integer, new_coeffs) == false
    new_mons = collect(monomials(new_quantization_constraints[k]));
    new_constraint = zero(R1)
    for l in 1:length(new_coeffs)
      if is_integer(new_coeffs[l]) == false
        new_constraint += reduce_numerator(new_coeffs[l]) * new_mons[l]
      end
    end
    push!(remaining_quantization_constraints, new_constraint)
  end
end
@req length(remaining_quantization_constraints) == 16 "Something must have changed!"





#################################################
# (11) Second and final round of solving the
# (11) quantization constraints.
#################################################

# Sccond round of solving the quantization constraints...
new_remaining_quantization_constraints = deepcopy(remaining_quantization_constraints);
offset = length(non_zero_quantization_conditions_matrix[1]) + length(non_zero_quantization_conditions_matrix);

for k in 1:length(new_remaining_quantization_constraints)

  println("---------------------------------------")
  println("Working on constraint $k (of $(length(new_remaining_quantization_constraints)))")
  println("---------------------------------------\n")
  println("$(new_remaining_quantization_constraints[k])")
  println("")

  # Analyze the k-th quantization constraint, to find the rational coefficient that we want to express in terms of other rational coefficients and integral quantities
  list_of_exps = collect(exponents(new_remaining_quantization_constraints[k]));
  list_of_exps = [[my_exp[l] for l in 1:length(non_zero_quantization_conditions_matrix[1])] for my_exp in list_of_exps];
  list_of_exps = [findfirst(l -> l == 1, my_exp) for my_exp in list_of_exps];
  list_of_exps = filter!(x -> x !== nothing, list_of_exps);

  # Initialize a value for good measure, then identify the true value...
  new_linear_combination = zero(R1)
  if length(list_of_exps) == 0

    # This might become tricky...
    println("Constraint among inserted integer variables encountered\n")
    my_coeffs = collect(coefficients(new_remaining_quantization_constraints[k]))
    if all(is_integer, my_coeffs)
      println("Can be ignored\n")
      new_remaining_quantization_constraints[k] = zero(R1)
      continue
    end

    # First, remove any integer multiple of an Ni, as those do not matter.
    new_coeffs = collect(coefficients(new_remaining_quantization_constraints[k]));
    new_mons = collect(monomials(new_remaining_quantization_constraints[k]));
    new_constraint = zero(R1)
    for l in 1:length(new_coeffs)
      if is_integer(new_coeffs[l]) == false
        new_constraint += reduce_numerator(new_coeffs[l]) * new_mons[l]
      end
    end
    println("Simplified constraint: $(new_constraint)\n")
    new_remaining_quantization_constraints[k] = new_constraint;

    # Extract key data about the constraint
    new_coeffs = collect(coefficients(new_constraint));
    new_mons = collect(monomials(new_constraint));
    denominator_list = [denominator(l) for l in new_coeffs]
    numerator_list = [numerator(l) for l in new_coeffs]
    
    # Bring all coefficients to the same denominator
    denominator_list_updated = [lcm(denominator_list) for l in 1:length(denominator_list)]
    numerator_list_updated = [numerator_list[l] * ZZ(denominator_list_updated[l] // denominator_list[l]) for l in 1:length(numerator_list)]
    println("New denominators: $denominator_list_updated")
    println("New numerators: $numerator_list_updated\n")

    # Try to solve the constraint. Alternatively, at least simplify the constraint...
    if findfirst(x -> (x == 1 || x == -1), numerator_list_updated) !== nothing

      # E.g. of the form 1/2 N1 + 1/2 N2 + 1/2 N3...
      println("We can solve this constraint\n")

      my_pos = findfirst(x -> (x == 1 || x == -1), numerator_list_updated)
      var_to_be_replaced = new_mons[my_pos];
      new_linear_combination = (-1) * sum(numerator_list_updated[l] * new_mons[l] for l in 1:length(numerator_list_updated))
      new_linear_combination += numerator_list_updated[my_pos] * new_mons[my_pos] + denominator_list_updated[1] * R1_gens[offset + k]
      new_linear_combination = numerator_list_updated[my_pos] * new_linear_combination
      my_replacement_dict_integral_quantities[var_to_be_replaced] = new_linear_combination;
      var_pos_to_be_replaced = findfirst(x -> x == var_to_be_replaced, R1_gens);

    elseif length(numerator_list_updated) == 1 && length(denominator_list_updated) == 1

      # E.g. of the form 3/17 * N1 or 13/25 * N1
      println("We can solve this constraint\n")

      var_to_be_replaced = new_mons[1]
      multies = [denominator_list_updated[1]]
      for l in 1:denominator_list_updated[1] - 1
        if numerator_list_updated[1] * l % denominator_list_updated[1] == 0
          push!(multies, l)
        end
      end
      multies = minimum(multies)
      new_linear_combination = multies * R1_gens[offset + k]
      my_replacement_dict_integral_quantities[var_to_be_replaced] = new_linear_combination;
      var_pos_to_be_replaced = findfirst(x -> x == var_to_be_replaced, R1_gens);

    else

      println("We cannot (yet) solve this constraint\n")
      continue

    end
  
  else

    var_pos_to_be_replaced = minimum(list_of_exps);
    var_to_be_replaced = R1_gens[var_pos_to_be_replaced];
    coeff_of_var_to_be_replaced = collect(coefficients(new_remaining_quantization_constraints[k]))[findfirst(l -> l == var_pos_to_be_replaced, list_of_exps)];
    Q_linear_replacement_combination = new_remaining_quantization_constraints[k] - coeff_of_var_to_be_replaced * var_to_be_replaced - R1_gens[offset + k];
    new_coeffs = [(-1) * l // coeff_of_var_to_be_replaced for l in collect(coefficients(Q_linear_replacement_combination))];
    new_mons = collect(monomials(Q_linear_replacement_combination));
    new_linear_combination = sum(new_coeffs[l] * new_mons[l] for l in 1:length(new_mons));
    my_replacement_dict[var_to_be_replaced] = new_linear_combination;

  end

  # Inform what we are doing
  println("Variable to be replaced: $(var_to_be_replaced)\n")
  println("Image: $(new_linear_combination)\n")

  # Create ring map to update the constraints with the relation in question.
  images = deepcopy(R1_gens);
  @req var_pos_to_be_replaced !== nothing "Inconsistency in determining position of variable to be replaced"
  if haskey(my_replacement_dict, var_to_be_replaced)
    images[var_pos_to_be_replaced] = my_replacement_dict[var_to_be_replaced];
  else
    images[var_pos_to_be_replaced] = my_replacement_dict_integral_quantities[var_to_be_replaced];
  end
  my_ring_map = hom(R1, R1, images);

  # Update the constraints
  @showprogress 1 "Updating list of constraints..." for l in 1:length(new_remaining_quantization_constraints)
    list_of_exps = collect(exponents(new_remaining_quantization_constraints[l]));
    list_of_exps = [list_of_exps[m][var_pos_to_be_replaced] for m in 1:length(list_of_exps)];
    monomial_position = findfirst(x -> x == 1, list_of_exps)
    if monomial_position !== nothing
      new_remaining_quantization_constraints[l] = my_ring_map(new_remaining_quantization_constraints[l])
    end
  end

  # Update the replacements
  println("Updating replacement dicts...\n")
  for (key, value) in my_replacement_dict
    my_replacement_dict[key] = my_ring_map(value)
  end
  for (key, value) in my_replacement_dict_integral_quantities
    my_replacement_dict[key] = my_ring_map(value)
  end

  # Display the updated constraint
  println("New constraint $k:\n$(new_remaining_quantization_constraints[k])\n")

end


# Ideally, my_replacement_dict should now hold the exact mappings that some of the rational coefficients ai have to undergo, so that the g4s satisfy the
# quantization condition.

my_replacement_dict = load("/home/i/1511.03209/my_replacement_dict.mrdi")
my_replacement_dict_integral_quantities = load("/home/i/1511.03209/my_replacement_dict_integral_quantities.mrdi")





#################################################
# (12) Final results
# (12) Final results
#################################################

# H^(2,2)(amb) is isomorphic to Q^2104

# Then identify all those G4 in H^(2,2)(amb) s.t. G4|_Y4 is vertical. This is a subspace of H^(2,2)(amb) of dimension 1712.
# The basis of which has been computed as the vector basis_of_vertical_g4_fluxes.

# Then we subject those G4s to the quantization condition. It turns out 224 of the rational coefficients must undergo a linear transformation.
# I.e. we have ai -> sum_{j in J}{lambda_j * a_j} + sum_{k in K}{mu_k * N_k}, where N_k in Z and a_j in Q are other rational coefficients.
# Consequently, the space of vertical G4-fluxes that satisfy the quantization condition is isomorphic to Q^(1488) \times N^{224}.

# Any vertical G4-flux A which satisfies the quantization condition can be expressed as follows:

# Pick 224 integer Ni.
# Pick 1488 rational numbers ai.
# Compute 224 rational numbers with the rules/linear transformation described in my_replacement_dict.
# Now consider sum_{i = 1}{1712}{ai * gi} with ai the rational numbers computed and gi the elements of the vector basis_of_vertical_g4_fluxes.

# Cool, huh!?!



save("/home/i/new_intersection_dict.mrdi", intersection_dict)
save("/home/i/new_special_intersection_dict.mrdi", special_intersection_dict)
save("/home/i/new_towards_quantization_conditions_matrix.mrdi", towards_quantization_conditions_matrix)
save("/home/i/basis_of_vertical_g4_fluxes.mrdi", basis_of_vertical_g4_fluxes)





#################################################
# (13) D3-tadpole
# (13) D3-tadpole
#################################################

# We want to compute G4 * G4. I believe it will be easiest to create a HUGE ring, such that we can represent a general vertical G4-flux,
# which is well quantized, as a polynomial in this ring.
# Then we compute the product of this polynomial with itself.
# It then remains to replace any quadrules of ambient space variables with their intersection number in Y4.
# Once we computed those numbers, we will get a quadratic form in the 1488 ais and 224 Nis.
# Since we know the Euler characteristic, it is then relatively easy to form the D3-tadpole constraint.
# If we can solve it, is a different matter of course...

# TO COME...



# blowup 207 renders the Tate polynomial inhomogeneous! Presumably, that is the first that I execute on top! But why???


