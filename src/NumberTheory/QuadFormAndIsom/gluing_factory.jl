###############################################################################
#
#  Accessors
#
###############################################################################

### Gluing factory

@doc raw"""
    ambient_modules(Fac::ZZLatGluingFactory) -> TorQuadModule, TorQuadModule

Return the two initial torsion quadratic modules from the factory `Fac`.
"""
ambient_modules(Fac::ZZLatGluingFactory) = Fac.ambient_modules

@doc raw"""
    local_classifying_groups(
      Fac::ZZLatGluingFactory,
    ) -> AutomorphismGroup{TorQuadModule}, AutomorphismGroup{TorQuadModule}

Return the groups used to do the computations of local gluings in the factory
`Fac`.
"""
local_classifying_groups(Fac::ZZLatGluingFactory) = Fac.local_classifying_groups

@doc raw"""
    conditions_modules(
      Fac::ZZLatGluingFactory,
    ) -> TorQuadModuleMap, TorQuadModuleMap

Return the embedding of the conditions modules in the ambient modules of
the factory `Fac`, computed via the conditions initiliazing `Fac.
"""
conditions_modules(Fac::ZZLatGluingFactory) = Fac.conditions_modules

@doc raw"""
    overlattice_parity(Fac::ZZLatGluingFactory) -> Symbol

Return the parity type of a primitive extension as setup in the construction
of the factory `Fac`.

It can be :even for even extensions, :odd for odd extensions and :both for
any integral extensions.
"""
overlattice_parity(Fac::ZZLatGluingFactory) = Fac.par

@doc raw"""
    possible_glue_order(Fac::ZZLatGluingFactor) -> Set{ZZRingElem}

Return the set of possible order for the glue groups, computed in the
initialization of the factory `Fac`.
"""
possible_glue_order(Fac::ZZLatGluingFactory) = Fac.glue_order

### Gluing

function _data(x::ZZLatGluing)
  return x.glue_map, x.inv_glue_map, x.glue_group_left, x.stabilizer_left, x.glue_group_right, x.stabilizer_right
end

function glue_map(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.glue_map
  elseif i == 2
    return x.inv_glue_map
  end
end

glue_maps(x::ZZLatGluing) = (glue_map(x, 1), glue_map(x, 2))

function glue_group(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.glue_group_left
  elseif i == 2
    return x.glue_group_right
  end
end

glue_groups(x::ZZLatGluing) = (glue_group(x, 1), glue_group(x, 2))

function stabilizer_glue_group(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.stabilizer_left
  elseif i == 2
    return x.stabilizer_right
  end
end

###############################################################################
#
#  Routine
#
###############################################################################

@doc raw"""
    test_overlattice(Fac::ZZLatGluingFactory, L::ZZLat) -> Bool

Return whether the given lattice satisfies the output conditions of
the factory `Fac`.

This means that:
- the parity of `L` is good,
- the genus of `L` is among the allowed ones (if any),
- the discriminant group of `L` is among the allowed ones (if any).
"""
function test_overlattice(Fac::ZZLatGluingFactory, L::ZZLat)
  par = overlattice_parity(Fac)
  if par === :even && !is_even(L)
    return false
  elseif par === :odd && is_even(L)
    return false
  end
  if isdefined(Fac, :genus_over)
    genus(L) in Fac.genus_over && return true
  end
  if isdefined(Fac, :form_over)
    M = gram_matrix_quadratic(first(normal_form(discriminant_group(L))))
    return M in Fac.form_over
  end
  return true
end

@doc raw"""
    is_trivial(Fac::ZZLatGluingFactory) -> Bool

Return whether the output conditions of the factory `Fac` are possible.
"""
function is_trivial(Fac::ZZLatGluingFactory)
  isempty(possible_glue_order(Fac)) && return true
  if isdefined(Fac, :genus_over)
    return isempty(Fac.genus_over)
  elseif isdefined(Fac, :form_over)
    return isempty(Fac.form_over)
  elseif isdefined(Fac, :glue_elementary_divisors)
    return isempty(Fac.glue_elementary_divisors)
  end
  return false
end

# The type ZZLatGluing is symmetrical, we just reverse the order of the fields
# to see it from right to left
function Oscar.inv(x::ZZLatGluing)
  phi, iphi, i1, j1, i2, j2 = _data(x)
  return ZZLatGluing(iphi, phi, i2, j2, i1, j1)
end

# Construct the torsion module D where we construct graphs of glue maps,
# which is just the orthogonal direct sum of q1 and q2. Compute also
# the embeddings O(q1), O(q2) -> O(D). Depending on par, q1 and q2, the
# module D is considered as a bilinear or quadratic module.
function _gluing_ambient(
    q1::TorQuadModule,
    q2::TorQuadModule,
    par::Symbol,
  )
  as_bilinear_module = par === :even ? false : true
  x = _direct_sum_with_embeddings_orthogonal_groups(q1, q2; as_bilinear_module)
  return ZZLatGluingAmbient(x...)
end

# Return all genera with given form, signature and parity conditions
# There are at most 2 of them
# par is either :even, or :odd, or :both
function _integer_genera(
  q::TorQuadModule,
  sign::NTuple{2, Int},
  par::Symbol,
)
  p, n = sign
  GKs = Set{ZZGenus}()
  if par != :odd
    ok, _G = is_genus_with_genus(q, (p, n); parity=2)
    ok && push!(GKs, _G)
  end

  if par != :even
    ok, _G = is_genus_with_genus(q, (p, n); parity=1)
    ok && push!(GKs, _G)
  end
  return GKs
end

@doc raw"""
    orthogonal_group_bilinear(T::TorQuadModule) -> AutomorphismGroup{TorQuadModule}

Return the orthogonal group for the bilinear form on ``T``, i.e. seeing ``T``
as a torsion bilinear form. If the modulus of the quadratic form on ``T`` is
equal to the modulus of its bilinear form, then this is the same as calling
`orthogonal_group(T)`.
"""
@attr AutomorphismGroup{TorQuadModule} function orthogonal_group_bilinear(T::TorQuadModule)
  return __orthogonal_group(T; as_bilinear_module=true)
end

function __as_finite_bilinear_module(
  q::TorQuadModule,
)
  n = modulus_bilinear_form(q)
  if n == modulus_quadratic_form(q)
    return q
  end
  
  qb =  torsion_quadratic_module(cover(q), relations(q); modulus=n, modulus_qf=n, gens=lift.(gens(q)))
  return qb
end

function __orthogonal_group(
  q::TorQuadModule;
  as_bilinear_module::Bool=false,
)
  if !as_bilinear_module || modulus_bilinear_form(q) == modulus_quadratic_form(q)
    return orthogonal_group(q)
  end
  qb = __as_finite_bilinear_module(q)
  return _orthogonal_group(q, _orthogonal_group_gens(qb); check=false)
end

###############################################################################
#
#  Initialization functions for ZZLatGluingFactory
#
###############################################################################

# Initialize the gluing factory with the conditions from the problem.
# By construction `Fac` should already know about the ambient modules and the
# wanted parity for the primitive extensions. If the local classifying groups
# have not been set, the function initialize them to be the orthogonal groups
# (maybe as finite bilinear modules, depending on the value of Fac.par) of the
# ambient modules.
function init_gluing_factory!(
  Fac::ZZLatGluingFactory;
  left_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  right_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_exponent::Union{Nothing, IntegerUnion}=nothing,
  glue_order::AbstractVector{T}=Int[],
  glue_elementary_divisors::Vector{Vector{ZZRingElem}}=Vector{ZZRingElem}[],
  form_over::Vector{TorQuadModule}=TorQuadModule[],
  genus_over::Vector{ZZGenus}=ZZGenus[],
) where T <: IntegerUnion
  @vprintln :ZZLatWithIsom 3 "Initialize gluing conditions"
  init_gluing_conditions!(
    Fac,
    left_glue_annihilator,
    right_glue_annihilator,
    glue_exponent,
    glue_order,
    glue_elementary_divisors,
    form_over,
  )
  assert_has_local_classifying_groups!(Fac)
  return nothing
end

# Using the glue annihilators and glue_exponent, the function sets up the
# "conditions modules", i.e. where to look for glue groups. From this, the
# fonction filters through the remaining initial conditions to determine the
# actual restrictions for the gluing computed in the factory. This list of
# actual conditions consists of:
# - a set of genera for the primitive extensions (optional: only if
#   genus_over is not empty),
# - a set of Gram matrix for the normal forms of the discriminant groups
#   of the primitive extensions (optional: only if form_over is not empty)
# - a set of elementary divisors for the glue groups (optional: only if
#   glue_elementary_divisors is not empty)
# - a set of orders for the glue groups.
#
# If genus_over (resp. form_over, glue_elementary_divisors) in input is not
# empty but it is after filtering (depending on the conditions modules and the
# input list glue_order), Fac.genus_over (resp. Fac.form_over,
# Fac.glue_elementary_divisors) is set to be the empty set and the following
# functions will abort (because of impossible conditions, Fac is "trivial").
#
# For Fac.glue_order, several things can happen:
# - if the initial conditions genus_over, form_over or glue_elementary_divisors
#   are not empty, then the function infers lists of possible orders from the
#   non-empty ones and intersect them. If moreover glue_order is not empty,
#   the previous list of infered orders is intersected with glue_order.
# - otherwise, if glue_order is non-empty, the functions filters which orders
#   could actually work in the context set up by the conditions modules.
# - otherwise, Fac.glue_order consists of all the divisors of the order of the
#   maximal abelian group embedding in both of the conditions modules.
# In the two first cases, if the final list of orders infered from the initial
# conditions is empty, Fac.glue_order is the empty set and the rest of the
# algorithm will abort as well.
function init_gluing_conditions!(
  Fac::ZZLatGluingFactory,
  left_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  right_glue_annihilator::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_exponent::Union{Nothing, IntegerUnion}=nothing,
  glue_order::AbstractVector{T}=ZZRingElem[],
  glue_elementary_divisors::Vector{Vector{ZZRingElem}}=Vector{ZZRingElem}[],
  form_over::Vector{TorQuadModule}=TorQuadModule[],
  genus_over::Vector{ZZGenus}=ZZGenus[],
) where T <: IntegerUnion

  # Remember if some conditions were initially set
  ignore_genus_over = isempty(genus_over)
  ignore_form_over = isempty(form_over)
  ignore_elem_divs = isempty(glue_elementary_divisors)
  empty_order_cond = isempty(glue_order)
  ignore_glue_order = ignore_genus_over && ignore_form_over && ignore_elem_divs && empty_order_cond

  qM, qN = ambient_modules(Fac)

  # Left glue group annihilated by the endomorphism left_glue_annihilator
  if !isnothing(left_glue_annihilator)
    @assert domain(left_glue_annihilator) === codomain(left_glue_annihilator) === qM
    VM, VMinqM = kernel(left_glue_annihilator)
  else
    VMinqM = id_hom(qM)
    VM = qM
  end

  # Right glue group annihilated by the endomorphism right_glue_annihilator
  if !isnothing(right_glue_annihilator)
    @assert domain(right_glue_annihilator) === codomain(right_glue_annihilator) === qN
    VN, VNinqN = kernel(right_glue_annihilator)
  else
    VNinqN = id_hom(qN)
    VN = qN
  end

  # If glue_exponent is a positive integer, the exponent of the glue groups
  # should divide it
  if !isnothing(glue_exponent) && glue_exponent > 0
    filter!(Base.Fix1(is_divisible_by, glue_exponent)∘last, elementary_divisors)
    VM, _VMinVM = torsion_subgroup(domain(VMinqM), e)
    VMinqM = compose(_VMinVM, VMinqM)
    VN, _VNinVN = torsion_subgroup(domain(VNinqN), e)
    VNinqN = compose(_VNinVN, VNinqN)
  end

  # This fixes the conditions modules where we do everything
  Fac.conditions_modules = (VMinqM, VNinqN)

  # For each primes of interest (i.e. the ones dividing the order of a
  # potential common subgroup, we prepare a dictionary which assign to
  # every potential list of elementary divisors, the local gluings
  elG = Hecke._maximal_common_subgroup_snf(abelian_group(VM), abelian_group(VN))
  pds = isempty(elG) ? ZZRingElem[] : prime_divisors(last(elG))
  ged = Dict{ZZRingElem, Dict{Vector{Int}, Vector{ZZLatGluing}}}(p => Dict{Vector{Int}, Vector{ZZLatGluing}}() for p in pds)
  Fac.glue_group_parent_snf = elG
  Fac.primes_of_interest = Set(pds)
  Fac.local_gluings_primary = ged

  _genus_over = Set{ZZGenus}()
  for G in genus_over
    !is_integral(G) && continue
    Fac.par === :even && !is_even(G) && continue
    Fac.par === :odd && is_even(G) && continue
    push!(_genus_over, G)
  end

  k1 = prod(elG; init=ZZ(1))
  _glue_order = Set(filter!(>(0), glue_order))
  # Since any glue group lies in the group with elementary divisor
  # elG, any good glue order condition must divide k1
  for o in _glue_order
    !is_divisible_by(k1, o) && delete!(_glue_order, o)
  end

  _glue_elem_divs = Set{Vector{ZZRingElem}}()
  # Since any glue group lies in the group with elementary divisor
  # elG, any good glue elementary divisor condition must be in smith
  # normal form and must be compatible with elG
  for v in glue_elementary_divisors
    if isempty(v)
      push!(_glue_elem_divs, v)
      continue
    end
    any(<=(0), v) && continue
    length(v) > length(elG) && continue
    flag = is_divisible_by(elG[end-length(v)+1], first(v))
    !flag && continue
    for i in 1:length(v)-1
      flag &= is_divisible_by(v[i+1], v[i])
      flag &= is_divisible_by(elG[end-i+1], v[end-i+1])
      !flag && break
    end
    if flag
      push!(_glue_elem_divs, v)
    end
  end

  _form_over = Set{QQMatrix}()
  # Any good primitive extension discriminant group condition must look
  # the discriminant group of an integral lattice with the correct parity
  for q in form_over
    !is_semi_regular(q) && continue
    !isone(modulus_bilinear_form(q)) && continue
    if Fac.par === :even
      modulus_quadratic_form(q) == 2 || continue
    elseif Fac.par === :odd
      modulus_quadratic_form(q) == 1 || continue
    else
      modulus_quadratic_form(q) in [1, 2] || continue
    end
    push!(_form_over, gram_matrix_quadratic(first(normal_form(q))))
  end

  # Each primitive extension genus condition and discriminant condition imposes
  # some condition on the order of a glue group. For this, said conditions must
  # be compatible with the initial discriminant groups of the problem: if they
  # are not, we discard them.
  # If they are, we compare the induced glue order conditions with the initial:
  #   - if there are no such initial conditions, we keep track of the induced
  #     ones;
  #   - if there are such initial conditions and the induced one is in the list,
  #     we keep track of it;
  #   - otherwise, we discard the primitive extension condition.
  k2 = order(qM)*order(qN)
  _tmp = Set{ZZRingElem}()
  for G in _genus_over
    bool, _k = divides(k2, numerator(abs(det(G))))
    if !bool # G can not occur as the genus of a primitive extension
      delete!(_genus_over, G)
      continue
    end
    bool, _k = is_square_with_sqrt(_k)
    if !bool # G can not occur as the genus of a primitive extension
      delete!(_genus_over, G)
      continue
    end
    # Either there are no order conditions on the glue groups, or _k is among
    # the possible ones, in which case we keep G. Otherwise we discard it
    if empty_order_cond || (_k in _glue_order)
      push!(_tmp, _k)
    else
      delete!(_genus_over, G)
      continue
    end
  end
  for q in _form_over
    bool, _k = divides(k2, order(q))
    if !bool # q can not occur as the discriminant group of a primitive extension
      delete!(_form_over, q)
      continue
    end
    bool, _k = is_square_with_sqrt(_k)
    if !bool # q can not occur as the discriminant group of a primitive extension
      delete!(_form_over, q)
      continue
    end
    # Either there are no order conditions on the glue groups, or _k is among
    # the possible ones, in which case we keep q. Otherwise we discard it
    if empty_order_cond || (_k in _glue_order)
      push!(_tmp, _k)
    else
      delete!(_form_over, q)
      continue
    end
  end

  if (!ignore_genus_over || !ignore_form_over)
    # _tmp is the new set of glue order conditions.
    if empty_order_cond
      union!(_glue_order, _tmp) # = _tmp since _glue_order is empty
    else
      intersect!(_glue_order, _tmp) # = _tmp since _tmp is contained in _glue_order
    end
  end

  # Each glue elementary divisors condition imposes again some condition on the
  # order of a glue group. We compare the induced glue order conditions with the
  # current set of glue order conditions:
  #   - if there are no such glue order conditions, we keep track of the
  #     induced ones;
  #   - if there are such glue order conditions and the induced one is in the
  #     list, we keep track of it;
  #   - otherwise, we discard the glue elementary divisors condition.
  if !ignore_elem_divs
    _tmp = Set{ZZRingElem}()
    for v in _glue_elem_divs
      _k = prod(v; init=ZZ(1))
      if (ignore_genus_over && ignore_form_over && empty_order_cond) || (_k in _glue_order)
        push!(_tmp, _k)
      else
        delete!(_glue_elem_divs, v)
        continue
      end
    end
    # Again _tmp is the new set of glue order conditions
    if ignore_genus_over && ignore_form_over && empty_order_cond
      union!(_glue_order, _tmp) # = _tmp since _glue_order is empty here
    else
      intersect!(_glue_order, _tmp) # = _tmp since _tmp is contained in _glue_order
    end
  end

  # If the conditions in input are all empty, meaning there are no restrictions
  # on glue order, we set the possible set of glue orders to be all the divisors
  # of k1
  if !ignore_glue_order
    Fac.glue_order = _glue_order
  else
    Fac.glue_order = Set(divisors(k1))
  end

  if !ignore_genus_over
    Fac.genus_over = _genus_over
  end

  if !ignore_form_over
    Fac.form_over = _form_over
  end

  if !ignore_elem_divs
    Fac.glue_elementary_divisors = _glue_elem_divs
  end
  return nothing
end

# We make sure that `Fac` has local classifying groups. These groups should
# contain whichever groups used for the classification of primitive extensions.
# In our case, the orthogonal groups of the initial discriminant groups are
# basically the largest we want to work with.
# Note: depending on the parity set up by the initial problem, we might need
# to consider the initial discriminant groups as torsion bilinear modules.
function assert_has_local_classifying_groups!(Fac::ZZLatGluingFactory)
  isdefined(Fac, :local_classifying_groups) && return nothing
  @vprintln :ZZLatWithIsom 3 "Compute local classifying groups"
  q1, q2 = ambient_modules(Fac)
  as_bilinear_module = (Fac.par !== :even)
  Oq1 = __orthogonal_group(q1; as_bilinear_module)
  Oq2 = __orthogonal_group(q2; as_bilinear_module)
  Fac.local_classifying_groups = (Oq1, Oq2)
  return nothing
end

###############################################################################
#
#  Local gluings
#
###############################################################################

# Return a complete set of representatives for the O-orbits of p-groups of W
# whose reverse sequence of p-valuations of elementary divisors is given by
# `subtype`. If `Fac` is decorated with a gluing context `Ctx` whose node
# `i` knows the result of this function, we fetch it from there.
# No checks are performed, so in fact here `O` should be one of the local
# classifying groups of `Fac`, and the data stored at the `i`th node of
# `Ctx` should have been computed with the same group `O`.
# If `Ctx` exists and has an `i`th node, the result of the function is
# automatically sent to the node if not already known.
function __subgroups_orbit_representatives_and_stabilizers_primary_subtype(
  Fac::ZZLatGluingFactory,
  Winq::TorQuadModuleMap,
  O::AutomorphismGroup{TorQuadModule},
  p::ZZRingElem,
  subtype::Vector{Int},
  i::Int = -1,
)
  W = domain(Winq)
  # Whether Ctx has the current module at the `i`th node
  bool = Fac.par === :even
  flag = isdefined(Fac, :Ctx) && (0 < i <= length(Fac.Ctx.modules)) && (Fac.Ctx.modules[i] === codomain(Winq))
  if flag && haskey(Fac.Ctx.orb_and_stab, (i, bool, p, subtype))
    # Retrieve known data
    subs = Fac.Ctx.orb_and_stab[(i, bool, p, subtype)]
  elseif first(subtype) == 1
    # The glue groups are elementary abelian p-groups of rank `length(subtype)`
    W = domain(Winq)
    _, WWinW = torsion_subgroup(W, p)
    WWinq = compose(WWinW, Winq)
    subs = Oscar.__subgroups_orbit_representatives_and_stabilizers_elementary(WWinq, O, p^(length(subtype)), p)
    if flag
      Fac.Ctx.orb_and_stab[(i, bool, p, subtype)] = subs
    end
  else
    # Non-elementary case
    subs = __subgroups_orbit_representatives_and_stabilizers_primary_subtype(Winq, O, p, subtype)
    if flag
      Fac.Ctx.orb_and_stab[(i, bool, p, subtype)] = subs
    end
  end
  return subs
end

# Compute all local gluings in `Fac` at the prime `p` and for
# glue groups whose reverse sequence of p-valuations of elementary
# divisors is given by `subtype`.
# Note: the factory `Fac` has some storage for such results. If they
# have already been computed, then we recover it from there and avoid
# redundant computations.
function _local_gluings_primary!(
  Fac::ZZLatGluingFactory,
  p::ZZRingElem,
  subtype::Vector{Int},
)
  # If the parity of the overlattice is not forced to be even, we treat
  # all torsion modules as bilinear modules
  as_bilinear_module = (Fac.par !== :even)
  ged = Fac.local_gluings_primary
  @assert haskey(ged, p)
  if haskey(ged[p], subtype)
    return ged[p][subtype]
  end
 
  loc_glue_p = ZZLatGluing[]
  flag1, flag2 = false, false
  i1, i2 = -1, -1
  if isdefined(Fac, :Ctx)
    # If there is a gluing context `Ctx` and the factory registered
    # a node identification with one of the initial discriminant group,
    # then we can fetch data there.
    # Note: should only be used with matching classifying groups
    vi = Fac.vertex_identification
    if vi[1] > 0
      flag1 = true
      i1 = vi[1]
    end
    if vi[2] > 0
      flag2 = true
      i2 = vi[2]
    end
  end

  O1, O2 = local_classifying_groups(Fac)
  V1inq1, V2inq2 = conditions_modules(Fac)
  V1 = domain(V1inq1)
  V2 = domain(V2inq2)

  W1, W1inV1 = primary_part(V1, p)
  W1inq1 = compose(W1inV1, V1inq1)
  W2, W2inV2 = primary_part(V2, p)
  W2inq2 = compose(W2inV2, V2inq2)

  subs1 = __subgroups_orbit_representatives_and_stabilizers_primary_subtype(Fac, W1inq1, O1, p, subtype, i1)
  subs2 = __subgroups_orbit_representatives_and_stabilizers_primary_subtype(Fac, W2inq2, O2, p, subtype, i2)

  for (H1inq1, stab1) in subs1
    H1 = domain(H1inq1)
    for (H2inq2, stab2) in subs2
      H2 = domain(H2inq2)
      ok, phi = is_anti_isometric_with_anti_isometry(H1, H2; as_bilinear_module)
      !ok && continue
      x = ZZLatGluing(phi, inv(phi), H1inq1, stab1, H2inq2, stab2)
      push!(loc_glue_p, x)
    end
  end
  ged[p][subtype] = loc_glue_p # Store known data
  return loc_glue_p
end

# Return the trivial gluing in the factory, i.e. the glue groups are trivial
function _trivial_gluing(
  Fac::ZZLatGluingFactory,
)
  q1, q2 = ambient_modules(Fac)
  s1, s2 = id_hom.(local_classifying_groups(Fac))
  z1, z1inq1 = sub(q1, TorQuadModuleElem[])
  z2, z2inq2 = sub(q2, TorQuadModuleElem[])
  phi = hom(z1, z2, zero_matrix(ZZ, 0, 0))
  iphi = hom(z2, z1, zero_matrix(ZZ, 0, 0))
  return ZZLatGluing(phi, iphi, z1inq1, s1, z2inq2, s2)
end

# Return all the glue maps from the factory `Fac` whose associated glue groups
# have order `o`.
function _local_glue_maps_ord(Fac::ZZLatGluingFactory, o::ZZRingElem)
  isone(o) && return ZZLatGluing[_trivial_gluing(Fac)]

  elG = Fac.glue_group_parent_snf
  _loc = Vector{ZZLatGluing}[]
  for (p, v) in factor(o)
    # Get all the p-types for sub-p-groups of elG of valuation v
    ptypes = Hecke._psubgroups_types(elG, p, v)
    loc_p = ZZLatGluing[]
    for subtype in ptypes
      # For each p-type, we get all potential local gluings along such
      # p-groups
      append!(loc_p, _local_gluings_primary!(Fac, p, subtype))
    end
    isempty(loc_p) && return ZZLatGluing[]
    push!(_loc, loc_p)
  end

  return merge_glue_maps(Fac, _loc)
end

# Return all the glue maps from the factory `Fac` whose associated glue groups
# have elementary divisors `v`.
function _local_glue_maps_eldiv(Fac::ZZLatGluingFactory, v::Vector{ZZRingElem})
  isempty(v) && return ZZLatGluing[_trivial_gluing(Fac)]

  pds = Fac.primes_of_interest
  _loc = Vector{ZZLatGluing}[]
  for p in pds
    # For each p, we retrieve the p-type from the given elementary divisors
    subtype = reverse!(Int[valuation(a, p) for a in v])
    filter!(!=(0), subtype)
    isempty(subtype) && continue
    loc_p = _local_gluings_primary!(Fac, p, subtype)
    isempty(loc_p) && return ZZLatGluing[]
    push!(_loc, loc_p)
  end

  return merge_glue_maps(Fac, _loc)
end

# Given a list of glue maps in `Fac`, each at pairwise distinct prime numbers,
# compute all the possible combination of global glue maps.
function merge_glue_maps(
  Fac::ZZLatGluingFactory,
  _loc::Vector{Vector{ZZLatGluing}},
)
  if isempty(_loc)
    return ZZLatGluing[]
  elseif isone(length(_loc))
    return only(_loc)
  end

  res = ZZLatGluing[]
  q1, q2 = ambient_modules(Fac)
  for x in Hecke.cartesian_product_iterator(_loc; inplace=true)
    gens1 = TorQuadModuleElem[]
    gens2 = TorQuadModuleElem[]
    # Merge generators for each prime
    for y in x
      _phi = glue_map(x, 1)
      H1 = domain(glue_group(y, 1))
      for a in gens(H1)
        push!(gens1, q1(lift(a)))
        push!(gens2, q2(lift(_phi(a))))
      end
    end
    H1, H1inq1 = sub(q1, gens1)
    H2, H2inq2 = sub(q2, gens2)
    # The orthogonal groups are direct sums of orthogonal groups
    # of each p-Sylow. So we get the global stabilizers by intersecting
    # the local stabilizers.
    j1 = stabilizer_glue_group(first(x), 1)
    j2 = stabilizer_glue_group(first(x), 2)
    for i in 2:length(_loc)
      _j1 = stabilizer_glue_group(x[i], 1)
      _j2 = stabilizer_glue_group(x[i], 2)
      stab1, jj, _ = intersect(domain(j1), domain(_j1))
      j1 = jj * j1
      stab2, jj, _ = intersect(domain(j2), domain(_j2))
      j2 = jj * j2
    end
    @hassert :ZZLatWithIsom 1 domain(j1) == first(stabilizer(codomain(j1), H1inq1))
    @hassert :ZZLatWithIsom 1 domain(j2) == first(stabilizer(codomain(j2), H2inq2))
    phi = hom(H1, H2, gens(H1), gens(H2))
    @hassert :ZZLatWithIsom 1 is_anti_isometry(phi; as_bilinear_module=(Fac.par !== :even))
    push!(res, ZZLatGluing(phi, inv(phi), H1inq1, j1, H2inq2, j2))
  end
  return res
end

# Return all the local glue maps from the factory `Fac`, following the initial
# conditions of the problem.
function local_glue_maps(
  Fac::ZZLatGluingFactory,
)
  res = ZZLatGluing[]
  @vprintln :ZZLatWithIsom 1 "Compute local glue maps"
  if isdefined(Fac, :glue_elementary_divisors)
    eldiv = Fac.glue_elementary_divisors
    for v in eldiv
      @vprintln :ZZLatWithIsom 2 "v = $v"
      append!(res, _local_glue_maps_eldiv(Fac, v))
    end
  else
    l = sort!(collect(possible_glue_order(Fac)))
    for o in l
      @vprintln :ZZLatWithIsom 2 "o = $o"
      append!(res, _local_glue_maps_ord(Fac, o))
    end
  end
  return res
end

@doc raw"""
    _local_glue_maps(
      q1::TorQuadModule,
      q2::TorQuadModule,
      parity::Symbol;
      Ctx::Union{Nothing, ZZLatGluingCtx}=nothing,
      vi::NTuple{2, Int}=(-1, -1),
      kwargs...,
    ) -> ZZLatGluingFactory, Vector{ZZLatGluing}

Return the local glue maps between `q1` and `q2`, up to the action of their
respective orthogonal groups.

The argument `parity` can be `:even`, `:odd` or `:both` depending on the kind
of glue maps one wants to compute.

In large scale computations, one can input a gluing context `Ctx` remembering
some orbit and stabilizers computations for certain `TorQuadModule` and for the
given parity. The other keyword `vi` is to refer to which nodes of `Ctx` the
groups `q1` and `q2` are attached (ignores if not in the range of `Ctx`).

The extra keyword arguments include conditions on the gluing. Currently:
  - `left_glue_annihilator::TorQuadModuleMap`: an endomorphism of `q1`;
    all glue groups inside `q1` are in its kernel. Set to be `nothing` by
    default.
  - `right_glue_annihilator::TorQuadModuleMap`: an endomorphism of `q2`;
    all glue groups inside `q2` are in its kernel. Set to be `nothing` by
    default.
  - `glue_exponent::IntegerUnion`: an integer; the exponent of the glue groups
    divide this value. Set to be `nothing` by default.
  - `glue_order::AbstractVector{IntegerUnion}`: a list of integers; the order
    of the glue groups is in this list if not empty. Set to be empty by
    default.
  - `glue_elementary_divisors::Vector{Vector{ZZRingElem}}`: a list of lists of
    integers; the elementary divisors of a glue group are in this list if not
    empty. Set to be empty by default.
  - `form_over::Vector{TorQuadModule}`: a list of torsion modules; the
    discriminant group of a primitive extension is in this list if not
    empty. Set to be empty by default.
  - `genus_over::Vector{ZZGenus}`: a list of integral genera; the genus
    of a primitive extension is in this list if not empty. Set to be empty
    by default.
These initial conditions are filtered and treated in the initialization of
a gluing factory `Fac` in which the local glue maps are computed.

The output consists of the gluing factory `Fac` setup by the initial
conditions and the list of associated local glue maps.
"""
function _local_glue_maps(
  q1::TorQuadModule,
  q2::TorQuadModule,
  parity::Symbol;
  Ctx=nothing,
  vi::NTuple{2, Int}=(-1, -1),
  kwargs...,
)
  Fac = ZZLatGluingFactory(q1, q2, parity)
  if !isnothing(Ctx)
    Fac.Ctx = Ctx
    Fac.vertex_identification = vi
  end
  init_gluing_factory!(Fac; kwargs...)
  # If the initial conditions give rise to no possible gluings, then
  # we return the empty list
  is_trivial(Fac) && return Fac, ZZLatGluing[]
  res = local_glue_maps(Fac)
  return Fac, res
end

###############################################################################
#
#  Operations on gluings
#
###############################################################################

# Split the orbit of the left group of x by the action of O1, with the unique
# glue map. Does not compute extra glue maps.
function _split_orbit_left_group(
  x::ZZLatGluing,
  O1::AutomorphismGroup{TorQuadModule},
)
  res = ZZLatGluing[]
  phi, _, H1inq1, s1, H2inq2, s2 = _data(x)
  @assert domain(O1) === codomain(H1inq1)

  H1, H2 = domain(H1inq1), domain(H2inq2)
  q1 = codomain(H1inq1)
  stab1 = domain(s1)
  Oq1 = codomain(s1)

  elq1 = elementary_divisors(q1)
  # Check if q1 is a free module over some finite ring, in which case we
  # transform Oq1 into a permutation group
  if isone(order(q1)) || elq1[1] == elq1[end]
    iso1 = isomorphism(PermGroup, Oq1)
  else
    iso1 = id_hom(Oq1)
  end

  @vprintln :ZZLatWithIsom 1 "Split orbit of glue group"
  splits = double_cosets(codomain(iso1), first(iso1(stab1)), first(iso1(O1)))
  @vprintln :ZZLatWithIsom 1 "    done: $(length(splits))"
  for _h in splits
    h = hom(iso1\(representative(_h)))
    ih = inv(h)
    # Image group
    _H1, _H1inq1 = sub(q1, elem_type(q1)[h(H1inq1(a)) for a in gens(H1)])
    # Pullback of glue map
    _phi = hom(_H1, H2, elem_type(H2)[phi(H1(lift(ih(_H1inq1(a))))) for a in gens(_H1)])
    _iphi = inv(_phi)
    # Push stabilizer
    _stab1, _ = Oscar._as_subgroup(Oq1, Oscar.GAPWrap.ConjugateSubgroup(GapObj(stab1), GapObj(Oq1(h))))
    # Intersect with the new classifying group O1
    _stab1, _, j1 = intersect(_stab1, O1)
    push!(res, ZZLatGluing(_phi, _iphi, _H1inq1, j1, H2inq2, s2))
  end
  return res
end

# Same as before but for the right. What we do: invert the gluing internally
# (at no cost), split on the left, and invert the outputs back.
function _split_orbit_right_group(
  x::ZZLatGluing,
  O2::AutomorphismGroup{TorQuadModule},
)
  return inv.(_split_orbit_left_group(inv(x), O2))
end

# For a gluing between subgroups of q1 and q2, and given an isometry
# f:q3 -> q1, pullack the glue map to q3 along f.
function _pullback_left(
  x::ZZLatGluing,
  f::TorQuadModuleMap,
)
  @vprintln :ZZLatWithIsom 1 "Pullback gluing"
  phi, iphi, i1, s1, i2, s2 = _data(x)
  @assert codomain(f) === codomain(i1)
  invf = inv(f)
  q0 = domain(f)
  # Transport the local classifying group to q0 (technically, we can allow something
  # else than O(q0) so pulling back the group is the safest thing to do)
  Oq0 = _orthogonal_group(q0, TorQuadModuleMap[f * hom(g) * invf for g in gens(codomain(s1))]; check=false)
  _, s0 = sub(Oq0, elem_type(Oq0)[Oq0(f * hom(s1(g)) * invf; check=false) for g in gens(domain(s1))])
  H0, i0 = sub(q0, TorQuadModuleElem[f\(i1(a)) for a in gens(domain(i1))])
  H0toH1 = hom(H0, domain(i1), TorQuadModuleElem[i1\(f(i0(a))) for a in gens(H0)])
  psi = H0toH1 * phi
  ipsi = inv(psi)
  return ZZLatGluing(psi, ipsi, i0, s0, i2, s2)
end

# Again, we invert, pullback and invert back.
function _pullback_right(
  x::ZZLatGluing,
  f::TorQuadModuleMap,
)
  return inv(_pullback_left(inv(x), f))
end

# Given a glue map between representatives of orbits of glue groups, we compute
# the representatives for the double cosets of glue maps between these glue groups,
# up to the action of the given classifying groups
function _all_glue_maps(
  x::ZZLatGluing,
)
  res = ZZLatGluing[]
  phi, iphi, H1inq1, s1, H2inq2, s2 = _data(x)
  q1 = codomain(H1inq1)
  q2 = codomain(H2inq2)
  H1 = domain(H1inq1)
  H2 = domain(H2inq2)
  stab1 = domain(s1)
  stab2 = domain(s2)

  OH1 = orthogonal_group(H1)
  imOH1 = elem_type(OH1)[OH1(restrict_automorphism(x, H1inq1); check=false) for x in gens(stab1)]
  act1 = hom(stab1, OH1, imOH1)
  im1, _ = image(act1)
  OH2 = orthogonal_group(H2)
  imOH2 = elem_type(OH2)[OH2(restrict_automorphism(x, H2inq2); check=false) for x in gens(stab2)]
  act2 = hom(stab2, OH2, imOH2)
  im2, _ = image(act2)

  elH2 = elementary_divisors(H2)
  if isone(order(H2)) || elH2[1] == elH2[end]
    iso2 = isomorphism(PermGroup, OH2)
  else
    iso2 = id_hom(OH2)
  end
  OH22 = first(iso2(OH2))
  im22 = first(iso2(im2))

  stab1gamma = elem_type(OH2)[OH2(iphi * hom(g) * phi; check=false) for g in gens(im1)]
  S1, _ = sub(OH2, stab1gamma)
  S2 = first(iso2(S1))
  for _g in double_cosets(OH22, S2, im22)
    g = hom(iso2\representative(_g))
    phig = compose(phi, g)
    iphig = compose(inv(g), iphi)
    push!(res, ZZLatGluing(phig, iphig, H1inq1, s1, H2inq2, s2))
  end
  return res
end

# Given a glue map, return a representative for the isometry class of the
# discriminant group of the associated primitive extension.
function _form_over(
  x::ZZLatGluing,
)
  phi = glue_map(x, 1)
  i1, i2 = glue_groups(x)
  H1, H2 = domain(i1), domain(i2)
  q1, q2 = codomain(i1), codomain(i2)
  D, (j1, j2) = direct_sum(q1, q2; cached=false, as_bilinear_module=true)
  H1inD = i1 * j1
  H2inD = i2 * j2
  _glue = Vector{QQFieldElem}[lift(H1inD(a)) + lift(H2inD(phi(a))) for a in gens(H1)]
  ext, _ = sub(D, D.(_glue))
  perp, _ = orthogonal_submodule(D, ext)
  disc = torsion_quadratic_module(cover(perp), cover(ext))
  return disc
end

function _overlattice_with_glue_stabilizer(
  x::ZZLatGluing,
  D::ZZLatGluingAmbient,
)
  phi, iphi, i1, s1, i2, s2 = _data(x)
  H1, H2 = domain(i1), domain(i2)
  q1, q2 = codomain(i1), codomain(i2)
  H1inD = i1 * D.j1
  H2inD = i2 * D.j2
  O1, O2 = orthogonal_group(H1), orthogonal_group(H2)
  im1 = elem_type(O1)[O1(restrict_automorphism(g, i1); check=false) for g in gens(domain(s1))]
  act1 = hom(domain(s1), O1, im1)
  im2 = elem_type(O2)[O2(restrict_automorphism(g, i2); check=false) for g in gens(domain(s2))]
  act2 = hom(domain(s2), O2, im2)
  V, extinD = _overlattice_with_graph(phi, H1inD, H2inD)
  @vprintln :ZZLatWithIsom 1 "Glue stabilizers"
  _disc, _stab = _glue_stabilizers(phi, act1, act2, D.k1, D.k2, extinD)
  qV = discriminant_group(V)
  OqV = orthogonal_group(qV)
  # qV and _disc are equal, but defined in julia with different base abelian group
  # so we cannot treat them as the same object. However we can identity them
  # with a canonical "identity map" phi2
  phi2 = hom(qV, _disc, elem_type(_disc)[_disc(lift(x)) for x in gens(qV)])
  @hassert :ZZLatWithIsom 1 is_isometry(phi2; as_bilinear_module=!is_even(V))
  iphi2 = inv(phi2)
  # Pullback the glue of the stabilizers along the canonical map phi2
  GV, _ = sub(OqV, elem_type(OqV)[OqV(phi2 * g * iphi2; check=false) for g in _stab])
  # This is where we require that q1, q2 are discriminant groups of some
  # lattices and phi indeed defines a primitive extension between the
  # respective relations
  M, T = _as_sublattices(V, relations(q1), relations(q2))
  return V, M, T, GV
end

###############################################################################
#
#  Primitive extensions
#
###############################################################################

### Unimodular case

function unimodular_primitive_extensions(
  M::ZZLat,
  N::ZZLat;
  right_action::Union{MatGroup{QQFieldElem, QQMatrix}, Nothing}=nothing,
  right_discriminant_action::Union{AutomorphismGroup{TorQuadModule}, Nothing}=nothing,
  left_action::Union{MatGroup{QQFieldElem, QQMatrix}, Nothing}=nothing,
  left_discriminant_action::Union{AutomorphismGroup{TorQuadModule}, Nothing}=nothing,
  first::Bool=false,
  exist_only::Bool=false,
  even::Bool=(is_even(M) && is_even(N)), #low priority
  parity::Symbol=(even ? :even : :both),
)
  @req is_integral(M) && is_integral(N) "Only available for integral lattices"

  if !isnothing(right_discriminant_action)
    GMbar = right_discriminant_action
  elseif !isnothing(right_action)
    GMbar, _ = image(discriminant_representation(M, right_action; full=false))
  elseif first || exist_only
    qM = discriminant_group(M)
    GMbar = Oscar._orthogonal_group(qM, TorQuadModuleMap[id_hom(qM)]; check=false)
  else
    GMbar, _ = image_in_Oq(M)
  end

  if !isnothing(left_discriminant_action)
    GNbar = left_discriminant_action
  elseif !isnothing(left_action)
    GNbar, _ = image(discriminant_representation(N, left_action; full=false))
  elseif first || exist_only
    qN = discriminant_group(N)
    GNbar = Oscar._orthogonal_group(qN, TorQuadModuleMap[id_hom(qN)]; check=false)
  else
    GNbar, _ = image_in_Oq(N)
  end

  if !is_even(M) || !is_even(N)
    (parity == :even) && NTuple{3, ZZLat}[]
  end
  return _unimodular_primitive_extensions(M, N, parity, GMbar, GNbar; first, exist_only)
end

function _unimodular_primitive_extensions(
  M::ZZLat,
  N::ZZLat,
  parity::Symbol,
  GM::AutomorphismGroup{TorQuadModule},
  GN::AutomorphismGroup{TorQuadModule};
  first::Bool=false,
  exist_only::Bool=false,
)
  results = NTuple{3, ZZLat}[]
  qM = discriminant_group(M)
  @assert qM === domain(GM)
  qN = discriminant_group(N)
  @assert qN === domain(GN)

  parity == :even && (!is_even(M) || !is_even(N)) && return false, results

  as_bilinear_module = (parity !== :even)
  ok, phi = is_anti_isometric_with_anti_isometry(qM, qN; as_bilinear_module)
  !ok && return false, results
  if exist_only
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      return true, results
    end
  end

  D, inj = direct_sum(qM, qN; cached=false, as_bilinear_module)
  qMinD, qNinD = inj
  if first
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      L, _ = _overlattice_with_graph(phi, qMinD, qNinD)
      M2, N2 = _as_sublattices(L, M, N)
      push!(results, (L, M2, N2))
      return true, results
    end
  end

  iphi = inv(phi)
  if parity == :even
    OqN = orthogonal_group(qN)
  else
    OqN = orthogonal_group_bilinear(qN)
  end

  genC = elem_type(OqN)[OqN(iphi * hom(g) * phi; check=false) for g in gens(GM)]
  GMphi, _ = sub(OqN, genC)
  elqN = elementary_divisors(qN)
  if isone(order(qN)) || elqN[1] == elqN[end]
    iso = isomorphism(PermGroup, OqN)
  else
    iso = id_hom(OqN)
  end
  reps = double_cosets(codomain(iso), Base.first(iso(GN)), Base.first(iso(GMphi)))
  for _g in reps
    g = iso\(representative(_g))
    phig = compose(phi, hom(g))
    L, _ = _overlattice_with_graph(phig, qMinD, qNinD)
    if parity == :even
      is_even(L) || continue
    elseif parity == :odd
      !is_even(L) || continue
    end
    exist_only && return true, results
    M2, N2 = _as_sublattices(L, M, N)
    push!(results, (L, M2, N2))
    first && return true, results
  end
  return length(results) > 0, results
end

function _primitive_extensions_coprime_left(
  M::ZZLat,
  N::ZZLat,
  parity::Symbol,
  GM::AutomorphismGroup{TorQuadModule},
  GN::AutomorphismGroup{TorQuadModule};
  first::Bool=false,
  exist_only::Bool=false,
)
  results = NTuple{3, ZZLat}[]
  qM = discriminant_group(M)
  @assert qM === domain(GM)
  qN = discriminant_group(N)
  @assert qN === domain(GN)

  pds = prime_divisors(order(qM))
  gensHN = TorQuadModuleElem[]
  for p in pds
    _H, j = primary_part(qN, p)
    append!(gensHN, j.(gens(_H)))
  end
  HN, HNinqN = sub(qN, gensHN)

  as_bilinear_module = (parity !== :even)
  ok, phi = is_anti_isometric_with_anti_isometry(qM, HN; as_bilinear_module)
  !ok && return false, results
  if exist_only
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      return true, results
    end
  end

  D, inj = direct_sum(qM, qN; cached=false, as_bilinear_module)
  qMinD, qNinD = inj
  HNinD = compose(HNinqN, qNinD)
  if first
    if parity != :odd || !is_even(M) || !is_even(N) || !is_anti_isometry(phi; as_bilinear_module=false)
      L, _ = _overlattice_with_graph(phi, qMinD, HNinD)
      M2, N2 = _as_sublattices(L, M, N)
      push!(results, (L, M2, N2))
      return true, results
    end
  end

  iphi = inv(phi)
  if parity == :even
    is_even(M) && is_even(N) || return results
    OHN = orthogonal_group(HN)
  else
    OHN = orthogonal_group_bilinear(HN)
  end
  
  imHN, resHN = restrict_automorphism_group(GN, HNinqN)
  genC = elem_type(OHN)[OHN(iphi * hom(g) * phi; check=false) for g in gens(GM)]
  GMphi, _ = sub(OHN, genC)
  elHN = elementary_divisors(HN)
  if isone(order(HN)) || elHN[1] == elHN[end]
    iso = isomorphism(PermGroup, OHN)
  else
    iso = id_hom(OHN)
  end
  reps = double_cosets(codomain(iso), Base.first(iso(imHN)), Base.first(iso(GMphi)))
  for _g in reps
    g = iso\(representative(_g))
    phig = compose(phi, hom(g))
    L, _ = _overlattice_with_graph(phig, qMinD, HNinD)
    if parity == :even
      is_even(L) || continue
    elseif parity == :odd
      !is_even(L) || continue
    end
    exist_only && return true, results
    M2, N2 = _as_sublattices(L, M, N)
    push!(results, (L, M2, N2))
    first && return true, results
  end
  return length(results) > 0, results
end

function _primitive_extensions_coprime_right(
  M::ZZLat,
  N::ZZLat,
  parity::Symbol,
  GM::AutomorphismGroup{TorQuadModule},
  GN::AutomorphismGroup{TorQuadModule};
  kwargs...
)
  r = _primitive_extensions_coprime_left(N, M, parity, GN, GM; kwargs...)
  t = [1,3,2]
  return eltype(r)[x[t] for x in r]
end

###############################################################################
#
#  Primitive embeddings
#
###############################################################################

### Unimodular top lattice

function _primitive_embeddings_in_unimodular(
  G1::ZZGenus,
  G2s::Vector{ZZGenus};
  check::Bool=true,
)
  @assert is_unimodular(G1)
  results = Vector{ZZLat, ZZLat, ZZLat}[]

  G2s = filter(<(rank(G1))∘rank, G2s)
  filter!(<(signature_pair(G1))∘signature_pair, G2s)
  if is_even(G1)
    filter!(is_even, G2s)
  end

  for G2 in G2s
    append!(results, _primitive_embeddings_in_unimodular_safe(G1, G2))
    GC.gc()
  end
  return results
end

function _primitive_embeddings_in_unimodular(
  G1::ZZGenus,
  G2::ZZGenus,
)
  @assert is_unimodular(G1)
  if is_even(G1) && !is_even(G2)
    return results
  elseif rank(G1) <= rank(G2)
    return results
  elseif !(signature(G2) < signature(G1))
    return results
  end
  return _primitive_embeddings_in_unimodular_safe(G1, G2)
end

function _primitive_embeddings_in_unimodular_safe(
  G1::ZZGenus,
  G2::ZZGenus,
)
  results = NTuple{3, ZZLat}[]
  if is_even(G1)
    parity = :even
    par = :even
  else
    parity = :odd
    par = :both
  end
  sign = signature_pair(G1) .- signature_pair(G2)
  q = rescale(discriminant_group(G2), -1; cached=false)
  GKs = _integer_genera(q, sign, par)
  isempty(GKs) && return results
 
  Ns = ZZLat[]
  for GK in GKs
    append!(Ns, representatives(GK))
  end
  Ms = representatives(G2)
  for M in Ms, N in Ns
    append!(results, unimodular_primitive_extensions(M, N; parity))
  end
  return results
end

### Coprime det

function _primitive_embeddings_coprime_det_safe(
  G1::ZZGenus,
  G2::ZZGenus,
)
  results = NTuple{3, ZZLat}[]
  if is_even(G1)
    parity = :even
    par = :even
  else
    parity = :odd
    par = :both
  end
  sign = signature_pair(G1) .- signature_pair(G2)
  q1 = discriminant_group(G1)
  q2 = discriminant_group(G2)
  q, _ = direct_sum(q1, rescale(q2, -1; cached=false); cached=false, as_bilinear_module=(par !== :even))
  GKs = _integer_genera(q, sign, par)
  isempty(GKs) && return results

  Ns = ZZLat[]
  for GK in GKs
    append!(Ns, representatives(GK))
  end
  Ms = representatives(G2)
  for M in Ms
    GM, _ = image_in_Oq(M)
    for N in Ns
      GN, _ = image_in_Oq(N)
      tmp = _primitive_extensions_coprime_left(M, N, parity, GM, GN)
      append!(results, tmp)
    end
  end
  return results
end

### General case

function _primitive_embeddings(
  G1s::Vector{ZZGenus},
  G2s::Vector{ZZGenus},
)
  res = NTuple{3, ZZLat}[]
  for G2 in G2s
    append!(res, _primitive_embeddings(G1s, G2))
  end
  return res
end

function _primitive_embeddings(
  G1s::Vector{ZZGenus},
  G2::ZZGenus,
)
  results = NTuple{3, ZZLat}[]
  G1s = filter(>(rank(G2))∘rank, G1s)
  filter!(>(signature_pair(G2))∘signature_pair, G1s)
  if !is_even(G2)
    filter!(!is_even, G1s)
  end
  q2 = discriminant_group(G2)
  Ctx = ZZLatGluingCtx()
  push!(Ctx.modules, q2)

  for G1 in G1s
    if is_unimodular(G1)
      append!(results, _primitive_embeddings_in_unimodular_safe(G1, G2))
    elseif isone(gcd(numerator(det(G1)), numerator(det(G2))))
      append!(results, _primitive_embeddings_coprime_det_safe(G1, G2))
    else
      append!(results, _primitive_embeddings_generic_safe(G1, G2; Ctx, vi=(-1, 1), q2))
    end
  end
  return results
end

function _primitive_embeddings(
  G1::ZZGenus,
  G2::ZZGenus,
)
  results = NTuple{3, ZZLat}[]
  if is_even(G1) && !is_even(G2)
    return results
  elseif rank(G1) <= rank(G2)
    return results
  elseif !(signature_pair(G2) < signature_pair(G1))
    return results
  end
  if is_unimodular(G1)
    return _primitive_embeddings_in_unimodular_safe(G1, G2)
  elseif isone(gcd(numerator(det(G1)), numerator(det(G2))))
    return _primitive_embeddings_coprime_det_safe(G1, G2)
  else
    return _primitive_embeddings_generic_safe(G1, G2)
  end
end

function _primitive_embeddings_generic_safe(
  G1::ZZGenus,
  G2::ZZGenus;
  Ctx=nothing,
  vi::NTuple{2, Int}=(-1, -1),
  q2::TorQuadModule=discriminant_group(G2),
  Ms::Vector{ZZLat}=ZZLat[],
)
  results = NTuple{3, ZZLat}[]
  R = rescale(representative(G1), -1; cached=false)
  U = hyperbolic_plane_lattice()
  if is_even(G1)
    parity = :even
    T, _ = direct_sum(R, U; cached=false)
  else
    parity = :both
    T, _ = direct_sum(R, U, U; cached=false)
  end
  q1n = discriminant_group(T)
  signK = signature_pair(G1) .- signature_pair(G2)

  Fac, lgm = _local_glue_maps(q2, q1n, parity; Ctx, vi=reverse(vi))
  isempty(lgm) && return results
  DKs = Dict{ZZGenus, Vector{ZZLat}}()
  
  for x in lgm
    qx = _form_over(x)
    qK = rescale(qx, -1; cached=false)
    GKs = _integer_genera(qK, signK, parity)
    isempty(GKs) && continue
    # At that point, we know that we have an extension
    for GK in GKs
      haskey(DKs, GK) && continue
      DKs[GK] = representatives(GK)
    end
    isempty(Ms) && append!(Ms, representatives(G2))
    for M in Ms
      qM = discriminant_group(M)
      D = _gluing_ambient(qM, q1n, parity)
      GM, _ = image_in_Oq(M)
      _ok, phiM = is_isometric_with_isometry(qM, q2)
      @assert _ok
      xM = _pullback_left(x, phiM)
      xMs = _split_orbit_left_group(xM, GM)
      for _x in xMs, y in _all_glue_maps(_x)
        V, M2, T2, GV = _overlattice_with_glue_stabilizer(y, D)
        resV = NTuple{3, ZZLat}[]
        qV = domain(GV)
        for GK in GKs, K in DKs[GK]
          _, pe = unimodular_primitive_extensions(V, K; right_discriminant_action=GV, parity)
          append!(resV, pe)
        end
        for (S, V2, W2) in resV
          T3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(T2), zero_matrix(QQ, rank(T2), degree(W2)-degree(T2))))
          L = lll(orthogonal_submodule(S, T3))
          genus(L) == G1 || continue
          M3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(M2), zero_matrix(QQ, rank(M2), degree(W2)-degree(M2))))
          @assert is_sublattice(L, M3)
          @assert is_primitive(L, M3)
          N = orthogonal_submodule(L, M3)
          bM = coordinates(basis_matrix(M3), L)
          bN = coordinates(basis_matrix(N), L)
          L = integer_lattice(; gram=gram_matrix(L))
          M3 = lattice_in_same_ambient_space(L, bM)
          N = lattice_in_same_ambient_space(L, bN)
          push!(results, (L, M3, N))
        end
      end
    end
  end
  return results
end

###############################################################################
#
#  Overlattices and embeddings
#
###############################################################################

function _overlattices_orbit_representatives_and_stabilizers(
  M::ZZLat;
  even::Bool=is_even(M),
  GM::AutomorphismGroup{TorQuadModule} = first(image_in_Oq(M)),
)
  even && !is_even(M) && return Tuple{ZZLat, AutomorphismGroup{TorQuadModule}}[]
  @assert domain(GM) === discriminant_group(M)
  res = Tuple{ZZLat, AutomorphismGroup{TorQuadModule}}[(M, GM)]
  o = abs(numerator(det(M)))
  if radical(o) == o
    return res
  end
  pds = prime_divisors(o)
  for p in pds, i in 1:length(res)
    L, GL = res[i]
    append!(res, _prime_power_overlattices_with_stabilizers(L, p, GL; even))
  end
  return res
end

function _prime_power_overlattices_with_stabilizers(
  M::ZZLat,
  p::IntegerUnion,
  GM::AutomorphismGroup{TorQuadModule};
  even::Bool=is_even(M),
)
  res = Tuple{ZZLat, AutomorphismGroup{TorQuadModule}}[]
  qM = discriminant_group(M)
  o = order(qM)
  v = valuation(o, p)
  if v <= 1
    return res
  end
  v = div(v, 2)
  if !even
    q = Hecke._as_finite_bilinear_module(qM)
  else
    q = qM
  end
  Gq = Oscar._orthogonal_group(q, matrix.(gens(GM)))
  _, H0inq = sub(q, elem_type(q)[])
  tmp = Tuple{TorQuadModule, TorQuadModuleMap, AutomorphismGroup{TorQuadModule}}[(q, H0inq, Gq)]
  for i in 1:v
    l = length(tmp)
    for j in 1:length(tmp)
      l -= 1
      (T, Hinq, G) = popfirst!(tmp)
      genq = q.(lift.(gens(domain(Hinq))))
      Tp, TpinT = torsion_subgroup(T, p)
      r = Tuple{TorQuadModuleMap, GAPGroupHomomorphism}[]
      for k in 1:valuation(order(Tp), p)
        append!(r, Oscar._stabilizer_isotropic_elementary(Tp, k))
      end
      is_empty(r) && continue
      OTp = codomain(last(first(r)))
      iso = isomorphism(PermGroup, OTp)
      GTp, restr = restrict_automorphism_group(G, TpinT; check=false)
      isoTp, _ = iso(GTp)
      for (_HinTp, j) in r
        _H = domain(_HinTp)
        for _h in double_cosets(codomain(iso), first(iso(domain(j))), isoTp)
          h = iso\(representative(_h))
          _, H2inq = sub(q, union!([q(lift(h(Tp(lift(a))))) for a in gens(_H)], genq))
          if valuation(last(elementary_divisors(domain(H2inq))), p) != i
            continue
          end
          flag = true
          if i != 1
            for (_, H3inq, _) in tmp[l+1:end]
              elementary_divisors(domain(H3inq)) != elementary_divisors(domain(H2inq)) && continue
              if first(is_conjugate_with_data(Gq, H2inq, H3inq))
                flag = false
                break
              end
            end
          end
          !flag && continue
          s, _ = stabilizer(Gq, H2inq)
          L = cover(domain(H2inq))
          mqL = even && is_even(L) ? 2 : 1
          qL = torsion_quadratic_module(dual(L), L; modulus=1, modulus_qf=mqL)
          s = Oscar._orthogonal_group(qL, [hom(qL, qL, [qL(lift(g(q(lift(a))))) for a in gens(qL)]) for g in gens(s)])
          push!(res, (L, s))
          push!(tmp, (qL, H2inq, s))
       end
      end
    end
  end
  return res
end

################################################################################
##
##  Temporary
##
################################################################################
#
## To connect to the LMFDB database:
## - fork https://github.com/thofma/LMFDB.jl
## - run `using LMFDB`
## - fetch the database: `db = LMFDB.LMFDBLite.LMFDBConnection()`
## - load some genera: z.B. `l6=LMFDB.genera(db; rank = 6, nplus=6, det = >=(1) & <=(div(300, 8, RoundUp)));`
#
#function _primitive_embeddings_lmfdb(
#  G1s::Vector{ZZGenus},
#  G2s::Vector{ZZGenus};
#  _save_path::Union{Nothing, String}=nothing,
#)
#  results = Tuple{String, String, ZZMatrix}[]
#  for G2 in G2s
#    @show G2
#    append!(results, _primitive_embeddings_lmfdb(G1s, G2; _save_path))
#    #_primitive_embeddings_lmfdb(G1s, G2, _save_path)
#    GC.gc()
#  end
#  return results
#  #return nothing
#end
#
#function _primitive_embeddings_lmfdb(
#  G1s::Vector{ZZGenus},
#  G2::ZZGenus;
#  _save_path::Union{Nothing, String}=nothing,
#)
#  results = Tuple{String, String, ZZMatrix}[]
#  G1s = filter(>(rank(G2))∘rank, G1s)
#  filter!(>(signature_pair(G2))∘signature_pair, G1s)
#  if !is_even(G2)
#    filter!(!is_even, G1s)
#  end
#  q2 = discriminant_group(G2)
#  Ctx = ZZLatGluingCtx()
#  push!(Ctx.modules, q2)
#
#  for G1 in G1s
#    @show G1
#    if is_unimodular(G1)
#      append!(results, _primitive_embeddings_in_unimodular_safe_lmfdb(G1, G2; _save_path))
#      #_primitive_embeddings_in_unimodular_safe_lmfdb(G1, G2, _save_path)
#    elseif isone(gcd(numerator(det(G1)), numerator(det(G2))))
#      append!(results, _primitive_embeddings_coprime_det_safe_lmfdb(G1, G2; _save_path))
#      #_primitive_embeddings_coprime_det_safe_lmfdb(G1, G2, _save_path)
#    else
#      append!(results, _primitive_embeddings_generic_safe_lmfdb(G1, G2; Ctx, vi=(-1, 1), q2, _save_path))
#      #_primitive_embeddings_generic_safe_lmfdb(G1, G2, _save_path; Ctx, vi=(-1, 1), q2)
#    end
#  end
#  return results
#  #return nothing
#end
#
#function save_data_lmfdb(
#  _save_path::String,
#  x::Tuple{String, String, ZZMatrix},
#)
#  str = x[1]*"|"*x[2]*"|["
#  v = x[3]
#  for j in 1:ncols(v)
#    k = v[1,j]
#    str *= "$k"
#    if j != ncols(v)
#      str *= ","
#    end
#  end
#  str *= "]\n"
#  f = open(_save_path, "a")
#  write(f, str)
#  close(f)
#end
#
#function save_data_lmfdb(
#  _save_path::String,
#  x::Tuple{String, String, Hecke.IntegerUnion, Vector{QQFieldElem}},
#)
#  str = x[1]*"|"*x[2]*"|$(x[3])|["
#  v = x[4]
#  for j in 1:length(v)
#    k = v[j]
#    str *= "$(numerator(k))"
#    if j != length(v)
#      str *= ","
#    end
#  end
#  str *= "]\n"
#  f = open(_save_path, "a")
#  write(f, str)
#  close(f)
#end
#
#function save_data_lmfdb(
#  _save_path::String,
#  x::Tuple{String, Int},
#)
#  str = x[1]*"|$(x[2])\n"
#  f = open(_save_path, "a")
#  write(f, str)
#  close(f)
#end
#
#function _primitive_embeddings_in_unimodular_safe_lmfdb(
#  G1::ZZGenus,
#  G2::ZZGenus;
#  _save_path::Union{Nothing, String}=nothing,
#)
#  results = Tuple{String, String, ZZMatrix}[]
#  if is_even(G1)
#    parity = :even
#    par = :even
#  else
#    parity = :odd
#    par = :both
#  end
#  sign = signature_pair(G1) .- signature_pair(G2)
#  q = rescale(discriminant_group(G2), -1; cached=false)
#  GKs = _integer_genera(q, sign, par)
#  isempty(GKs) && return results
#  #isempty(GKs) && return nothing
#
#  Ns = ZZLat[]
#  for GK in GKs
#    append!(Ns, representatives(GK))
#  end
#  Ms = get_attribute(G2, :representatives)
#  for M in Ms
#    label_bottom = get_attribute(M, :lmfdb_label)
#    for N in Ns
#      _, tmp = unimodular_primitive_extensions(M, N; parity)
#      @show length(tmp)
#      for (T, M2, N2) in tmp
#        @assert gram_matrix(M) == gram_matrix(M2)
#        genus(T) != G1 && continue #Should never happen but still
#        for T2 in get_attribute(G1, :representatives)
#          if is_definite(G1)
#            ok, f_top = is_isometric_with_isometry(T, T2)
#          elseif length(get_attribute(G1, :representatives)) == 1
#            ok = true
#          else
#            ok = is_isometric(T, T2)
#          end
#          !ok && continue
#          label_top = get_attribute(T2, :lmfdb_label)
#          if !is_definite(G1)
#            v = zero_matrix(ZZ, 0, 0)
#          else
#            v = map_entries(ZZ, coordinates(basis_matrix(N2), T)*f_top)
#          end
#          __x = (label_top, label_bottom, v)
#          @show __x
#          push!(results, __x)
#          if !isnothing(_save_path)
#            save_data_lmfdb(_save_path, __x)
#          end
#        end
#      end
#    end
#  end
#  return results
#  #return nothing
#end
#
#### Coprime det
#
#function _primitive_embeddings_coprime_det_safe_lmfdb(
#  G1::ZZGenus,
#  G2::ZZGenus;
#  _save_path::Union{Nothing, String}=nothing,
#)
#  results = Tuple{String, String, ZZMatrix}[]
#  if is_even(G1)
#    parity = :even
#    par = :even
#  else
#    parity = :odd
#    par = :both
#  end
#  sign = signature_pair(G1) .- signature_pair(G2)
#  q1 = discriminant_group(G1)
#  q2 = discriminant_group(G2)
#  q, _ = direct_sum(q1, rescale(q2, -1; cached=false); cached=false, as_bilinear_module=(par !== :even))
#  GKs = _integer_genera(q, sign, par)
#  isempty(GKs) && return results
#  #isempty(GKs) && return nothing
#
#  Ns = ZZLat[]
#  for GK in GKs
#    append!(Ns, representatives(GK))
#  end
#  Ms = get_attribute(G2, :representatives)
#  for M in Ms
#    GM, _ = image_in_Oq(M)
#    label_bottom = get_attribute(M, :lmfdb_label)
#    for N in Ns
#      GN, _ = image_in_Oq(N)
#      _, tmp = _primitive_extensions_coprime_left(M, N, parity, GM, GN)
#      @show length(tmp)
#      for (T, M2, N2) in tmp
#        @assert gram_matrix(M) == gram_matrix(M2)
#        genus(T) != G1 && continue #Should never happen but still
#        for T2 in get_attribute(G1, :representatives)
#          if is_definite(G1)
#            ok, f_top = is_isometric_with_isometry(T, T2)
#          elseif length(get_attribute(G1, :representatives)) == 1
#            ok = true
#          else
#            ok = is_isometric(T, T2)
#          end
#          !ok && continue
#          label_top = get_attribute(T2, :lmfdb_label)
#          if !is_definite(G1)
#            v = zero_matrix(ZZ, 0, 0)
#          else
#            v = map_entries(ZZ, coordinates(basis_matrix(N2), T)*f_top)
#          end
#          __x = (label_top, label_bottom, v)
#          @show __x
#          push!(results, __x)
#          if !isnothing(_save_path)
#            save_data_lmfdb(_save_path, __x)
#          end
#        end
#      end
#    end
#  end
#  return results
#  #return nothing
#end
#
#
#function _primitive_embeddings_generic_safe_lmfdb(
#  G1::ZZGenus,
#  G2::ZZGenus;
#  _save_path::Union{Nothing, String}=nothing,
#  Ctx=nothing,
#  vi::NTuple{2, Int}=(-1, -1),
#  q2::TorQuadModule=discriminant_group(G2),
#  Ms::Vector{ZZLat}=ZZLat[],
#)
#  results = Tuple{String, String, ZZMatrix}[]
#  R = rescale(representative(G1), -1; cached=false)
#  U = hyperbolic_plane_lattice()
#  if is_even(G1)
#    parity = :even
#    T, _ = direct_sum(R, U; cached=false)
#  else
#    parity = :both
#    T, _ = direct_sum(R, U, U; cached=false)
#  end
#  q1n = discriminant_group(T)
#  signK = signature_pair(G1) .- signature_pair(G2)
#
#  Fac, lgm = _local_glue_maps(q2, q1n, parity; Ctx, vi=reverse(vi))
#  isempty(lgm) && return results
#  #isempty(lgm) && return nothing
#  @show length(lgm)
#  DKs = Dict{ZZGenus, Vector{ZZLat}}()
#
#  for x in lgm
#    qx = _form_over(x)
#    qK = rescale(qx, -1; cached=false)
#    GKs = _integer_genera(qK, signK, parity)
#    isempty(GKs) && continue
#    # At that point, we know that we have an extension
#    for GK in GKs
#      haskey(DKs, GK) && continue
#      DKs[GK] = representatives(GK)
#    end
#    isempty(Ms) && append!(Ms, get_attribute(G2, :representatives))
#    for M in Ms
#      label_bottom = get_attribute(M, :lmfdb_label)
#      qM = discriminant_group(M)
#      D = _gluing_ambient(qM, q1n, parity)
#      GM, _ = image_in_Oq(M)
#      _ok, phiM = is_isometric_with_isometry(qM, q2)
#      @assert _ok
#      xM = _pullback_left(x, phiM)
#      xMs = _split_orbit_left_group(xM, GM)
#      for _x in xMs, y in _all_glue_maps(_x)
#        V, M2, T2, GV = _overlattice(y, D)
#        resV = NTuple{3, ZZLat}[]
#        qV = domain(GV)
#        for GK in GKs, K in DKs[GK]
#          _, pe = unimodular_primitive_extensions(V, K; right_discriminant_action=GV, parity)
#          append!(resV, pe)
#        end
#        @show length(resV)
#        for (S, V2, W2) in resV
#          T3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(T2), zero_matrix(QQ, rank(T2), degree(W2)-degree(T2))))
#          L = lll(orthogonal_submodule(S, T3))
#          genus(L) == G1 || continue
#          M3 = lattice_in_same_ambient_space(S, hcat(basis_matrix(M2), zero_matrix(QQ, rank(M2), degree(W2)-degree(M2))))
#          @assert is_sublattice(L, M3)
#          @assert is_primitive(L, M3)
#          N = orthogonal_submodule(L, M3)
#          bM = coordinates(basis_matrix(M3), L)
#          bN = coordinates(basis_matrix(N), L)
#          L = integer_lattice(; gram=gram_matrix(L))
#          M3 = lattice_in_same_ambient_space(L, bM)
#          N = lattice_in_same_ambient_space(L, bN)
#          for L2 in get_attribute(G1, :representatives)
#            if is_definite(G1)
#              ok, f_top = is_isometric_with_isometry(L, L2)
#            elseif length(get_attribute(G1, :representatives)) == 1
#              ok = true
#            else
#              ok = is_isometric(L, L2)
#            end
#            !ok && continue
#            label_top = get_attribute(L2, :lmfdb_label)
#            if is_definite(G1)
#              v = map_entries(ZZ, coordinates(bN*f_top, L2))
#            else
#              v = zero_matrix(ZZ, 0, 0)
#            end
#            __x = (label_top, label_bottom, v)
#            @show __x
#            push!(results, __x)
#            if !isnothing(_save_path)
#              save_data_lmfdb(_save_path, __x)
#            end
#          end
#        end
#      end
#    end
#  end
#  return results
#end
#
#function _prime_index_overlattices_lmfdb(
#  Gs::Vector{ZZGenus};
#  _save_path::Union{Nothing, String}=nothing,
#)
#  res = Tuple{String, String, Int, Vector{QQFieldElem}}[]
#  for G in Gs
#    q = discriminant_group(G)
#    o = order(q)
#    if isone(o)
#      continue
#    elseif radical(o) == o
#      continue
#    end
#    if is_even(G)
#      q = Hecke._as_finite_bilinear_module(q)
#    end
#    pds = prime_divisors(last(elementary_divisors(q)))
#    D = Dict{eltype(pds), Tuple{TorQuadModule, Vector{Tuple{TorQuadModuleMap, GAPGroupHomomorphism}}}}()
#    for p in pds
#      qp, _ = torsion_subgroup(q, p)
#      r = Oscar._stabilizer_isotropic_elementary(qp, 1)
#      D[p] = (qp, r)
#    end
#    for M in get_attribute(G, :representatives)
#      GM, _ = image_in_Oq(M)
#      label_bottom = get_attribute(M, :lmfdb_label)
#      T = domain(GM)
#      if is_even(G)
#        T = Hecke._as_finite_bilinear_module(T)
#        GM = Oscar._orthogonal_group(T, matrix.(gens(GM)))
#      end
#      for p in pds
#        qp, tmp = D[p]
#        is_empty(tmp) && continue
#        Tp, TpinT = torsion_subgroup(T, p)
#        OTp = orthogonal_group_bilinear(Tp)
#        iso = isomorphism(PermGroup, OTp)
#        Gp, _ = restrict_automorphism_group(GM, TpinT)
#        isoGp, _ = iso(Gp)
#        ok, phi = is_isometric_with_isometry(Tp, qp)
#        @assert ok
#        iphi = inv(phi)
#        for (Hinqp, j) in tmp
#          _H, _HinTp = sub(Tp, elem_type(Tp)[phi\(qp(lift(a))) for a in gens(domain(Hinqp))])
#          _, jTp = sub(codomain(iso), elem_type(codomain(iso))[iso(OTp(phi * hom(j(g)) * iphi; check=false)) for g in gens(domain(j))])
#          for _h in double_cosets(codomain(iso), domain(jTp), isoGp)
#            h = iso\(representative(_h))
#            _Hh, _ = sub(Tp, [h(Tp(lift(a))) for a in gens(_H)])
#            L = cover(_Hh)
#            _v = p*lift(first(gens(_Hh)))
#            G2 = genus(L)
#            k = findfirst(isequal(G2), Gs)
#            @assert !isnothing(k)
#            flag=false
#            if length(get_attribute(Gs[k], :representatives)) == 1
#              flag=true
#              label_top = get_attribute(only(get_attribute(Gs[k], :representatives)), :lmfdb_label)
#            else
#              for L2 in get_attribute(Gs[k], :representatives)
#                !is_isometric(L, L2) && continue
#                flag=true
#                label_top = get_attribute(L2, :lmfdb_label)
#              end
#            end
#            __x = (label_bottom, label_top, p, _v)
#            push!(res, __x)
#            if !isnothing(_save_path)
#              save_data_lmfdb(_save_path, __x)
#            end
#            @assert flag
#          end
#        end
#      end
#    end
#  end
#  return res
#end
#
#function _embeddings_in_K3_lattice(
#  Gs::Vector{ZZGenus};
#  _save_path::Union{Nothing, String}=nothing,
#)
#  L = rescale(k3_lattice(), -1)
#  set_attribute!(L, :lmfdb_label, "")
#  GK3 = genus(L)
#  set_attribute!(GK3, :representatives, ZZLat[L])
#  res = Tuple{String, Int}[]
#  for G in Gs
#    !is_even(G) && continue
#    p, n = signature_pair(G)
#    n in [1,2] || continue
#    @show G
#    r = _primitive_embeddings_in_unimodular_safe_lmfdb(GK3, G)
#    u = unique!([rr[2] for rr in r])
#    for s in u
#      i = count(rr -> rr[2] == s, r)
#      push!(res, (s, i))
#      if !isnothing(_save_path)
#        save_data_lmfdb(_save_path, (s, i))
#      end
#    end
#  end
#  return res
#end
#
#function _treat_case(
#  G1s::Vector{ZZGenus},
#  G2s::Vector{ZZGenus};
#  __save_path::Union{String, Nothing}=nothing,
#)
#  if !isnothing(__save_path)
#    _save_path = joinpath(__save_path, "prim_emb.txt")
#    touch(_save_path)
#    prim_emb = _primitive_embeddings_lmfdb(G1s, G2s; save_path)
#    _save_path = joinpath(__save_path, "overlat.txt")
#    touch(_save_path)
#    overlat = _prime_index_overlattices_lmfdb(G2s; save_path)
#    _save_path = joinpath(__save_path, "k3_embed.txt")
#    touch(_save_path)
#    k3_emb = _embeddings_in_K3_lattice(G2s; _save_path)
#  else
#    prim_emb = _primitive_embeddings_lmfdb(G1s, G2s)
#    overlat = _prime_index_overlattices_lmfdb(G2s)
#    k3_emb = _embeddings_in_K3_lattice(G2s)
#  end
#  return prim_emb, overlat, k3_emb
#end
