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

# Return all the data stored in the fields of `x`
function _data(x::ZZLatGluing)
  return x.glue_map, x.inv_glue_map, x.glue_group_left, x.stabilizer_left, x.glue_group_right, x.stabilizer_right
end

# Return the glue map (i = 1) or its inverse (i = 2) for the associated
# gluing `x`
function glue_map(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.glue_map
  elseif i == 2
    return x.inv_glue_map
  end
end

# Return the glue map of `x` and its inverse
glue_maps(x::ZZLatGluing) = (glue_map(x, 1), glue_map(x, 2))

# Return the left glue group (i = 1) or the right glue group (i = 2)
# associated to `x`
function glue_group(x::ZZLatGluing, i::Int)
  if isone(i)
    return x.glue_group_left
  elseif i == 2
    return x.glue_group_right
  end
end

# Return the glue groups of `x`
glue_groups(x::ZZLatGluing) = (glue_group(x, 1), glue_group(x, 2))

# Return the stabilizer of the left glue group (i = 1) or of the
# right glue group (i = 2) via their embedding in the associated
# classifying group of `x`
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
# to see it from right to left. Mathematically, we just inverse the gluing
# and see it from right to left
function Oscar.inv(x::ZZLatGluing)
  phi, iphi, i1, j1, i2, j2 = _data(x)
  return ZZLatGluing(iphi, phi, i2, j2, i1, j1)
end

# Construct the torsion module D where we construct graphs of glue maps,
# which is just the orthogonal direct sum of q1 and q2. Compute also
# the embeddings O(q1), O(q2) -> O(D). Depending on par, q1 and q2, the
# module D is considered as a bilinear or quadratic module
function _gluing_ambient(
    q1::TorQuadModule,
    q2::TorQuadModule;
    as_bilinear_module::Bool=false,
  )
  x = _direct_sum_with_embeddings_orthogonal_groups(q1, q2; as_bilinear_module)
  return ZZLatGluingAmbient(x...)
end

# Same as before but we input a gluing, which remembers two discriminant
# modules where we glue
function _gluing_ambient(
  x::ZZLatGluing;
  as_bilinear_module::Bool=false,
)
  q1, q2 = codomain.(glue_groups(x))
  return _gluing_ambient(q1, q2; as_bilinear_module)
end

# Return all genera with given form, signature and parity conditions
# There are at most 2 of them
# `par` is either :even, or :odd, or :both
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

# TODO: EXPORT ???
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

# Forget about the quadratic form on `q`, if any. The difference with the Hecke
# function is that we also record the same set of generators
# TODO: Modify the Hecke version to specify whether to keep the same
# generators and delete this version
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

# If `q` is bilinear or `as_bilinear_module` is `true`, return the orthogonal
# group of the associated torsion bilinear form, otherwise the orthogonal group
# of the torsion quadratic form
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
  glue_annihilator_left::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_annihilator_right::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_exponent::Union{Nothing, IntegerUnion}=nothing,
  glue_order::AbstractVector{T}=ZZRingElem[],
  glue_elementary_divisors::AbstractVector=Vector{ZZRingElem}[],
  form_over::Vector{TorQuadModule}=TorQuadModule[],
  genus_over::Vector{ZZGenus}=ZZGenus[],
) where T <: IntegerUnion
  @vprintln :ZZLatWithIsom 3 "Initialize gluing conditions"
  init_gluing_conditions!(
    Fac,
    glue_annihilator_left,
    glue_annihilator_right,
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
# - a set of Gram matrices for the normal forms of the discriminant groups
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
  glue_annihilator_left::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_annihilator_right::Union{Nothing, TorQuadModuleMap}=nothing,
  glue_exponent::Union{Nothing, IntegerUnion}=nothing,
  glue_order::AbstractVector{T}=ZZRingElem[],
  glue_elementary_divisors::AbstractVector=Vector{ZZRingElem}[],
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

  # Left glue group annihilated by the endomorphism glue_annihilator_left
  if !isnothing(glue_annihilator_left)
    @assert domain(glue_annihilator_left) === codomain(glue_annihilator_left) === qM
    VM, VMinqM = kernel(glue_annihilator_left)
  else
    VMinqM = id_hom(qM)
    VM = qM
  end

  # Right glue group annihilated by the endomorphism glue_annihilator_right
  if !isnothing(glue_annihilator_right)
    @assert domain(glue_annihilator_right) === codomain(glue_annihilator_right) === qN
    VN, VNinqN = kernel(glue_annihilator_right)
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
  # potential common subgroup of the conditions modules) we prepare a
  # dictionary which assign a list of gluings to every potential list of
  # elementary divisors for the p-part of a glue group
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
  _glue_order = Set(filter(>(0), glue_order))
  # Since any glue group lies in the group with elementary divisors
  # elG, any good glue order condition must divide k1
  for o in _glue_order
    !is_divisible_by(k1, o) && delete!(_glue_order, o)
  end

  _glue_elem_divs = Set{Vector{ZZRingElem}}()
  # Since any glue group lies in the group with elementary divisors
  # elG, any good glue elementary divisors condition must be in smith
  # normal form and must be compatible with elG
  for _v in glue_elementary_divisors
    if is_empty(_v)
      push!(_glue_elem_divs, ZZRingElem[])
      continue
    end
    _v isa AbstractVector || continue
    eltype(_v) <: RationalUnion || continue
    all(is_integer, _v) || continue
    v = map(ZZ, _v)
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
  _form_over = Set{TorQuadModule}()
  # Any good primitive extension discriminant group condition must look like
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
    push!(_form_over, q)
  end

  # Each primitive extension genus condition and discriminant condition imposes
  # some condition on the order of a glue group. For this, said condition must
  # be compatible with the initial discriminant groups of the problem: if they
  # are not, we discard them.
  # If they are, we compare the induced glue order conditions with the initial
  # ones:
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

  union!(_form_over, discriminant_group.(_genus_over))

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

  # If we allowed any parity for the primitive extensions but the initial
  # conditions only allow for a certain parity, we update Fac.par to the
  # correct value
  if Fac.par === :both
    _t = unique!(modulus_quadratic_form.(_form_over))
    if length(_t) == 1
      Fac.par = isone(unique(_t)) ? :odd : even
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

  # Save the set of potential genera if any specified in input
  if !ignore_genus_over
    Fac.genus_over = _genus_over
  end

  # Save the set of potential gram matrices of the normal forms of the
  # discriminant form of a primitive extension if any specified in input
  if !ignore_form_over
    Fac.form_over = Set(gram_matrix_quadratic.(first.(normal_form.(_form_over))))
  end

  # Save the set of potential elementary divisors of a glue group if
  # any specified in input
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

# Return a complete set of representatives for the O-orbits of p-subgroups of W
# whose reverse sequence of p-valuations of elementary divisors is given by
# `subtype`
#
# If `Fac` is decorated with a gluing context `Ctx` whose node
# `(i, bool, p, subtype)` exists, we fetch the output of this function from
# there. Here `bool` mention whether we consider W as a quadratic module (in
# some large scale computations, this may vary so we have to be careful that
# the group `O` might differ for the same other inputs, depending on whether
# we see W as a bilinear or quadratic modules)
#
# No checks are performed, so in fact here `O` should be one of the local
# classifying groups of `Fac`, it should match the value of `bool` (in the
# sense above) and the data stored at the `(i, bool, p, subtype)` node of
# `Ctx` should have been computed with this exact group `O`.
# If `Ctx` exists but the node `(i, bool, p, subtype)` does not exist, the
# function creates the node and automatically store the output there
function __subgroups_orbit_representatives_and_stabilizers_primary_subtype(
  Fac::ZZLatGluingFactory,
  Winq::TorQuadModuleMap,
  O::AutomorphismGroup{TorQuadModule},
  p::ZZRingElem,
  subtype::Vector{Int},
  i::Int = -1,
)
  W = domain(Winq)
  bool = Fac.par === :even
  if isdefined(Fac, :Ctx) && haskey(Fac.Ctx.orb_and_stab, (i, bool, p, subtype))
    # Retrieve known data
    subs = Fac.Ctx.orb_and_stab[(i, bool, p, subtype)]
  elseif first(subtype) == 1
    # The glue groups are elementary abelian p-groups of rank `length(subtype)`
    W = domain(Winq)
    _, WWinW = torsion_subgroup(W, p)
    WWinq = compose(WWinW, Winq)
    subs = Oscar.__subgroups_orbit_representatives_and_stabilizers_elementary(WWinq, O, p^(length(subtype)), p)
    if isdefined(Fac, :Ctx)
      Fac.Ctx.orb_and_stab[(i, bool, p, subtype)] = subs
    end
  else
    # Non-elementary case
    subs = __subgroups_orbit_representatives_and_stabilizers_primary_subtype(Winq, O, p, subtype)
    if isdefined(Fac, :Ctx)
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
    # then we can fetch data there
    # Note: should only be used with matching classifying groups and parity
    # conditions
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

  # Preparation to compute glue maps; we want to avoid some redundant computations
  _data1 = Array{GAPGroupHomomorphism}(undef, 1, length(subs1))
  _data2 = Array{GAPGroupHomomorphism}(undef, 1, length(subs2))

  for i in 1:length(subs1)
    H1inq1, s1 = subs1[i]
    H1 = domain(H1inq1)
    OH1 = __orthogonal_group(H1; as_bilinear_module)
    imOH1 = elem_type(OH1)[OH1(restrict_automorphism(x, H1inq1); check=false) for x in gens(domain(s1))]
    act1 = hom(domain(s1), OH1, imOH1)
    _, j1 = image(act1)
    _data1[i] = j1
  end

  for j in 1:length(subs2)
    H2inq2, s2 = subs2[j]
    H2 = domain(H2inq2)
    OH2 = __orthogonal_group(H2; as_bilinear_module)
    imOH2 = elem_type(OH2)[OH2(restrict_automorphism(x, H2inq2); check=false) for x in gens(domain(s2))]
    act2 = hom(domain(s2), OH2, imOH2)
    _, j2 = image(act2)
    _data2[j] = j2
  end

  # compute all glue maps for each pair of representatives of
  # orbits of potential glue p-groups, up to the actions
  # of O(q1) and O(q2)
  for j in 1:length(subs2)
    H2inq2, s2 = subs2[j]
    j2 = _data2[j]
    H2 = domain(H2inq2)
    im2 = domain(j2)
    OH2 = codomain(j2)
    elH2 = elementary_divisors(H2)
    if isone(order(H2)) || elH2[1] == elH2[end]
      iso2 = isomorphism(PermGroup, OH2)
    else
      iso2 = id_hom(OH2)
    end
    OH22 = first(iso2(OH2))
    im22 = first(iso2(im2))
    for i in 1:length(subs1)
      H1inq1, s1 = subs1[i]
      j1 = _data1[i]
      im1 = domain(j1)
      H1 = domain(H1inq1)
      ok, phi = is_anti_isometric_with_anti_isometry(H1, H2; as_bilinear_module)
      !ok && continue
      iphi = inv(phi)
      stab1gamma = elem_type(OH2)[OH2(iphi * hom(g) * phi; check=false) for g in gens(im1)]
      S1, _ = sub(OH2, stab1gamma)
      S2 = first(iso2(S1))
      for _g in double_cosets(OH22, im22, S2)
        g = hom(iso2\representative(_g))
        psi = compose(phi, g)
        x = ZZLatGluing(psi, inv(psi), H1inq1, s1, H2inq2, s2)
        push!(loc_glue_p, x)
      end
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
      _phi = glue_map(y, 1)
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
    local_glue_maps(
      q1::TorQuadModule,
      q2::TorQuadModule;
      parity::Symbol=:default,
      Ctx::Union{Nothing, ZZLatGluingCtx}=nothing,
      vi::NTuple{2, Int}=(-1, -1),
      kwargs...,
    ) -> ZZLatGluingFactory, Vector{ZZLatGluing}

Return the local glue maps between `q1` and `q2`, up to the action of their
respective orthogonal groups.

The keyword argument `parity` can be `:default`, `:even`, `:odd` or `:both`
depending on the kind of glue maps one wants to compute. If set to `:default`,
then the algorithm decides depending on the quadratic moduli of `q1` and `q2`

In large scale computations, one can input a gluing context `Ctx` remembering
some orbit and stabilizers computations for certain `TorQuadModule` and for the
given parity. The other keyword `vi` is to refer to which nodes of `Ctx` the
groups `q1` and `q2` are attached (ignores if not in the range of `Ctx`).

The extra keyword arguments include conditions on the gluing. Currently:
  - `glue_annihilator_left::TorQuadModuleMap`: an endomorphism of `q1`;
    all glue groups inside `q1` are in its kernel. Set to be `nothing` by
    default.
  - `glue_annihilator_right::TorQuadModuleMap`: an endomorphism of `q2`;
    all glue groups inside `q2` are in its kernel. Set to be `nothing` by
    default.
  - `glue_exponent::IntegerUnion`: an integer; the exponent of the glue groups
    divides this value. Set to be `nothing` by default.
  - `glue_order::AbstractVector{IntegerUnion}`: a list of integers; the order
    of the glue groups is in this list if not empty. Set to be empty by
    default.
  - `glue_elementary_divisors::AbstractVector`: a list of lists of integers;
    the elementary divisors of a glue group are in this list if not empty.
    Set to be empty by default.
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
function local_glue_maps(
  q1::TorQuadModule,
  q2::TorQuadModule;
  parity::Symbol=:default,
  Ctx::Union{ZZLatGluingCtx, Nothing}=nothing,
  vi::NTuple{2, Int}=(-1, -1),
  kwargs...,
)
  Fac = gluing_factory(q1, q2; parity, Ctx, vi)
  init_gluing_factory!(Fac; kwargs...)
  # If the initial conditions give rise to no possible gluings, then
  # we return the empty list
  is_trivial(Fac) && return Fac, ZZLatGluing[]
  res = local_glue_maps(Fac)
  return Fac, res
end

# Construct the gluing factory attached to `q1` and `q2`, with the given
# `parity`. Possibility to attach a gluing context `Ctx` and vertices
# identifiers `vi`
function gluing_factory(
  q1::TorQuadModule,
  q2::TorQuadModule;
  parity::Symbol=:default,
  Ctx::Union{ZZLatGluingCtx, Nothing}=nothing,
  vi::NTuple{2, Int}=(-1, -1),
)
  Fac = ZZLatGluingFactory(q1, q2)
  if parity === :default
    if modulus_quadratic_form(q1) == modulus_quadratic_form(q2) == 2
      Fac.par = :even
    else
      Fac.par = :odd
    end
  else
    @req parity in [:even, :odd, :both] "Non-supported parity input: $(parity)"
    Fac.par = parity
  end

  if !isnothing(Ctx)
    Fac.Ctx = Ctx
    Fac.vertex_identification = vi
  end
  return Fac
end

###############################################################################
#
#  Operations on gluings
#
###############################################################################

# Split the orbit of x on the left under the action of O1. We first split the
# orbit of the glue group and then for each new glue group, we split the
# double coset of the (pullback of the) glue map.
# If `as_bilinear_module == true`, consider all TorQuadModule as bilinear
# modules.
function _split_orbit_left(
  x::ZZLatGluing,
  O1::AutomorphismGroup{TorQuadModule};
  as_bilinear_module::Bool=false,
)
  res = ZZLatGluing[]
  phi, _, H1inq1, s1, H2inq2, s2 = _data(x)
  @assert domain(O1) === codomain(H1inq1)

  H1, H2 = domain(H1inq1), domain(H2inq2)
  q1 = codomain(H1inq1)
  stab1 = domain(s1)
  Oq1 = codomain(s1)
  if !is_subgroup(O1, Oq1)[1] # TODO: try to find why this can happen
    O1, _ = sub(Oq1, elem_type(Oq1)[Oq1(matrix(g); check=false) for g in gens(O1)])
  end

  # Do some preparations; to split the double coset of glue map, we have
  # to take into account the left action of `s2`
  OH2 = __orthogonal_group(H2; as_bilinear_module)
  imOH2 = elem_type(OH2)[OH2(restrict_automorphism(x, H2inq2); check=false) for x in gens(domain(s2))]
  act2 = hom(domain(s2), OH2, imOH2)
  im2, _ = image(act2)

  # Check if q1 is a free module over some finite ring, in which case we
  # transform Oq1 into a permutation group
  elq1 = elementary_divisors(q1)
  if isone(order(q1)) || elq1[1] == elq1[end]
    iso1 = isomorphism(PermGroup, Oq1)
  else
    iso1 = id_hom(Oq1)
  end

  @vprint :ZZLatWithIsom 1 "Split orbit of glue group"
  splits = double_cosets(codomain(iso1), first(iso1(stab1)), first(iso1(O1)))
  @vprintln :ZZLatWithIsom 1 "                        done: $(length(splits))"
  for _h in splits
    h = hom(iso1\(representative(_h)))
    ih = inv(h)
    # Image group
    _H1, _H1inq1 = sub(q1, elem_type(q1)[h(H1inq1(a)) for a in gens(H1)])
    # Pullback of glue map
    _phi = hom(_H1, H2, elem_type(H2)[phi(H1(lift(ih(_H1inq1(a))))) for a in gens(_H1)])
    _iphi = inv(_phi)
    @hassert :ZZLatWithIsom 1 is_anti_isometry(_phi; as_bilinear_module)
    # Conjugate stabilizer
    _stab1, _ = Oscar._as_subgroup(Oq1, Oscar.GAPWrap.ConjugateSubgroup(GapObj(stab1), GapObj(Oq1(h))))
    # Intersect with the new classifying group O1
    __stab1, _, j1 = intersect(_stab1, O1)

    # We now split the double coset of the glue map
    _OH1 = __orthogonal_group(_H1; as_bilinear_module)
    _imOH1 = elem_type(_OH1)[_OH1(restrict_automorphism(x, _H1inq1); check=false) for x in gens(_stab1)]
    act1 = hom(_stab1, _OH1, _imOH1)
    _ims1, _ = image(act1)
    _imO1, _ = act1(__stab1)

    # Same as before, if possible, use permutation groups (seems faster
    # for double coset computations)
    elH1 = elementary_divisors(H1)
    if isone(order(H1)) || elH1[1] == elH1[end]
      _iso1 = isomorphism(PermGroup, _ims1)
    else
      _iso1 = id_hom(_ims1)
    end
    # This is how s2 acts on the right glue group via the initial glue map _phi
    # Then we intersect with the old clasisfying group stab1, get a group L1
    # then we compute a transversal for the double cosets L1\_ims1/_imO1 to
    # find double cosets of glue maps up to the action of s2 and O1
    _im2phi, _ = sub(_OH1, elem_type(_OH1)[_OH1(_phi * hom(g) * _iphi; check=false) for g in gens(im2)])
    L1, _ = intersect(_im2phi, _ims1)
    for _g in double_cosets(codomain(_iso1), first(_iso1(L1)), first(_iso1(_imO1)))
      psi =  hom(_iso1\(representative(_g))) * _phi
      ipsi = inv(psi)
      push!(res, ZZLatGluing(psi, ipsi, _H1inq1, j1, H2inq2, s2))
    end
  end
  return res
end

# Same as before but for the right. What we do: invert the gluing internally
# (at no cost), split on the left, and invert the outputs back.
# If `as_bilinear_module == true`, consider all TorQuadModule as bilinear
# modules.
function _split_orbit_right(
  x::ZZLatGluing,
  O2::AutomorphismGroup{TorQuadModule};
  as_bilinear_module::Bool=false,
)
  return inv.(_split_orbit_left(inv(x), O2; as_bilinear_module))
end

# For a gluing between subgroups of q1 and q2, and given an isometry
# f:q3 -> q1, pullack the glue map to q3 along f.
# If `as_bilinear_module == true`, consider all TorQuadModule as bilinear
# modules.
function _pullback_left(
  x::ZZLatGluing,
  f::TorQuadModuleMap;
  as_bilinear_module::Bool=false,
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
  @hassert :ZZLatWithIsom 1 domain(s0) == first(stabilizer(Oq0, i0))
  H0toH1 = hom(H0, domain(i1), TorQuadModuleElem[i1\(f(i0(a))) for a in gens(H0)])
  psi = H0toH1 * phi
  @hassert :ZZLatWithIsom 1 is_anti_isometry(psi; as_bilinear_module)
  ipsi = inv(psi)
  return ZZLatGluing(psi, ipsi, i0, s0, i2, s2)
end

# Again, we invert, pullback and invert back.
# If `as_bilinear_module == true`, consider all TorQuadModule as bilinear
# modules.
function _pullback_right(
  x::ZZLatGluing,
  f::TorQuadModuleMap;
  as_bilinear_module::Bool=false,
)
  return inv(_pullback_left(inv(x), f; as_bilinear_module))
end

# Given a glue map, return a representative for the isometry class of the
# discriminant group of the associated primitive extension.
function _form_over(
  x::ZZLatGluing,
  parity::Symbol,
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
  mq = parity === :even ? 2 : 1
  disc = torsion_quadratic_module(cover(perp), cover(ext); modulus=1, modulus_qf=mq)
  return disc
end

# If `as_bilinear_module == true`, consider all TorQuadModule as bilinear
# modules.
function _overlattice_with_glue_stabilizer(
  x::ZZLatGluing,
  D::ZZLatGluingAmbient;
  as_bilinear_module::Bool=false,
)
  return __overlattice(x, D; as_bilinear_module, glue_stab=true)
end

function _overlattice(
  x::ZZLatGluing,
  D::ZZLatGluingAmbient;
  as_bilinear_module::Bool=false,
)
  return __overlattice(x, D; as_bilinear_module, glue_stab=false)[1:3]
end

function _overlattice(
  x::ZZLatGluing;
  as_bilinear_module::Bool=false,
)
  D = _gluing_ambient(x; as_bilinear_module)
  return __overlattice(x, D; as_bilinear_module, glue_stab=false)[1:3]
end

function __overlattice(
  x::ZZLatGluing,
  D::ZZLatGluingAmbient;
  as_bilinear_module::Bool=false,
  glue_stab::Bool=false,
)
  phi, iphi, i1, s1, i2, s2 = _data(x)
  H1, H2 = domain(i1), domain(i2)
  q1, q2 = codomain(i1), codomain(i2)
  H1inD = i1 * D.j1
  H2inD = i2 * D.j2
  O1, O2 = __orthogonal_group(H1; as_bilinear_module), __orthogonal_group(H2; as_bilinear_module)
  im1 = elem_type(O1)[O1(restrict_automorphism(g, i1); check=false) for g in gens(domain(s1))]
  act1 = hom(domain(s1), O1, im1)
  im2 = elem_type(O2)[O2(restrict_automorphism(g, i2); check=false) for g in gens(domain(s2))]
  act2 = hom(domain(s2), O2, im2)
  V, extinD = _overlattice_with_graph(phi, H1inD, H2inD)
  if glue_stab
    @vprintln :ZZLatWithIsom 1 "Glue stabilizers"
    _disc, _stab = _glue_stabilizers(phi, act1, act2, D.k1, D.k2, extinD)
    qV = discriminant_group(V)
    OqV = __orthogonal_group(qV; as_bilinear_module)
    # qV and _disc are equal, but defined in julia with different base abelian group
    # so we cannot treat them as the same object. However we can identity them
    # with a canonical "identity map" phi2
    phi2 = hom(qV, _disc, elem_type(_disc)[_disc(lift(x)) for x in gens(qV)])
    @hassert :ZZLatWithIsom 1 is_isometry(phi2; as_bilinear_module=!is_even(V))
    iphi2 = inv(phi2)
    # Pullback the glue of the stabilizers along the canonical map phi2
    GV, _ = sub(OqV, elem_type(OqV)[OqV(phi2 * g * iphi2; check=false) for g in _stab])
  else
    GV = O1
  end
  # This is where we require that q1, q2 are discriminant groups of some
  # lattices and phi indeed defines a primitive extension between the
  # respective relations
  M, T = _as_sublattices(V, relations(q1), relations(q2))
  return V, M, T, GV
end

# Given a glue map between representatives of orbits of glue groups, we compute
# the representatives for the double cosets of glue maps between these glue groups,
# up to the action of the given classifying groups.
# If `as_bilinear_module == true`, consider all TorQuadModule as bilinear
# modules.
function _all_glue_maps(
  x::ZZLatGluing;
  as_bilinear_module::Bool=false,
)
  res = ZZLatGluing[]
  phi, iphi, H1inq1, s1, H2inq2, s2 = _data(x)
  q1 = codomain(H1inq1)
  q2 = codomain(H2inq2)
  H1 = domain(H1inq1)
  H2 = domain(H2inq2)
  stab1 = domain(s1)
  stab2 = domain(s2)

  OH1 = __orthogonal_group(H1; as_bilinear_module)
  imOH1 = elem_type(OH1)[OH1(restrict_automorphism(x, H1inq1); check=false) for x in gens(stab1)]
  act1 = hom(stab1, OH1, imOH1)
  im1, _ = image(act1)
  OH2 = __orthogonal_group(H2; as_bilinear_module)
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
