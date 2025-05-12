###################################################################################################################
#
# Constructors
#
###################################################################################################################

@doc raw"""
    elliptic_surface(generic_fiber::EllipticCurve,
                     euler_characteristic::Int,
                     mwl_gens::Vector{<:EllipticCurvePoint}=EllipticCurvePoint[];
                     is_basis::Bool=true)
                     -> EllipticSurface

Return the relatively minimal elliptic surface with generic fiber ``E/k(t)``.

This is also known as the Kodaira-Néron model of ``E``.

Input:
- `generic_fiber` -- an elliptic curve over a function field
- `euler_characteristic` -- the Euler characteristic of the Kodaira-Néron model of ``E``.
- `mwl_gens` -- a vector of rational points of the generic fiber
- `is_basis` -- if set to `false` compute a reduced basis from `mwl_gens`

# Examples
```jldoctest
julia> Qt, t = polynomial_ring(QQ, :t);

julia> Qtf = fraction_field(Qt);

julia> E = elliptic_curve(Qtf, [0,0,0,0,t^5*(t-1)^2]);

julia> X3 = elliptic_surface(E, 2)
Elliptic surface
  over rational field
with generic fiber
  -x^3 + y^2 - t^7 + 2*t^6 - t^5

```
"""
function elliptic_surface(generic_fiber::EllipticCurve{BaseField},
                          euler_characteristic::Int,
                          mwl_gens::Vector{<:EllipticCurvePoint}=EllipticCurvePoint[];
                          resolution_strategy::Symbol=:iterative,
                          is_basis::Bool=true) where {
                          BaseField <: FracFieldElem{<:PolyRingElem{<:FieldElem}}}
  @req all(parent(i)==generic_fiber for i in mwl_gens) "not a vector of points on $(generic_fiber)"
  S = EllipticSurface(generic_fiber, euler_characteristic, mwl_gens; resolution_strategy)
  if is_basis
    return S
  end
  update_mwl_basis!(S, mwl_gens)
  return S
end

@doc raw"""
    kodaira_neron_model(E::EllipticCurve) -> EllipticSurface
    
Return the Kodaira-Neron model of the elliptic curve `E`.
"""
kodaira_neron_model(E::EllipticCurve) = elliptic_surface(E)

@doc raw"""
    elliptic_surface(g::MPolyRingElem, P::Vector{<:RingElem})

Transform a bivariate polynomial `g` of the form `y^2 - Q(x)` with `Q(x)` of
degree at most ``4`` to Weierstrass form, apply Tate's algorithm and 
return the corresponding relatively minimal elliptic surface 
as well as the coordinate transformation.
"""
function elliptic_surface(
    g::MPolyRingElem, P::Vector{<:RingElem}; 
    minimize::Bool=true, resolution_strategy::Symbol=:iterative
  )
  R = parent(g)
  (x, y) = gens(R)
  P = base_ring(R).(P)
  g2, phi2 = transform_to_weierstrass(g, x, y, P);
  Y2, phi1 = _elliptic_surface_with_trafo(g2; minimize)
  return Y2, phi2 * phi1  
end

@doc raw"""
    fibration_on_weierstrass_model(X::EllipticSurface)
    
Return the elliptic fibration ``W \to \mathbb{P}^1`` where ``W`` is the Weierstrass model of ``X``.
"""
function fibration_on_weierstrass_model(X::EllipticSurface)
  if !isdefined(X, :fibration_weierstrass_model)
    weierstrass_model(X) # trigger caching
  end
  return X.fibration_weierstrass_model
end

@doc raw"""
    fibration(X::EllipticSurface)
    
Return the elliptic fibration ``X \to \mathbb{P}^1``.
"""
function fibration(X::EllipticSurface)
  if !isdefined(X, :fibration)
    X.fibration = compose(weierstrass_contraction(X), fibration_on_weierstrass_model(X))
  end
  return X.fibration
end

###################################################################################################################
#
# Basic attributes, properties
#
###################################################################################################################
  
# Unlocks Scheme functionality
function underlying_scheme(X::EllipticSurface)
  S = X
  if isdefined(S,:Y)
    return S.Y
  end
  # trigger the computation
  weierstrass_contraction(S)
  return underlying_scheme(S)
end

base_ring(X::EllipticSurface) = coefficient_ring(base_ring(base_field(generic_fiber(X))))

@doc raw"""
    generic_fiber(X::EllipticSurface) -> EllipticCurve

Return the generic fiber as an elliptic curve.
"""
generic_fiber(X::EllipticSurface) = X.E

@doc raw"""
    weierstrass_chart(X::EllipticSurface)

Return the Weierstrass chart of ``X`` on its `weierstrass_model`.
"""
weierstrass_chart(X::EllipticSurface) = weierstrass_model(X)[1][1][1]

@doc raw"""
    euler_characteristic(X::EllipticSurface) -> Int

Return the Euler characteristic ``\chi(\mathcal{O}_X)``.
"""
euler_characteristic(X::EllipticSurface) = X.euler_characteristic
 

########################################################################################################
#
# Printing
#
########################################################################################################
  
function Base.show(io::IO, X::EllipticSurface)
  io = pretty(io)
  if is_terse(io)
    print(io, "Elliptic surface")
  else
    E = generic_fiber(X)
    print(io, "Elliptic surface with generic fiber ", equation(E))
  end
end

function Base.show(io::IO, ::MIME"text/plain", X::EllipticSurface)
  io = pretty(io)
  println(io, "Elliptic surface")
  println(io, Indent(), "over ", Lowercase(), base_ring(X))
  println(io, Dedent(), "with generic fiber")
  print(io, Indent(), Lowercase(), equation(generic_fiber(X)), Dedent())
  if isdefined(X, :Y)
    println(io)
    println(io, "and relatively minimal model")
    print(io, Indent(), Lowercase(), X.Y, Dedent())
  end
  print(io, Dedent())
end

########################################################################################################
#
# Updating the Mordell-Weil generators.
#
########################################################################################################
@doc raw"""
    set_mordell_weil_basis!(X::EllipticSurface, mwl_basis::Vector{EllipticCurvePoint})

Set a basis for the Mordell-Weil sublattice of ``X`` or at least of a sublattice.

This invalidates previous computations depending on the generators of the
Mordell Weil lattice such as the `algebraic_lattice`. Use with care.

The points in `mwl_basis` must be linearly independent.
"""
function set_mordell_weil_basis!(X::EllipticSurface, mwl_basis::Vector{<:EllipticCurvePoint})
  @req all(parent(P) == generic_fiber(X) for P in mwl_basis) "points must lie on the generic fiber"
  X.MWL = mwl_basis
  # clear old computations
  if has_attribute(X, :algebraic_lattice)
    delete!(X.__attrs, :algebraic_lattice)
  end
  if has_attribute(X, :mordell_weil_sublattice)
    delete!(X.__attrs, :mordell_weil_sublattice)
  end
end
  
@doc raw"""
    _compute_mwl_basis(X::EllipticSurface, mwl_gens::Vector{<:EllipticCurvePoint}) -> ZZLat, Vector{<:EllipticCurvePoint}

Return a tuple `(M, B)` where  `B` is an LLL-reduced basis of the sublattice `M` of the
Mordell-Weil lattice of ``X`` generated by `mwl_gens`.
"""
function _compute_mwl_basis(X::EllipticSurface, mwl_gens::Vector{<:EllipticCurvePoint})
  # it would be good to have the height pairing implemented
  basis,tors, SX = _algebraic_lattice(X, mwl_gens)
  basisTriv, GTriv = trivial_lattice(X)
  r = length(basisTriv)
  l = length(mwl_gens)
  V = ambient_space(SX)
  rk = rank(V)
  G = ZZ.(gram_matrix(V))
  # project away from the trivial lattice
  pr_mwl = orthogonal_projection(V,basis_matrix(SX)[1:r, :])
  BMWL = pr_mwl.matrix[r+1:end, :]
  GB = gram_matrix(V,BMWL)
  @assert rank(GB) == rk-r
  _, u = hnf_with_transform(ZZ.(denominator(GB) * GB))
  B = u[1:rk-r,:] * BMWL

  MWL = lll(lattice(V, B, isbasis=false))
  u = solve(BMWL, basis_matrix(MWL); side=:left)
  u = ZZ.(u)
  mwl_basis = [sum(u[i,j] * mwl_gens[j] for j in 1:length(mwl_gens)) for i in 1:nrows(u)]
  return MWL, mwl_basis
end

@doc raw"""
    update_mwl_basis!(X::EllipticSurface, mwl_gens::Vector{<:EllipticCurvePoint})

Compute a reduced basis of the sublattice of the Mordell-Weil lattice spanned
by `mwl_gens` and set these as the new generators of the Mordell-Weil lattice of
``X``.
"""
function update_mwl_basis!(X::EllipticSurface, mwl_gens::Vector{<:EllipticCurvePoint})
  mwl, mwl_basis = _compute_mwl_basis(X, mwl_gens)
  set_mordell_weil_basis!(X, mwl_basis)
end

@doc raw"""
    algebraic_lattice_primitive_closure(X::EllipticSurface, p) -> Vector{<:EllipticCurvePoint}

Return sections ``P_1,\dots P_n`` of the generic fiber, such that together with
the generators of the algebraic lattice ``A``, they generate
```math
\frac{1}{p} A \cap N
```
where ``N`` is the numerical lattice of ``X``.

The algorithm proceeds by computing division points in the Mordell-Weil subgroup of `X`
and using information coming from the discriminant group of the algebraic lattice
to do so.
"""
algebraic_lattice_primitive_closure(X::EllipticSurface, p) = algebraic_lattice_primitive_closure(X, ZZ(p))

function algebraic_lattice_primitive_closure(X::EllipticSurface, p::ZZRingElem)
  S = X
  L = algebraic_lattice(S)[3]
  @req is_even(L) "not implemented"
  Ld  = intersect(dual(L) , (1//p * L))
  D = torsion_quadratic_module(Ld, L, modulus = 1, modulus_qf=2)
  candidates = [x for x in D if !iszero(x) && iszero(quadratic_product(x))]
  t = length(trivial_lattice(S)[1])
  r = rank(L)
  cc = [(x->mod(x,p)).(p*lift(c)[t+1:end]) for c in candidates]
  unique!(cc)
  cc = [c for c in cc if !iszero(c)]
  pts = [division_points(sum(v[i]*S.MWL[i] for i in 1:(r-t)),p) for v in cc]
  return [i[1] for i in pts if length(i)>0]
end

function algebraic_lattice_primitive_closure!(X::EllipticSurface, prime)
  pts = algebraic_lattice_primitive_closure(X, prime)
  update_mwl_basis!(X, vcat(pts, X.MWL))
  return pts
end

@doc raw"""
    algebraic_lattice_primitive_closure!(X::EllipticSurface)

Compute the primitive closure of the algebraic lattice of ``X`` inside its
numerical lattice and update the generators of its Mordell--Weil group accordingly.

The algorithm works by computing suitable division points in its Mordell Weil group.
"""
function algebraic_lattice_primitive_closure!(X::EllipticSurface)
  S = X
  L = algebraic_lattice(S)[3]
  for p in prime_divisors(ZZ(det(L)))
    while true
      pts = algebraic_lattice_primitive_closure!(S, p)
      if length(pts)==0
        break
      end
    end
  end
  return S
end


###################################################################################################################
#
# Kodaira-Néron Model
#
###################################################################################################################



@doc raw"""
    weierstrass_chart_on_minimal_model(X::EllipticSurface)

Return an affine chart ``U`` of ``X`` which is isomorphic to the `weierstrass_chart` 
of ``X`` on its `weierstrass_model`, but with all singular fibers removed.

More precisely, the affine coordinates of ``U`` are ``(x,y,t)`` and the chart is 
constructed as the vanishing locus of
```math
y^2 + a_1(t) xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6
```
minus the reducible singular fibers.
"""
weierstrass_chart_on_minimal_model(X::EllipticSurface) = X[1][1]

@doc raw"""
    weierstrass_model(X::EllipticSurface) -> CoveredScheme, CoveredClosedEmbedding

Return the Weierstrass model ``W`` of ``X`` and the inclusion in
its ambient projective bundle
```math
S\subseteq \mathbb{P}( \mathcal{O}_{\mathbb{P}^1}(-2s) \oplus \mathcal{O}_{\mathbb{P}^1}(-3s) \oplus \mathcal{O}_{\mathbb{P}^1}).
```
"""
function weierstrass_model(X::EllipticSurface)
  if isdefined(X, :Weierstrassmodel)
    return X.Weierstrassmodel, X.inc_Weierstrass
  end

  s = euler_characteristic(X)
  E = generic_fiber(X)

  kt = base_ring(base_field(E))
  k = coefficient_ring(kt)

  IP1 = projective_space(k, 1)
  c = standard_covering(IP1)
  # rename the variables on the affine charts
  # to a more readable version
  if k isa FqField
    OO(c[1]).data.S = [:t]
    OO(c[2]).data.S = [:s]
  else
    OO(c[1]).S = [:t]
    OO(c[2]).S = [:s]
  end

  O0 = twisting_sheaf(IP1, 0)
  O4 = twisting_sheaf(IP1, -2*s)
  O6 = twisting_sheaf(IP1, -3*s)

  bundleE = direct_sum([O0, O4, O6])

  P_proj = projectivization(bundleE, var_names=[:z, :x, :y])
  P = covered_scheme(P_proj)
  pr = covered_projection_to_base(P_proj)
  @assert has_decomposition_info(default_covering(P))

  # Create the singular Weierstrass model S of the elliptic K3 surface X
  a = a_invariants(E)
  U = affine_charts(P)[1]  # the standard Weierstrass chart
  (x, y, t) = gens(OO(U))
  @assert all(denominator(i)==1 for i in a)
  a = [numerator(a)(t) for a in a]
  (a1,a2,a3,a4,a6) = a
  ft = y^2  + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x+a6)
  I = IdealSheaf(P, U, [ft]; check=false)

  inc_S = CoveredClosedEmbedding(P, I)
  Scov = domain(inc_S)  # The ADE singular elliptic K3 surface
  X.Weierstrasschart = Scov[1][1]
  X.fibration_weierstrass_model = compose(inc_S, pr)

  X.Weierstrassmodel = Scov
  X.inc_Weierstrass = inc_S

  set_attribute!(Scov, :is_irreducible=>true)
  set_attribute!(Scov, :is_reduced=>true)
  set_attribute!(Scov, :is_integral=>true)
  set_attribute!(Scov, :is_equidimensional=>true)
  return Scov, inc_S
end
  
@doc raw"""
    weierstrass_contraction(X::EllipticSurface) -> SchemeMor

Return the contraction morphism of ``X`` to its Weierstrass model.

This triggers the computation of the `underlying_scheme` of ``X``
as a blowup from its Weierstrass model. It may take a few minutes.
"""
function weierstrass_contraction(X::EllipticSurface)
  algorithm = X.resolution_strategy
  if algorithm == :iterative
    return weierstrass_contraction_iterative(X)
  elseif algorithm == :simultaneous
    return weierstrass_contraction_simultaneous(X)
  else
    error("algorithm not recognized")
  end
end

@doc raw"""
    _separate_singularities!(X::EllipticSurface) -> Covering

Create a covering of the ambient projective bundle ``P``
of the Weierstrass model ``W`` of ``X`` such that each chart
(of ``X``) contains at most one singular point of ``W``.
Append this covering to the list of coverings of ``X`` and return it.
"""
function _separate_singularities!(X::EllipticSurface)
  S, inc_S = weierstrass_model(X)
  P = codomain(inc_S)

  I_sing = ideal_sheaf_of_singular_locus(S)
  I_sing_P = SimplifiedIdealSheaf(pushforward(inc_S)(I_sing))

  # Refine the covering over the reducible singular fibers
  # to make sure that there is only a single singular point in each chart
  refined_charts = AbsAffineScheme[]
  U = P[1][1]  # the weierstrass_chart
  IsingU = I_sing_P(U)::MPolyIdeal
  if isone(IsingU)
    # we want one smooth weierstrass chart
    push!(refined_charts, U)
    set_attribute!(U, :is_smooth => true)
  else
    # there is at most one singularity in every fiber
    # project the singular locus to an affine chart of P1
    disc = gens(eliminate(IsingU, coordinates(U)[1:2]))[1]
    # The t-coordinates of the reducible fibers
    redfib = [f[1] for f in factor(disc)]
    # One chart with all reducible fibers taken out
    UU = PrincipalOpenSubset(U, redfib)
    set_attribute!(UU, :is_smooth => true)
    push!(refined_charts, UU)
    if length(redfib)==1
      # We need to recreate U as a PrincipalOpenSubset of itself here 
      # in order to maintain the correct tree-structure for refinements.
      # In any Covering no patch is allowed to be an ancestor of another. 
      push!(refined_charts, PrincipalOpenSubset(U, one(OO(U))))
    else
      for i in 1:length(redfib)
        # We take out all but the i-th singular fiber
        r = copy(redfib)
        g = r[i]
        deleteat!(r, i)
        Uref = PrincipalOpenSubset(U, r)
        push!(refined_charts, Uref)
      end
    end
  end

  # Create a chart which contains the fiber over s=0
  # and no other reducible singular fibers
  # these are visible in the charts that we have already
  # i.e. we add the fiber at s=0 and remove all other singular fibers
  V = P[1][4]
  IsingV = I_sing_P(V)
  if isone(IsingV)
    push!(refined_charts, V)
  else
    # reducible singular fibers
    disc = gens(eliminate(IsingV, coordinates(V)[1:2]))[1]
    (x,y,s) = coordinates(V)
    b, d = divides(disc, s)
    if b
      disc = d
    end
    redfib = [f[1] for f in factor(disc)]
    if length(redfib)> 0
      push!(refined_charts, PrincipalOpenSubset(V, redfib))
    else
      push!(refined_charts, V)
    end
  end

  # no extra singularities in the X = 1 chart
  # therefore we just exclude all the singularities visible here
  for W in [P[1][2],P[1][5]]
    Ising = I_sing_P(W)
    if isone(Ising)
      push!(refined_charts, W)
      continue
    end
     (z,y,s_or_t) = coordinates(W)
    # reducible singular fibers
    local disc = gens(eliminate(Ising, [z, s_or_t]))[1]
    local redfib = [p for (p,e) in factor(disc)]
    push!(refined_charts, PrincipalOpenSubset(W, redfib))
  end

  # no extra singularities on the the zero section
  # This is the Y = 1 chart
  # therefore we just exclude all the singularities visible here
  for W in [P[1][3],P[1][6]]
    local Ising = I_sing_P(W)
    if isone(Ising)
      push!(refined_charts, W)
      continue
    end
    local (z,x,s_or_t) = coordinates(W)
    # reducible singular fibers
    local disc = gens(eliminate(Ising, [x, s_or_t]))[1]
    local redfib = [p for (p,e) in factor(disc)]
    push!(refined_charts, PrincipalOpenSubset(W, redfib))
  end


  Cref = Covering(refined_charts)
  inherit_gluings!(Cref, P[1])
  push!(P.coverings, Cref)
  @assert has_decomposition_info(default_covering(P))
  inherit_decomposition_info!(P, Cref)
  @assert has_decomposition_info(Cref)
  # Now we have an extra covering where each chart just contains a single singularity

  @assert scheme(I_sing) === S
  @assert scheme(I_sing_P) === P
  return Cref
end


function weierstrass_contraction_simultaneous(Y::EllipticSurface)
  if isdefined(Y, :blowup)
    return Y.blowup
  end
  S, inc_S = weierstrass_model(Y)
  @assert has_attribute(S, :is_equidimensional) && get_attribute(S, :is_equidimensional) === true
  
  X0 = codomain(inc_S)
  Y0 = S
  set_attribute!(Y0, :is_reduced=>true)
  set_attribute!(Y0, :is_irreducible=>true)
  set_attribute!(Y0, :is_equidimensional=>true)
  inc_Y0 = inc_S
  I_sing_Y0 = AbsIdealSheaf[ideal_sheaf_of_singular_locus(Y0)]
  #I_sing_Y0 = AbsIdealSheaf[simplify(ideal_sheaf_of_singular_locus(Y0))]
  I_sing_X0 = simplify.(pushforward(inc_Y0).(I_sing_Y0))

  # Prepare a covering which has a permanent weierstrass chart
  U0 = X0[1][1]
  disc = numerator(discriminant(generic_fiber(Y)))
  U = PrincipalOpenSubset(U0, evaluate(disc, gens(OO(U0))[3]))
  _find_chart(U, default_covering(X0))
  Cref = Covering(vcat([U], affine_charts(X0)))
  inherit_gluings!(Cref, X0[1])
  push!(X0.coverings, Cref)
  @assert has_decomposition_info(default_covering(X0))
  inherit_decomposition_info!(X0, Cref)
  @assert has_decomposition_info(Cref)

  ambient_exceptionals = EffectiveCartierDivisor[]
  varnames = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j,:k,:l,:m,:n,:o,:p,:q,:r,:u,:v,:w]
  projectionsX = BlowupMorphism[]
  projectionsY = AbsCoveredSchemeMorphism[]
  count = 0

  @vprint :EllipticSurface 2 "Blowing up Weierstrass model simultaneously in all singular points\n"
  while true
    count = count+1
    @vprint :EllipticSurface 1 "blowup number: $(count)\n"
    @vprint :EllipticSurface 1 "number of ideal sheaves to be blown up: $(length(I_sing_X0))\n"
    if length(I_sing_X0)==0
      # stop if smooth
      break
    end
    # make sure we have a Weierstrass chart kept.
    if count == 1
      cov = Cref
    else
      cov = simplified_covering(X0)
    end
    # take the first ideal sheaf and blow it up
    J = SimplifiedIdealSheaf(I_sing_X0[1])
    pr_X1 = blow_up(J, covering=cov, var_name=varnames[1+mod(count, length(varnames))])

    X1 = domain(pr_X1)
    @vprint :EllipticSurface 1 "$(X1)\n"
    E1 = exceptional_divisor(pr_X1)

    @vprint :EllipticSurface 2 "computing strict transforms\n"
    # compute the exceptional divisors
    ambient_exceptionals = EffectiveCartierDivisor[strict_transform(pr_X1, e) for e in ambient_exceptionals]
    # move the divisors coming originally from S up to the next chart
    push!(ambient_exceptionals, E1)

    Y1, inc_Y1, pr_Y1 = strict_transform(pr_X1, inc_Y0)
    # Speed up the computation of singular loci
    set_attribute!(Y1, :is_irreducible=> true)
    set_attribute!(Y1, :is_reduced=>true)
    set_attribute!(Y1, :is_integral=>true)
    set_attribute!(Y1, :is_equidimensional=>true)

    # transform the singular loci
    I_sing_X0 = AbsIdealSheaf[pullback(pr_X1, J) for J in I_sing_X0[2:end]]

    # Add eventual new components
    @vprint :EllipticSurface 2 "computing singular locus\n"
    I_sing_new = ideal_sheaf_of_singular_locus(Y1; focus=pullback(inc_Y1, ideal_sheaf(E1)))
    #I_sing_new = pushforward(inc_Y1, I_sing_new) + ideal_sheaf(E1) # new components only along the exc. set
    I_sing_new = simplify(pushforward(inc_Y1, I_sing_new))

    @vprint :EllipticSurface 2 "decomposing singular locus\n"
    !is_one(I_sing_new) && push!(I_sing_X0, I_sing_new)

    push!(projectionsX, pr_X1)
    push!(projectionsY, pr_Y1)
    simplify!(Y1)

    # set up for the next iteration
    Y0 = Y1
    inc_Y0 = inc_Y1
    X0 = X1
    # Speed up the computation of singular loci
    set_attribute!(Y0, :is_irreducible=> true)
    set_attribute!(Y0, :is_reduced=>true)
    set_attribute!(Y0, :is_integral=>true)
    set_attribute!(Y0, :is_equidimensional=>true)
    set_attribute!(X0, :is_irreducible=> true)
    set_attribute!(X0, :is_reduced=>true)
    set_attribute!(X0, :is_integral=>true)
  end
  Y.Y = Y0
  Y.blowups = projectionsY

  # We need to rewrap the last maps so that the domain is really Y
  last_pr = pop!(projectionsY)
  last_pr_wrap = CoveredSchemeMorphism(Y, codomain(last_pr), covering_morphism(last_pr))

  push!(projectionsY, last_pr_wrap)
  Y.ambient_blowups = projectionsX

  Y.ambient_exceptionals = ambient_exceptionals
  piY = CompositeCoveredSchemeMorphism(reverse(projectionsY))
  Y.blowup = piY

  inc_Y0_wrap = CoveredClosedEmbedding(Y, codomain(inc_Y0), covering_morphism(inc_Y0), check=false)
  Y.inc_Y = inc_Y0_wrap

  set_attribute!(Y, :is_irreducible=> true)
  set_attribute!(Y, :is_reduced=>true)
  set_attribute!(Y, :is_integral=>true)
  return piY
end

     
function weierstrass_contraction_iterative(Y::EllipticSurface)
  if isdefined(Y, :blowup)
    return Y.blowup
  end
  S, inc_S = weierstrass_model(Y)
  @assert has_attribute(S, :is_equidimensional) && get_attribute(S, :is_equidimensional) === true
  Crefined = _separate_singularities!(Y)
  # Blow up singular points (one at a time) until smooth
  # and compute the strict transforms of the `divisors`
  # collect the exceptional divisors
  # blowup ambient spaces: X0 → X   ⊂
  # blowup pi: Y  → (S singular weierstrass model)
  #
  # initialization for the while loop
  X0 = codomain(inc_S)
  Y0 = S
  set_attribute!(Y0, :is_reduced=>true)
  set_attribute!(Y0, :is_irreducible=>true)
  set_attribute!(Y0, :is_equidimensional=>true)
  inc_Y0 = inc_S
  I_sing_Y0 = maximal_associated_points(ideal_sheaf_of_singular_locus(Y0))::Vector{<:AbsIdealSheaf}
  I_sing_X0 = pushforward(inc_Y0).(I_sing_Y0)


  ambient_exceptionals = EffectiveCartierDivisor[]
  varnames = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j,:k,:l,:m,:n,:o,:p,:q,:r,:u,:v,:w]
  projectionsX = BlowupMorphism[]
  projectionsY = AbsCoveredSchemeMorphism[]
  count = 0

  @vprint :EllipticSurface 2 "Blowing up Weierstrass model\n"
  @vprint :EllipticSurface 2 "in $(Crefined)\n"
  while true
    count = count+1
    @vprint :EllipticSurface 1 "blowup number: $(count)\n"
    @vprint :EllipticSurface 1 "number of singular points: $(length(I_sing_X0))\n"
    if length(I_sing_X0)==0
      # stop if smooth
      break
    end
    # make sure there is only one singular point per chart
    if count == 1
      cov = Crefined
    else
      # the following leads to difficult bugs
      cov = simplified_covering(X0)
    end
    # take the first singular point and blow it up
    J = SimplifiedIdealSheaf(I_sing_X0[1])
    pr_X1 = blow_up(J, covering=cov, var_name=varnames[1+mod(count, length(varnames))])

    # Set the attribute so that the strict_transform does some extra work
    #isomorphism_on_open_subset(pr_X1)

    X1 = domain(pr_X1)
    @vprint :EllipticSurface 1 "$(X1)\n"
    E1 = exceptional_divisor(pr_X1)

    @vprint :EllipticSurface 2 "computing strict transforms\n"
    # compute the exceptional divisors
    ambient_exceptionals = EffectiveCartierDivisor[strict_transform(pr_X1, e) for e in ambient_exceptionals]
    # move the divisors coming originally from S up to the next chart
    push!(ambient_exceptionals, E1)

    Y1, inc_Y1, pr_Y1 = strict_transform(pr_X1, inc_Y0)
    # Speed up the computation of singular loci
    set_attribute!(Y1, :is_irreducible=> true)
    set_attribute!(Y1, :is_reduced=>true)
    set_attribute!(Y1, :is_integral=>true)
    set_attribute!(Y1, :is_equidimensional=>true)

    # transform the singular loci
    I_sing_X0 = AbsIdealSheaf[pullback(pr_X1, J) for J in I_sing_X0[2:end]]

    # Add eventual new components
    @vprint :EllipticSurface 2 "computing singular locus\n"
    I_sing_new = ideal_sheaf_of_singular_locus(Y1; focus=pullback(inc_Y1, ideal_sheaf(E1)))
    #I_sing_new = pushforward(inc_Y1, I_sing_new) + ideal_sheaf(E1) # new components only along the exc. set
    I_sing_new = pushforward(inc_Y1, I_sing_new)

    @vprint :EllipticSurface 2 "decomposing singular locus\n"
    I_sing_X0 = vcat(I_sing_X0, maximal_associated_points(I_sing_new))

    push!(projectionsX, pr_X1)
    push!(projectionsY, pr_Y1)
    simplify!(Y1)

    # set up for the next iteration
    Y0 = Y1
    inc_Y0 = inc_Y1
    X0 = X1
    # Speed up the computation of singular loci
    set_attribute!(Y0, :is_irreducible=> true)
    set_attribute!(Y0, :is_reduced=>true)
    set_attribute!(Y0, :is_integral=>true)
    set_attribute!(Y0, :is_equidimensional=>true)
    set_attribute!(X0, :is_irreducible=> true)
    set_attribute!(X0, :is_reduced=>true)
    set_attribute!(X0, :is_integral=>true)
  end
  Y.Y = Y0
  Y.blowups = projectionsY

  # We need to rewrap the last maps so that the domain is really Y
  last_pr = pop!(projectionsY)
  last_pr_wrap = CoveredSchemeMorphism(Y, codomain(last_pr), covering_morphism(last_pr))
  #set_attribute!(last_pr_wrap, :isomorphism_on_open_subset, get_attribute(last_pr, :isomorphism_on_open_subset))

  push!(projectionsY, last_pr_wrap)
  Y.ambient_blowups = projectionsX

  Y.ambient_exceptionals = ambient_exceptionals
  piY = CompositeCoveredSchemeMorphism(reverse(projectionsY))
  Y.blowup = piY

  inc_Y0_wrap = CoveredClosedEmbedding(Y, codomain(inc_Y0), covering_morphism(inc_Y0), check=false)
  Y.inc_Y = inc_Y0_wrap

  set_attribute!(Y, :is_irreducible=> true)
  set_attribute!(Y, :is_reduced=>true)
  set_attribute!(Y, :is_integral=>true)
  return piY
end

###################################################################################################################
#
# Intersection theory, Neron-Severi-lattice, Mordell-Weil lattice
#
###################################################################################################################

@doc raw"""
    algebraic_lattice(X) -> Vector{AbsWeilDivisor}, ZZLat

Return the sublattice ``L`` of the numerical lattice spanned by fiber components,
torsion sections and the sections provided at the construction of ``X``.

The first return value is a list ``B`` of vectors corresponding to the standard basis of the ambient space of `L`.
The second consists of generators ``T`` for the torsion part of the Mordell-Weil group. 
Together ``B`` and ``L`` generate the algebraic lattice. However, they are linearly dependent.
The third return value is the lattice ``L``.

!!! warning 
    The ordering of the fiber components is in general not canonical and may change in between sessions. 
"""
@attr Any function algebraic_lattice(X::EllipticSurface)
  return _algebraic_lattice(X,X.MWL)
end

function _algebraic_lattice(X::EllipticSurface, mwl_basis::Vector{<:EllipticCurvePoint})
  basisTriv, GTriv = trivial_lattice(X)
  r = length(basisTriv)
  l = length(mwl_basis)
  sections = [section(X, i) for i in mwl_basis]
  n = l+r
  GA = zero_matrix(ZZ, n, n)
  GA[1:r,1:r] = GTriv
  GA[r+1:n,r+1:n] = -euler_characteristic(X)*identity_matrix(ZZ, l)
  basisA = vcat(basisTriv, sections)
  @vprint :EllipticSurface 2 "computing intersection numbers\n"
  for i in 1:n
      @vprint :EllipticSurface 3 "\nrow $(i): \n"
      for j in max(i + 1, r + 1):n
      if i!=2 && i <= r
        I = components(basisA[i])[1]
        J = components(basisA[j])[1]
        if isone(I+J)
          ij = 0
        else
          ij = 1
        end
      else
        @vprint :EllipticSurface 4 "$(j) "
        ij = intersect(basisA[i],basisA[j])
      end
      GA[i,j] = ij
      GA[j,i] = GA[i,j]
    end
  end
  GA_QQ = change_base_ring(QQ,GA)
  # primitive closure of the trivial lattice comes from torsion sections
  tors = [section(X, h) for h in mordell_weil_torsion(X)]
  torsV = QQMatrix[]
  for T in tors
    @vprint :EllipticSurface 2 "computing basis representation of torsion point $(T)\n"
    vT = zero_matrix(QQ, 1, n)
    for i in 1:r
      if i== 2
        vT[1,i] = intersect(basisA[i], T)
      else
        @assert length(components(basisTriv[i])) == 1
        I = sum(components(basisA[i]))
        J = components(T)[1]
        if !isone(I+J)
          @assert i!=2 # O does not meet any torsion section
          vT[1,i] = 1
        end
      end
    end
    for i in r+1:n
      vT[1,i] = intersect(T,basisA[i])
    end
    push!(torsV, solve(GA_QQ, vT; side=:left))
  end
  gen_tors = zip(tors, torsV)
  push!(torsV, identity_matrix(QQ,n))
  V = quadratic_space(QQ,GA)
  L = lattice(V, reduce(vcat,torsV), isbasis=false)
  return basisA, collect(gen_tors), L
end

@doc raw"""
    mordell_weil_sublattice(X::EllipticSurface) -> Vector{EllipticCurvePoint}, ZZLat

Return the (sublattice) of the Mordell-Weil lattice of ``X``  spanned
by the sections of ``X`` supplied at its construction.

The Mordell-Weil lattice is represented in the same vector space as the
algebraic lattice (with quadratic form rescaled by ``-1``).
"""
@attr ZZLat function mordell_weil_sublattice(X::EllipticSurface)
  S = X
  NS = algebraic_lattice(S)[3]
  t = length(trivial_lattice(S)[1])
  trivNS = basis_matrix(NS)[1:t,:]
  R = basis_matrix(NS)[t+1:end,:]
  V = ambient_space(NS)
  P = orthogonal_projection(V, trivNS)
  mwl = rescale(lattice(V,R*P.matrix),-1)
  return mwl
end

@doc raw"""
    mordell_weil_torsion(X::EllipticSurface) -> Vector{EllipticCurvePoint}

Return the torsion part of the Mordell-Weil group of the generic fiber of ``X``.
"""
@attr Vector{EllipticCurvePoint} function mordell_weil_torsion(X::EllipticSurface)
  S = X
  E = generic_fiber(S)
  O = E([0,1,0])
  N = trivial_lattice(S)[2]
  tors = EllipticCurvePoint[]
  d = det(N)
  for p in prime_divisors(d)
    if valuation(d, p) == 1
      continue
    end
    r = 1
    i = 0
    dp = typeof(O)[]
    while true
      i = i+1
      @vprint :EllipticSurface 2 "computing $(p^i)-torsion"
      dp = division_points(O, p^i)
      if length(dp) == r
        break
      end
      r = length(dp)
    end
    for pt in dp
      if pt != O
        push!(tors, pt)
      end
    end
  end
  return tors
end

@doc raw"""
    section(X::EllipticSurface, P::EllipticCurvePoint) -> EllipticSurfaceSection

Given a rational point ``P\in E(C)`` of the generic fiber ``E/C`` of ``\pi\colon X \to C``,
return its closure in ``X`` as a `AbsWeilDivisor`.
"""
function section(X::EllipticSurface, P::EllipticCurvePoint)
  if iszero(P[1])&&iszero(P[3])
    return zero_section(X)
  end
  return EllipticSurfaceSection(X, P)
end


function _section_on_weierstrass_ambient_space(X::EllipticSurface, P::EllipticCurvePoint)
  S0,incS0 = weierstrass_model(X)
  X0 = codomain(incS0)
  if P[3] == 0
    # zero section
    V = X0[1][3]
    (z,x,t) = coordinates(V)
    return IdealSheaf(X0, V, [x,z])
  end
  U = X0[1][1]
  (x,y,t) = coordinates(U)
  b = P
  return ideal_sheaf(X0,U,[OO(U)(i) for i in [x*denominator(b[1])(t)-numerator(b[1])(t),y*denominator(b[2])(t)-numerator(b[2])(t)]]; check=false)
end


@doc raw"""
    zero_section(X::EllipticSurface) -> AbsWeilDivisor

Return the zero section of the relatively minimal elliptic
fibration ``\pi\colon X \to C``.
"""
@attr EllipticSurfaceSection zero_section(X::EllipticSurface) = EllipticSurfaceSection(X, generic_fiber(X)([0,1,0]))


@doc raw"""
    basis_representation(X::EllipticSurface, D::AbsWeilDivisor)

Return the vector representing the numerical class of ``D``
with respect to the basis of the ambient space of `algebraic_lattice(X)`.
"""
function basis_representation(X::EllipticSurface, D::AbsWeilDivisor)
  basis_ambient,_, NS = algebraic_lattice(X)
  G = gram_matrix(ambient_space(NS))
  n = length(basis_ambient)
  v = zeros(ZZRingElem, n)
  @vprint :EllipticSurface 3 "computing basis representation of $D\n"
  kk = base_ring(X)
  if iszero(characteristic(kk)) && has_attribute(X, :good_reduction_map)
    X_red_raw, bc = raw_good_reduction(X)
    red_dict = IdDict{AbsWeilDivisor, AbsWeilDivisor}(D=>_reduce_as_prime_divisor(bc, D) for D in basis_ambient)
    D_red = _reduce_as_prime_divisor(bc, D)
    for (i, E) in enumerate(basis_ambient)
      @vprintln :EllipticSurface 4 "intersecting in positive characteristic with $(i): $(basis_ambient[i])"
      v[i] = intersect(red_dict[E], D_red)
    end
  else
    for i in 1:n
      @vprintln :EllipticSurface 4 "intersecting with $(i): $(basis_ambient[i])"

      v[i] = intersect(basis_ambient[i], D)
    end
  end
  @vprint :EllipticSurface 3 "done computing basis representation\n"
  return v*inv(G)
end


###################################################################################################################
#
# Fibers and trivial lattice
#
###################################################################################################################

@doc raw"""
    trivial_lattice(X::EllipticSurface) -> Vector{AbsWeilDivisor}, ZZMatrix

Return a basis for the trivial lattice as well as its Gram matrix.

The trivial lattice is the sublattice of the numerical lattice spanned by fiber components and
the zero section of ``X``.

!!! warning 
    The ordering of the basis is in general not canonical and may change in between sessions. 
"""
function trivial_lattice(X::EllipticSurface)
  T = _trivial_lattice(X)[1:2]
  return T
end

@doc raw"""
    _trivial_lattice(X::EllipticSurface; reducible_singular_fibers_in_PP1=_reducible_fibers_disc(S))
  
Internal function. Returns a list consisting of:
- basis of the trivial lattice
- gram matrix
- fiber_components without multiplicities

The keyword argument `reducible_singular_fibers_in_PP1` must be a list of vectors of length `2` over 
the base field representing the points in projective space over which there are reducible fibers.
Specify it to force this ordering of the basis vectors of the ambient space of the `algebraic_lattice`
"""
@attr Any function _trivial_lattice(X::EllipticSurface; reducible_singular_fibers_in_PP1=_reducible_fibers_disc(X))
  S = X
  O = zero_section(S)
  pt0, F = fiber(S)
  set_attribute!(components(O)[1], :_self_intersection, -euler_characteristic(S))
  basisT = [F, O]
  grams = [ZZ[0 1;1 -euler_characteristic(S)]]
  sing = reducible_singular_fibers_in_PP1
  f = [[pt, fiber_components(S,pt)] for  pt in sing]
  fiber_componentsS = []
  for (pt, ft) in f
    @vprint :EllipticSurface 2 "normalizing fiber: over $pt \n"
    Ft0 = standardize_fiber(S, ft)
    @vprint :EllipticSurface 2 "$(Ft0[1]) \n"
    append!(basisT , Ft0[3][2:end])
    push!(grams,Ft0[4][2:end,2:end])
    push!(fiber_componentsS, vcat([pt], collect(Ft0)))
  end
  G = block_diagonal_matrix(grams)
  # make way for some more pretty printing
  for (pt,root_type,_,comp) in fiber_componentsS
    for (i,I) in enumerate(comp)
      name = string(root_type[1], root_type[2])
      set_attribute!(components(I)[1], :name, string("Component ", name, "_", i-1," of fiber over ", Tuple(pt)))
      set_attribute!(components(I)[1], :_self_intersection, -2)
    end
  end
  return basisT, G, fiber_componentsS
end

function _reducible_fibers_disc(X::EllipticSurface; sort::Bool=true)
  E = generic_fiber(X)
  j = j_invariant(E)
  d = numerator(discriminant(E))
  kt = parent(d)
  k = coefficient_ring(kt)
  sing = Vector{elem_type(k)}[]
  for (p,v) in factor(d)
    if v == 1
      continue
    end
    r = k.(roots(p))
    if length(r) == 0
      error("not all reducible fibers are visible over $(base_ring(S))")
    end
    @assert length(r) ==1
    rt = r[1]
    if v == 2
      # not a type II fiber
      if j!=0
        push!(sing, [rt,k(1)])
      end
    end
    if v > 2
      push!(sing, [rt,k(1)])
    end
  end
  if sort
    sort!(sing, by=x->by_total_order(x[1]))
  end
  # fiber over infinity is always last (if it is there)
  if degree(d) <= 12*euler_characteristic(X) - 2
    push!(sing, k.([1, 0]))
  end
  return sing
end

function by_total_order(x::FqFieldElem) 
  return [lift(ZZ, i) for i in absolute_coordinates(x)]
end
  
by_total_order(x::QQFieldElem) = x
  
function by_total_order(x::NumFieldElem)
  return absolute_coordinates(x)
end 


@doc raw"""
    reducible_fibers(X::EllipticSurface)

Return the reducible fibers of ``X``.

The output format is the following:
A list `[F1, ..., Fn]` where each entry `Fi` represents a reducible fiber.

The list ``F`` has the following entries:
- A point ``P \in \mathbb{P}^{1}`` such that ``F = \pi^{-1}(P)``;
- The ADE-type of the fiber;
- The fiber ``F`` as a Weil divisor, including its multiplicities;
- The irreducible components of the fiber. The first component intersects the zero section;
- Their intersection matrix.
"""
function reducible_fibers(X::EllipticSurface)
  return _trivial_lattice(X)[3]
end


@doc raw"""
    standardize_fiber(X::EllipticSurface, f::Vector{<:AbsWeilDivisor})

Internal method. Used to prepare for [`reducible_fibers`](@ref).
`f` must be the list of the components of the reducible fiber `F`.
Output a list of tuples with each tuple as follows
- the root type of ``F``, e.g. `(:A, 3)`
- the class of ``F`` as a divisor with the appropriate multiplicities
- the irreducible components `[F0,...Fn]` of `F` sorted such that the first entry `F0` is the one intersecting the zero section. The others are sorted in some standard way
- gram matrix of the intersection of [F0,...,Fn], it is an extended ADE-lattice.
"""
function standardize_fiber(X::EllipticSurface, f::Vector{<:AbsWeilDivisor})
  S = X
  @hassert :EllipticSurface 2 all(is_prime(i) for i in f)
  f = copy(f)
  O = components(zero_section(S))[1]
  local f0
  for (i,D) in enumerate(f)
    if !isone(O+components(D)[1])
      f0 = D
      deleteat!(f,i)
      break
    end
  end
  r = length(f)
  G = -2*identity_matrix(ZZ, r)
  @vprintln :EllipticSurface 2 "computing intersections:"
  for i in 1:r
    @vprint :EllipticSurface 3 "\nrow $(i): \n"
    for j in 1:i-1
      @vprint :EllipticSurface 4 "$(j) "
      # we know the intersections are 0 or 1, so we can replace the line below by a shortcut.
      # G[i, j] = G[j, i] = intersect(f[i], f[j])
      # In the examples treated, this led to roughly a factor 3 in speed and memory consumption. 
      if isone(components(f[i])[1]+components(f[j])[1])
        G[i,j] = 0
      else
        G[i,j] = 1
      end
      G[j,i] = G[i,j]
    end
  end
  L = integer_lattice(gram=G)
  rt,_ = root_lattice_recognition(L)
  @assert length(rt)==1
  rt = rt[1]
  R = root_lattice(rt[1], rt[2])
  b, I = _is_equal_up_to_permutation_with_permutation(G, -gram_matrix(R))
  @assert b
  gensF = vcat([f0], f[I])
  Gext, v = extended_ade(rt[1],rt[2])
  Fdiv = sum(v[i]*gensF[i] for i in 1:length(gensF))
  return rt, Fdiv, gensF, Gext
end

@doc raw"""
    fiber_cartier(X::EllipticSurface, P::Vector = ZZ.([0,1])) -> EffectiveCartierDivisor

Return the fiber of ``\pi\colon X \to C`` over ``P\in C`` as a Cartier divisor.
"""
function fiber_cartier(X::EllipticSurface, P::Vector = ZZ.([0,1]))
  S = X
  S0,_ = weierstrass_model(S)
  underlying_scheme(S) # cache stuff
  D = IdDict{AbsAffineScheme, RingElem}()
  k = base_ring(S0)
  P = k.(P)

  if P[1]!=0 && P[2]!=0
    t0 = P[1] * inv(P[2])
    for i in 1:3
      U = S0[1][i]
      (_,_,t) = coordinates(U)
      D[U] = t-t0
    end
    s0 = inv(t0)
    for i in 4:6
      U = S0[1][i]
      (_,_,s) = coordinates(U)
      D[U] = s - s0
    end
  elseif P[1] != 0
    # it is the fiber at [1:0]
    for i in 4:6
      U = S0[1][i]
      (_,_,s) = coordinates(U)
      D[U] = s
    end
    for i in 1:3
      U = S0[1][i]
      (_,_,t) = coordinates(U)
      D[U] = one(parent(t))
    end
  elseif P[2] != 0
    # the fiber at [0:1]
    for i in 1:3
      U = S0[1][i]
      (_,_,t) = coordinates(U)
      D[U] = t
    end
    for i in 4:6
      U = S0[1][i]
      (_,_,s) = coordinates(U)
      D[U] = one(parent(s))
    end
  else
    error("[0,0] is not a point in projective space")
  end
  F = EffectiveCartierDivisor(S0, D, trivializing_covering=S0[1], check=false)
  return pullback(S.blowup)(F)
end

@doc raw"""
    fiber_components(X::EllipticSurface, P) -> Vector{<:AbsWeilDivisor}

Return the fiber components of the fiber over the point ``P \in C``.
"""
function fiber_components(X::EllipticSurface, P; algorithm=:exceptional_divisors)
  S = X
  @vprint :EllipticSurface 2 "computing fiber components over $(P)\n"
  P = base_ring(S).(P)
  W = codomain(S.inc_Weierstrass)
  Fcart = fiber_cartier(S, P)
  if isone(P[2])
    U = default_covering(W)[1]
    (x,y,t) = coordinates(U)
    F = PrimeIdealSheafFromChart(W, U, ideal(t - P[1]))
  elseif isone(P[1])
    U = default_covering(W)[4]
    (x,y,s) = coordinates(U)
    F = PrimeIdealSheafFromChart(W, U, ideal(s - P[2]))
  end 
  FF = ideal_sheaf(Fcart)
  EE = exceptional_divisors(S)
  EP = filter(E->issubset(FF, E), EE)
  for bl in S.ambient_blowups
    F = strict_transform(bl, F)
  end
  F = pullback(S.inc_Y, F)
  F = weil_divisor(F, ZZ)
  fiber_components = [weil_divisor(E, ZZ; check=false) for E in EP]
  push!(fiber_components, F)
  return fiber_components
end
  
@attr Vector{<:AbsIdealSheaf} function exceptional_divisors(X::EllipticSurface)
  S = X
  PP = AbsIdealSheaf[]
  @vprintln :EllipticSurface 2 "computing exceptional divisors"
  # If we have resolution_strategy=:simultaneous, then the following is a non-trivial preprocessing step.
  ambient_pts = reduce(vcat, maximal_associated_points.(ideal_sheaf.(S.ambient_exceptionals)))
  for I in ambient_pts
    @vprintln :EllipticSurface 4 "decomposing divisor "
    mp = maximal_associated_points(pullback(S.inc_Y, I); use_decomposition_info=true)
    append!(PP, mp)
  end 
  @vprintln :EllipticSurface 3 "done"
  return PP
end
  
function fiber(X::EllipticSurface)
  b, pt, F = irreducible_fiber(X)
  if b
    W = weil_divisor(F)
    # we manually set the self-intersection
    set_attribute!(components(W)[1], :_self_intersection, 0)
  else
    # all fibers are reducible pick the one over [0 : 1]
    k = base_ring(X)
    pt = [k(0),k(1)]
    f = fiber_components(X, pt)
    fiber_type, W, componentsW, gramW = standardize_fiber(X, f)
  end
  set_attribute!(W, :name=> "Fiber over ($(pt[1]), $(pt[2]))")
  return pt, W
end


@doc raw"""
    irreducible_fiber(X::EllipticSurface) -> Bool, Point, EffectiveCartierDivisor

Return an irreducible fiber as a cartier divisor and whether it exists.

The return value is a triple `(b, pt, F)` where

- `b` is `true` if an irreducible fiber exists over the base field of `S`
- `pt` the base point of the fiber
- `F` the irreducible fiber which projects to `pt`
"""
function irreducible_fiber(X::EllipticSurface)
  S = X
  W = weierstrass_model(S)
  d = numerator(discriminant(generic_fiber(S)))
  kt = parent(d)
  k = coefficient_ring(kt)
  r = [k.(roots(i[1])) for i in factor(d) if i[2]>=2]
  sing = reduce(append!,r, init=[])
  pt = k.([0,0]) # initialize
  found = false
  if is_finite(k)
    for i in k
      if !(i in sing)  # true if the fiber over [i,1] is irreducible
        pt = k.([i,1])
        found = true
        break
      end
    end
    if !found && (degree(d) >= 12*euler_characteristic(S) - 1)  # irreducible at infinity?
      pt = k.([1, 0])
      found = true
    end
  else
    i = k(0)
    while true
      i = i+1
      if !(i in sing)
        pt = k.([i,1])
        found = true
        break
      end
    end
  end
  F = fiber_cartier(S, pt)
  return found, pt, F
end
