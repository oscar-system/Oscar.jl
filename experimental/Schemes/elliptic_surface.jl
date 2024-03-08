export elliptic_surface, trivial_lattice, weierstrass_model, weierstrass_chart, algebraic_lattice, zero_section, section, weierstrass_contraction, fiber_components, generic_fiber, reducible_fibers, fibration_type, mordell_weil_lattice, elliptic_parameter, set_mordell_weil_basis!, EllipticSurface, weierstrass_chart_on_minimal_model, transform_to_weierstrass

@doc raw"""
    EllipticSurface{BaseField<:Field, BaseCurveFieldType} <: AbsCoveredScheme{BaseField}

The type of a relatively minimal elliptic surface.

A genus $1$-fibration is a proper map
$\pi \colon X \to C$ to a curve $C$ whose fibers are curves of
(arithmetic) genus $1$.
The fibration is relatively minimal if its fibers do not contain any ``(-1)``-curves.
We call the fibration elliptic if it comes equipped with a section.
This turns the generic fiber of $\pi$ into an elliptic curve $E/k(C)$ where
$k(C)$ is the function field of the curve $C$.

For now functionality is restricted to $C = \mathbb{P}^1$.
"""
@attributes mutable struct EllipticSurface{BaseField<:Field, BaseCurveFieldType} <: AbsCoveredSurface{BaseField}
  Y::CoveredScheme{BaseField}  # the underlying_scheme
  E::EllipticCurve{BaseCurveFieldType}
  MWL::Vector{EllipticCurvePoint{BaseCurveFieldType}} # basis for the mordell weil group
  MWLtors::Vector{EllipticCurvePoint{BaseCurveFieldType}} # torsion sections
  Weierstrasschart::AbsAffineScheme
  Weierstrassmodel::CoveredScheme
  inc_Weierstrass::CoveredClosedEmbedding # inclusion of the weierstrass chart in its ambient projective bundle
  inc_Y::CoveredClosedEmbedding # inclusion of Y in its ambient blown up projective bundle
  euler_characteristic::Int
  # the following are temporary until we have a dedicated type for
  # iterated blow ups
  blowup::AbsCoveredSchemeMorphism
  blowups::Vector{<:AbsCoveredSchemeMorphism}
  # exceptionals not used for now
  ambient_blowups::Vector{<:BlowupMorphism}
  ambient_exceptionals::Vector{<:EffectiveCartierDivisor}

  function EllipticSurface(generic_fiber::EllipticCurve{F}, euler_characteristic::Int, mwl_basis::Vector{<:EllipticCurvePoint}) where F
    B = typeof(coefficient_ring(base_ring(base_field(generic_fiber))))
    S = new{B,F}()
    S.E = generic_fiber
    S.MWL = mwl_basis
    S.euler_characteristic = euler_characteristic
    set_attribute!(S, :is_irreducible=>true)
    set_attribute!(S, :is_reduced=>true)
    set_attribute!(S, :is_integral=>true)
    return S
  end

end

base_ring(X::EllipticSurface) = coefficient_ring(base_ring(base_field(generic_fiber(X))))

@doc raw"""
    set_mordell_weil_basis!(X::EllipticSurface, mwl_basis::Vector{EllipticCurvePoint})

Set generators for the Mordell-Weil lattice of ``X`` or at least of a sublattice.

This invalidates previous computations depending on the generators of the
mordell weil lattice such as the `algebraic_lattice`. Use with care.
"""
function set_mordell_weil_basis!(X::EllipticSurface, mwl_basis::Vector{<:EllipticCurvePoint})
  @req all(parent(P) == generic_fiber(X) for P in mwl_basis) "points must lie on the generic fiber"
  X.MWL = mwl_basis
  # clear old computations
  if has_attribute(X, :algebraic_lattice)
    deleteat!(X.__attrs, :algebraic_lattice)
  end
  if has_attribute(X, :mordell_weil_lattice)
    deleteat!(X.__attrs, :mordell_weil_lattice)
  end
end

@doc raw"""
    elliptic_surface(generic_fiber::EllipticCurve, euler_characteristic::Int, mwl_basis::Vector{<:EllipticCurvePoint}=EllipticCurvePoint[]) -> EllipticSurface

Return the relatively minimal elliptic surface with generic fiber ``E/k(t)``.

This is also known as the Kodaira-Néron model of ``E``.

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

julia> Base.show(stdout, X3)
Elliptic surface with generic fiber -x^3 + y^2 - t^7 + 2*t^6 - t^5
```
"""
function elliptic_surface(generic_fiber::EllipticCurve{BaseField}, euler_characteristic::Int,
                          mwl_basis::Vector{<:EllipticCurvePoint}=EllipticCurvePoint[]) where {
                          BaseField <: FracFieldElem{<:PolyRingElem{<:FieldElem}}}
  @req all(parent(i)==generic_fiber for i in mwl_basis) "not a vector of points on $(generic_fiber)"
  S = EllipticSurface(generic_fiber, euler_characteristic, mwl_basis)
  return S
end

kodaira_neron_model(E::EllipticCurve) = elliptic_surface(E)

function underlying_scheme(S::EllipticSurface)
  if isdefined(S,:Y)
    return S.Y
  end
  # trigger the computation
  weierstrass_contraction(S)
  return underlying_scheme(S)
end

@doc raw"""
    generic_fiber(S::EllipticSurface) -> EllipticCurve

Return the generic fiber as an elliptic curve.
"""
generic_fiber(S::EllipticSurface) = S.E

@doc raw"""
    weierstrass_chart(X::EllipticSurface)

Return the Weierstrass chart of ``X`` on its `weierstrass_model`.
"""
weierstrass_chart(X::EllipticSurface) = weierstrass_model(X)[1][1][1]

@doc raw"""
    weierstrass_chart_on_minimal_model(X::EllipticSurface)

Return an affine chart ``U`` of ``X`` which is isomorphic to the `weierstrass_chart` 
of ``X`` on its `weierstrass_model`, but with all singular fibers removed.

More precisely, the affine coordinates of ``U`` are ``(x,y,t)`` and the chart is 
constructed as the vanishing locus of
``y^2 + a_1(t) xy + a_3 y = x^3 + a_2 x^2 + a_4 x + a_6``
minus the reducible singular fibers.
"""
weierstrass_chart_on_minimal_model(X::EllipticSurface) = X[1][1]

@doc raw"""
    euler_characteristic(X::EllipticSurface) -> Int

Return $\chi(\mathcal{O}_X)$.
"""
euler_characteristic(X::EllipticSurface) = X.euler_characteristic

@doc raw"""
    algebraic_lattice(X) -> Vector{WeilDivisor}, ZZLat

Return the sublattice `L` of ``Num(X)`` spanned by fiber components,
torsion sections and the sections provided at the construction of ``X``.

The first return value is the basis of the ambient space of `L`.
The second consists of additional generators for `L` coming from torsion sections.
The third is ``L``.
"""
@attr function algebraic_lattice(X)
  mwl_basis = X.MWL
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
  GAinv = inv(change_base_ring(QQ,GA))
  # primitive closure of the trivial lattice comes from torsion sections
  tors = [section(X, h) for h in mordell_weil_torsion(X)]
  torsV = QQMatrix[]
  for T in tors
    @vprint :EllipticSurface 2 "computing basis representation of $(T)\n"
    vT = zero_matrix(ZZ, 1, n)
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
    push!(torsV, vT*GAinv) # not necessarily integral .... todo
  end
  gen_tors = zip(tors, torsV)
  push!(torsV, identity_matrix(QQ,n))
  V = quadratic_space(QQ,GA)
  L = lattice(V, reduce(vcat,torsV), isbasis=false)
  return basisA, collect(gen_tors), L
end

@doc raw"""
    mordell_weil_lattice(S::EllipticSurface) -> Vector{EllipticCurvePoint}, ZZLat

Return the (sublattice) of the Mordell-Weil lattice of ``S``  spanned
by the sections of ``S`` supplied at its construction.
"""
@attr ZZLat function mordell_weil_lattice(S::EllipticSurface)
  NS = algebraic_lattice(S)[3]
  t = length(trivial_lattice(S)[1])
  trivNS = basis_matrix(NS)[1:t,:]
  R = basis_matrix(NS)[t+1:end,:]
  V = ambient_space(NS)
  P = orthogonal_projection(V, trivNS)
  mwl = rescale(P(NS),-1)
  #Todo basis
  return lll(mwl)
end

@doc raw"""
    mordell_weil_torsion(S::EllipticSurface) -> Vector{EllipticCurvePoint}

Return the torsion part of the Mordell-Weil group of the generic fiber of ``S``.
"""
@attr function mordell_weil_torsion(S::EllipticSurface)
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
    while true
      i = i+1
      @vprint :EllipticSurface 2 "computing $(p^i)-torsion"
      global dp = division_points(O, p^i)
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

function Base.show(io::IO, S::EllipticSurface)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Elliptic surface")
  else
    E = generic_fiber(S)
    print(io, "Elliptic surface with generic fiber ", equation(E))
  end
end

function Base.show(io::IO, ::MIME"text/plain", S::EllipticSurface)
  io = pretty(io)
  println(io, "Elliptic surface")
  println(io, Indent(), "over ", Lowercase(), base_ring(S))
  println(io, Dedent(), "with generic fiber")
  print(io, Indent(), Lowercase(), equation(generic_fiber(S)), Dedent())
  if isdefined(S, :Y)
    println(io)
    println(io, "and relatively minimal model")
    print(io, Indent(), Lowercase(), S.Y, Dedent())
  end
  print(io, Dedent())
end

@doc raw"""
    weierstrass_model(X::EllipticSurface) -> CoveredScheme, CoveredClosedEmbedding

Return the Weierstrass model ``S`` of ``X`` and the inclusion in
its ambient projective bundle
$$S\subseteq \mathbb{P}( \mathcal{O}_{\mathbb{P}^1}(-2s) \oplus \mathcal{O}_{\mathbb{P}^1}(-3s) \oplus \mathcal{O}_{\mathbb{P}^1}).$$
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

  P_proj = projectivization(bundleE, var_names=["z", "x", "y"])
  P = covered_scheme(P_proj)
  @assert has_decomposition_info(default_covering(P))

  # Create the singular Weierstrass model S of the elliptic K3 surface X
  a = a_invariants(E)
  U = affine_charts(P)[1]  # the standard Weierstrass chart
  (x, y, t) = gens(OO(U))
  @assert all(denominator(i)==1 for i in a)
  a = [numerator(a)(t) for a in a]
  (a1,a2,a3,a4,a6) = a
  ft = y^2  + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x+a6)
  I = IdealSheaf(P, U, [ft])

  inc_S = CoveredClosedEmbedding(P, I)
  Scov = domain(inc_S)  # The ADE singular elliptic K3 surface
  X.Weierstrasschart = Scov[1][1]

  X.Weierstrassmodel = Scov
  X.inc_Weierstrass = inc_S

  set_attribute!(Scov, :is_irreducible=>true)
  set_attribute!(Scov, :is_reduced=>true)
  set_attribute!(Scov, :is_integral=>true)
  return Scov, inc_S
end

@doc raw"""
    _separate_singularities!(X::EllipticSurface) -> Covering

Create a covering of the ambient projective bundle $P$
of the Weierstrass model $S$ of $X$ such that each chart
(of $X$) contains at most one singular point of $S$.
Append this covering to the list of coverings of $X$ and return it.
"""
function _separate_singularities!(X::EllipticSurface)
  S, inc_S = weierstrass_model(X)
  P = codomain(inc_S)

  I_sing = ideal_sheaf_of_singular_locus(S)
  I_sing_P = radical(pushforward(inc_S)(I_sing))

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
    local Ising = I_sing_P(W)
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

@doc raw"""
    weierstrass_contraction(X::EllipticSurface) -> SchemeMor

Return the contraction morphism of ``X`` to its Weierstrass model.

This triggers the computation of the `underlying_scheme` of ``X``
as a blowup from its Weierstrass model. It may take a few minutes.
"""
function weierstrass_contraction(X::EllipticSurface)
  Y = X
  if isdefined(Y, :blowup)
    return Y.blowup
  end
  S, inc_S = weierstrass_model(Y)
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
  inc_Y0 = inc_S


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
    @vprint :EllipticSurface 2 "computing singular locus\n"
    I_sing_Y0 = ideal_sheaf_of_singular_locus(Y0)
    @vprint :EllipticSurface 2 "decomposing singular locus\n"
    I_sing_Y0 = maximal_associated_points(I_sing_Y0)
    I_sing_X0 = pushforward(inc_Y0).(I_sing_Y0)
    @vprint :EllipticSurface 1 "number of singular points: $(length(I_sing_Y0))\n"
    if length(I_sing_Y0)==0
      # stop if smooth
      break
    end
    # make sure there is only one singular point per chart
    if count == 1
      cov = Crefined
    else
      # the following leads to difficult bugs
      #=
      cov0 = simplified_covering(X0)
      cov1 = _separate_disjoint_components(I_sing_X0, covering=cov0)
      cov = _one_patch_per_component(cov1, I_sing_X0)
      push!(X0.coverings, cov)
      @assert has_decomposition_info(default_covering(X0))
      inherit_decomposition_info!(X0, cov)
      @assert has_decomposition_info(cov)
      =#
      cov = simplified_covering(X0)
      #inherit_decomposition_info!(cov, X0)
    end
    # take the first singular point and blow it up
    J = radical(I_sing_X0[1]) # radical to have small number of generators
    pr_X1 = blow_up(J, covering=cov, var_name=varnames[1+mod(count, length(varnames))])

    # Set the attribute so that the strict_transform does some extra work
    isomorphism_on_open_subset(pr_X1)

    X1 = domain(pr_X1)
    @vprint :EllipticSurface 1 "$(X1)\n"
    E1 = exceptional_divisor(pr_X1)

    @vprint :EllipticSurface 2 "computing strict transforms\n"
    # compute the exceptional divisors
    ambient_exceptionals = EffectiveCartierDivisor[strict_transform(pr_X1, e) for e in ambient_exceptionals]
    # move the divisors coming originally from S up to the next chart
    push!(ambient_exceptionals, E1)

    Y1, inc_Y1, pr_Y1 = strict_transform(pr_X1, inc_Y0)

    push!(projectionsX, pr_X1)
    push!(projectionsY, pr_Y1)
    simplify!(Y1)

    # set up for the next iteration
    Y0 = Y1
    inc_Y0 = inc_Y1
    X0 = X1
    set_attribute!(Y0, :is_irreducible=> true)
    set_attribute!(Y0, :is_reduced=>true)
    set_attribute!(Y0, :is_integral=>true)
    set_attribute!(X0, :is_irreducible=> true)
    set_attribute!(X0, :is_reduced=>true)
    set_attribute!(X0, :is_integral=>true)
  end
  Y.Y = Y0
  Y.blowups = projectionsY

  # We need to rewrap the last maps so that the domain is really Y
  last_pr = pop!(projectionsY)
  last_pr_wrap = CoveredSchemeMorphism(Y, codomain(last_pr), covering_morphism(last_pr))
  set_attribute!(last_pr_wrap, :isomorphism_on_open_subset, get_attribute(last_pr, :isomorphism_on_open_subset))

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

#  global divisors0 = [strict_transform(pr_X1, e) for e in divisors0]
#exceptionals_res = [pullback(inc_Y0)(e) for e in exceptionals]
@doc raw"""
    _trivial_lattice(S::EllipticSurface)

Internal function. Returns a list consisting of:
- basis of the trivial lattice
- gram matrix
- fiber_components without multiplicities
"""
@attr function _trivial_lattice(S::EllipticSurface)
  #=
  inc_Y = S.inc_Y
  X = codomain(inc_Y)
  exceptionals_res = [pullback(inc_Y0)(e) for e in exceptionals]
  ExWeil = WeilDivisor.(exceptional_res)
  tmp = []
  ExWeil = reduce(append!, [components(i) for i in ExWeil], init= tmp)
  =#
  W = weierstrass_model(S)
  d = numerator(discriminant(generic_fiber(S)))
  j = j_invariant(generic_fiber(S))
  kt = parent(d)
  k = coefficient_ring(kt)
  # find the reducible fibers
  sing = elem_type(k)[]
  for (p,v) in factor(d)
    if v == 1
      continue
    end
    r = k.(roots(p))
    if length(r) == 0
      error("not all reducible fibers are visible over $(base_ring(X))")
    end
    @assert length(r) ==1
    rt = r[1]
    if v == 2
      # not a type II fiber
      if j!=0
        push!(sing, rt)
      end
    end
    if v > 2
      push!(sing, rt)
    end
  end
  # the reducible fibers are over the points in sing
  # and possibly the point at infinity
  f = [[k.([i,1]), fiber_components(S,[i,k(1)])] for  i in sing]
  if degree(d) <= 12*euler_characteristic(S) - 2
    pt = k.([1, 0])
    push!(f, [pt, fiber_components(S, pt)])
  end

  O = zero_section(S)
  pt0, F = fiber(S)
  set_attribute!(components(O)[1], :_self_intersection, -euler_characteristic(S))


  basisT = [F, O]

  grams = [ZZ[0 1;1 -euler_characteristic(S)]]
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

@doc raw"""
    trivial_lattice(X::EllipticSurface) -> Vector{WeilDivisor}, ZZMatrix

Return a basis for the trivial lattice as well as its gram matrix.

The trivial lattice is the lattice spanned by fiber components and
the zero section of $X$.
"""
function trivial_lattice(X::EllipticSurface)
  T = _trivial_lattice(X)[1:2]
  return T
end

@doc raw"""
    reducible_fibers(S::EllipticSurface)

Return the reducible fibers of $S$.

The output format is the following:
A list [F1, ..., Fn] where each entry Fi represents a reducible fiber.

The list $F$ has the following entries:
- A point $P \in \mathbb{P}^{1}$ such that $F = \pi^{-1}(P)$;
- The ADE-type of the fiber;
- The fiber $F$ as a Weil divisor, including its multiplicities;
- The irreducible components of the fiber. The first component intersects the zero section;
- Their intersection matrix.
"""
function reducible_fibers(S::EllipticSurface)
  return _trivial_lattice(S)[3]
end


@doc raw"""
    standardize_fiber(S::EllipticSurface, f::Vector{<:WeilDivisor})

Internal method. Used to prepare for [`reducible_fibers`](@ref).
`f` must be the list of the components of the reducible fiber `F`.
Output a list of tuples with each tuple as follows
- the root type of ``F``, e.g. `(:A, 3)`
- the class of ``F`` as a divisor with the appropriate multiplicities
- the irreducible components `[F0,...Fn]` of `F` sorted such that the first entry `F0` is the one intersecting the zero section. The others are sorted in some standard way
- gram matrix of the intersection of [F0,...,Fn], it is an extended ADE-lattice.
"""
function standardize_fiber(S::EllipticSurface, f::Vector{<:WeilDivisor})
  @req all(is_prime(i) for i in f) "not a vector of prime divisors"
  f = copy(f)
  O = components(zero_section(S))[1]
  for (i,D) in enumerate(f)
    if !isone(O+components(D)[1])
      global f0 = D
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
      # we know the intersections are 0 or 1
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
    fiber_cartier(S::EllipticSurface, P::Vector = ZZ.([0,1])) -> EffectiveCartierDivisor

Return the fiber of $\pi\colon X \to C$ over $P\in C$ as a Cartier divisor.
"""
function fiber_cartier(S::EllipticSurface, P::Vector = ZZ.([0,1]))
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
    fiber_components(S::EllipticSurface, P) -> Vector{<:WeilDivisor}

Return the fiber components of the fiber over the point $P \in C$.
"""
function fiber_components(S::EllipticSurface, P)
  @vprint :EllipticSurface 2 "computing fiber components over $(P)\n"
  F = fiber_cartier(S, P)
  @vprint :EllipticSurface 2 "decomposing fiber   "
  comp = maximal_associated_points(ideal_sheaf(F))
  @vprint :EllipticSurface 2 "done decomposing fiber\n"
  return [weil_divisor(c, check=false) for c in comp]
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
    irreducible_fiber(S::EllipticSurface) -> Bool, Point, EffectiveCartierDivisor

Return an irreducible fiber as a cartier divisor and whether it exists.

The return value is a triple `(b, pt, F)` where

- `b` is `true` if an irreducible fiber exists over the base field of `S`
- `pt` the base point of the fiber
- `F` the irreducible fiber which projects to `pt`
"""
function irreducible_fiber(S::EllipticSurface)
  W = weierstrass_model(S)
  d = numerator(discriminant(generic_fiber(S)))
  kt = parent(d)
  k = coefficient_ring(kt)
  r = [k.(roots(i[1])) for i in factor(d) if i[2]>=2]
  sing = reduce(append!,r, init=[])
  pt = k.([0,0]) # initialize
  found = false
  if degree(d) >= 12*euler_characteristic(S) - 1  # irreducible at infinity?
    pt = k.([1, 0])
    found = true
  else
    if is_finite(k)
      for i in k
        if !(i in sing)  # true if the fiber over [i,1] is irreducible
          pt = k.([i,1])
          found = true
          break
        end
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
  end
  F = fiber_cartier(S, pt)
  return found, pt, F
end

@doc raw"""
    section(X::EllipticSurface, P::EllipticCurvePoint)

Given a rational point $P\in E(C)$ of the generic fiber $E/C$ of $\pi\colon X \to C$,
return its closure in $X$ as a `WeilDivisor`.
"""
function section(X::EllipticSurface, P::EllipticCurvePoint)
  if iszero(P[1])&&iszero(P[3])
    return zero_section(X)
  end
  return _section(X, P)
end

function _section(X::EllipticSurface, P::EllipticCurvePoint)
  @vprint :EllipticSurface 3 "Computing a section from a point on the generic fiber\n"
  S0,incS0 = weierstrass_model(X)
  X0 = codomain(incS0)
  if P[3] == 0
    # zero section
    V = X0[1][3]
    (z,x,t) = coordinates(V)
    PX = IdealSheaf(X0, V, [x,z])
  else
    U = X0[1][1]
    (x,y,t) = coordinates(U)
    b = P
    PX = ideal_sheaf(X0,U,[OO(U)(i) for i in [x*denominator(b[1])(t)-numerator(b[1])(t),y*denominator(b[2])(t)-numerator(b[2])(t)]])
  end
  for f in X.ambient_blowups
    PX = strict_transform(f , PX)
  end
  PY = pullback(X.inc_Y, PX)
  set_attribute!(PY, :name, string("section: (",P[1]," : ",P[2]," : ",P[3],")"))
  set_attribute!(PY, :_self_intersection, -euler_characteristic(X))
  W =  WeilDivisor(PY, check=false)
  set_attribute!(W, :is_prime=>true)
  return W
end

@doc raw"""
    zero_section(S::EllipticSurface) -> WeilDivisor

Return the zero section of the relatively minimal elliptic
fibration \pi\colon X \to C$.
"""
@attr zero_section(S::EllipticSurface) = _section(S, generic_fiber(S)([0,1,0]))

################################################################################
#
# Some linear systems on elliptic surfaces
#
################################################################################

@doc raw"""
    _prop217(E::EllipticCurve, P::EllipticCurvePoint, k)

Compute a basis for the linear system
``|O + P + kF|``
on the  minimal elliptic (K3) surface defined by E.
Here F is the class of a fiber O the zero section
and P any non-torsion section.

The return value is a list of pairs ``(a(t),b(t))``

```jldoctest
julia> kt,t = polynomial_ring(GF(29),:t);

julia> ktfield = fraction_field(kt);

julia> bk = [((17*t^4 + 23*t^3 + 18*t^2 + 2*t + 6, 8*t^5 + 2*t^4 + 6*t^3 + 25*t^2 + 24*t + 5 )),
             ((17*t^6 + 3*t^5 + 16*t^4 + 4*t^3 + 13*t^2 + 6*t + 5)//(t^2 + 12*t + 7), (4*t^8 + 19*t^7 + 14*t^6 + 18*t^5 + 27*t^4 + 13*t^3 + 9*t^2 + 14*t + 12)//(t^3 + 18*t^2 + 21*t + 13) ),
             ((17*t^6 + 10*t^5 + 24*t^4 + 15*t^3 + 22*t^2 + 27*t + 5)//(t^2 + 16*t + 6), (20*t^8 + 24*t^7 + 22*t^6 + 12*t^5 + 21*t^4 + 21*t^3 + 9*t^2 + 21*t + 12)//(t^3 + 24*t^2 + 18*t + 19) ),
             ((17*t^8 + 21*t^7 + 20*t^5 + 24*t^4 + 21*t^3 + 4*t^2 + 9*t + 13)//(t^4 + 17*t^3 + 12*t^2 + 28*t + 28), (23*t^11 + 25*t^10 + 8*t^9 + 7*t^8 + 28*t^7 + 16*t^6 + 7*t^5 + 23*t^4 + 9*t^3 + 27*t^2 + 13*t + 13)//(t^6 + 11*t^5 + 14*t^4 + 13*t^3 + 6*t^2 + 18*t + 12) )];

julia> E = elliptic_curve(ktfield,[3*t^8+24*t^7+22*t^6+15*t^5+28*t^4+20*t^3+16*t^2+26*t+16, 24*t^12+27*t^11+28*t^10+8*t^9+6*t^8+16*t^7+2*t^6+10*t^5+3*t^4+22*t^3+27*t^2+10*t+3]);

julia> bk = [E(collect(i)) for i in bk];

julia> Oscar._prop217(E,bk[2],2)
5-element Vector{Tuple{FqPolyRingElem, FqPolyRingElem}}:
 (t^2 + 12*t + 7, 0)
 (t^3 + 8*t + 3, 0)
 (t^4 + 23*t + 2, 0)
 (25*t + 22, 1)
 (12*t + 28, t)

julia> Oscar._prop217(E,bk[1],1)
2-element Vector{Tuple{FqPolyRingElem, FqPolyRingElem}}:
 (1, 0)
 (t, 0)
```
"""
function _prop217(E::EllipticCurve, P::EllipticCurvePoint, k)
  @req !iszero(P[3]) "P must not be torsion" # seems like we cannot check this
  xn = numerator(P[1])
  xd = denominator(P[1])
  yn = numerator(P[2])
  yd = denominator(P[2])
  OP = divexact(max(degree(xd), degree(xn) - 4), 2)
  dega = k + 2*OP
  degb = k + 2*OP - 2 - divexact(degree(xd), 2) #?
  base = base_field(E)
  Bt = base_ring(base)
  B = coefficient_ring(Bt)

  R,ab = polynomial_ring(base,vcat([Symbol(:a,i) for i in 0:dega],[Symbol(:b,i) for i in 0:degb]),cached=false)
  Rt, t1 = polynomial_ring(R,:t)
  a = reduce(+,(ab[i+1]*t1^i for i in 0:dega), init=zero(Rt))
  b = reduce(+,(ab[2+dega+j]*t1^j for j in 0:degb), init=zero(Rt))
  c = a*xn(t1) - b*yn(t1)
  r = mod(c, xd(t1))
  # setup the linear equations for coefficients of r to vanish
  # and for the degree of c to be bounded above by
  # k + 2*OP + 4 + degree(xd)
  eq1 = collect(coefficients(r))
  eq2 = [coeff(c,i) for i in (k + 2*OP + 4 + degree(xd) + 1):degree(c)]
  eqns = vcat(eq1, eq2)

  # collect the equations as a matrix
  cc = [[coeff(j, abi) for abi in ab] for j in eqns]
  M = matrix(B, length(eqns), length(ab), reduce(vcat,cc, init=elem_type(base)[]))
  # @assert M == matrix(base, cc) # does not work if length(eqns)==0
  K = kernel(M; side = :right)
  kerdim = ncols(K)
  result = Tuple{elem_type(Bt),elem_type(Bt)}[]
  t = gen(Bt)
  for j in 1:kerdim
    aa = reduce(+, (K[i+1,j]*t^i for i in 0:dega), init=zero(Bt))
    bb = reduce(+, (K[dega+i+2,j]*t^i for i in 0:degb), init=zero(Bt))
    push!(result, (aa, bb))
  end
  # confirm the computation
  @assert kerdim == 2*k + OP # prediced by Riemann-Roch
  for (a,b) in result
    @assert mod(a*xn - b*yn, xd) == 0
    @assert degree(a) <= k + 2*OP
    @assert degree(b) <= k + 2*OP - 2 - 1//2*degree(xd)
    @assert degree(a*xn - b*yn) <= k + 2*OP + 4 + degree(xd)
  end
  return result
end

function iszero(P::EllipticCurvePoint)
  return iszero(P[1]) && isone(P[2]) && iszero(P[3])
end

@doc raw"""
    linear_system(X::EllipticSurface, P::EllipticCurvePoint, k::Int64) -> LinearSystem

Compute the linear system ``|O + P + k F|`` on the elliptic surface ``X``.
Here ``F`` is the class of the fiber over ``[0:1]``, ``O`` the zero section
and ``P`` any section given as a point on the generic fiber.

The linear system is represented in terms of the Weierstrass coordinates.
"""
function linear_system(X::EllipticSurface, P::EllipticCurvePoint, k::Int64)
  euler_characteristic(X) == 2 || error("linear system implemented only for elliptic K3s")
  #FS = function_field(weierstrass_model(X)[1])
  FS = function_field(X)
  U = weierstrass_chart_on_minimal_model(X)
  (x,y,t) = ambient_coordinates(U)

  sections = elem_type(FS)[]
  if iszero(P[3])
    append!(sections, [FS(t)^(i-k) for i in 0:k])
    append!(sections, [FS(t)^(i-k)*FS(x) for i in 0:k-4])
  else
    xn = numerator(P[1])
    xd = denominator(P[1])
    yn = numerator(P[2])
    yd = denominator(P[2])

    I = saturated_ideal(defining_ideal(U))
    IP = ideal([x*xd(t)-xn(t),y*yd(t)-yn(t)])
    issubset(I, IP) || error("P does not define a point on the Weierstrasschart")

    @assert gcd(xn, xd)==1
    @assert gcd(yn, yd)==1
    ab = _prop217(generic_fiber(X), P, k)
    d = divexact(yd, xd)(t)
    den = t^k*(x*xd(t) - xn(t))
    for (a,b) in ab
      c = divexact(b*yn - a*xn, xd)
      num = a(t)*x+b(t)*d*y + c(t)
      push!(sections, FS(num//den))
    end
  end
  return sections
end

@doc raw"""
    two_neighbor_step(X::EllipticSurface, F1::Vector{QQFieldElem})

Given an isotropic nef divisor ``F1`` with ``F1.F = 2``,
compute the linear system ``|F1|`` and return the corresponding generic fiber
as a double cover `C` of the projective line branched over four points.

Input:
``F1`` is represented as a vector in the `algebraic_lattice(X)`

Output:
A tuple `(C, (x1, y1, t1))` defined as follows.
- `C` is given by a polynomial `y1^2 - q(x1)` in `k(t)[x1,y1]` with `q` of degree 3 or 4.
- (x1,y1,t1) are expressed as rational functions in terms of the weierstrass coordinates `(x,y,t)`.
"""
function two_neighbor_step(X::EllipticSurface, F::Vector{QQFieldElem})
  E = generic_fiber(X)
  basisNS, tors, NS = algebraic_lattice(X)
  V = ambient_space(NS)
  @req inner_product(V, F, F)==0 "not an isotropic divisor"
  @req euler_characteristic(X) == 2 "not a K3 surface"
  F0 = zeros(QQ,degree(NS)); F0[1]=1

  @req inner_product(V, F, F0) == 2 "not a 2-neighbor"

  D1, D, P, l, c = horizontal_decomposition(X, F)
  u = _elliptic_parameter(X, D1, D, P, l, c)
  @assert scheme(parent(u)) === X
  pr = weierstrass_contraction(X)
  WX, _ = weierstrass_model(X)
  # The following is a cheating version of the command u = pushforward(pr)(u)
  u = function_field(WX)(u[weierstrass_chart_on_minimal_model(X)])
  @assert scheme(parent(u)) === weierstrass_model(X)[1]

  # Helper function
  my_const(u::MPolyRingElem) = is_zero(u) ? zero(coefficient_ring(parent(u))) : first(coefficients(u))

  # transform to a quartic y'^2 = q(x)
  if iszero(P[3])  #  P = O
    eqn1, phi1 = _elliptic_parameter_conversion(X, u, case=:case1)
    eqn2, phi2 = _normalize_hyperelliptic_curve(eqn1)
#   function phi_func(x)
#     y = phi1(x)
#     n = numerator(y)
#     d = denominator(y)
#     return phi2(n)//phi2(d)
#   end
#   phi = MapFromFunc(domain(phi1), codomain(phi2), phi_func)
#   # TODO: Verify that the construction below also works and replace by that, eventually.
#   phi_alt = compose(phi1, extend_domain_to_fraction_field(phi2))
#   @assert phi.(gens(domain(phi))) == phi_alt.(gens(domain(phi)))
   phi = compose(phi1, extend_domain_to_fraction_field(phi2))
  elseif iszero(2*P) # P is a 2-torsion section
    eqn1, phi1 = _elliptic_parameter_conversion(X, u, case=:case3)
    #eqn1, phi1 = _conversion_case_3(X, u)
    (x2, y2) = gens(parent(eqn1))

    # Make sure the coefficient of y² is one (or a square) so that 
    # completing the square works. 
    c = my_const(coeff(eqn1, [x2, y2], [0, 2]))::AbstractAlgebra.Generic.FracFieldElem
    eqn1 = inv(unit(factor(c)))*eqn1

    eqn2, phi2 = _normalize_hyperelliptic_curve(eqn1)
    phi = compose(phi1, extend_domain_to_fraction_field(phi2))
  else  # P has infinite order
    eqn1, phi1 = _elliptic_parameter_conversion(X, u, case=:case2)
    #eqn1, phi1 = _conversion_case_2(X, u)
    (x2, y2) = gens(parent(eqn1))
    
    # Make sure the coefficient of y² is one (or a square) so that 
    # completing the square works. 
    c = my_const(coeff(eqn1, [x2, y2], [0, 2]))::AbstractAlgebra.Generic.FracFieldElem
    eqn1 = inv(unit(factor(c)))*eqn1

    eqn2, phi2 = _normalize_hyperelliptic_curve(eqn1)
    phi = compose(phi1, extend_domain_to_fraction_field(phi2))
  end

  return eqn2, phi
end

@doc raw"""
    horizontal_decomposition(X::EllipticSurface, L::Vector{QQFieldElem}) -> WeilDivisor, EllipticCurvePoint

Given a divisor ``L`` as a vector in the `algebraic_lattice(X)`
find a linearly equivalent divisor ``(n-1) O + P + V = D ~ L`` where
``O`` is the zero section, ``P`` is any section and ``V`` is vertical.

Returns a tuple `(D1, D, P, l, c)` where `D` and `P` are as above and
``D <= D1 = (n-1)O + P + n_1F_1 + ... n_k F_k`` with ``l = n_1 + ... n_k`` minimal
and the `F_i` are some other fibers.
The rational function `c=c(t)` has divisor of zeros and poles``
(c) = lF - n_0F_1 + ... n_k F_k``
"""
function horizontal_decomposition(X::EllipticSurface, F::Vector{QQFieldElem})
  E = generic_fiber(X)
  basisNS, tors, NS = algebraic_lattice(X)
  V = ambient_space(NS)
  @req inner_product(V, F, F)==0 "not an isotropic divisor"
  @req euler_characteristic(X) == 2 "not a K3 surface"
  # how to give an ample divisor automagically in general?
  # @req is_nef(X, F) "F is not nef"
  l = F[1]
  rk_triv = nrows(trivial_lattice(X)[2])
  n = rank(NS)
  @assert degree(NS) == rank(NS)
  P0 = sum([ZZ(F[i])*X.MWL[i-rk_triv] for i in (rk_triv+1):n], init = E([0,1,0]))
  P0_div = section(X, P0)
  @vprint :EllipticSurface 2 "Computing basis representation of $(P0)\n"
  p0 = basis_representation(X, P0_div) # this could be done from theory alone
  F1 = F - p0  # should be contained in the QQ-trivial-lattice
  F2 = F1
  if all(isone(denominator(i)) for i in F1)
    # no torsion
    P = P0
    P_div = P0_div
  else
    found = false
    for (i,(T, tor)) in enumerate(tors)
      d = F2 - _vec(tor)
      if all(isone(denominator(i)) for i in d)
        found = true
        T0 = mordell_weil_torsion(X)[i]
        P = P0 + T0
        break
      end
    end
    @assert found
    P_div = section(X, P)
    p = basis_representation(X, P_div)
    F2 = F1 - p
    @assert all( isone(denominator(i)) for i in F2)
  end
  D = P_div
  D = D + ZZ(F2[2])*zero_section(X)
  D1 = D
  F2 = ZZ.(F2); F2[2] = 0
  l = F2[1] # number of fibers that we need
  # find the fiber components meeting O necessary
  F3 = F2
  (_,_,t) = ambient_coordinates(weierstrass_chart_on_minimal_model(X))
  c = t^0
  for (pt, rt, fiber, comp, gram) in reducible_fibers(X)
    Fib0 = comp[1]
    f0 = zeros(QQFieldElem, length(basisNS))
    for i in 1:length(basisNS)
      if !isone(components(Fib0)[1]+components(basisNS[i])[1])
        if length(comp)==2 && 2<i<=rk_triv
          f0[i] = 2
        else
          f0[i] = 1
        end
      end
    end
    f0 = f0 * inv(gram_matrix(ambient_space(NS)))
    @assert inner_product(ambient_space(NS), f0,f0) == -2
    nonzero = [i for i in 3:rk_triv if f0[i]!=0]
    if pt[2]==0 # at infinity
      t0 = t
    else
      t0 = t//(t*pt[2]-pt[1])
    end
    while any(F3[i]<0 for i in nonzero)
      F3 = F3 - f0
      D = D + Fib0
      D1 = D1 + fiber
      c = c*t0
      @assert F3[1]>=0
    end
  end
  pt, _ = fiber(X)
  if pt[2]==0 # at infinity
    t0 = t
  else
    t0 = t//(t*pt[2]-pt[1])
  end
  c = c*t0^ZZ(F3[1])
  D = D + F3[1]*basisNS[1]
  D1 = D1 + F3[1]*basisNS[1]
  F4 = copy(F3); F4[1]=0
  @assert all(F4[i]>=0 for i in 1:length(basisNS))
  D = D + sum(ZZ(F4[i])*basisNS[i] for i in 1:length(basisNS))
  @assert D<=D1
  l = Int(l)
  return D1, D, P, l, c
end

@doc raw"""
  elliptic_parameter(X::EllipticSurface, F::Vector{QQFieldElem}) -> LinearSystem

Return the elliptic parameter ``u`` of the divisor class `F`. 

The input `F` must be given with respect to the basis of
`algebraic_lattice(X)` and be an isotropic nef divisor. 
This method assumes that $X$ is a K3 surface.
"""
function elliptic_parameter(X::EllipticSurface, F::Vector{QQFieldElem})
  D1, D, P, l, c = horizontal_decomposition(X, F)
  return _elliptic_parameter(X, D1, D, P, l, c)
end

@doc raw"""
    _elliptic_parameter(X::EllipticSurface, D::WeilDivisor, l, c)

Compute the linear system of ``D = (n-1) O + P + V``.
where V is vertical and `l` is the coefficient of the fiber class.
Assumes `D` nef and `D^2=0`.
Typically ``D`` is the output of `horizontal_decomposition`.
"""
function _elliptic_parameter(X::EllipticSurface, D1::WeilDivisor, D::WeilDivisor, P::EllipticCurvePoint, l::Int, c)
  S, piS = weierstrass_model(X);
  piX = weierstrass_contraction(X)
  c = function_field(X)(c)
  L = [i*c for i in linear_system(X, P, l)];
  LonX = linear_system(L, D1, check=false);

  LsubF, Tmat = subsystem(LonX, D);

  LsubFonS = [sum(Tmat[i,j]*L[j] for j in 1:ncols(Tmat)) for i in 1:nrows(Tmat)]

  @assert length(LsubFonS)==2

  u2 = LsubFonS[2]//LsubFonS[1]
  return u2
end


@doc raw"""
  extended_ade(ADE::Symbol, n::Int)

Return the dual intersection matrix of an extended ade Dynkin diagram
as well as the isotropic vector (with positive coefficients in the roots).
"""
function extended_ade(ADE::Symbol, n::Int)
  R = change_base_ring(ZZ,gram_matrix(root_lattice(ADE,n)))
  G = block_diagonal_matrix([ZZ[2;],R])
  if ADE == :E && n == 8
    G[1,n] = -1
    G[n,1] = -1
  end
  if ADE == :E && n == 7
    G[1,2] = -1
    G[2,1] = -1
  end
  if ADE == :E && n == 6
    G[1,n+1] = -1
    G[n+1,1] = -1
  end
  if ADE == :A && n > 0
    G[1,2] = -1
    G[2,1] = -1
    G[1,n+1] = -1
    G[n+1,1] = -1
  end
  if ADE == :A && n ==1 0
    G[1,2]= -2
    G[2,1] = -2
  end
  if ADE == :D
    G[1,n] = -1
    G[n,1] = -1
  end
  @assert rank(G) == n
  return -G, kernel(G; side = :left)
end

function basis_representation(X::EllipticSurface, D::WeilDivisor)
  basis_ambient,_, NS = algebraic_lattice(X)
  G = gram_matrix(ambient_space(NS))
  n = length(basis_ambient)
  v = zeros(ZZRingElem, n)
  @vprint :EllipticSurface 3 "computing basis representation of $D\n"
  for i in 1:n
    v[i] = intersect(basis_ambient[i], D)
  end
  @vprint :EllipticSurface 3 "done computing basis representation\n"
  return v*inv(G)
end

################################################################################
#
# patches for Oscar
#
################################################################################


########################################################################
# Internal functionality for Weierstrass transformation 
########################################################################


@doc raw"""
    _normalize_hyperelliptic_curve(g::MPolyRingElem, parent=nothing)

Transform ``a(x)y^2 + b(x)y - h(x)`` in ``K(t)[x,y]`` to ``y'^2 - h(x')``
"""
function _normalize_hyperelliptic_curve(g::MPolyRingElem; parent::Union{MPolyRing, Nothing}=parent(g))
  R = Oscar.parent(g)
  @assert ngens(R) == 2 "polynomial must be bivariate"
  F = fraction_field(R)
  kt = coefficient_ring(R)
  (x, y) = gens(R)

  # Prepare the output ring
  if parent===nothing
    R1, (x1, y1) = R, gens(R)
  else
    R1 = parent
    @assert coefficient_ring(R1) == coefficient_ring(R) "coefficient ring of output is incompatible with input"
    (x1, y1) = gens(R1)
  end

  # Get the coefficients of g as a univariate polynomial in y
  ktx, X = polynomial_ring(kt, :X, cached=false)
  ktxy, Y = polynomial_ring(ktx, :y, cached=false)

  # Maps to transform to univariate polynomials in y
  split_map_R = hom(R, ktxy, [ktxy(X), Y])
  split_map_R1 = hom(R1, ktxy, [ktxy(X), Y])
  G = split_map_R(g)
  @assert degree(G) == 2 "polynomial must be of degree 2 in its second variable"

  #complete the square
  h, b, a = collect(coefficients(G))
  h = -h
  u = unit(factor(a))
  a = inv(u)*a
  b = inv(u)*b
  success, sqa = is_square_with_sqrt(a)
  @assert success "leading coefficient as univariate polynomial in the second variable must be a square"

  F1 = fraction_field(R1)
  psi = hom(R1, F, F.([x, (2*evaluate(a, x)*y + evaluate(b, x))//(2*evaluate(sqa, x))]))
  conv = MapFromFunc(ktx, R1, f->evaluate(f, x1))
  (a1, b1, sqa1) = conv.([a, b, sqa])
  phi = hom(R, F1, F1.([x1, (2*sqa1*y1-b1)//(2*a1)]))
  phiF = MapFromFunc(F, F1, x-> phi(numerator(x))//phi(denominator(x)))
  # the inverse map if wanted
  # psiF = MapFromFunc(F1, F, x-> psi(numerator(x))//psi(denominator(x)))
  # @assert all(phiF(psiF(F1(i)))==i for i in gens(R1))

  # absorb squares into y1
  g1 = numerator(phi(g))
  G1 = split_map_R1(g1)
  ff = factor(first(coefficients(G1)))
  c = prod([p^div(i, 2) for (p, i) in ff], init=one(ktx))
  #d = sqrt(my_coeff(g1, y1, 2))
  d = last(coefficients(split_map_R1(g1)))
  success, d = is_square_with_sqrt(d)
  @assert success "leading coefficient must be a square"

  phi1 = hom(R1, F1, [F1(x1), F1(evaluate(c, x1), evaluate(d, x1))*y1])
  phiF1 = MapFromFunc(F1, F1, x-> phi1(numerator(x))//phi1(denominator(x)))
  phi2 = compose(phi, phiF1)
  g2 = numerator(phi1(g1))
  #c = my_coeff(g2, y1, 2)
  c = last(coefficients(split_map_R1(g2)))
  g2 = divexact(g2, evaluate(c, x1))
  @assert degree(g2, gen(parent, 1)) <= 4 "degree in the first variable is too high"
  @assert degree(g2, gen(parent, 1)) >= 3 "degree in the first variable is too low"
  return g2, phi2
end


@doc raw"""
    elliptic_surface(g::MPolyRingElem, x::MPolyRingElem, y::MPolyRingElem, P::Vector{<:RingElem})

Transform a bivariate polynomial `g` of the form `y^2 - Q(x)` with `Q(x)` of
degree at most ``4`` to Weierstrass form and return the corresponding
elliptic surface as well as the coordinate transformation.
"""
function elliptic_surface(g::MPolyRingElem, P::Vector{<:RingElem})
  R = parent(g)
  (x, y) = gens(R)
  P = base_ring(R).(P)
  g2, phi2 = transform_to_weierstrass(g, x, y, P);
  E2 = elliptic_curve(g2, x, y);
  t2 = gen(base_field(E2));
  Y2 = elliptic_surface(E2, 2)

  W = weierstrass_chart(Y2)
  RW = ambient_coordinate_ring(W); FRW = fraction_field(RW);
  (x, y, t) = gens(RW) # Output of phi2 lives in (k(t))(x, y), so we need to convert
  psi = hom(base_ring(codomain(phi2)), FRW, f->evaluate(f, t), FRW.([x,y])); # a helper function
  psi_ext = extend_domain_to_fraction_field(psi)
  return Y2, phi2 * psi_ext
end

@doc raw"""
    transform_to_weierstrass(g::MPolyRingElem, x::MPolyRingElem, y::MPolyRingElem, P::Vector{<:RingElem})

Transform a bivariate polynomial `g` of the form `y^2 - Q(x)` with `Q(x)` of degree ``≤ 4``
to Weierstrass form. This returns a pair `(f, trans)` where `trans` is an endomorphism of the 
`fraction_field` of `parent(g)` and `f` is the transform. The input `P` must be a rational point 
on the curve defined by `g`, i.e. `g(P) == 0`.
"""
function transform_to_weierstrass(g::MPolyRingElem, x::MPolyRingElem, y::MPolyRingElem, P::Vector{<:RingElem})
  R = parent(g)
  F = fraction_field(R)

  # In case of variables in the wrong order, switch and transform the result.
  if x == R[2] && y == R[1]
    switch = hom(R, R, reverse(gens(R)))
    g_trans, trans = transform_to_weierstrass(switch(g), y, x, reverse(P))
    new_trans = MapFromFunc(F, F, f->begin
                                switch_num = switch(numerator(f))
                                switch_den = switch(denominator(f))
                                interm_res = trans(F(switch_num))//trans(F(switch(den)))
                                num = numerator(interm_res)
                                den = denominator(interm_res)
                                switch(num)//switch(den)
                            end
                           )
    return switch(g_trans), new_trans
  end
  @assert ngens(R) == 2 "input polynomial must be bivariate"
  @assert x in gens(R) "second argument must be a variable of the parent of the first"
  @assert y in gens(R) "third argument must be a variable of the parent of the first"

  kk = coefficient_ring(R)
  kkx, X = polynomial_ring(kk, :x, cached=false)
  kkxy, Y = polynomial_ring(kkx, :y, cached=false)

  imgs = [kkxy(X), Y]
  split_map = hom(R, kkxy, imgs)

  G = split_map(g)
  @assert degree(G) == 2 "input polynomial must be of degree 2 in y"
  @assert all(h->degree(h)<=4, coefficients(G)) "input polynomial must be of degree <= 4 in x"
  @assert iszero(coefficients(G)[1]) "coefficient of linear term in y must be zero"
  @assert isone(coefficients(G)[2]) "leading coefficient in y must be one"

  length(P) == 2 || error("need precisely two point coordinates")
  (px, py) = P
  #    assert g.subs({x:px,y:py})==0
  @assert iszero(evaluate(g, P)) "point does not lie on the hypersurface"
  gx = -evaluate(g, [X + px, zero(X)])
  coeff_gx = collect(coefficients(gx))
  A = coeff(gx, 4)
  B = coeff(gx, 3)
  C = coeff(gx, 2)
  D = coeff(gx, 1)
  E = coeff(gx, 0)
  #E, D, C, B, A = coeff_gx
  if !iszero(E)
    b = py
    a4, a3, a2, a1, a0 = A,B,C,D,E
    A = b
    B = a1//(2*b)
    C = (4*a2*b^2-a1^2)//(8*b^3)
    D = -2*b

    x1 = x//y
    y1 = (A*y^2+B*x*y+C*x^2+D*x^3)//y^2
    x1 = x1+px

    # TODO: The following are needed for the inverse. To be added eventually.
    # x2 = (y-(A+B*x+C*x^2))//(D*x^2)
    # y2 = x2//x
    # x2 = evaluate(x2, [x-px, y])
    # y2 = evaluate(y2, [x-px, y])

    # @assert x == evaluate(x1, [x2, y2])
    # @assert y == evaluate(y1, [x2, y2])
  else
    # TODO compute the inverse transformation (x2,y2)
    x1 = 1//x
    y1 = y//x^2
    g1 = numerator(evaluate(g, [x1, y1]))
    c = coeff(g1, [x], [3])
    x1 = evaluate(x1, [-x//c, y//c])
    y1 = evaluate(y1, [-x//c, y//c])
    x1 = x1+px
    #@assert x == evaluate(x1, [x2, y2])
    #@assert y == evaluate(y1, [x2, y2])
  end
  @assert F === parent(x1) "something is wrong with caching of fraction fields"
  # TODO: eventually add the inverse.
  trans = MapFromFunc(F, F, f->evaluate(numerator(f), [x1, y1])//evaluate(denominator(f), [x1, y1]))
  f_trans = trans(F(g))
  fac = [a[1] for a in factor(numerator(f_trans)) if isone(a[2]) && _is_in_weierstrass_form(a[1])]
  isone(length(fac)) || error("transform to weierstrass form did not succeed")

  # normalize the output
  result = first(fac)
  result = inv(first(coefficients(coeff(result, gens(parent(result)), [3, 0]))))*result

  return result, trans
end

function _is_in_weierstrass_form(f::MPolyRingElem)
  R = parent(f)
  @req ngens(R) == 2 "polynomial must be bivariate"
  # Helper function
  my_const(u::MPolyRingElem) = is_zero(u) ? zero(coefficient_ring(parent(u))) : first(coefficients(u))

  (x, y) = gens(R)
  f = -inv(my_const(coeff(f, [x, y], [0, 2]))) * f
  isone(-coeff(f, [x, y], [0, 2])) || return false
  isone(coeff(f, [x, y], [3, 0])) || return false
  
  a6 = coeff(f, [x,y], [0,0])
  a4 = coeff(f, [x,y], [1,0])
  a2 = coeff(f, [x,y], [2,0])
  a3 = -coeff(f, [x,y], [0,1])
  a1 = -coeff(f, [x,y], [1,1])
  a_invars = [my_const(i) for i in [a1,a2,a3,a4,a6]]
  (a1,a2,a3,a4,a6) = a_invars
  return f == (-(y^2 + a1*x*y + a3*y) + (x^3 + a2*x^2 + a4*x + a6))
end

function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:MPolyRingElem}, a::Vector{T}) where {T<:RingElem}
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end

function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:PolyRingElem}, a::RingElem)
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end

function extend_domain_to_fraction_field(phi::Map{<:MPolyRing, <:Ring})
  ext_dom = fraction_field(domain(phi))
  return MapFromFunc(ext_dom, codomain(phi), x->phi(numerator(x))*inv(phi(denominator(x))))
end

########################################################################
# The three conversions from Section 39.1 in 
#   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface" 
# pp. 44--45.
########################################################################

function _elliptic_parameter_conversion(X::EllipticSurface, u::VarietyFunctionFieldElem; 
    case::Symbol=:case1, names=[:x, :y, :t]
  )
  @req variety(parent(u)) === weierstrass_model(X)[1] "function field element must live on the weierstrass model of the first argument"
  @req length(names) == 3 "need 3 variable names x, y, t"
  U = weierstrass_chart(X)
  R = ambient_coordinate_ring(U)
  x, y, t = gens(R)
  loc_eqn = first(gens(modulus(OO(U))))
  E = generic_fiber(X)::EllipticCurve
  f = equation(E)
  kk = base_ring(X)
  kkt_frac_XY = parent(f)::MPolyRing
  (xx, yy) = gens(kkt_frac_XY)
  kkt_frac = coefficient_ring(kkt_frac_XY)::AbstractAlgebra.Generic.FracField
  kkt = base_ring(kkt_frac)::PolyRing
  T = first(gens(kkt))

# kk = base_ring(U)
# kkt, T = polynomial_ring(kk, :T, cached=false)
# kkt_frac = fraction_field(kkt)
# kkt_frac_XY, (xx, yy) = polynomial_ring(kkt_frac, [:X, :Y], cached=false)
  R_to_kkt_frac_XY = hom(R, kkt_frac_XY, [xx, yy, kkt_frac_XY(T)])

  f_loc = first(gens(modulus(OO(U))))
  @assert f == R_to_kkt_frac_XY(f_loc) && _is_in_weierstrass_form(f) "local equation is not in Weierstrass form"
  a = a_invariants(E)

  u_loc = u[U]::AbstractAlgebra.Generic.FracFieldElem # the representative on the Weierstrass chart

  # Set up the ambient_coordinate_ring of the new Weierstrass-chart
  kkt2, t2 = polynomial_ring(kk, names[3], cached=false)
  kkt2_frac = fraction_field(kkt2)
  S, (x2, y2) = polynomial_ring(kkt2_frac, names[1:2], cached=false)
  FS = fraction_field(S)

  # Helper function
  my_const(u::MPolyRingElem) = is_zero(u) ? zero(coefficient_ring(parent(u))) : first(coefficients(u))

  # We verify the assumptions made on p. 44 of
  #   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface"
  # for the first case considered there.
  @assert all(x->isone(denominator(x)), a) "local equation does not have the correct form"
  a = numerator.(a)
  @assert iszero(a[1]) "local equation does not have the correct form"
  @assert degree(a[2]) <= 4 "local equation does not have the correct form"
  @assert iszero(a[3]) "local equation does not have the correct form"
  @assert degree(a[4]) <= 8 "local equation does not have the correct form"
  @assert degree(a[5]) <= 12 "local equation does not have the correct form" # This is really a₆ in the notation of the paper, a₅ does not exist.
  # reduce fraction
  u_frac = R_to_kkt_frac_XY(numerator(u_loc))//R_to_kkt_frac_XY(denominator(u_loc))
  u_num = numerator(u_frac)
  u_den = denominator(u_frac)
  if case == :case1
    # D = 2O
    u_poly = u_num*inv(u_den) # Will throw if the latter is not a unit
    # Extract a(t) and b(t) as in the notation of the paper
    a_t = my_const(coeff(u_poly, [xx, yy], [0, 0]))
    b_t = my_const(coeff(u_poly, [xx, yy], [1, 0]))

    a_t = evaluate(a_t, x2)
    b_t = evaluate(b_t, x2)
    phi = hom(R, FS, FS.([(t2 - a_t)//b_t, y2, x2]))
    f_trans = phi(f_loc)
    return numerator(f_trans), phi
  elseif case == :case2
    # D = O + P
    @assert degree(u_num, 2) == 1 && degree(u_num, 1) <= 1 "numerator does not have the correct degree"
    @assert degree(u_den, 1) == 1 && degree(u_den, 2) == 0 "denominator does not have the correct degree"

    # We expect a form as on p. 44, l. -4
    denom_unit = my_const(coeff(u_den, [xx, yy], [1, 0]))
    x0 = -inv(denom_unit)*my_const(coeff(u_den, [xx, yy], [0, 0]))
    b_t = inv(denom_unit)*my_const(coeff(u_num, [xx, yy], [0, 1]))
    u_num = u_num - denom_unit * b_t * yy
    a_t = inv(denom_unit)*my_const(coeff(u_num, [xx, yy], [1, 0]))
    u_num = u_num - denom_unit * a_t * (xx - x0)
    @assert is_constant(u_num) "numerator is not in the correct form"
    y0 = my_const(coeff(u_num, [xx, yy], [0, 0])) * inv(denom_unit * b_t)

    @assert a_t + b_t*(yy + y0)//(xx - x0) == u_frac "decomposition failed"
    # We have 
    #
    #   y ↦ (u - a_t) * (x - x₀) / b_t - y₀ = (t₂ - a_t(x₂)) * (y₂ - x₀(x₂)) / b_t(x₂) - y₀(x₂)
    #   x ↦ y₂
    #   t ↦ x₂
    phi = hom(R, FS, FS.([y2, (t2 - evaluate(a_t, x2)) * (y2 - evaluate(x0, x2)) // evaluate(b_t, x2) - evaluate(y0, x2), x2]))
    f_trans = phi(f_loc)
    eqn1 = numerator(f_trans)
    # According to 
    #   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface" 
    # p. 45, l. 1 we expect the following cancelation to be possible:
    divisor_num = evaluate(numerator(x0), x2)
    divisor_den = evaluate(denominator(x0), x2)
    divisor = divisor_den * y2 - divisor_num
    success, eqn1 = divides(eqn1, divisor) # This division must only be possible in the ring K(x2)[y2].
                                           # Hence, multiplying by the denominator `divisor_den` is 
                                           # merely an educated guess.
    @assert success "division failed"
    return eqn1, phi
  elseif case == :case3
    # D = O + T

    @assert u_den == xx "elliptic parameter was not brought to the correct form"
    @assert degree(u_num, 1) <= 1 && degree(u_num, 2) <= 1 "numerator does not have the correct degrees"
    a_t = my_const(coeff(u_num, [xx, yy], [1, 0]))
    b_t = my_const(coeff(u_num, [xx, yy], [0, 1]))

    # New Weierstrass equation is of the form 
    #
    #   x^2 = h(t, u)
    #
    # so y₂ = x, x₂ = t, and t₂ = u.
    #
    # We have u = a_t + b_t * y/x ⇒ y = (u - a_t) * x / b_t = (t₂ - a_t(x₂)) * y₂ / b_t(x₂)
    phi = hom(R, FS, FS.([y2, (t2 - evaluate(a_t, x2)) * y2 // evaluate(b_t, x2), x2]))
    f_trans = phi(f_loc)
    eqn1 = numerator(f_trans)
    # According to 
    #   A. Kumar: "Elliptic Fibrations on a generic Jacobian Kummer surface" 
    # p. 45, l. 15 we expect the following cancelation to be possible:
    success, eqn1 = divides(eqn1, y2)
    @assert success "equation did not come out in the anticipated form"
    return eqn1, phi
  else 
    error("case not recognized")
  end
end
