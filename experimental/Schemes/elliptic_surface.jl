export elliptic_surface, trivial_lattice, weierstrass_model, weierstrass_chart, algebraic_lattice, zero_section, section, relatively_minimal_model, fiber_components

@attributes mutable struct EllipticSurface{BaseField<:Field, BaseCurveFieldType} <: AbsCoveredScheme{BaseField}
  Y::CoveredScheme{BaseField}
  E::EllCrv{BaseCurveFieldType}
  MWL::Vector{EllCrvPt{BaseCurveFieldType}}
  MWLtors::Vector{EllCrvPt{BaseCurveFieldType}}
  Weierstrasschart::AbsSpec
  Weierstrassmodel::CoveredScheme
  inc_Weierstrass
  inc_Y
  bundle_number::Int
  blowup
  blowups
  exceptionals
  ambient_blowups
  ambient_exceptionals

  function EllipticSurface(generic_fiber::EllCrv{F}, bundle_number::Int) where F
    B = typeof(coefficient_ring(base_ring(base_field(generic_fiber))))
    S = new{B,F}()
    S.E = generic_fiber
    S.bundle_number = bundle_number
    return S
  end

end

elliptic_surface(generic_fiber::EllCrv, s::Int) = EllipticSurface(generic_fiber, s)

function elliptic_surface(generic_fiber::EllCrv, s::Int, mwl_gens::Vector{<:EllCrvPt})
  @req all(parent(i)==S.E for i in mwl_gens) "not a vector of points on $(generic_fiber)"
  S = elliptic_surface(generic_fiber, s)
  S.MWL = mwl_gens
end

function underlying_scheme(S::EllipticSurface)
  if isdefined(S,:X)
    return S.Y
  end
  return relatively_minimal_model(generic_fiber(S), S.bundle_number)
end

generic_fiber(S::EllipticSurface) = S.E
weierstrass_chart(S::EllipticSurface) = S.Weierstrasschart

@attr function algebraic_lattice(S)
  isdefined(S, :MWL) || error("no generators for the Mordell-Weil group available")
  return algebraic_lattice(S, S.MWL)
end

function algebraic_lattice(S::EllipticSurface, mwl_gens::Vector{<:EllCrvPt})
  basis, _, G = trivial_lattice(S)
  l = length(mwl_gens)
  r = length(basis)
  sections = [section(S, i) for i in mwl_gens]
  n = l+r
  GA = zero_matrix(ZZ, n, n)
  GA[1:r,1:r] = G
  GA[r+1:n,r+1:n] = -2*identity_matrix(ZZ, l)
  gensA = vcat(basis, sections)
  @vprint :ellipticK3 2 "computing intersection numbers"
  for i in 1:n
    @vprint :ellipticK3 2 "\nrow $(i): \n"
    for j in max(i + 1, r + 1):n
      @vprint :ellipticK3 2 "$(j) "
      GA[i,j] = intersect(gensA[i],gensA[j])
      GA[j,i] = GA[i,j]
    end
  end
  @assert rank(GA) == n "todo: treat torsion sections" # need to adapt mordell_weil then too
  return gensA, integer_lattice(GA)
end

@attr ZZLat function mordell_weil_lattice(S::EllipticSurface)
  NS = algebraic_lattice(S)
  t = length(trivial_lattice(S)[1])
  trivNS = basis_matrix(NS)[1,:t]
  V = ambient_space(NS)
  P = orthogonal_complement(V, trivNS)
  mwl = rescale(lattice(V, basis_matrix(NS)*P, is_basis=false),-1)
  return lll(mwl)
end

@attr ZZLat function mordell_weil_torsion(S::EllipticSurface)
  error("not implemented")
end

function Base.show(io::IOContext, S::EllipticSurface)
  io = pretty(io)
  if get(io, :supercompact, false)
    print(io, "Elliptic surface")
  else
    print(io, "Elliptic surface with generic fiber ", equation(E))
  end
end

function Base.show(io::IO, ::MIME"text/plain", S::EllipticSurface)
  io = pretty(io)
  println(io, "Elliptic surface")
  println(io, Indent(), "with generic fiber")
  print(io, Indent(), Lowercase(), generic_fiber(S), Dedent())
  if isdefined(S, :X)
    println(io, "")
    print(io, "with relative minimal model ", Lowercase(), S.Y)
  end
  print(io, Dedent())
end

@doc raw"""

Return the elliptic surface defined by $E$ in the projectivized bundle
P = O(-2s) + O(-3s) + O(1) over P^1 and the corresponding inclusion
S \to P
"""
function weierstrass_model(S::EllipticSurface)
  if isdefined(S, :Weierstrassmodel)
    return S.Weierstrassmodel, S.inc_Weierstrass
  end

  s = S.bundle_number
  E = generic_fiber(S)

  kt = base_ring(base_field(E))
  k = coefficient_ring(kt)
  delta = factor(discriminant(E), kt).fac
  reducible_singular_fibers = [p for p in keys(delta) if delta[p]>1]

  IP1 = projective_space(k, 1)
  c = standard_covering(IP1)
  # rename the variables on the affine charts
  # to a more readable version
  OO(c[1]).S = [:t]
  OO(c[2]).S = [:s]


  O0 = twisting_sheaf(IP1, 0)
  O4 = twisting_sheaf(IP1, -2*s)
  O6 = twisting_sheaf(IP1, -3*s)

  bundleE = direct_sum([O0, O4, O6])

  X_proj = projectivization(bundleE, var_names=["z", "x", "y"])
  X = covered_scheme(X_proj)

  # Create the singular Weierstrass model S of the elliptic K3 surface
  a = a_invars(E)
  U = affine_charts(X)[1]  # the standard Weierstrass chart
  (x, y, t) = gens(OO(U))
  @assert all(denominator(i)==1 for i in a)
  a = [numerator(a)(t) for a in a]
  (a1,a2,a3,a4,a6) = a
  ft = y^2  + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x+a6)
  I = IdealSheaf(X, U, [ft])

  inc_S = CoveredClosedEmbedding(X, I)
  Scov = domain(inc_S)  # The ADE singular elliptic K3 surface
  S.Weierstrasschart = Scov[1][1]

  S.Weierstrassmodel = Scov
  S.inc_Weierstrass = inc_S

  return Scov, inc_S
end


# todo factor into two functions?
# refine covering for blowup
# resolve singularities
function _separate_singularities!(E::EllipticSurface)
  S, inc_S = weierstrass_model(E)
  X = codomain(inc_S)

  I_sing = ideal_sheaf_of_singular_locus(S)
  I_sing_X = radical(pushforward(inc_S)(I_sing))



  # Refine the covering over the reducible singular fibers
  # to make sure that there is only a single singular point in each chart
  refined_charts = AbsSpec[]
  U = X[1][1]  # the weierstrass_chart
  IsingU = I_sing_X(U)
  if isone(IsingU)
    push!(refined_charts, U)
  else
    # there is at most one singularity in every fiber
    # project the singular locus to an affine chart of P1
    disc = gens(eliminate(IsingU, coordinates(U)[1:2]))[1]
    redfib = [f[1] for f in factor(disc)]
    for i in 1:length(redfib)
      r = copy(redfib)
      deleteat!(r, i)
      push!(refined_charts, PrincipalOpenSubset(U, r))
    end
  end

  # Create a chart which contains the fiber over s=0
  # and no other reducible singular fibers
  # these are visible in the charts that we have already
  # i.e. we add the fiber at s=0 and remove all other singular fibers
  V = X[1][4]
  IsingV = I_sing_X(V)
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
    push!(refined_charts, PrincipalOpenSubset(V, redfib))
  end

  # no extra singularities in the X = 1 chart
  # therefore we just exclude all the singularities visible here
  for W in [X[1][2],X[1][5]]
    local Ising = I_sing_X(W)
    if isone(Ising)
      push!(refined_charts, W)
      continue
    end
     (z,y,s_or_t) = coordinates(W)
    # reducible singular fibers
    local disc = gens(eliminate(Ising, [z, s_or_t]))[1]
    local redfib = collect(keys(factor(disc).fac))
    push!(refined_charts, PrincipalOpenSubset(W, redfib))
  end

  # no extra singularities on the the zero section
  # This is the Y = 1 chart
  # therefore we just exclude all the singularities visible here
  for W in [X[1][3],X[1][6]]
    local Ising = I_sing_X(W)
    if isone(Ising)
      push!(refined_charts, W)
      continue
    end
    local (z,x,s_or_t) = coordinates(W)
    # reducible singular fibers
    local disc = gens(eliminate(Ising, [x, s_or_t]))[1]
    local redfib = collect(keys(factor(disc).fac))
    push!(refined_charts, PrincipalOpenSubset(W, redfib))
  end


  Cref = Covering(refined_charts)
  inherit_glueings!(Cref, X[1])
  push!(X.coverings, Cref)
  # Now we have an extra covering where each chart just contains a single singularity

  @assert scheme(I_sing) === S
  @assert scheme(I_sing_X) === X
  return Cref
end


function relatively_minimal_model(E::EllipticSurface)
  if isdefined(E, :blowups)
    return E.Y, E.blowup
  end
  S, inc_S = weierstrass_model(E)
  Crefined = _separate_singularities!(E)
  # Blow up singular points (one at a time) until smooth
  # and compute the strict transforms of the `divisors`
  # collect the exceptional divisors
  # blowup ambient spaces: X0 → X   ⊂
  # blowup pi: (K3 = Y0)  → (S singular weierstrass model)
  #
  # initialization for the while loop
  X0 = codomain(inc_S)
  Y0 = S
  inc_Y0 = inc_S

  exceptionals = []
  varnames = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j,:k,:l]
  projectionsX = []
  projectionsY = []
  count = 0

  @vprint :ellipticK3 2 "Blowing up Weierstrass model\n"
  @vprint :ellipticK3 2 "in $(Crefined)\n"
  while true
    count = count+1
    @vprint :ellipticK3 1 "blowup number: $(count)\n"
    @vprint :ellipticK3 2 "computing singular locus\n"
    I_sing_Y0 = ideal_sheaf_of_singular_locus(Y0)
    @vprint :ellipticK3 2 "decomposing singular locus\n"
    I_sing_Y0 = maximal_associated_points(I_sing_Y0)
    @vprint :ellipticK3 1 "number of singular points: $(length(I_sing_Y0))\n"
    if length(I_sing_Y0)==0
      # stop if smooth
      break
    end
    # take the first singular point and blow it up
    I_sing_X0_1 = radical(pushforward(inc_Y0)(I_sing_Y0[1]))
    if count == 1
      cov = Crefined
    else
      cov = simplified_covering(X0)
    end
    pr_X1 = blow_up(I_sing_X0_1, covering=cov, var_name=varnames[mod(count, length(varnames))])
    X1 = domain(pr_X1)
    @vprint :ellipticK3 1 "$(X1)\n"
    E1 = exceptional_divisor(pr_X1)

    @vprint :ellipticK3 2 "computing strict transforms\n"
    # compute the exceptional divisors
    exceptionals = [strict_transform(pr_X1, e) for e in exceptionals]
    # move the divisors coming originally from S up to the next chart
    push!(exceptionals, E1)

    Y1, inc_Y1, pr_Y1 = strict_transform(pr_X1, inc_Y0)

    push!(projectionsX, pr_X1)
    push!(projectionsY, pr_Y1)
    simplify!(Y1)

    # set up for the next iteration
    Y0 = Y1
    inc_Y0 = inc_Y1
    X0 = X1
  end
  E.Y = Y0
  E.blowups = projectionsY
  E.ambient_blowups = projectionsX
  E.ambient_exceptionals = exceptionals
  piY = reduce(*, reverse(projectionsY))
  E.blowup = piY
  E.inc_Y = inc_Y0
  return Y0, piY
end

#  global divisors0 = [strict_transform(pr_X1, e) for e in divisors0]
#exceptionals_res = [pullback(inc_Y0)(e) for e in exceptionals]

@attr function trivial_lattice(S::EllipticSurface)
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
  kt = parent(d)
  k = coefficient_ring(kt)
  r = [k.(roots(i[1])) for i in factor(d) if i[2]>=2]
  sing = reduce(append!,r, init=[])
  f = []
  f = [[k.([i,1]), fiber_components(S,[i,k(1)])] for  i in sing]
  if degree(d) <= 12*S.bundle_number - 2
    pt = k.([1, 0])
    push!(f, [pt, fiber_components(S, pt)])
  end
  O = zero_section(S)
  F = fiber_components(S, k.([23,1]))
  @assert length(F) == 1 "todo: find an irreducible fiber"
  basisT = [F[1], O]
  @assert S.bundle_number == 2

  grams = [ZZ[0 1;1 -2]]
  # TODO: the -2 selfintersection is probably some K3 artefact
  fiber_components_meeting_O = []
  for (pt, ft) in f
    @vprint :ellipticK3 2 "normalizing fiber: "
    rt, f0, f1, G = standardize_fiber(S, ft)
    @vprint :ellipticK3 2 "$rt \n"
    append!(basisT , f1)
    push!(grams,G)
    push!(fiber_components_meeting_O, (pt, rt, f0))
  end
  G = block_diagonal_matrix(grams)

  return basisT, fiber_components_meeting_O, G
end

function standardize_fiber(S::EllipticSurface, f::Vector{<:WeilDivisor})
  f = copy(f)
  O = zero_section(S)
  for (i,D) in enumerate(f)
    if intersect(O,D)==1
      global f0 = D
      deleteat!(f,i)
      break
    end
  end
  r = length(f)
  G = -2*identity_matrix(ZZ, r)
  @vprint :ellipticK3 2 "computing intersection numbers:"
  for i in 1:r
    @vprint :ellipticK3 2 "\nrow $(i): \n"
    for j in 1:i-1
      @vprint :ellipticK3 2 "$(j) "
      G[i,j] = intersect(f[i],f[j])
      G[j,i] = G[i,j]
    end
  end
  L = integer_lattice(gram=G)
  rt,_ = root_lattice_recognition(L)
  @assert length(rt)==1
  rt = rt[1]
  R = root_lattice(rt[1], rt[2])
  b, I = is_isomorphic_with_permutation(G, -gram_matrix(R))
  @assert b
  return rt, f0, f[I], G[I,I]
end


function fiber_cartier(S::EllipticSurface, P::Vector = ZZ.([0,1]))
  S0,_ = weierstrass_model(S)
  _ = relatively_minimal_model(S) # cache stuff
  D = IdDict{AbsSpec, RingElem}()
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
  F = EffectiveCartierDivisor(S0, D, trivializing_covering=S0[1])
  return pullback(S.blowup)(F)
end

function fiber_components(S::EllipticSurface, P=[0,1])
  @vprint :ellipticK3 2 "computing fiber components over $(P)\n"
  F = fiber_cartier(S, P)
  @vprint :ellipticK3 2 "decomposing fiber\n"
  comp = maximal_associated_points(ideal_sheaf(F))
  return [weil_divisor(c) for c in comp]
end

function section(S::EllipticSurface, P::EllCrvPt)
  S0,incS0 = weierstrass_model(S)
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
  for f in S.ambient_blowups
    PX = strict_transform(f , PX)
  end
  PY = pullback(S.inc_Y, PX)
  return WeilDivisor(PY)
end

@attr zero_section(S::EllipticSurface) = section(S, generic_fiber(S)([0,1,0]))


function is_isomorphic_with_map(G1::Graph, G2::Graph)
  f12 = Polymake.graph.find_node_permutation(G1.pm_graph, G2.pm_graph)
  if isnothing(f12)
    return false, Vector{Int}()
  end
  return true, Polymake.to_one_based_indexing(f12)
end

function graph(G::MatElem)
  n = nrows(G)
  g = Graph{Undirected}(n)
  for i in 1:n
    # small hack to single out the fiber
    if G[i,i]==0
      add_edge!(g,i,i)
    end
    for j in 1:i-1
      if G[i,j] == 1
        add_edge!(g,i,j)
      end
    end
  end
  return g
end

function is_isomorphic_with_permutation(A1::MatElem, A2::MatElem)
  b, T = is_isomorphic_with_map(graph(A1),graph(A2))
  @assert b || A1[T,T] == A2
  return b, T
end
