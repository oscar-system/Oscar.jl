  # TODO:
  # Kill all standard Specs
  # shortcut to get the coefficient(s) of a weil divisor
  # intersect Weil with Weil on a surface
  #=
  julia> in_linear_system(linsys[1], D,check=false)
false
 =#
 # subsystem typo |D + P| should be |D - n P|
 # is_subsistem disable check in order_on_divisor

  add_verbosity_scope(:ellipticK3)
  set_verbosity_level(:ellipticK3, 2)
  k = GF(29)
  IP1 = projective_space(k, ["s", "t"])

  O0 = twisting_sheaf(IP1, 0)
  O4 = twisting_sheaf(IP1, -4)
  O6 = twisting_sheaf(IP1, -6)

  E = direct_sum([O0, O4, O6])

  X_proj = projectivization(E, var_names=["z", "x", "y"])

  # remove some charts not needed to cover S
  X = covered_scheme(X_proj)
  C = Covering([X[1][i] for i in [1,3,4,6]])
  for U in patches(C)
    for V in patches(C)
      add_glueing!(C, glueings(X[1])[(U, V)])
    end
  end


  # The generic fiber of the elliptic fibration
  # as an elliptic curve over k(t)
  kt, t = polynomial_ring(k, :t)
  kP1 = fraction_field(kt)
  t = gen(kP1)
  E = EllipticCurve(kP1, [(21*t^7+6*t^6+11*t^4),(21*t^10+15*t^9+17*t^7+18*t^6+t^5)])
  # A basis for the Mordell-Weil group of E
  mwl_basis = [E(collect(i)) for i in [(7*t^3 + 24*t^2 + 9, 9*t^5 + 18*t^4 + 12*t^3 + 8*t^2 + 2),
 (13*t^3 + 9*t + 24, 2*t^5 + 7*t^4 + 6*t^3 + 16*t^2 + 16*t + 22),
 (t^3 + 24*t^2 + 22*t + 5, 10*t^5 + 6*t^4 + 6*t^3 + 8*t^2 + 14*t + 3),
 ((17*t^5 + 14*t^4 + 28*t^2 + 2*t + 1)//(t^2 + 11*t + 23),
  (t^8 + 19*t^7 + 2*t^6 + t^5 + 26*t^4 + 2*t^3 + 26*t + 28)//(t^3 + 2*t^2 + 11*t + 25)),
 ((11*t^5 + 22*t^4 + 23*t^3 + 14*t^2 + 17*t + 16)//(t^2 + 24*t + 28),
  (22*t^8 + 6*t^7 + 21*t^6 + 24*t^5 + 4*t^4 + 23*t^3 + 25*t^2 + 15*t + 6)//(t^3 + 7*t^2 + 26*t + 17))]]

  # The section P used for the fibration hop  Fiber_new = O + P + Vertical
  P = -mwl_basis[1] - mwl_basis[5]

  # Create the singular Weierstrass model S of the elliptic K3 surface
  a = a_invars(E)
  U = affine_charts(X)[1]  # the standard Weierstrass chart
  (x, y, t) = gens(OO(U))
  @assert all(denominator(i)==1 for i in a)
  a = [numerator(a)(t) for a in a]
  (a1,a2,a3,a4,a6) = a
  ft = y^2  + a1*x*y + a3*y - (x^3 + a2*x^2 + a4*x+a6)
  I = IdealSheaf(X, U, [ft])

  inc_S = Oscar.CoveredClosedEmbedding(X, I)
  @test I === image_ideal(inc_S)
  S = domain(inc_S)  # The ADE singular elliptic K3 surface

  I_sing = Oscar.ideal_sheaf_of_singular_locus(S)
  I_sing_X = radical(pushforward(inc_S)(I_sing))

  @test scheme(I_sing) === S
  @test scheme(I_sing_X) === X

  # The mwl_basis as ideal sheaves
  UX = X[1][1]
  sections = [ideal_sheaf(X,X[1][1],[OO(UX)(i) for i in [x*denominator(b[1])(t)-numerator(b[1])(t),y*denominator(b[2])(t)-numerator(b[2])(t)]]) for b in mwl_basis]
  PonX = ideal_sheaf(X,X[1][1],[OO(UX)(i) for i in [x*denominator(P[1])(t)-numerator(P[1])(t),y*denominator(P[2])(t)-numerator(P[2])(t)]])

  #X = CoveredScheme(C)
  (x,y,t) = coordinates(UX)
  fiber = IdealSheaf(X, UX, [t-1,ft]) # we want a smooth fiber
  (z,x,t) = coordinates(X[1][3])
  zero_section = IdealSheaf(X, X[1][3], [x,z])

  # components of the singular fibers visible on S
  (x,y,t) = coordinates(UX)
  A1_0 = IdealSheaf(X, UX, [t, ft])
  (x,y,s) = ambient_coordinates(X[1][4])
  E8_0 = IdealSheaf(X,X[1][4], OO(X[1][4]).([s,gens((modulus(OO(S[1][4]))))[1]]))


  # S = K3 --> X ambient space
  # all relevant divisors as ideal sheaves on X
  divisors = vcat([PonX, zero_section, fiber, A1_0, E8_0], sections)


  # initialization for the while loop
  X0 = X
  Y0 = S
  inc_Y0 = inc_S
  # blow up points (one at a time) until smooth
  exceptionals = []
  divisors0 = divisors
  count = 0
  varnames = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j,:k,:l]
  projectionsX = []
  projectionsY = []
  while true
    global count = count+1
    @vprint :ellipticK3 1 "blowup number: $(count)\n"
    @vprint :ellipticK3 2 "computing singular locus\n"
    I_sing_Y0 = Oscar.ideal_sheaf_of_singular_locus(Y0)
    @vprint :ellipticK3 2 "decomposing singular locus\n"
    I_sing_Y0 = Oscar.maximal_associated_points(I_sing_Y0)
    @vprint :ellipticK3 1 "number of singular points: $(length(I_sing_Y0))\n"
    if length(I_sing_Y0)==0
      # stop if smooth
      break
    end
    # take the first singular point and blow it up
    I_sing_X0_1 = radical(pushforward(inc_Y0)(I_sing_Y0[1]))
    prX1 = blow_up(I_sing_X0_1, covering=Oscar.simplified_covering(X0),
                  var_name=varnames[count])
    X1 = domain(prX1)
    E1 = exceptional_divisor(prX1)

    @vprint :ellipticK3 2 "computing strict transforms\n"
    # compute the exceptional divisors
    global exceptionals = [strict_transform(prX1, e) for e in exceptionals]
    # move the divisors coming originally from S up to the next chart
    global divisors0 = [strict_transform(prX1, e) for e in divisors0]
    push!(exceptionals, E1)

    Y1, inc_Y1, pr_Y1 = strict_transform(prX1, inc_Y0)

    #=
    @vprint :ellipticK3 "computing strict transforms of previous exceptional divisors\n"
    global weil_divs = Any[saturation(pullback(pr_Y1)(e), ideal_sheaf(pullback(inc_Y1)(E1))) for e in weil_divs]
    @vprint :ellipticK3 "computing exceptional weil-divisor\n"
    W1 = [radical(i) for i in keys(coefficient_dict(weil_divisor(pullback(inc_Y1)(E1))))]
    #append!(weil_divs, W1)
    @vprint :ellipticK3 "done\n"
    =#

    push!(projectionsX, prX1)
    push!(projectionsY, pr_Y1)
    simplify!(Y1)

    # set up for the next iteration
    global Y0 = Y1
    global inc_Y0 = inc_Y1
    global X0 = X1
  end
  # Restrict the exceptional divisors:
  # X > S <- Y0 < X0
  @vprint :ellipticK3 2 "Pull back divisors to Y0\n"
  # Pulling back the divisors to Y0 may introduce multiplicities
  exceptionals_res = [pullback(inc_Y0)(e) for e in exceptionals]
  divisors_res = [pullback(inc_Y0)(e) for e in divisors0]

  # Compute the intersection matrix of the exceptional divisors:
  Ex = exceptionals_res
  @vprint :ellipticK3 2 "Exceptional Cartier to Weil divisors\n"
  ExWeil = weil_divisor.(Ex)   # too slow to be a test
  @vprint :ellipticK3 2 "done\n"
  n = length(Ex)
  A = zero_matrix(ZZ, n, n)
  for i in 1:n
    A[i, i] = -ZZ(2)
    for j in i+1:n
     @assert length(components(ExWeil[i]))==1
     ci = coefficient_dict(ExWeil[i])[components(ExWeil[i])[1]]
     cj = coefficient_dict(ExWeil[j])[components(ExWeil[j])[1]]
     A[i, j] = divexact(integral(intersect(ExWeil[i], Ex[j])), ci*cj)
     A[j, i] = A[i,j]
   end
  end
  # The exceptional Weil divisors have multiplicities
  # this comes from the intersection with Y0
  # Compute the reduced one
  tmp = [collect(keys(coefficient_dict(i))) for i in ExWeil]
  @assert all(length(i)==1 for i in tmp)
  ExWeil = [i[1] for i in tmp]
  I = [2,9,4,6,5,8,7,3,1]; A[I,I]
  ExWeil = ExWeil[I]

  # f: Y -- > S the total blowup
  f = projectionsY[1]
  for g in projectionsY[2:end]
    global f = g*f
  end
  (x, y, t) = coordinates(S[1][1])
  fstar = pullback(f[Y0[1][1]])
  # |D| = | P + O + 1 F|
  linsys = prop217(Y0, E, P, 1, fstar(x), fstar(y), fstar(t), fstar(t)-1)
  PonY, OonY, FonY = [WeilDivisor(i, ZZ,check=false) for i in divisors_res[1:3]]
  D = OonY + PonY + FonY
  L = linear_system(linsys, D, check=false)
