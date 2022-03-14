  kk = QQ
  n=3
  #R, x = PolynomialRing(kk, ["x$i" for i in 0:n])
  #S, _ = grade(R, [1 for i in 0:n])
  Pn = projective_space(kk, n, var_name="x")
  R, _ = PolynomialRing(kk, ["a$(i)$(j)$(k)" for i in 0:n for j in i:n for k in j:n])
  S, a = grade(R, [1 for i in 1:length(gens(R))])
  PW = ProjectiveScheme(S)
  PW_cov = as_covered_scheme(PW)
  U000 = patches(PW_cov)[1]

  function to_ijk_index(l::Int, n::Int)
    table = [(i,j,k) for i in 0:n for j in i:n for k in j:n]
    return table[l]
  end

  function to_linear_index( ijk::Tuple{Int, Int, Int}, n::Int)
    table = [(i,j,k) for i in 0:n for j in i:n for k in j:n]
    for l in 1:length(table)
      table[l] == ijk && return l
    end
    @show ijk
    @show table
    return -1
  end

  a = gens(base_ring(OO(U000)))
  pushfirst!(a, one(base_ring(OO(U000))))

  PolyType = elem_type(base_ring(OO(U000)))
  f_0ii = PolyType[3*a[to_linear_index((0,i,i), n)] - a[to_linear_index((0, 0, i), n)]^2 for i in 1:n]
  f_0ij = PolyType[3*a[to_linear_index((0,i,j), n)] - a[to_linear_index((0, 0, i), n)]*a[to_linear_index((0,0,j), n)] for i in 1:n for j in i+1:n]
  f_iii = PolyType[9*a[to_linear_index((i,i,i), n)] - a[to_linear_index((0, 0, i), n)]*a[to_linear_index((0,i,i), n)] for i in 1:n]
#  f_iij = PolyType[3*b[to_linear_index((i,i,j), n)] - b[to_linear_index((0, i, i), n)]*b[to_linear_index((0,0,j), n)] for i in 1:n for j in i+1:n]
#  f_ijj = PolyType[3*b[to_linear_index((i,j,j), n)] - b[to_linear_index((0, i, i), n)]*b[to_linear_index((0,0,j), n)] for i in 1:n for j in i+1:n]
  f_iij = PolyType[3*a[to_linear_index((i,i,j), n)] - a[to_linear_index((0, i, i), n)]*a[to_linear_index((0,0,j), n)] for i in 1:n for j in i+1:n]
  f_ijj = PolyType[3*a[to_linear_index((i,j,j), n)] - a[to_linear_index((0, j, j), n)]*a[to_linear_index((0,0,i), n)] for i in 1:n for j in i+1:n]
  f_ijk = PolyType[3*a[to_linear_index((i,j,k), n)] - a[to_linear_index((0, i, j), n)]*a[to_linear_index((0,0,k), n)] for i in 1:n for j in i+1:n for k in j+1:n]

  f = f_0ii
  f = vcat(f, f_0ij)
  f = vcat(f, f_iii)
  f = vcat(f, f_iij)
  f = vcat(f, f_ijj)
  f = vcat(f, f_ijk)
  @show length(f)
  var_names_V1 = Vector{Symbol}(vcat(
                                     [Symbol("b0$i$i") for i in 1:n],
                                     [Symbol("b0$i$j") for i in 1:n for j in i+1:n],
                                     [Symbol("b$i$i$i") for i in 1:n],
                                     [Symbol("b$i$i$j") for i in 1:n for j in i+1:n],
                                     [Symbol("b$i$j$j") for i in 1:n for j in i+1:n],
                                     [Symbol("b$i$j$k") for i in 1:n for j in i+1:n for k in j+1:n]
                                    ))
  index_table = vcat(
                     [(0,i, i) for i in 1:n],
                     [(0, i, j) for i in 1:n for j in i+1:n],
                     [(i,i,i) for i in 1:n],
                     [(i,i,j) for i in 1:n for j in i+1:n],
                     [(i,j,j) for i in 1:n for j in i+1:n],
                     [(i,j,k) for i in 1:n for j in i+1:n for k in j+1:n]
                    )

  function b_lookup(a) 
    @show a
    for i in 1:length(index_table)
      a == index_table[i] && return i
    end
    return -1
  end

  IB0 = IdealSheaf(PW_cov, f)

#
#  PolyType = elem_type(S)
#  f_0ii = PolyType[3*a[to_linear_index((0,i,i), n)]*a[1] - a[to_linear_index((0, 0, i), n)]^2 for i in 1:n]
#  f_0ij = PolyType[3*a[to_linear_index((0,i,j), n)]*a[1] - a[to_linear_index((0, 0, i), n)]*a[to_linear_index((0,0,j), n)] for i in 1:n for j in i+1:n]
#  f_iii = PolyType[9*a[to_linear_index((i,i,i), n)]*a[1] - a[to_linear_index((0, 0, i), n)]*a[to_linear_index((0,i,i), n)] for i in 1:n]
#  f_iij = PolyType[3*a[to_linear_index((i,i,j), n)]*a[1] - a[to_linear_index((0, i, i), n)]*a[to_linear_index((0,0,j), n)] for i in 1:n for j in i+1:n]
#  f_ijk = PolyType[3*a[to_linear_index((i,j,k), n)]*a[1] - a[to_linear_index((0, i, j), n)]*a[to_linear_index((0,0,k), n)] for i in 1:n for j in i+1:n for k in j+1:n]
#
#  f = f_0ii
#  f = vcat(f, f_0ij)
#  f = vcat(f, f_iii)
#  f = vcat(f, f_iij)
#  f = vcat(f, f_ijk)
#
#  IB0_alt = IdealSheaf(PW, f)
#  return
  
  Z = subscheme(IB0)
#  CZ = default_covering(Z)
#  for i in 1:npatches(CZ)
#    U = CZ[i]
#    @show i
#    @show U
#    a = is_equidimensional_and_smooth(U)
#    @show a
#    if a
#      @show as_smooth_local_complete_intersection(U)
#    end
#  end
# 
#  D = as_smooth_local_complete_intersection(IB0, verbose=true)
#  U = collect(keys(D))
#  using ProfileView
#  results = []
#  for i in 1:length(U)
#    @show i
#    @show D[U[i]][2]
#    #push!(results, blow_up(U[i], D[U[i]][2]))
#  end
# 
#  line_cond = PolyType[]
#  push!(line_cond, b[to_linear_index((0,0,1), n)]^2*b[to_linear_index((0,1,1), n)]^2 + 18*b[to_linear_index((0,0,1), n)]*b[to_linear_index((1,1,1), n)] - 4*b[to_linear_index((0,0,1), n)]^3*b[to_linear_index((1,1,1), n)] - 27*b[to_linear_index((1,1,1),n)])

#as_smooth_lci_of_cod(Z[1][2], 7, verbose=true)
#D = as_smooth_local_complete_intersection(PW_cov[1][2], IB0[PW_cov[1][2]], verbose=true)
#(P, PC, p, E) = blow_up(PW_cov[1][2], collect(values(D))[1][2], is_regular_sequence=true)
#IB0_prep = as_smooth_lci(IB0, verbose=true, check=false)
(P, PC, p, E) = blow_up(PW_cov[1][1], IB0[PW_cov[1][1]], var_names=var_names_V1, is_regular_sequence=true)
#(P, PC, p, E) = blow_up(PW_cov[1][2], collect(values(D))[1][2])
U = PC[1][1]
#as_smooth_lci_of_cod(U, 6, verbose=true)

#l = _non_degeneration_cover(U, jacobi_matrix(gens(modulus(OO(U)))), 6, check=false, verbose=true);
_, CS = simplify!(PC)


S = homogeneous_coordinate_ring(P)
b = gens(S)

g_iii =  elem_type(S)[3*b[b_lookup((i,i,i))] - 2*a[to_linear_index((0,0,i), n)]*b[b_lookup((0,i,i))] for i in 1:n]
g_iij =  elem_type(S)[3*b[b_lookup((i,i,j))] - a[to_linear_index((0,0,i), n)] * b[b_lookup((0,i,j))] for i in 1:n for j in i+1:n]
g_ijj =  elem_type(S)[3*b[b_lookup((i,j,j))] - a[to_linear_index((0,0,j), n)] * b[b_lookup((0,i,j))] for i in 1:n for j in i+1:n] # correct?
g_ijk =  elem_type(S)[3*b[b_lookup((i,j,k))] - a[to_linear_index((0,0,i), n)] * b[b_lookup((0,j,k))] - a[to_linear_index((0,0,j), n)] * b[b_lookup((0,i,k))] for i in 1:n for j in i+1:n for k in j+1:n]

g = vcat(g_iii, g_iij, g_ijj, g_ijk, g_ijk)
IB1 = IdealSheaf(P, CS, g)

B1 = subscheme(IB1)
@show is_equidimensional_and_smooth(B1)

P2, PC2, p2, E2 = blow_up(IB1)

_, CS2 = simplify!(PC2)
