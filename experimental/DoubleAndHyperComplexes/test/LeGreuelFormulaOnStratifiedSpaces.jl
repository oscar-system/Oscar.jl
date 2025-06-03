########################################################################
# This file contains the code to compute Examples 5.1 and 5.2 from 
# arXiv:2411.02682. 
#
# Above that it serves to test the strand complexes for multigraded 
# rings and modules and the Eagon-Northcott complexes.
########################################################################

@testset "Example 5.1" begin
  # We create a ring with the variables of the ambient space and the parameters.
  P, a = QQ[vcat(["a_{$(i), $j}" for i in 1:2 for j in 1:2], [:x, :y, :z, :w])...]
  a = a[1:4]
  (x, y, z, w) = gens(P)[5:end]

  # We create a multigraded ring for IP^1 x IP^1 over that space.
  S, _ = P[:s, :t, :u, :v]
  S, (s, t, u, v) = grade(S, [[1, 0], [1, 0], [0, 1], [0, 1]])

  # We set up graded modules to represent the required bundles.
  Z2 = grading_group(S)
  S2 = graded_free_module(S, [zero(Z2) for _ in 1:2])
  S2x2 = tensor_product(S2, S2);
  tm = Oscar.tensor_pure_function(S2x2)

  A = S[a[1] a[2] a[3] a[4]]
  enc = Oscar.EagonNorthcottComplex(A; F=S2x2);

  T1 = s*S2[1] + t*S2[2] # tautological bundle of the first factor
  T2 = u*S2[1] + v*S2[2] # tautological bundle of the second factor
  # the sub-bundle serving as the Nash-bundle for the generic determinantal variety
  T, inc = sub(S2x2, [tm(T1, S2[1]), tm(T1, S2[2]), tm(S2[1], T2), tm(S2[2], T2)])

  # We compute the induced Eagon-Northcott complex on the sub-bundle
  rest = Oscar.InducedENC(enc, T);
  # and a free resolution thereof.
  res = Oscar.CartanEilenbergResolution(rest);
  tot_res = total_complex(res);
  tot_res_simp = simplify(tot_res);

  # equations for the Nash transform
  S1 = graded_free_module(S, [zero(Z2)])
  M1 = S[s x y; t z w]
  I1 = ideal(S, minors(M1, 2))
  M2 = S[u v; x y; z w]
  I2 = ideal(S, minors(M2, 2))

  J, inc = (I1 + I2)*S1
  OX = cokernel(inc)
  OX_res, _ = free_resolution(Oscar.SimpleFreeResolution, OX);

  # restriction of the resolution of the Nash bundle
  all_res = tensor_product(OX_res, tot_res_simp);
  tot = simplify(total_complex(all_res));

  tot_0 = Oscar.StrandComplex(tot, zero(Z2))

  # To compute the direct image we need to know where to truncate.
  bnd = Oscar._regularity_bound(tot, 0:10)

  k = Int(bnd[1]) - 1
  irr = ideal(S, [s^k, t^k])

  k = Int(bnd[2]) - 1
  irr = irr*ideal(S, [u^k, v^k]) 

  # the generic direct image complex
  coh = Oscar._derived_pushforward(tot, gens(irr));
  
  @test [ngens(coh.original_complex[i]) for i in -4:10] == [0, 169, 1817, 8412, 22848, 40967, 51448, 46849, 31338, 15027, 4749, 830, 54, 0, 0]
  
  # Test for the cohomology spectral sequence in `experimental/Schemes/SpectralSequences.jl`
  css = Oscar.CohomologySpectralSequence(S, tot);
  cssp = css[4];
  for i in 0:9
    for j in 0:-1:-2
      if i == 0 && j == 0
        @test !is_zero(cssp[i, j])
      else
        @test is_zero(cssp[i, j])
      end
    end
  end

  #= The code below takes too long for the CI.
  # It is kept here for the purpose of reproducibility of the paper's results. 
  
  @test [ngens(coh[i]) for i in 11:-1:-2] == [0, 0, 0, 0, 0, 0, 0, 1, 9, 16, 9, 1, 0, 0]
  
  # We substitute for a specific function.
  R, (x, y, z, w) = QQ[:x, :y, :z, :w]
  k = 5
  l = [-k*x^(k-1), 0, 0, 1]
  subs = hom(P, R, vcat(R.(l), gens(R)))

  tot_0_subs, _ = change_base_ring(subs, tot_0);
  subs_coh, _ = change_base_ring(subs, coh);
  @time H0, _ = simplify(homology(subs_coh, 0)[1]);
  @test vector_space_dimension(H0) == 4

  R, (x, y, z, w) = QQ[:x, :y, :z, :w]
  k = 5
  l = [1, 0, 0, 1]
  subs = hom(P, R, vcat(R.(l), gens(R)))
  subs_coh, _ = change_base_ring(subs, coh);
  @test all(is_zero(homology(subs_coh, i)[1]) for i in 0:5)
  
  =#
end

#= The following tests are disabled to cut down on cost for CI.
# 
# We keep them here to test against regression. As this is part 
# of the data associated to a paper, we will eventually move it 
# to Zenodo in order to adhere to the MarDi principles.
@testset "Example 5.2" begin
  # All steps similar to the above.
  P, v = QQ[vcat(["a_{$(i), $j}" for i in 1:2 for j in 1:2], ["b_{$(i), $j}" for i in 1:2 for j in 1:2], ["c"], [:x, :y, :z, :w])...]
  a = v[1:4]
  b = v[5:8]
  c = v[9]
  (x, y, z, w) = gens(P)[10:end]
  S, _ = P[:s, :t, :u, :v]
  S, (s, t, u, v) = grade(S, [[1, 0], [1, 0], [0, 1], [0, 1]])

  Z2 = grading_group(S)
  S2 = graded_free_module(S, [zero(Z2) for _ in 1:2])
  S2x2 = tensor_product(S2, S2);
  tm = Oscar.tensor_pure_function(S2x2)

  A = S[a[1] a[2] a[3] a[4]; b[1] b[2] b[3] b[4]]
  enc = Oscar.EagonNorthcottComplex(A; F=S2x2);

  T1 = s*S2[1] + t*S2[2]
  T2 = u*S2[1] + v*S2[2]
  T, inc = sub(S2x2, [tm(T1, S2[1]), tm(T1, S2[2]), tm(S2[1], T2), tm(S2[2], T2)])

  rest = Oscar.InducedENC(enc, T);
  res = Oscar.CartanEilenbergResolution(rest);
  tot_res = total_complex(res);
  tot_res_simp = simplify(tot_res);

  S1 = graded_free_module(S, [zero(Z2)])
  C = c*S1[1]
  Kc = koszul_complex(Oscar.KoszulComplex, C);
  K = tensor_product(Kc, tot_res_simp);

  tot = simplify(total_complex(K));

  S1 = graded_free_module(S, [zero(Z2)])
  M1 = S[s x y; t z w]
  I1 = ideal(S, minors(M1, 2))
  M2 = S[u v; x y; z w]
  I2 = ideal(S, minors(M2, 2))

  J, inc = (I1 + I2)*S1
  OX = cokernel(inc)
  OX_res, _ = free_resolution(Oscar.SimpleFreeResolution, OX);
  all_res = tensor_product(OX_res, tot);
  tot = simplify(total_complex(all_res));

  bnd = Oscar._regularity_bound(tot, 0:10)

  k = Int(bnd[1]) - 1
  irr = ideal(S, [s^k, t^k])

  k = Int(bnd[2]) - 1
  irr = irr*ideal(S, [u^k, v^k]) 

  coh = Oscar._derived_pushforward(tot, gens(irr));

  R, (x, y, z, w) = QQ[:x, :y, :z, :w]
  k = 5
  c = [y - z]
  a = [0, 1, -1, 0]
  b = [-k*x^(k-1), 0, 0, 1]
  subs = hom(P, R, vcat(R.(a), R.(b), R.(c), gens(R)))
  subs_coh, _ = change_base_ring(subs, coh);
  @time H0, _ = simplify(homology(subs_coh, 0)[1]);
  @test vector_space_dimension(H0) == 6

  R, (x, y, z, w) = QQ[:x, :y, :z, :w]
  k = 5
  c = [y - z]
  a = [0, 1, -1, 0]
  b = [-1, 0, 0, 1]
  subs = hom(P, R, vcat(R.(a), R.(b), R.(c), gens(R)))
  subs_coh, _ = change_base_ring(subs, coh);
  @time H0, _ = simplify(homology(subs_coh, 0)[1]);
  @test vector_space_dimension(H0) == 2
end
=#

