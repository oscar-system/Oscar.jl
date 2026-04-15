###################
# (1) Sample spaces
###################

# 2d sample spaces
two_dim_sample_varieties = [
hirzebruch_surface(NormalToricVariety, -3),
hirzebruch_surface(NormalToricVariety, -2),
hirzebruch_surface(NormalToricVariety, -1),
projective_space(NormalToricVariety, 2),
hirzebruch_surface(NormalToricVariety, 1),
hirzebruch_surface(NormalToricVariety, 2),
hirzebruch_surface(NormalToricVariety, 3),
del_pezzo_surface(NormalToricVariety, 1),
del_pezzo_surface(NormalToricVariety, 2),
del_pezzo_surface(NormalToricVariety, 3)]

# Friedemann's beloved toric variety
P1 = convex_hull([0 0;2 2;1 3])
P2 = convex_hull([0 0;2 0;1 2])
P3 = convex_hull([3 0;1 1])
P = minkowski_sum(minkowski_sum(P1, P2), P3)
friedemann_variety = normal_toric_variety(P)

# 3d sample spaces
P1 = projective_space(NormalToricVariety, 1)
three_dim_sample_varieties = [
hirzebruch_surface(NormalToricVariety, -3) * P1,
hirzebruch_surface(NormalToricVariety, -2) * P1,
hirzebruch_surface(NormalToricVariety, -1) * P1,
projective_space(NormalToricVariety, 2) * P1,
hirzebruch_surface(NormalToricVariety, 1) * P1,
hirzebruch_surface(NormalToricVariety, 2) * P1,
hirzebruch_surface(NormalToricVariety, 3) * P1,
del_pezzo_surface(NormalToricVariety, 1) * P1,
del_pezzo_surface(NormalToricVariety, 2) * P1,
del_pezzo_surface(NormalToricVariety, 3) * P1,
projective_space(NormalToricVariety, 3),
#friedemann_variety # currently (November 11, 2025) computationally too challenging for the following scan
]

# compile complete list of sample varieties
sample_spaces = [P1, two_dim_sample_varieties..., three_dim_sample_varieties...]



###################
# (2) Comparer
###################

function all_cohomologies_with_local_cohomology(l::ToricLineBundle)
  g = divisor_class(toric_divisor_class(l))
  cplx = Oscar.cohomology_model(Oscar.ToricCtx(toric_variety(l)), g)
  return [ZZ(rank(cplx[-k])) for k in 0:dim(toric_variety(l))]
end

function martins_comparer(l::ToricLineBundle)
  return all_cohomologies(l) == all_cohomologies_with_local_cohomology(l)
end



###################
# (3) Scanner
###################

R = 1
function run_exhaustive_tester()
  for (i, v) in enumerate(sample_spaces)
    r = rank(picard_group(v))
    all_degrees = Iterators.product((-R:R for _ in 1:r)...)
    for my_deg_tuple in all_degrees
      my_deg = collect(my_deg_tuple)
      l = toric_line_bundle(v, my_deg)
      @req martins_comparer(l) "Found counterexample for min_k computation for variety at position $i and line bundle degree $my_deg"
    end
  end
end

run_exhaustive_tester()
