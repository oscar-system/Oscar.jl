export BlowupMorphism
export center, projection, exceptional_divisor

mutable struct BlowupMorphism{
                              CodomainType<:AbsCoveredScheme
                             }
  projective_bundle::CoveredProjectiveScheme 
  codomain::CodomainType   # i.A. ein CoveredScheme
  center::IdealSheaf      # on codomain
  projection::AbsCoveredSchemeMorphism
  domain::AbsCoveredScheme # i.A. ein CoveredScheme
  exceptional_divisor::WeilDivisor

  function BlowupMorphism(
      IP::CoveredProjectiveScheme,
      I::IdealSheaf
    )
    X = base_scheme(IP)
    X === scheme(I) || error("ideal sheaf not compatible with blown up variety")
    return new{typeof(X)}(IP, X, I)
  end
end

function domain(p::BlowupMorphism)
  if !isdefined(p, :domain)
    p.domain = covered_scheme(p.projective_bundle)
  end
  return p.domain
end

codomain(p::BlowupMorphism) = p.codomain
center(p::BlowupMorphism) = p.center

function projection(p::BlowupMorphism)
  if !isdefined(p, :projection)
    p.projection = covered_projection_to_base(p.projective_bundle)
  end
  return p.projection
end

# TODO: Find better name!
covered_projective_scheme(p::BlowupMorphism) = p.projective_bundle

function exceptional_divisor(p::BlowupMorphism)
  if !isdefined(p, :exceptional_divisor)
    Y = domain(p)
    pr = projection(p)
    pr_cov = covering_morphism(pr)
    IC = center(p)
    ideals_on_patches = IdDict{AbsSpec, Ideal}()
    mult_flag = false # To compute the multiplicity only once 
    multiplicity = one(ZZ)
    for U in domain(pr_cov)
      pr_U = pr_cov[U]
      V = codomain(pr_U)
      J = ideal(OO(U), pullback(pr_U).(gens(IC(V))))
      J_rad = radical(J)
      # TODO: Pick the "simplest" chart for the following computation.
      # We could also first use elimpart for simplification.
      if !mult_flag && !(J == J_rad)
        R = base_ring(OO(U))
        S = MPolyComplementOfPrimeIdeal(saturated_ideal(J_rad))
        R_loc, loc_map = localization(R, S)
        F = FreeMod(R_loc, 1)
        JF, _ = sub(F, [g*F[1] for g in gens(saturated_ideal(J))])
        M, _ = quo(F, JF)
        multiplicity = length(M)
        mult_flag = true
      end
      ideals_on_patches[U] = J_rad
    end
    IE = IdealSheaf(Y, ideals_on_patches)
    p.exceptional_divisor = WeilDivisor(Y, ZZ, IdDict{IdealSheaf, elem_type(ZZ)}(IE => multiplicity))
  end
  return p.exceptional_divisor
end

### Auxiliary function which were not there; to be moved

@attr function radical(I::MPolyQuoLocalizedIdeal)
  L = base_ring(I)
  J = pre_image_ideal(I)
  J_rad = radical(J)
  return ideal(L, [g for g in L.(gens(J_rad)) if !iszero(g)])
end

@attr function radical(I::MPolyLocalizedIdeal)
  L = base_ring(I)
  J = pre_saturated_ideal(I)
  J_rad = radical(J)
  return ideal(L, [g for g in L.(gens(J_rad))])
end

@attr function dim(I::MPolyQuoLocalizedIdeal)
  return dim(pre_image_ideal(I))
end

@attr function dim(I::MPolyLocalizedIdeal)
  return dim(saturated_ideal(I))
end
