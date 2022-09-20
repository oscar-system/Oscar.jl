export ToricSpec

export affine_normal_toric_variety

export ToricCoveredScheme

export normal_toric_variety

@attributes mutable struct ToricSpec{BRT, RT, SpecType<:Spec} <: AbsSpec{BRT, RT}
  X::SpecType
  antv::AffineNormalToricVariety
  dual_cone::Cone
  hb::fmpz_mat # Hilbert basis of the dual cone corresponding to the variables

  function ToricSpec(antv::AffineNormalToricVariety; 
      R::MPolyRing=base_ring(toric_ideal(antv))
    )
    C = cone(antv)
    Cdual = polarize(C)
    hb = matrix(ZZ, hilbert_basis(Cdual))
    I = toric_ideal(R, hb)
    X = Spec(R, I)

    return new{typeof(base_ring(X)), typeof(OO(X)), typeof(X)}(X, antv, Cdual, hb)
  end
end

underlying_scheme(X::ToricSpec) = X.X
affine_normal_toric_variety(X::ToricSpec) = X.antv
#antv(X::ToricSpec) = affine_normal_toric_variety(X)
antv = affine_normal_toric_variety
cone(X::ToricSpec) = cone(antv(X))
dual_cone(X::ToricSpec) = X.dual_cone
hilbert_basis(X::ToricSpec) = X.hb

@Markdown.doc """
    torus_inclusions(X::ToricSpec)

For an affine toric variety ``X`` this returns a list `l` 
containing the inclusions ``Tʳ⁽ⁱ⁾ ↪ X`` of the different 
tori. 
"""
function torus_inclusions(X::ToricSpec)::Vector{<:AbsSpecMor}
  #TODO: Fill in 
end

@Markdown.doc """
    torus_action(X::ToricSpec)

For an affine toric variety ``X`` with a dense open torus ``T``
this returns a quintuple of morphisms `(pT, pX, incX, mult)` 
consisting of 

 * the projection ``T × X → T`` of the product with the torus ``T`` to ``T``
 * the projection ``T × X → X``
 * the inclusion ``X ↪ T × X`` taking ``x`` to ``(1, x)``
 * the group action ``T × X → X``.
"""
function torus_action(X::ToricSpec)::AbsSpecMor
  #TODO: Fill in 
end

function Base.show(io::IO, X::ToricSpec) 
  print(io, "Spec of a toric variety with cone spanned by $(rays(cone(X)))")
end

@attributes mutable struct ToricCoveredScheme{BRT, 
                                              CoveredSchemeType<:CoveredScheme
                                             } <: AbsCoveredScheme{BRT}
  X::CoveredSchemeType
  ntv::NormalToricVariety

  function ToricCoveredScheme(ntv::NormalToricVariety)
    F = fan(ntv)
    C = maximal_cones(F)

    antv = affine_open_covering(ntv)
    rings = [PolynomialRing(coefficient_ring(antv[k]), ["x_$(i)_$(k)" for i in 1:ngens(base_ring(toric_ideal(antv[k])))], cached=false)[1] for k in 1:length(antv)]
    patch_list = [ToricSpec(U, R=R) for (U, R) in zip(antv, rings)]
    cov = Covering(patch_list)

    for i in 1:(length(patch_list)-1)
      X = patch_list[i]
      CX = cone(X)
      CXdual = dual_cone(X)
      for j in i+1:length(patch_list)
        Y = patch_list[j]
        CY = cone(Y)
        CYdual = dual_cone(Y)

        ### Now glue X and Y 
        # the cone of the intersection
        facet = intersect(CX, CY)
        # harvest the dual cone
        dual_facet = polarize(facet)
        dim(facet) == dim(CX) - 1 || continue

        # extract the candidates for the localization element
        candidates = dual_facet.pm_cone.HILBERT_BASIS_GENERATORS[2]
        v = candidates[1,:]
        for j in 1:nrows(candidates)
          v = candidates[j,:]
          if contains(CXdual, v) && contains(CYdual, -v)
            break
          elseif contains(CXdual, -v) && contains(CYdual, v)
            v = -v
            break
          end
        end
        (contains(CXdual, v) && contains(CYdual, -v)) || error("no element found for localization")

        # We have found it and it is v. 
        # Express v (resp. -v) in the hilbert bases of the two patches 
        # and localize at the corresponding elements. 
        vmat = matrix(ZZ.(collect(v)))
        AX = transpose(hilbert_basis(X))
        Id = identity_matrix(ZZ, ncols(AX))
        sol = solve_mixed(AX, vmat, Id, zero(vmat))

        fX = prod([x^k for (x, k) in zip(gens(ambient_ring(X)), sol)])
        U = PrincipalOpenSubset(X, OO(X)(fX))
        
        AY = transpose(hilbert_basis(Y))
        Id = identity_matrix(ZZ, ncols(AY))
        sol = solve_mixed(AY, -vmat, Id, zero(vmat))

        fY = prod([x^(k>0 ? 1 : 0) for (x, k) in zip(gens(ambient_ring(Y)), sol)])
        V = PrincipalOpenSubset(Y, OO(Y)(fY))
        
        # set up the glueing isomorphisms 
        # We assume that the v element must be part of the Hilbert basis 
        # already. 
        l = 0
        for j in 1:ncols(AX)
          if vmat == AX[:, j] 
            l = j
            break
          end
        end
        Idext = identity_matrix(ZZ, ncols(AX))
        Idext[l,l] = 0
        img_gens = [solve_mixed(AX, AY[:, k], Idext, zero(MatrixSpace(ZZ, ncols(Idext), 1))) for k in 1:ncols(AY)]
        fres = hom(OO(V), OO(U), [prod([(k >= 0 ? x^k : inv(x)^(-k)) for (x, k) in zip(gens(OO(U)), w)]) for w in img_gens])
        
        l = 0
        for j in 1:ncols(AY)
          if -vmat == AY[:, j] 
            l = j
            break
          end
        end
        Idext = identity_matrix(ZZ, ncols(AY))
        Idext[l,l] = 0
        img_gens = [solve_mixed(AY, AX[:, k], Idext, zero(MatrixSpace(ZZ, ncols(Idext), 1))) for k in 1:ncols(AX)]
        gres = hom(OO(U), OO(V), [prod([(k >= 0 ? x^k : inv(x)^(-k)) for (x, k) in zip(gens(OO(V)), w)]) for w in img_gens])

        set_attribute!(gres, :inverse, fres)
        set_attribute!(fres, :inverse, gres)
        f = SpecMor(U, V, fres)
        g = SpecMor(V, U, gres)
        set_attribute!(g, :inverse, f)
        set_attribute!(f, :inverse, g)

        G = SimpleGlueing(X, Y, f, g)
        add_glueing!(cov, G)
      end
    end

    X = CoveredScheme(cov)
    return new{typeof(base_ring(X)), typeof(X)}(X, ntv)
  end
end

underlying_scheme(X::ToricCoveredScheme) = X.X
normal_toric_variety(X::ToricCoveredScheme) = X.ntv
ntv = normal_toric_variety
fan(X::ToricCoveredScheme) = fan(normal_toric_variety(X))
