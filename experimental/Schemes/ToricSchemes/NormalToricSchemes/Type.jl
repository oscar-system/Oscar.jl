export ToricCoveredScheme


@attributes mutable struct ToricCoveredScheme{BRT, CoveredSchemeType<:CoveredScheme} <: AbsCoveredScheme{BRT}
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

        # Now glue X and Y
        
        # Step 1: Find the cone of the intersection of X and Y.
        facet = intersect(CX, CY)
        # harvest the dual cone
        dual_facet = polarize(facet)
        if !(dim(facet) == dim(CX) - 1)
          continue
        end

        # Step 2: Extract the candidates for the localization element.
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
        
        # Step 3: We express the loclization element (it is v, resp. -v) in the Hilbert basis of the two patches.
        # -> Then localize at this element.
        vmat = matrix(ZZ.(collect(v)))
        
        AX = transpose(hilbert_basis(X))
        Id = identity_matrix(ZZ, ncols(AX))
        sol = solve_mixed(AX, vmat, Id, zero(vmat))
        fX = prod([x^k for (x, k) in zip(gens(ambient_coordinate_ring(X)), sol)])
        U = PrincipalOpenSubset(X, OO(X)(fX))
        
        AY = transpose(hilbert_basis(Y))
        Id = identity_matrix(ZZ, ncols(AY))
        sol = solve_mixed(AY, -vmat, Id, zero(vmat))
        fY = prod([x^(k>0 ? 1 : 0) for (x, k) in zip(gens(ambient_coordinate_ring(Y)), sol)])
        V = PrincipalOpenSubset(Y, OO(Y)(fY))
        
        # Step 4: Set up the glueing isomorphisms.
        # (Assumption: v is part of the Hilbert basis.)
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
    
    # TODO: Improve the gluing (lazy gluing) or try to use the Hasse diagram.
    # TOOD: For now, we conjecture, that the composition of the computed glueings is sufficient to deduce all glueings.
    fill_transitions!(cov)
    X = CoveredScheme(cov)
    return new{typeof(base_ring(X)), typeof(X)}(X, ntv)
  end
end

function Base.show(io::IO, X::ToricCoveredScheme)
  print(io, "Scheme of a toric variety with fan spanned by $(rays(fan(X)))")
end
