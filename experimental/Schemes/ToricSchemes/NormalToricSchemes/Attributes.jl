@doc Markdown.doc"""
    underlying_scheme(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this returns
the underlying scheme. In other words, by applying
this method, you obtain a scheme that has forgotten
its toric origin.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]

julia> underlying_scheme(toric_scheme)
covered scheme with 3 affine patches in its default covering
```
"""
function underlying_scheme(Z::ToricCoveredScheme)
  if !isdefined(Z, :X)
    F = fan(normal_toric_variety(Z))
    C = maximal_cones(F)

    antv = affine_open_covering(normal_toric_variety(Z))
    var_names = [Symbol.(["x_$(i)_$(k)" for i in 1:ngens(base_ring(toric_ideal(antv[k])))]) for k in 1:length(antv)]
    #rings = [polynomial_ring(coefficient_ring(antv[k]), ["x_$(i)_$(k)" for i in 1:ngens(base_ring(toric_ideal(antv[k])))], cached=false)[1] for k in 1:length(antv)]
    patch_list = [ToricSpec(U, var_names=n) for (U, n) in zip(antv, var_names)]
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
        img_gens = [solve_mixed(AX, AY[:, k], Idext, zero(matrix_space(ZZ, ncols(Idext), 1))) for k in 1:ncols(AY)]
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
        img_gens = [solve_mixed(AY, AX[:, k], Idext, zero(matrix_space(ZZ, ncols(Idext), 1))) for k in 1:ncols(AX)]
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
    # TODO: For now, we conjecture, that the composition of the computed glueings is sufficient to deduce all glueings.
    fill_transitions!(cov)
    X = CoveredScheme(cov)
    Z.X = X
  end
  return Z.X
end


@doc Markdown.doc"""
    normal_toric_variety(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this returns
the underlying normal toric variety.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]

julia> normal_toric_variety(toric_scheme)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor
```
"""
normal_toric_variety(X::ToricCoveredScheme) = X.ntv


@doc Markdown.doc"""
    fan(X::ToricCoveredScheme)

For a toric covered scheme ``X``, this returns
the fan of the underlying toric variety.

# Examples
```jldoctest
julia> P2 = projective_space(NormalToricVariety, 2)
Normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> toric_scheme = ToricCoveredScheme(P2)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0], [0, 1], [-1, -1]]

julia> fan(toric_scheme)
Polyhedral fan in ambient dimension 2
```
"""
fan(X::ToricCoveredScheme) = fan(normal_toric_variety(X))

