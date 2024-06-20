struct PhylogeneticRing
  ring::Ring
  gens::Dict
end

function phylogenetic_ring(varindices::Vector{Tuple{Vararg{Int64}}}; var_name::VarName="p", F::Field=QQ)
  varnames = ["$var_name[$(join(x, ", "))]" for x in varindices] 

  S, s = polynomial_ring(F, varnames)
  p = Dict([varindices[i] => s[i] for i in 1:length(varindices)])

  return PhylogeneticRing(S, p)
end

function Base.show(io::IO, R::PhylogeneticRing)

  coeffs = base_ring(R.ring)
  k = ngens(R.ring)
 
  print(io, "Phylogenetic ring over $(coeffs) in $(k) variables", "\n", chop(string(gens(R.ring)), head = 16, tail = 1))
end

function ring(R::PhylogeneticRing)
  R.ring
end

function parameterization(pm::PhylogeneticModel; var_name::VarName="p", F::Field=QQ)

  p = probability_map(pm);
  p = compute_equivalent_classes(p)
  parametrization = p.parametrization
  indices = collect(keys(parametrization))

  R = phylogenetic_ring(indices, var_name=var_name, F=F)
  S = probability_ring(pm)
  S, = polynomial_ring(F, vcat([string(x) for x in gens(S)]))

  hom(Oscar.ring(R), S, reduce(vcat, [change_coefficient_ring(F, parametrization[k]) for k in indices]))

end

function parameterization(pm::GroupBasedPhylogeneticModel; F::Field=QQ, coordinates::Symbol=:fourier)

  if coordinates == :probabilities
     return parameterization(phylogenetic_model(pm), var_name="p", F=F)

  elseif coordinates == :fourier
      q = fourier_map(pm);
      q = compute_equivalent_classes(q)
      parametrization = q.parametrization
      indices = collect(keys(parametrization))
  
      R = phylogenetic_ring(indices, var_name="q", F=F)
      S = fourier_ring(pm)
      S, = polynomial_ring(F, vcat([string(x) for x in gens(S)]))

      return hom(Oscar.ring(R), S, reduce(vcat, [change_coefficient_ring(F, parametrization[k]) for k in indices]))
  end

end

function vanishing_ideal(pm::PhylogeneticModel; F::Field=QQ, algorithm::Symbol=:eliminate)
  parametrization = parameterization(pm, F=F, coordinates=::probabilities)

  if algorithm == :kernel
      invariants = kernel(parametrization)
  else
      S = domain(parametrization)
      R = codomain(parametrization)

      elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat([string(x) for x in gens(R)], [string(x) for x in gens(S)]))
      inject_R = hom(R, elim_ring, elim_gens[1:ngens(R)])
      inject_S = hom(S, elim_ring, elim_gens[ngens(R)+1:ngens(elim_ring)])

      I = ideal([inject_S(s) - inject_R(parametrization(s)) for s in gens(S)])
     
      if algorithm == :eliminate
          invariants = eliminate(I, elim_gens[1:ngens(R)])
      elseif algorithm == :f4
          # F = GF(32003) #1206.6940
          invariants = ideal(groebner_basis_f4(I, eliminate=ngens(R)))
      # elseif algorithm == :4ti2
      # elseif algorithm == :MultigradedImplicitization
      end

      project_S = hom(elim_ring, S, z -> z, [repeat([1], ngens(R)); gens(S)])
      invariants = project_S(invariants)
  end
  invariants
end

function vanishing_ideal(pm::GroupBasedPhylogeneticModel; F::Field=QQ, coordinates::Symbol=:fourier, algorithm::Symbol = :eliminate)
  parametrization = parameterization(pm, F=F, coordinates=coordinates)

  if algorithm == :kernel
      invariants = kernel(parametrization)
  else
      S = domain(parametrization)
      R = codomain(parametrization)

      elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat([string(x) for x in gens(R)], [string(x) for x in gens(S)]))
      inject_R = hom(R, elim_ring, elim_gens[1:ngens(R)])
      inject_S = hom(S, elim_ring, elim_gens[ngens(R)+1:ngens(elim_ring)])

      I = ideal([inject_S(s) - inject_R(parametrization(s)) for s in gens(S)])
     
      if algorithm == :eliminate
          invariants = eliminate(I, elim_gens[1:ngens(R)])
      elseif algorithm == :f4
          # F = GF(32003) #1206.6940
          invariants = ideal(groebner_basis_f4(I, eliminate=ngens(R)))
      # elseif algorithm == :4ti2
      # elseif algorithm == :MultigradedImplicitization
      end

      project_S = hom(elim_ring, S, z -> z, [repeat([1], ngens(R)); gens(S)])
      invariants = project_S(invariants)
  end
  invariants
end
