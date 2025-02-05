struct PhylogeneticRing
  ring::Ring
  gens::Dict
end

# seems as though this is the same as what is going on for the other rings and we can probably have
# an abstract type like stats ring or something? not sure about the name yet
function phylogenetic_ring(F::Field, varindices::Vector{Tuple{Vararg{Int64}}}; var_name::VarName="p")
  varnames = ["$var_name[$(join(x, ", "))]" for x in varindices]
  S, s = polynomial_ring(F, varnames)
  p = Dict([varindices[i] => s[i] for i in 1:length(varindices)])

  return PhylogeneticRing(S, p)
end

phylogenetic_ring(varindices::Vector{Tuple{Vararg{Int64}}}; var_name::VarName="p") = phylogenetic_model(QQ, varindices; var_name=var_name)
ring(R::PhylogeneticRing) = R.ring
base_ring(R::PhylogeneticRing) = base_ring(ring(R))

function Base.show(io::IO, R::PhylogeneticRing)
  coeffs = base_ring(R)
  k = ngens(ring(R))
  print(io, "Phylogenetic ring over $(coeffs) in $(k) variables", "\n", chop(string(gens(ring(R))), head = 16, tail = 1))
end

################################################################################
# Parametrizations
function parameterization(F::Field, pm::PhylogeneticModel; var_name::VarName="p")
  p = probability_map(pm);
  p = compute_equivalent_classes(p)
  parametrization = p.parametrization
  indices = collect(keys(parametrization))

  R = phylogenetic_ring(F, indices, var_name=var_name)
  S = probability_ring(pm)
  S, = polynomial_ring(F, vcat([string(x) for x in gens(S)]))

  hom(ring(R), S, reduce(vcat, [change_coefficient_ring(F, parametrization[k]) for k in indices]))
end

parameterization(pm::PhylogeneticModel; var_name::VarName="p") = parameterization(QQ, pm; var_name=var_name)

function parameterization(F::Field, pm::GroupBasedPhylogeneticModel, coordinates::Symbol=:fourier)
  if coordinates == :probabilities
    return parameterization(F, phylogenetic_model(pm); var_name="p")

  elseif coordinates == :fourier
    q = fourier_map(pm);
    q = compute_equivalent_classes(q)
    parametrization = q.parametrization
    indices = collect(keys(parametrization))

    R = phylogenetic_ring(F, indices, var_name="q")
    S = fourier_ring(pm)
    S, = polynomial_ring(F, vcat([string(x) for x in gens(S)]))

    return hom(Oscar.ring(R), S, reduce(vcat, [change_coefficient_ring(F, parametrization[k]) for k in indices]))
  else
    error("Couldn't recognize coordinates type to be used in parameterization")
  end
end

################################################################################
# Vanishing Ideal

const PhyloModelUnion = Union{PhylogeneticModel, GroupBasedPhylogeneticModel}

function vanishing_ideal(F::Field, pm::PhyloModelUnion, param::MPolyAnyMap;
                         algorithm::Symbol = :eliminate)
  if algorithm == :kernel
    invariants = kernel(param)
  else
    S = domain(param)
    R = codomain(param)

    #elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat([string(x) for x in gens(R)], [string(x) for x in gens(S)]))
    elim_ring, elim_gens = polynomial_ring(coefficient_ring(R), vcat(R.S, S.S))
    inject_R = hom(R, elim_ring, elim_gens[1:ngens(R)])
    inject_S = hom(S, elim_ring, elim_gens[ngens(R)+1:ngens(elim_ring)])

    I = ideal([inject_S(s) - inject_R(param(s)) for s in gens(S)])

    if algorithm == :eliminate
      invariants = eliminate(I, elim_gens[1:ngens(R)])
    elseif algorithm == :f4
      # F = GF(32003) #1206.6940
      invariants = ideal(groebner_basis_f4(I, eliminate=ngens(R)))

    elseif algorithm == :m4ti2 # symbols cannot start with numbers
      @req is_binomial(I) "Ideal must be a binomial ideal to use 4ti2 algorithm"
      # we can probably avoid constructing the ideal in the future for this case
      lat_gens = map(e -> e[2] - e[1],
                     collect.(map(i -> exponents(i; ordering=lex(elim_ring)),
                                  gens(I))))
      markov_basis = Oscar.markov4ti2(matrix(ZZ, lat_gens))
      invariants = eliminate(binomial_exponents_to_ideal(elim_ring, markov_basis),
                             elim_gens[1:ngens(R)])
      # elseif algorithm == :MultigradedImplicitization
    end

    project_S = hom(elim_ring, S, z -> z, [repeat([1], ngens(R)); gens(S)])
    invariants = project_S(invariants)
  end
  invariants
end

vanishing_ideal(F::Field, pm::PhylogeneticModel;
                algorithm::Symbol=:eliminate) = vanishing_ideal(F, pm, parameterization(F, pm); algorithm=algorithm)

vanishing_ideal(pm::PhylogeneticModel;
                algorithm::Symbol=:eliminate) = vanishing_ideal(QQ, pm; algorithm=algorithm)

vanishing_ideal(pm::GroupBasedPhylogeneticModel;
                coordinates::Symbol=:fourier,
                algorithm::Symbol=:eliminate) = vanishing_ideal(QQ, pm, parameterization(QQ, pm, coordinates); algorithm=algorithm)
