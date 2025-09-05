
#phylogenetic_ring(varindices::Vector{Tuple{Vararg{Int64}}}; var_name::VarName="p") = phylogenetic_ring(QQ, varindices; var_name=var_name)
#ring(R::PhylogeneticRing) = R.ring
#base_ring(R::PhylogeneticRing) = base_ring(ring(R))

# function Base.show(io::IO, R::PhylogeneticRing)
#   coeffs = base_ring(R)
#   k = ngens(ring(R))
#   print(io, "Phylogenetic ring over $(coeffs) in $(k) variables", "\n", chop(string(gens(ring(R))), head = 16, tail = 1))
# end
# 
################################################################################
# Parametrizations
# function parametrization(pm::PhylogeneticModel; var_name::VarName="p")
#   # The field should be a parameter of the model?
#   F = QQ
#   p = probability_map(pm);
#   p = compute_equivalent_classes(p)
#   parametrization = p.parametrization
#   indices = collect(keys(parametrization))
# 
#   R = phylogenetic_ring(F, indices, var_name=var_name)
#   S = parameter_ring(pm)
#   S, = polynomial_ring(F, vcat([string(x) for x in gens(S)]))
# 
#   hom(ring(R), S, reduce(vcat, [change_coefficient_ring(F, parametrization[k]) for k in indices]))
# end
# 
#parameterization(pm::PhylogeneticModel; var_name::VarName="p") = parametrization(QQ, pm; var_name=var_name)

function parametrization(F::Field, pm::GroupBasedPhylogeneticModel, coordinates::Symbol=:fourier)
  if coordinates == :probabilities
    return parametrization(F, phylogenetic_model(pm); var_name="p")

  elseif coordinates == :fourier
    q = fourier_map(pm);
    q = compute_equivalent_classes(q)
    parametrization = q.parametrization
    indices = collect(keys(parametrization))

    R = phylogenetic_ring(F, indices, var_name="q")
    S = fourier_ring(pm)
    S, = polynomial_ring(F, vcat([string(x) for x in gens(S)]))

    return hom(ring(R), S, reduce(vcat, [change_coefficient_ring(F, parametrization[k]) for k in indices]))
  else
    error("Couldn't recognize coordinates type to be used in parametrization")
  end
end
