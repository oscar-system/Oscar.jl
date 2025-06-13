# This file implements explicit models of complex reflection groups.
#
# References:
#
# * Lehrer, G. I., & Taylor, D. E. (2009). Unitary reflection groups (Vol. 20, p. viii). Cambridge University Press, Cambridge.
# 
# * Thiel, U. (2014). On restricted rational Cherednik algebras. TU Kaiserslautern.
#
# * Marin, I., & Michel, J. (2010). Automorphisms of complex reflection groups. Represent. Theory, 14, 747–788.
#
# Ulrich Thiel, 2023 

function complex_reflection_group(G::ComplexReflectionGroupType, model::Symbol=:default)
  
  # this will be the list of matrix groups corresponding to the components of G
  component_groups = MatrixGroup[]

  # list of models of the components
  modellist = []

  for C in components(G)

    t = C.type[1]

    # Get default model for C
    if model == :default
      if isa(t,Int)
        Cmodel = :Magma
      else
        (m,p,n) = t
        if m == 1 && p == 1
          Cmodel = :CHEVIE
        else
          Cmodel = :Magma
        end
      end
    else
      Cmodel = model
    end
    
    if Cmodel == :LT
      matgrp = complex_reflection_group_LT(t)
    elseif Cmodel == :Magma
      matgrp = complex_reflection_group_Magma(t)
    elseif Cmodel == :CHEVIE
      matgrp = complex_reflection_group_CHEVIE(t)
    else
      error("Specified model not found")
    end

    # set attributes that are already known from type
    set_attribute!(matgrp, :order, order(C))
    set_attribute!(matgrp, :is_complex_reflection_group, true)
    set_attribute!(matgrp, :complex_reflection_group_type, C)
    set_attribute!(matgrp, :complex_reflection_group_model, [Cmodel])
    set_attribute!(matgrp, :is_irreducible, true)

    # add to list
    push!(component_groups, matgrp)
    push!(modellist, model)
  end

  if length(component_groups) == 1
    # treat this as a special case because direct_product([G]) will have type direct
    # product which looks weird for a single group.
    return component_groups[1]
  else
    matgrp = direct_product(component_groups)

    set_attribute!(matgrp, :order, order(G))
    set_attribute!(matgrp, :is_complex_reflection_group, true)
    set_attribute!(matgrp, :complex_reflection_group_type, G)
    set_attribute!(matgrp, :complex_reflection_group_model, modellist)
    set_attribute!(matgrp, :is_irreducible, false)

    return matgrp
  end

end

# Convenience constructors
complex_reflection_group(i::Int, model::Symbol=:default) = complex_reflection_group(ComplexReflectionGroupType(i), model)

complex_reflection_group(m::Int, p::Int, n::Int, model::Symbol=:default) = complex_reflection_group(ComplexReflectionGroupType(m,p,n), model)

complex_reflection_group(X::Vector, model::Symbol=:default) = complex_reflection_group(ComplexReflectionGroupType(X), model)


############################################################################################
# Getter functions
############################################################################################
function complex_reflection_group_type(G::MatrixGroup)
  if has_attribute(G, :complex_reflection_group_type)
    return get_attribute(G, :complex_reflection_group_type)
  end
  return nothing
  # this should be upgraded later to work with a general matrix group (identifying the
  # type from scratch is not so easy though)
end

function complex_reflection_group_model(G::MatrixGroup)
  if has_attribute(G, :complex_reflection_group_model)
    return get_attribute(G, :complex_reflection_group_model)
  end
  return nothing
end

############################################################################################
# The dual of a complex reflection group, i.e. the dual of the natural representation of
# the matrix group. This is of course a general operation of matrix group but this function
# here (applied to a complex reflection group) sets several attributes 
############################################################################################
function complex_reflection_group_dual(W::MatrixGroup)

  if !is_complex_reflection_group(W)
    error("Group is not a complex reflection group")
  end
  WD = matrix_group([transpose(matrix(w^-1)) for w in gens(W)])

  set_attribute!(WD, :order, order(W))
  set_attribute!(WD, :is_complex_reflection_group, true)
  set_attribute!(WD, :complex_reflection_group_type, complex_reflection_group_type(W))

  return WD

end


###########################################################################################
# Checking if a matrix group is a complex reflection group.
###########################################################################################
function is_complex_reflection_group(G::MatrixGroup{T}) where T <: QQAlgFieldElem
  
  # First, check if we already know that G is a complex reflection group
  if has_attribute(G, :is_complex_reflection_group)
    return get_attribute(G, :is_complex_reflection_group)
  end
  
  # Now, check if the (fixed) generating set happens to consist of reflections.
  if all(is_complex_reflection, gens(G))
    set_attribute!(G, :is_complex_reflection_group, true)
    return true
  end

  # Last attempt: determine all reflections and see if they generate the group (this is
  # then equivalent to being a complex reflection group).
  refls = collect(complex_reflections(G))
  H,f = sub(G, [G(matrix(w)) for w in refls])
  if H == G
    set_attribute!(G, :is_complex_reflection_group, true)
    return true
  end

  return false
end

###########################################################################################
# Cartan matrix
###########################################################################################
function complex_reflection_group_cartan_matrix(W::MatrixGroup)

  if !is_complex_reflection_group(W)
    throw(ArgumentError("Group is not a complex reflection group"))
  end

  # We collect roots and coroots of the generators of W
  roots = []
  coroots = []

  for g in gens(W)
    b,g_data = is_complex_reflection_with_data(g)
    push!(roots, root(g_data))
    push!(coroots, coroot(g_data))
  end

  K = base_ring(W)
  n = length(roots)
  C = matrix(K,n,n,[ canonical_pairing(coroots[j], roots[i]) for i=1:n for j=1:n ]) 

  return C
end
