function pivots(M::MatElem)
  rref!(M)
  pivots = Vector{Int}[]

  for i in 1:nrows(M)
    for j in 1:ncols(M)
      if !iszero(M[i,j])
        push!(pivots, [i,j])
        break
      end
    end
  end

  return pivots
end

# compute the homogeneity space of an ideal and return it in row reduced echelon form
function homogeneity_space(I::MPolyIdeal)
  n = ngens(base_ring(I))
  homogeneity_eqs = zero_matrix(ZZ, 1, n)

  for f in gens(I)
    exps = [exponent_vector(f,i) for i in 1:length(f)]
    eqs = matrix(ZZ, vcat([exps[1] - exps[i] for i in 2:length(exps)]))
    homogeneity_eqs = vcat(homogeneity_eqs, eqs)
  end
  
  return rref(transpose(nullspace(homogeneity_eqs)[2]))[2]
end

# finds the maximal grading of the kernel of a polynomial map phi
# by computing the homogeneity space of the elimination ideal
function max_grading(phi::MPolyAnyMap)
  codom = codomain(phi)
  dom = domain(phi)
  elim_ring, z = polynomial_ring(QQ, vcat(symbols(codom), symbols(dom)))
  lift_codom = hom(codom, elim_ring, gens(elim_ring)[1:ngens(codom)])
  lift_dom = hom(dom, elim_ring, z[ngens(codom) + 1:ngens(elim_ring)])

  return homogeneity_space(ideal([lift_dom(x) - lift_codom(phi(x)) for x in gens(dom)]))
end

# compute the (transpose) of the jacobian
function jacobian(phi::MPolyAnyMap)
  return matrix(codomain(phi), [[derivative(phi(j), i) for j in gens(domain(phi))] for i in gens(codomain(phi))])
end

# compute the jacobian and evalauate it at the point pt
function jacobian_at_point(phi::MPolyAnyMap, pt::fpFieldElem)
  return matrix(codomain(phi), [[evaluate(derivative(phi(j), i), pt) for j in gens(domain(phi))] for i in gens(codomain(phi))])
end

# compute the jacobian at a random point with parameters sampled from a finite field
function jacobian_at_rand_point(phi::MPolyAnyMap; char::UInt=UInt(32003))
  K = fpField(char)
  
  # remake codomain over finite field
  R = codomain(phi)
  fq_R, fq_R_gens = polynomial_ring(K, ngens(R))

  # randomly sample parameter values in the finite field
  pt = rand(K, ngens(fq_R))

  return matrix(K, [[evaluate(derivative(fq_R(phi(j)), i), pt) for j in gens(domain(phi))] for i in fq_R_gens])
end

# find the indices of the variable support of all monomials in mon_basis
function component_support(mon_basis::Vector{<: MPolyDecRingElem})
  supp = var_index.(unique!(reduce(vcat, vars.(mon_basis))))
end

function component_of_kernel(deg::FinGenAbGroupElem,
                             phi::MPolyAnyMap,
                             prev_gens::Vector{MPolyDecRingElem},
                             jac::MatElem)
  if isempty(prev_gens)
    mon_basis = monomial_basis(domain(phi), deg)
  else
    # computes the basis for domain in degree deg but removes all monomials 
    # which correspond to previously computed relations G
    gen_shifts = reduce(vcat, [
      [g*b for b in monomial_basis(domain(phi), deg - degree(g))]
      for g in prev_gens
        ])
    mons = unique!(reduce(vcat, [collect(monomials(f)) for f in gen_shifts];
                          init=MPolyDecRingElem[]))
    isempty(mons) && return MPolyDecRingElem[]
    
    coeffs = matrix(QQ, [[coeff(f, mon) for f in gen_shifts] for mon in mons])
    bad_monomials = [mons[p[1]] for p in pivots(coeffs)]

    mon_basis = [m for m in monomial_basis(domain(phi), deg) if !(m in bad_monomials)]
  end

  length(mon_basis) < 2 && return MPolyDecRingElem[]

  # find the indices of all variables involved in this component
  supp = component_support(mon_basis)
  
  # check if the corresponding submatrix drops rank
  # if it does not, then this component cannot contain any generator of the kernel
  if rank(jac[:, supp]) == length(supp)
    return MPolyDecRingElem[]
  end

  image_polys = [phi(m) for m in mon_basis]
  mons = unique!(reduce(vcat, [collect(monomials(f)) for f in image_polys]))
  coeffs = matrix(QQ, [[coeff(f, mon) for f in image_polys] for mon in mons])
  (r, K) = nullspace(coeffs)
  return mon_basis * K
end


function compute_component(mon_basis::Vector{<:MPolyDecRingElem}, phi::MPolyAnyMap)
  image_polys = [phi(m) for m in mon_basis]
  mons = unique!(reduce(vcat, [collect(monomials(f)) for f in image_polys]))
  coeffs = matrix(QQ, [[coeff(f, mon) for f in image_polys] for mon in mons])
  (r, K) = nullspace(coeffs)
  return mon_basis * K
end

function filter_component(deg::FinGenAbGroupElem, mon_basis::Vector{<: MPolyDecRingElem}, jac::MatElem)
  # if the basis only has 1 element, then there are no generators by assumption
  length(mon_basis) < 2 && return [true, deg]

  # find the indices of all variables involved in this component
  supp = component_support(mon_basis)

  # check if the corresponding submatrix drops rank
  # if it does not, then this component cannot contain any generator of the kernel
  if rank(jac[:, supp]) == length(supp)
    return [true, deg]
  end

  return [false, deg]
end


function components_of_kernel(d::Int, phi::MPolyAnyMap;
                              wp::Union{OscarWorkerPool, Nothing}=nothing)
  # Compute a maximal grading on the domain
  A = Oscar.max_grading(phi)[:, ngens(codomain(phi)) + 1:end]
  
  # grade the domain
  graded_dom = grade(domain(phi), A)[1]
  
  # grade the domain by the standard grading as well and compute the jacobian at a random point over a finite field
  total_deg_dom = grade(domain(phi), [1 for i in 1:ngens(domain(phi))])[1]
  graded_phi = hom(graded_dom, codomain(phi), [phi(x) for x in gens(domain(phi))])
  jac = Oscar.jacobian_at_rand_point(graded_phi)

  # create a dictionary to store the generators by their degree
  gens_dict = Dict{FinGenAbGroupElem, Vector{<:MPolyDecRingElem}}()

  for i in 1:d
    all_mons = graded_dom.(monomial_basis(total_deg_dom, [i]))
    all_degs = degree.(all_mons)
    mon_bases = Dict{FinGenAbGroupElem, Vector{<:MPolyDecRingElem}}()

    # avoid redundant degree computation
    for (d, m) in zip(all_degs, all_mons)
      if haskey(mon_bases, d)
        push!(mon_bases[d], m)
      else
        mon_bases[d] = [m]
      end
    end

    # find the previous generators 
    if isempty(gens_dict)
      prev_gens = MPolyDecRingElem[]
    else
      prev_gens = reduce(vcat, values(gens_dict))
    end

    # first filter out all easy cases
    # this could be improved to do some load-balancing
    if isnothing(wp)
      filter_results = pmap(filter_component,
                            [deg for deg in keys(mon_bases)],
                            [mon_bases[deg] for deg in keys(mon_bases)],
                            [jac for _ in keys(mon_bases)]; distributed=false)
    else
      filter_results = pmap(filter_component, wp,
                            [deg for deg in keys(mon_bases)],
                            [mon_bases[deg] for deg in keys(mon_bases)],
                            [jac for _ in keys(mon_bases)])
    end

    remain_degs = [result[2] for result in filter_results if !result[1]]

    # now we compute all of the remaining cases which we cannot filter with the Jacobian of phi
    # this could also be improved to do some load-balancing
    if !isempty(remain_degs)
      if isnothing(wp)
        results = pmap(compute_component, map(deg -> mon_bases[deg], remain_degs),
                       [graded_phi for _ in remain_degs];
                       distributed=false)
      else
        results = pmap(compute_component, wp, map(deg -> mon_bases[deg], remain_degs),
                       [graded_phi for _ in remain_degs])
      end
    else
      results =[]
    end

    merge!(gens_dict, Dict(zip(remain_degs, results)))
  end

  return gens_dict
end

