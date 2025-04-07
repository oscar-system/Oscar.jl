# Some generic helper functions which I could not find in Julia or OSCAR
# Determine the pivots by computing the RREF and then looking for the first nonzero entry of each row
function pivots(M::MatElem)
  (r, R) = rref(M)
  pivots = []

  for i in 1:nrows(R)
      for j in 1:ncols(R)

          if R[i,j] != 0
              append!(pivots, [[i,j]])
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
  elim_ring, z = polynomial_ring(QQ, vcat([string(i) for i in symbols(codom)], [string(i) for i in symbols(dom)]))
  lift_codom = hom(codom, elim_ring, gens(elim_ring)[1:ngens(codom)])
  lift_dom = hom(dom, elim_ring, gens(elim_ring)[ngens(codom) + 1:ngens(elim_ring)])

  return homogeneity_space(ideal([lift_dom(x) - lift_codom(phi(x)) for x in gens(dom)]))
end

# compute the (transpose) of the jacobian
function jacobian(phi)
  return matrix(codomain(phi), [[derivative(phi(j), i) for j in gens(domain(phi))] for i in gens(codomain(phi))])
end

# compute the jacobian and evalauate it at the point pt
function jacobian_at_point(phi, pt)
  return matrix(codomain(phi), [[evaluate(derivative(phi(j), i), pt) for j in gens(domain(phi))] for i in gens(codomain(phi))])
end

# compute the jacobian at a random point with parameters sampled from a finite field
function jacobian_at_rand_point(phi::MPolyAnyMap; char::Int=32003)
  K = GF(char) #TODO: implement rand for fpField so we can replace GF
  
  # remake codomain over finite field
  R = codomain(phi)
  fq_R = polynomial_ring(K, string.(gens(R)))[1]

  # randomly sample parameter values in the finite field
  pt = [rand(K) for _ in 1:ngens(fq_R)]

  return matrix(K, [[evaluate(derivative(fq_R(phi(j)), i), pt) for j in gens(domain(phi))] for i in gens(fq_R)])
end


# find the indices of the variable support of all monomials in mon_basis
function component_support(mon_basis)
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
    gen_shifts = reduce(vcat, [[g*b for b in monomial_basis(domain, deg - degree(g))] for g in prev_gens])
    mons = unique!(reduce(vcat, [collect(monomials(f)) for f in gen_shifts]))
    
    coeffs = matrix(QQ, [[coeff(f, mon) for f in M] for mon in mons])
    bad_monomials = [mons[p[1]] for p in pivots(coeffs)]

    mon_basis = [m for m in monomial_basis(domain(phi), deg) if !(m in bad_monomials)]
  end

  if length(mon_basis) < 2
    return MPolyDecRingElem[]
  end

  # find the indices of all variables involved in this component
  supp = component_support(mon_basis)
  
  # check if the corresponding submatrix drops rank
  # if it does not, then this component cannot contain any generator of the kernel
  if rank(jac[:, supp]) == length(supp)
    return MPolyDecRingElem[]
  end

  image_polys = [phi(m) for m in mon_basis]
  mons = unique!(reduce(vcat, [collect(monomials(f)) for f in image_polys]))
  coeffs = matrix(QQ, [[coeff(f, mon) for f in M] for mon in mons])
  (r, K) = nullspace(coeffs)
  return mon_basis * K
end

struct KernelComponentTask <: ParallelTask
  deg::FinGenAbGroupElem
  phi::MPolyAnyMap
  prev_gens::Vector{<:MPolyDecRingElem}
  jac::MatElem
end

function _compute(kct::KernelComponentTask)
  return component_of_kernel(kct.deg, kct.phi, kct.prev_gens, kct.jac)
end

function components_of_kernel(d::Int, phi::MPolyAnyMap; parallel::Bool=false)
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
    all_degs = unique!(degree.(all_mons))

    if isempty(gens_dict)
      prev_gens = MPolyDecRingElem[]
    else
      prev_gens = reduce(vcat, values(gens_dict))
    end
    if !isempty(all_degs)
      results = parallel_all([KernelComponentTask(deg, graded_phi, prev_gens, jac) for deg in all_degs])
    else
      results =[]
    end

    merge!(gens_dict, Dict(zip(all_degs, results)))
  end

  return gens_dict
end
