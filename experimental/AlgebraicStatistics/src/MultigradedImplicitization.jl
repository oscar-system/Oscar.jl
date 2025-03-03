# Some generic helper functions which I could not find in Julia or OSCAR
# Determine the pivots by computing the RREF and then looking for the first nonzero entry of each row
function pivots(M)
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


# given a vector of monomials mons and a row vector M of polynomials finds a matrix coeffs such that mons*coeffs = M
# this replaces the function coefficients which we use in macaulay2
@everywhere function mons_and_coeffs(mons, M)
    coeffs = [[coeff(f, mons[j]) for f in M] for j in 1:length(mons)]

    return matrix(QQ, coeffs)
end


# compute the homogeneity space of an ideal and return it in row reduced echelon form
function homogeneity_space(I)

  F = gens(I)
  n = length(gens(parent(F[1])))
  homogeneity_eqs = matrix(ZZ, [[0 for i in 1:n]])

  for f in F

    exps = [exponent_vector(f,i) for i in 1:length(f)]
    eqs = matrix(ZZ, vcat([exps[1]-exps[i] for i in 2:length(exps)]))
    homogeneity_eqs = vcat(homogeneity_eqs, eqs)
  end
  
  return rref(transpose(nullspace(homogeneity_eqs)[2]))[2]
end


# finds the maximal grading of the kernel of a polynomial map phi
# by computing the homogeneity space of the elimination ideal
function max_grading(phi)
  codom = codomain(phi)
  dom = domain(phi)
  elim_ring, z = polynomial_ring(QQ, vcat([string(i) for i in symbols(codom)], [string(i) for i in symbols(dom)]))

  lift_codom = hom(codom, elim_ring, gens(elim_ring)[1:ngens(codom)])
  lift_dom = hom(dom, elim_ring, gens(elim_ring)[ngens(codom)+1:ngens(elim_ring)])

  return homogeneity_space(ideal([lift_dom(x) - lift_codom(phi(x)) for x in gens(dom)]))
end


# computes the basis for domain in degree deg but removes all monomials 
# which correspond to previously computed relations G
@everywhere function find_basis_in_degree(domain, deg, prev_gens) 
  if length(prev_gens) == 0
      return monomial_basis(domain, deg)
  end
  gen_shifts = reduce(vcat, [[g*b for b in monomial_basis(domain, deg - degree(g))] for g in prev_gens])
  mons = unique!(reduce(vcat, [collect(monomials(f)) for f in gen_shifts]))
  coeffs = mons_and_coeffs(mons, gen_shifts)
  bad_monomials = [mons[p[1]] for p in pivots(coeffs)]

  return [m for m in monomial_basis(domain, deg) if !(m in bad_monomials)]
end


@everywhere function component_of_kernel(deg, phi, prev_gens)

  mon_basis = find_basis_in_degree(domain(phi), deg, prev_gens)

  if length(mon_basis) < 2
    return MPolyDecRingElem[]
  end

  image_polys = [phi(m) for m in mon_basis]
  mons = unique!(reduce(vcat, [collect(monomials(f)) for f in image_polys]))
  coeffs = mons_and_coeffs(mons, image_polys)
  (r, K) = nullspace(coeffs)
  return mon_basis*K
end


function jacobian(phi)
  return matrix(codomain(phi), [[derivative(phi(j), i) for j in gens(domain(phi))] for i in gens(codomain(phi))])
end


function jacobian_at_point(phi, pt)
  return matrix(codomain(phi), [[evaluate(derivative(phi(j), i), pt) for j in gens(domain(phi))] for i in gens(codomain(phi))])
end


function components_of_kernel(d, phi)
  
  Oscar.put_params(channels, phi)
  A = max_grading(phi)[:, (ngens(codomain(phi))+1):end]

  #
  #A = vcat(matrix(ZZ, [[1 for i in 1:ncols(A)]]), A)
  graded_dom = grade(domain(phi), A)[1]
  grA = grading_group(graded_dom)
  Oscar.put_params(channels, grA)
  Oscar.put_params(channels, graded_dom)
      
  total_deg_dom = grade(domain(phi), [1 for i in 1:ngens(domain(phi))])[1]
  phi = hom(graded_dom, codomain(phi), [phi(x) for x in gens(domain(phi))])
  jac = jacobian_at_point(phi, [rand(KK) for i in codomain()])

  gens_dict = Dict{FinGenAbGroupElem, Vector{<:MPolyDecRingElem}}()

  for i in 1:d
    all_mons = [graded_dom(m) for m in monomial_basis(total_deg_dom, [i])]
    all_degs = unique!([degree(m) for m in all_mons])

    if length(collect(values(gens_dict))) == 0
        prev_gens = MPolyDecRingElem[]
    else
        prev_gens = reduce(vcat, values(gens_dict))
    end
    
    results = pmap(component_of_kernel,
                    all_degs,
                    [phi for _ in all_degs] ,
                    [prev_gens for _ in all_degs])
    merge!(gens_dict, Dict(zip(all_degs, results)))
  end

  return gens_dict
end
