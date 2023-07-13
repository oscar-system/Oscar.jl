mutable struct AffineMatrixGroup_GlorSl #not a scheme! As soon as coordinate matrix is defined, the parent ring is defined.
  is_special::Bool
  dimension::Int64
  coordinate_matrix::MatSpaceElem
  det_coordinate_matrix::MPolyRingElem
  ring::RingType
  inv_superring::MPolyRing

  function MatrixGroupScheme(sl::Bool, dim::Int64)
    z = new()
    z.is_special = sl
    z.dimension = dim
    return z
  end
end

function coordinate_matrix(G::AffineMatrixGroup_GlorSl)
  K = G.ring
  n = G.dimension
  F,z,x,a = mixed_pol_ring(G)
  #this doesn't ensure that the determinant is non zero. Include that later somehow. 
  M = matrix([z[i,j] for i in 1:n, j in 1:n])
  G.coordinate_matrix = M
  return M
end

function det_coordinate_matrix(G::AffineMatrixGroup_GlorSl)
  if typeof(G.coordinate_matrix) == MatSpaceElem #not sure about this  
    #^ this is to check if the function coordinate_matrix is already called. We don't repeat computation
    det = det(G.coordinate_matrix)
    G.det_coordinate_matrix = det
  else
    M = coordinate_matrix(G)
    det = det(M)
    G.det_coordinate_matrix = det
  end
  return det
end

#returns the polynomial ring needed in the example in page 197 of Derksen 
#F,z,x,a has base ring as G.ring, and variables z[i,,j] (nxn), a[i] as the basis of vector space V
#where V is the vector space of degree "degree" monomials in variables x[1] to x[n]
function mixed_pol_ring(degree::Int64, G::AffineMatrixGroup_GlorSl)
  n = G.dimension
  m = binomial(n + degree - 1, n - 1)
  R = base_ring(parent(V[1]))
  F,z,x,a = PolynomialRing(R, "z" => (1:n, 1:n), "x"=>(1:n), "a" => (1:m))
  return F,z,x,a
end

#do we cache these values?
#mu_star(degre::Int64, G::AffineMatrixGroup_GlorSl)
#computes the image of the basis of degree "degree" under action of G. 
#This returns a dictionary of where the a_i's go under action of the group. 
#The a_i's belong to a really large polynomial ring unfortunately. 

#maybe this is too much computation in one step - this returns mu_star of every element
#of degree basis of degree = "degree". 
function mu_star(degree::Int64, G::AffineMatrixGroup_GlorSl)
  n = G.dimension
  m = binomial(n + degree - 1, n - 1)
  F,z,x,a = mixed_pol_ring(degree,G)
  function new_vars(H::AffineMatrixGroup_GlorSl)
    n = H.dimension
    D = Vector{MPolyElem}(undef,0)
    v = F()
    for j in 1:n
      for i in 1:n
        v = v + x[i]*z[i,j]
      end
      push!(D,v)
      v = F()
    end
    return D
  end
  D = new_vars(G)
  function degree_basis__(d:Int64)
    v = F(1)
    gens = gens(F)[(n^2)1:(n^2)n] #just the x's
    C = zero_matrix(Int64,n,n)
    for i in 1:n
      C[i,i] = -1 #find a better way to write this TODO
    end
    d = [0 for i in 1:n]
    A = [1 for i in 1:n]
    b = [d]
    P = Polyhedron((C,d),(A,b))
    L = lattice_points(P)
    W = Vector{MPolyRingElem}(undef,0)
    for l in L
      for i in 1:n
        v = v*gens[i]^l[i]
      end
    push!(W,v)
    v = F(1)
    end
    return W
  end
  W = degree_basis(degree)
  g = F(1)
  h = F()
  for j in 1:m
    V = exponent_vector(W[j],1)
    for i in 1:n
      if V[n^2 + m + i] != 0
        g = g*D[i]^(V[n^2 + i])
      end
    end
    h += a[j]*g
  end
  H = collect(monimials(h))
  #add all the coefficients of monomials with same degree basis element
  true = 1
  for j in 1:m
    V = exponent_vector(W[j],1) 
    h = F()
    B = Dict([])
    for mon in H
      V2 = exponent_vector(mon,1)
      for i in 1:n 
        if V[n^2 + i] != V2[n^2 + i]
          true = 0
          break
        end
      end
      if true == 1
        h += mon
      end
      true = 1
    end
    k = h // W[j]
    push!(B, a[j] => (W[j],k)) #This dict appears somewhere else make sure the change is made TODO 
  end
  return B
end


#In the end, we implement Algorithm 2.5.8 in sturmfels to compute the invariant ring. 
#Within this alg we use Cohen's omega process to compute the reynold's operator of different monomials. 
function primary_invariants(G::AffineMatrixGroup_GlorSl)
  t = 1
  Q = Vector{MPolyElem}(undef,0)
  n = G.dimension
  R, y = PolynomialRing(G.ring, "y" => 1:n) #our actual variables are y_is #should we cache this ring?? 
  G.inv_superring = R
  #but we will do our mu_star and reynolds computation in terms of a_is from mixed_pol_ring
  while true
    t += 1
    D = mu_star(t, G)
    for basis_elem in degree_basis(R, t) 
      w = reynolds__(basis_elem, D, G, R) 
      if!(w in ideal(R, Q)) 
        push!(Q,w)
      end
      if n == length(Q)
        return Q
      end
    end
  end
end

#omega process
function omegap(p::Int64, G::AffineMatrixGroup_GlorSl, f::MPolyElem)
  det = (G.det_coordinate_matrix)^p
  for (coeff,monomial) in (coefficients(det),monomials(det))
    exp_vect = exponent_vector(monomial, 1)
    x = f
    for i in 1:length(exp_vect)
      for j in 1:exp_vect[i]
      x = derivative(x, gens(mixedpolring)[i]) 
      if x == 0
        break
      end
      end
    end
    h += coeff*x
  end
  return h
end

#this returns reynolds(elem/(detZ)^p). How to return just reynolds(elem)??
#unfinished
function reynolds__(elem::MPolyElem, D::Dict{Any,Any}, G::AffineMatrixGroup_GlorSl, R::MPolyRing)
  V = exponent_vector(elem,1)
  W = collect(keys(D))
  for (elem2) in W
    if exponent_vector(elem2,1)[n^2 + 1: n^2 + n] == V 
      mu_st = get(D, (elem2))[2] #this is in terms of a_is
    end
  end
  #find out degree of mu_star
  d = total_degree(mu_star)
  if !(divides(d,n)[1])
    return R()
  else 
    p = divexact(d,n)
  end
  mixedpolring = parent(mu_st)
  det = (det_coordinate_matrix(G))^p
  h = omegap(p, G, mu_st)
  #h is the numerator
  #we need to divide by c_{p,n}
  cpn = omegap(p, G, G.det_coordinate_matrix)
  elem_by_detzp = h//cpn
  #this is in terms of a_is. Make it in terms of yi_s. 
  result = R()
  x = R(1)
  for monomial in monomials(elem_by_detzp)
    (coeff,factors) = factorisation!(monomial)
    for (factor, deg) in factors
      W = get(D,factor)[1]
      x = x*monom_in_R(W,R)
    end
    result += coeff*x
  end
  return x
end

function monom_in_R(elem::MPolyElem, R::MPolyRing)
  gens = gens(R)
  x = R(1)
  V = exponent_vector(elem, 1)
  for i in 1:length(V)
    if V[i] != 0
      for j in 1:V[i]
        x = x*gens[i - n^2]
      end
    end
  end
  return x
end

#need this to convert elements of type ai to type yi
function factorisation!(elem::MPolyElem)
  coeffs = collect(coefficients(elem))
  length(coeffs) == 1 || @error("can only factorise monomials")
  R = elem.parent
  gens = gens(R)
  factors = Vector{Tuple{MPolyElem, Int64}}
  v = exponent_vector(monomial,1)
  for i in 1:length(v)
    push!(factors, (gens[i], v[i]))
  end
  return coeffs[1], factors
end

#used to compute primary invariants
function degree_basis(R::MPolyRing, t::Int64)
  v = R(1)
  gens = gens(R)
  n = ngens(R)
  C = zero_matrix(Int64,n,n)
  for i in 1:n
    C[i,i] = -1 #find a better way to write this TODO
  end
  d = [0 for i in 1:n]
  A = [1 for i in 1:n]
  b = [d]
  P = Polyhedron((C,d),(A,b))
  L = lattice_points(P)
  W = Vector{MPolyRingElem}(undef,0)
  for l in L
    for i in 1:n
      v = v*gens[i]^l[i]
    end
  push!(W,v)
  v = R(1)
  end
  return W
end

#the user will only use all this to compute the invariant ring. 

