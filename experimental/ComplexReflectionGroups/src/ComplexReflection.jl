

###########################################################################################
# Structure for complex reflections
###########################################################################################
struct ComplexReflection{T <: QQAlgFieldElem}
  base_ring::QQAlgField #this should be the parent of T
  base_space::AbstractAlgebra.Generic.FreeModule{T}
  matrix::MatElem{T}
  root::AbstractAlgebra.Generic.FreeModuleElem{T}
  root_line::AbstractAlgebra.Generic.Submodule{T}
  root_line_inclusion::AbstractAlgebra.Generic.ModuleHomomorphism{T}
  coroot::AbstractAlgebra.Generic.FreeModuleElem{T}
  coroot_form::AbstractAlgebra.Generic.ModuleHomomorphism{T}
  hyperplane::AbstractAlgebra.Generic.Submodule{T}
  hyperplane_inclusion::AbstractAlgebra.Generic.ModuleHomomorphism{T}
  hyperplane_basis::Vector{AbstractAlgebra.Generic.FreeModuleElem{T}}
  eigenvalue::T
  order::Int
  is_unitary::Bool
end

# Printing
function Base.show(io::IO, ::MIME"text/plain", w::ComplexReflection)
  print(io, "Complex reflection of order ", w.order)
end

# Getter functions
function base_ring(w::ComplexReflection)
  return w.base_ring
end

function base_space(w::ComplexReflection)
  return w.base_space
end

function matrix(w::ComplexReflection)
  return w.matrix
end

function root(w::ComplexReflection)
  return w.root
end

function root_line(w::ComplexReflection)
  return w.root_line
end

function root_line_inclusion(w::ComplexReflection)
  return w.root_line_inclusion
end

function coroot(w::ComplexReflection)
  return w.coroot
end

function coroot_form(w::ComplexReflection)
  return w.coroot_form
end

function hyperplane(w::ComplexReflection)
  return w.hyperplane
end

function hyperplane_inclusion(w::ComplexReflection)
  return w.hyperplane_inclusion
end

function hyperplane_basis(w::ComplexReflection)
  return w.hyperplane_basis
end

function eigenvalue(w::ComplexReflection)
  return w.eigenvalue
end

function order(w::ComplexReflection)
  return w.order
end

function is_unitary(w::ComplexReflection)
  return w.is_unitary
end

###########################################################################################
# Construction of complex reflection from root and coroot
###########################################################################################
function complex_reflection(root::AbstractAlgebra.Generic.FreeModuleElem{T}, coroot::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: QQAlgFieldElem

  V = parent(root)
  K = base_ring(V)
  n = dim(V)
  basis = gens(V)

  if root == 0
    throw(ArgumentError("root needs to be nonzero."))
  end

  if coroot == 0
    throw(ArgumentError("coroot needs to be nonzero."))
  end

  zeta = 1 - canonical_pairing(root,coroot)

  if zeta == 1
    throw(ArgumentError("root lies in kernel of coroot."))
  end

  b,d = is_root_of_unity_with_data(zeta)

  if !b
    throw(ArgumentError("Eigenvalue is not a root of unity."))
  end

  w = matrix([ basis[i] - coroot[i]*root for i=1:n ])

  coroot_form = linear_form(coroot)

  H, Hincl = kernel(coroot_form)

  Hbasis = [ Hincl(h) for h in gens(H) ]

  root_line, root_line_inclusion = sub(V, [root])

  return ComplexReflection(K, V, w, root, root_line, root_line_inclusion, coroot, coroot_form, H, Hincl, Hbasis, zeta, d, is_unitary(w))

end

###########################################################################################
# Construction of unitary reflection from root
###########################################################################################
function unitary_reflection(root::AbstractAlgebra.Generic.FreeModuleElem{T}, zeta::T) where T <: QQAlgFieldElem

  V = parent(root)
  K = base_ring(V)
  n = dim(V)

  if root == 0
    throw(ArgumentError("root needs to be nonzero."))
  end

  conj = complex_conjugation(K)

  coeff = (1-zeta)//scalar_product(root,root)

  coroot = coeff*V([ conj(root[i]) for i=1:n ])

  return complex_reflection(root,coroot)

end


function unitary_reflection(root::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: QQAlgFieldElem

  V = parent(root)
  K = base_ring(V)

  zeta = K(-1)
  
  return unitary_reflection(root,zeta)

end

###########################################################################################
# Construction of a complex reflection from matrix
###########################################################################################
function is_complex_reflection_with_data(w::MatElem{T}; debug::Bool=false) where T <: QQAlgFieldElem
  
  if !is_square(w)
    return false, nothing
  end

  n = number_of_rows(w)
  K = base_ring(w)
  I = identity_matrix(K, n)
  V = vector_space(K,n)
  A = I-w
  f = module_homomorphism(V,V,A)

  H, Hincl = kernel(f)

  if dim(H) != n-1
    if debug
      println("Fixed space is not a hyperplane")
    end
    return false, nothing
  end

  # Now that we know that the fixed space of w is a hyperplane, the only remaining thing we
  # need to check is whether w is of finite order and that w is actually an element of GL,
  # i.e. invertible. We first do the following necessary check because we have to compute it
  # for the data anyways.

  R, Rincl = image(f)

  # If H \cap R != {0} then w is not of finite order: if w is of finite order, then
  # image(1-w) is the orthogonal complement of ker(1-w) with respect to a w-invariant scalar
  # product, which exists since w is of finite order.
  
  HcapR, HcapRincl = intersect(H,R)

  if dim(HcapR) > 0
    if debug
      println("Matrix is not of finite order")
    end
    return false, nothing
  end

  # It could still happen that w is not of finite order, e.g. w = matrix(QQ,2,2,[2 0; 0 1]).
  # We now do the actual test whether w is of finite order.

  # If H \cap R = {0}, then R must be 1-dimensional and we take as root a (the) basis
  # vector.
  alpha = Rincl(gens(R)[1])

  # Now, we need to determine the scalar zeta by which w acts on alpha (we do not yet know
  # whether w is a root of unity (which is equivalent to w being of finite order).
  alpha_w = alpha*w

  # We need a non-zero entry of alpha for this
  i=1
  while alpha[i] == 0
    i += 1
  end
  zeta = alpha_w[i]//alpha[i]

  # w is not invertible if and only if zeta == 0.
  # This happens for example for w = matrix(QQ,2,2,[1 0; 0 0])
  if zeta == 0
    if debug
      println("Matrix is not invertible")
    end
    return false, nothing
  end

  # zeta needs to be of finite order (so that w is of finite order), i.e. zeta needs
  # to be a root of unity.
  b, d = is_root_of_unity_with_data(zeta)

  if b == false
    if debug
      println("Matrix is not of finite order")
    end
    return false, nothing
  end

  # Now, w is indeed a reflection. Next, determine the coroot.
  # We have h*(I-w) = 0 for any h \in H. This means that h*v = 0 for any column v of I-w.
  # We thus take L_H to be a non-zero column of I-w, considered as a row vector.
  # Normalization yields the coroot of w with respect to alpha.
  i=1
  while is_zero(A[:,i])
    i += 1
  end
  L_H = V(A[:,i])
  alpha_check = (1-zeta)//canonical_pairing(L_H,alpha) * L_H

  w_data = ComplexReflection(K, V, w, alpha, R, Rincl, alpha_check, linear_form(alpha_check), H, Hincl, [ Hincl(h) for h in gens(H) ], zeta, d, is_unitary(w))

  return true, w_data

end

function is_complex_reflection(w::MatElem{T}; debug::Bool=false) where T <: QQAlgFieldElem

  b,data = is_complex_reflection_with_data(w; debug=debug)
  return b

end

function complex_reflection(w::MatElem{T}) where T <: QQAlgFieldElem

  b,w_data = is_complex_reflection_with_data(w)

  if !b
    throw(ArgumentError("Matrix is not a reflection"))
  end

  return w_data

end

###########################################################################################
# Construction of a complex reflection from matrix group element
###########################################################################################
function is_complex_reflection_with_data(w::MatrixGroupElem{T}) where T <: QQAlgFieldElem

  G = parent(w)

  # If G has reflections assigned (see function below), try to find w in the list and 
  # return the data.
  if has_attribute(G, :complex_reflections)
    refls = collect(get_attribute(G, :complex_reflections))
    i = findfirst(g->matrix(g)==matrix(w), refls)
    if i !== nothing
      w_data = refls[i]
      return true, w_data
    end
  end

  # Otherwise compute
  return is_complex_reflection_with_data(matrix(w))
end

function is_complex_reflection(w::MatrixGroupElem{T}) where T <: QQAlgFieldElem
  b, w_data = is_complex_reflection_with_data(w)
  return b
end

function complex_reflection(w::MatrixGroupElem{T}) where T <: QQAlgFieldElem

  b,w_data = is_complex_reflection_with_data(w)

  if !b
    throw(ArgumentError("Group element is not a complex reflection"))
  end

  return w_data

end

###########################################################################################
# Determining all reflections of a matrix group.
###########################################################################################
function complex_reflections(G::MatrixGroup{T}) where T <: QQAlgFieldElem
   
  if has_attribute(G, :complex_reflections)
    return get_attribute(G, :complex_reflections)
  end

  refls = Set(ComplexReflection[])

  # This is not efficient yet: we should loop only over conjugacy classes
  for g in G
    b,g_data = is_complex_reflection_with_data(g)
    if b 
      push!(refls, g_data)
    end
  end

  set_attribute!(G, :complex_reflections, refls)

  return refls

end
