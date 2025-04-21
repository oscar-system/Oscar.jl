# I would like to credit Dario Mathiä who produced an initial version of the following code

import Chevie: complex_reflection_group, charinfo, CharTable, CyclotomicNumbers
export wreath_macs

# Tools

#gives indices of a subarray
function subarray_indices(sub,arr)
  res = Int[]
  for i in 1:length(sub)
    push!(res,findfirst(x->x==sub[i],arr))
  end
  return res
end

# Computes the b-invariant of a partition
function b1(lambda::Partition)
  if length(lambda)==0
    return 0
  end
  return sum(i->(i-1)*lambda[i],1:length(lambda))
end

# Computes the b-invariant of a multipartition
function btot(lbb::Multipartition)
  return b1(map(lambda -> sum(lambda)[1], lbb)) + length(lbb)*sum(i -> b1(lbb[i]),1:length(lbb))
end

function beta_number_to_partition(beta::Vector{Integer})
  lb=Integer[]
  lbeta=length(beta)
  for j in 1:lbeta
    nbholes=beta[j]-beta[1]-(j-1)
    if nbholes >= 1
    pushfirst!(lb, nbholes)
    end
  end
  return partition(lb)
end

function inv_perm(wperm::PermGroupElem, r::Integer)
   res=Integer[]
   list_wperm=Vector(wperm,r)
   for i in 1:r
      append!(res,indexin(i,list_wperm)[1])
   end
   return perm(res)
end

# Residue notation is in terms of rows and columns

function core(coroot::Vector{Integer},r::Integer)
  beta=Integer[]
  m=minimum(coroot)
  M=maximum(coroot)
  for k in m:M
    for j in 1:r
       if k <= coroot[j]
          append!(beta, k*r+(j-1))
       end
    end
  end
  beta=sort!(beta)
  return beta_number_to_partition(beta)
end

#Inspired by Sagemath's algorithm to obtain a partition from core and quotient

function tau_om(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Integer})
   r=length(lbb)
   w0=perm([r-i for i in 0:(r-1)])
   wperm=inv_perm(wperm,r)
   lbb_perm=multipartition(permuted(lbb.mp,wperm*w0))
   gamma=core([-x for x in coroot],r)
   lg=length(gamma)
   k=r*maximum(length(lbb_perm[i]) for i in 1:r) + lg
   v=[[gamma[i]-(i-1) for i in 1:lg]; [-i for i in lg:(k-1)]]
   w=[[x for x in v if mod((x-i),r) == 0] for i in 1:r]
   new_w=[]
   for i in 1:r
     lw=length(w[i])
     lq=length(lbb_perm[i])
     new_w=[new_w;[w[i][j] + r*lbb_perm[i][j] for j in 1:lq]]
     new_w=[new_w;[w[i][j] for j in lq+1:lw]]
   end
   sort!(new_w,rev=true)
   new_w=[new_w[i]+(i-1) for i in 1:length(new_w)]
   filter!(x-> x!=0,new_w)
   return partition(new_w)
end

#Coroots are given in the canonical basis (e_i) of Z^I.
#Note also that \alpha^{\vee}_i=e_{i-1}-e_i.
#The coroot lattice is thus the lattice of null sum elements in Z^I.

function chev_to_osc(charTirr::Matrix{Cyc{Int64}}, r::Integer, K::AbsSimpleNumField)
  l=size(charTirr)[1]
  res=AbsSimpleNumFieldElem[]
  a=gen(K)
  for j in 1:l
    for i in 1:l
      if r > 2
        coeff=CyclotomicNumbers.coefficients(charTirr[i,j])
        temp=K(0)
        for k in 1:length(coeff)
          temp+=coeff[k]*a^(k-1)
        end
        push!(res,temp)
      else
        push!(res,K(charTirr[i,j]))
      end
    end
  end
  return reshape(res,l,l)
end

#Reoders the CharTable from Chevie
function reorder(charTirr::Matrix{Cyc{Int64}}, r::Integer, n::Integer, Modules::Vector{Vector{Vector{Vector{Int64}}}}, K::AbsSimpleNumField)
  mps=multipartitions(n,r)
  new_ord=Int64[]
    for lbb in mps
      i=findfirst(x-> x==lbb,Modules)
      push!(new_ord,i)
    end
  return chev_to_osc(charTirr,r,K)[new_ord,:]
end

#Tools from representation theory

# Computes the fake degree of a multipartition cf. Ste89 Thm 5.3 and Prop. 3.3.2 Haiman cdm
function fake_deg(lbb::multipartition, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}, var::AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem})
  r=size(lbb)[1]
  n=sum(lbb)[1]
  res=Q(1)
  for p in 1:n
    res=res*(1-var^(r*p))
  end
  for k in 1:r
    for i in 1:length(lbb[k])
      for j in 1:lbb[k][i]
        hookfactor=(Q(1)-var^(r*(1+lbb[k][i]+length(findall(row->row >= j,lbb[k]))-j-i)))
        res=res//hookfactor
      end
    end
  end
  res=res*var^btot(lbb)
  return res
end

# computes C_Delta defined in PhD relation (1.45)
function C_Delta(r::Integer, n::Integer, K::AbsSimpleNumField, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}, var::AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem})
  G=complex_reflection_group(r,1,n)
  Modules=charinfo(G).charparams
  charTable=CharTable(G)
  charTirr=charTable.irr
  Modules=[multipartition(map(x->partition(x),lbb[1])) for lbb in Modules]
  l=length(Modules)
  charTirr=reorder(charTirr,r,n,Modules,K)
  mult=multipartitions(n,r)
  rows=zero_matrix(Q,l,0)
  for i in 1:l
    col=zero_matrix(Q,l,1)
    for j in 1:l
      v=matrix(K,l,1,map(k->charTirr[i,k]*charTirr[j,k],1:l))
      x=solve(matrix(K,transpose(charTirr)),v,side=:right)
      col=col+fake_deg(mult[j],Q,var)*x
    end
    rows=hcat(rows,col)
  end
  return rows
end

# computes the character of the simple representations of the Cherednik algebra.
function C_L(r::Integer, n::Integer, K::AbsSimpleNumField, Q::AbstractAlgebra.Generic.FracField{AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem}}, var::AbstractAlgebra.Generic.MPoly{AbsSimpleNumFieldElem})
  mps=multipartitions(n,r)
  l=length(mps)
  C_D=C_Delta(r,n,K,Q,var)
  diag=[]
  for i in 1:l
    for j in 1:l
      if i==j
        push!(diag,var^btot(mps[i])//fake_deg(mps[i],Q,var))
      else
        push!(diag,0)
      end
    end
  end
  diag=matrix(Q,l,l,diag)
  return diag * C_D
end

function bigger_ord(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Integer})
  r=length(lbb)
  n=sum(lbb)
  mps=multipartitions(n,r)
  lbbquot=tau_om(lbb,wperm,coroot)
  res=map(x->x.mp,filter(x-> dominates(tau_om(x,wperm,coroot),lbbquot),mps))
  return res
end

function smaller_ord(lbb::Multipartition, wperm::PermGroupElem, coroot::Vector{Integer})
  r=length(lbb)
  n=sum(lbb)
  mps=multipartitions(n,r)
  lbbquot=tau_om(lbb,wperm,coroot)
  res=map(x->x.mp,filter(x-> dominates(lbbquot,tau_om(x,wperm,coroot)),mps))
  return res
end

function wreath_macs(n::Integer, r::Integer, wperm::PermGroupElem, coroot::Vector{Integer})
  K,_ = cyclotomic_field(r)
  R, (q,t) = polynomial_ring(K, [:q,:t])
  Q = fraction_field(R)

  mps=multipartitions(n,r)
  l=length(mps)

  c_L=C_L(r,n,K,Q,t)
  c_L_q=map_entries(f->numerator(f)(0,q)//denominator(f)(0,q),c_L)
  c_L_tinv=map_entries(f->numerator(f)(0,1//t)//denominator(f)(0,1//t),c_L)

  rows=[]
  for lbb in mps
    smallers=smaller_ord(lbb, wperm, coroot) #tinv
    biggers=bigger_ord(lbb, wperm, coroot) #q
    smaller_indices=sort!(subarray_indices(smallers,mps))
    bigger_indices=sort!(subarray_indices(biggers,mps))
    sub_smaller_tinv=vcat(map(i->c_L_q[i:i,:],smaller_indices)...)
    sub_bigger_q=vcat(map(i->c_L_tinv[i:i,:],bigger_indices)...)
    M=vcat(sub_smaller_tinv,sub_bigger_q)
    B=kernel(M, side=:left)
    rows=vcat(rows,B[1,1:length(smaller_indices)]*sub_smaller_tinv)
  end
  c_L_qt = matrix(Q,l,l,rows)

  triv=[partition([n])]
  for i in 2:r
    push!(triv,partition([]))
  end
  triv=multipartition(triv)
  index_triv=findfirst(x-> x==triv,mps)
  diag=[]
  for i in 1:l
    for j in 1:l
      if i==j
        push!(diag,1//c_L_qt[i,index_triv])
      else
        push!(diag,Q(0))
      end
    end
  end
  diag=matrix(Q,l,l,diag)
  c_L_qt_H=diag*c_L_qt
  println(multipartitions(n,r))
  return c_L_qt_H
end


#Example:
n = 1
r = 3
wperm=@perm r (3,1,2)
println(multipartitions(n,r))
println(wreath_macs(n,r,wperm,[1,-1,0]))
