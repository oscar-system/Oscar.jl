
function can_compute(fac::HCTensorProductMapFactory, c::AbsHyperComplex, p::Int, I::Tuple)
  i = collect(I)
  # decompose I into its individual indices
  j = Vector{Vector{Int}}()
  k = 0
  for f in fac.factors
    push!(j, i[k+1:k+dim(f)])
    k = k + dim(f)
  end
  d = length.(j)
  l = findfirst(l->sum(d[1:l]; init=0)>=p, 1:length(j))
  l === nothing && error("mapping direction out of bounds")
  return can_compute_map(fac.factors[l], p - sum(d[1:l-1]; init=0), Tuple(j[l]))
end

