Qx, x = QQ["x"]
k, = number_field(x^4 - 13*x^2 + 16)
candidates = []
for F in abelian_normal_extensions(k, [2], ZZ(10)^13)
  K, = absolute_simple_field(number_field(F))
  if !is_tamely_ramified(K)
    continue
  end
  if !is_isomorphic(galois_group(K)[1], quaternion_group(8))
    continue
  end
  push!(candidates, K)
end
# output
┌ Warning: Assignment to `K` in soft scope is ambiguous because a global variable by the same name exists: `K` will be treated as a new local. Disambiguate by using `local K` to suppress this warning or `global K` to assign to the existing global variable.
└ @ none:2
