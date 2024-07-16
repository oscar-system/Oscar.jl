for i in 1:8, j in 1:number_of_small_groups(i)
  G = small_group(i, j)
  if is_abelian(G)
    continue
  end
  QG = QQ[G]
  ZG = integral_group_ring(QG)
  println(i, " ", j, " ", describe(G), ": ", locally_free_class_group(ZG))
end

# output
┌ Warning: Assignment to `G` in soft scope is ambiguous because a global variable by the same name exists: `G` will be treated as a new local. Disambiguate by using `local G` to suppress this warning or `global G` to assign to the existing global variable.
└ @ none:2
6 1 S3: Z/1
8 3 D8: Z/1
8 4 Q8: Z/2
