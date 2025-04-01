@testset "MultigradedImplicitization" begin
  R, (x11, x12, x21, x22, y11, y12, y21, y22) = QQ["x11", "x12", "x21", "x22", "y11", "y12", "y21", "y22"]

  S, (s1, s2, t1, t2, m1, m2, n1, n2) = QQ["s1", "s2", "t1", "t2", "m1", "m2", "n1", "n2"]

  phi = hom(R, S, [s1*t1, s1*t2, s2*t1, s2*t2, m1*n1, m1*n2, m2*n1, m2*n2])

  dicts = components_of_kernel(2, phi)
end
