function _johnson_solid(::Val{9})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 -mre*(sre5+1)//4 1//2;
       1//2 -mre*(sre5+1)//4 -1//2;
       -1//2 -mre*(sre5+1)//4 1//2;
       -1//2 -mre*(sre5+1)//4 -1//2;
       (1+sre5)//4 mre*(sre5-1)//4 1//2;
       (1+sre5)//4 mre*(sre5-1)//4 -1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -1//2;
       0 mre 1//2;
       0 mre -1//2;
       0 0 (1+mre*(sre5-1))//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J9: Elongated pentagonal pyramid")
  return res
end

function _johnson_solid(::Val{10})
  Qx, x = QQ["x"]
  NF, qr8 = number_field(x^4 - 8, "a")
  ENF, qre8 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  sre2 = (qre8^2)//2
  V = [0 0 (2*sre2 + qre8)//4;
       1//2 1//2 qre8//4;
       1//2 -1//2 qre8//4;
       -1//2 1//2 qre8//4;
       -1//2 -1//2 qre8//4;
       0 sre2//2 -qre8//4;
       0 -sre2//2 -qre8//4;
       sre2//2 0 -qre8//4;
       -sre2//2 0 -qre8//4]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J10: Gyroelongated square pyramid")
  return res
end

function _johnson_solid(::Val{13})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 -mre*(sre5+1)//4 0;
       -1//2 -mre*(sre5+1)//4 0;
       (1+sre5)//4 mre*(sre5-1)//4 0;
       -(1+sre5)//4 mre*(sre5-1)//4 0;
       0 mre 0;
       0 0 mre*(sre5-1)//2;
       0 0 -mre*(sre5-1)//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J13: Pentagonal bipyramid")
  return res
end

function _johnson_solid(::Val{16})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 -mre*(sre5+1)//4 1//2;
       1//2 -mre*(sre5+1)//4 -1//2;
       -1//2 -mre*(sre5+1)//4 1//2;
       -1//2 -mre*(sre5+1)//4 -1//2;
       (1+sre5)//4 mre*(sre5-1)//4 1//2;
       (1+sre5)//4 mre*(sre5-1)//4 -1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -1//2;
       0 mre 1//2;
       0 mre -1//2;
       0 0 (1+mre*(sre5-1))//2;
       0 0 -(1+mre*(sre5-1))//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J16: Elongated pentagonal bipyramid")
  return res
end

function _johnson_solid(::Val{17})
  Qx, x = QQ["x"]
  NF, qr8 = number_field(x^4 - 8, "a")
  ENF, qre8 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  sre2 = (qre8^2)//2
  V = [0 0 (2*sre2 + qre8)//4;
       0 0 -(2*sre2 + qre8)//4;
       1//2 1//2 qre8//4;
       1//2 -1//2 qre8//4;
       -1//2 1//2 qre8//4;
       -1//2 -1//2 qre8//4;
       0 sre2//2 -qre8//4;
       0 -sre2//2 -qre8//4;
       sre2//2 0 -qre8//4;
       -sre2//2 0 -qre8//4]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J17: Gyroelongated square bipyramid")
  return res
end

function _johnson_solid(::Val{18})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 sre3//2 1//2;
       1//2 sre3//2 -1//2;
       1//2 -sre3//2 1//2;
       1//2 -sre3//2 -1//2;
       -1//2 sre3//2 1//2;
       -1//2 sre3//2 -1//2;
       -1//2 -sre3//2 1//2;
       -1//2 -sre3//2 -1//2;
       1 0 1//2;
       1 0 -1//2;
       -1 0 1//2;
       -1 0 -1//2;
       1//2 -sre3//6 (3+2*sre2*sre3)//6;
       -1//2 -sre3//6 (3+2*sre2*sre3)//6;
       0 sre3//3 (3+2*sre2*sre3)//6]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J18: Elongated triangular cupola")
  return res
end

function _johnson_solid(::Val{20})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5-1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5-1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5-1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5-1))//2;
       0 mre (1+mre*(sre5-1))//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J20: Elongated pentagonal cupola")
  return res
end

function _johnson_solid(::Val{21})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       0 mre (1+mre*(sre5+1))//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       0 -mre*(sre5+1)//2 (1+2*mre)//2;
       (1+sre5)//4 mre*(sre5+3)//4 (1+2*mre)//2;
       -(1+sre5)//4 mre*(sre5+3)//4 (1+2*mre)//2;
       (3+sre5)//4 -mre//2 (1+2*mre)//2;
       -(3+sre5)//4 -mre//2 (1+2*mre)//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J21: Elongated pentagonal rotunda")
  return res
end

function _johnson_solid(::Val{22})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  sr2, sr3 = srv
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (sr3-1), "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre2 = EMF(sr2)
  sre3 = EMF(sr3)
  V = [1//2 -sre3//6 (2*sre2*sre3+3*mre)//6;
       -1//2 -sre3//6 (2*sre2*sre3+3*mre)//6;
       0 sre3//3 (2*sre2*sre3+3*mre)//6;
       1//2 sre3//2 mre//2;
       1//2 -sre3//2 mre//2;
       -1//2 sre3//2 mre//2;
       -1//2 -sre3//2 mre//2;
       1 0 mre//2;
       -1 0 mre//2;
       sre3//2 1//2 -mre//2;
       sre3//2 -1//2 -mre//2;
       -sre3//2 1//2 -mre//2;
       -sre3//2 -1//2 -mre//2;
       0 1 -mre//2;
       0 -1 -mre//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J22: Gyroelongated triangular cupola")
  return res
end

function _johnson_solid(::Val{23})
  NF, sr2 = quadratic_field(2)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (sr2+2), "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-2 - 2*sr2 + (sr2+2)*mr)//8, "b")
  ELF, lre = Hecke.embedded_field(LF, real_embeddings(LF)[6])
  sre2 = ELF(sr2)
  mre = ELF(mr)
  V = [1//2 1//2 (sre2+2*lre)//2;
       1//2 -1//2 (sre2+2*lre)//2;
       -1//2 1//2 (sre2+2*lre)//2;
       -1//2 -1//2 (sre2+2*lre)//2;
       1//2 (1+sre2)//2 lre;
       1//2 -(1+sre2)//2 lre;
       -1//2 (1+sre2)//2 lre;
       -1//2 -(1+sre2)//2 lre;
       (1+sre2)//2 1//2 lre;
       (1+sre2)//2 -1//2 lre;
       -(1+sre2)//2 1//2 lre;
       -(1+sre2)//2 -1//2 lre;
       0 mre//sre2 -lre;
       0 -mre//sre2 -lre;
       mre//sre2 0 -lre;
       -mre//sre2 0 -lre;
       mre//2 mre//2 -lre;
       mre//2 -mre//2 -lre;
       -mre//2 mre//2 -lre;
       -mre//2 -mre//2 -lre]
  res = convex_hull(ELF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J23: Gyroelongated square cupola")
  return res
end

function _johnson_solid(::Val{24})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-4 - 2*sr5 + (3*sr5+5)*mr)//8, "b")
  ELF, lre = Hecke.embedded_field(LF, real_embeddings(LF)[4])
  sre5 = ELF(sr5)
  mre = ELF(mr)
  V = [1//2 -mre*(sre5+1)//4 mre*(sre5-1)//2+lre;
       -1//2 -mre*(sre5+1)//4 mre*(sre5-1)//2+lre;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5-1)//2+lre;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5-1)//2+lre;
       0 mre mre*(sre5-1)//2+lre;
       1//2 mre*(sre5+5)//4 lre;
       1//2 -mre*(sre5+5)//4 lre;
       -1//2 mre*(sre5+5)//4 lre;
       -1//2 -mre*(sre5+5)//4 lre;
       (3+sre5)//4 sre5*mre//2 lre;
       (3+sre5)//4 -sre5*mre//2 lre;
       -(3+sre5)//4 sre5*mre//2 lre;
       -(3+sre5)//4 -sre5*mre//2 lre;
       (1+sre5)//2 0 lre;
       -(1+sre5)//2 0 lre;
       mre*(sre5+5)//4 1//2 -lre;
       mre*(sre5+5)//4 -1//2 -lre;
       -mre*(sre5+5)//4 1//2 -lre;
       -mre*(sre5+5)//4 -1//2 -lre;
       sre5*mre//2 (3+sre5)//4 -lre;
       sre5*mre//2 -(3+sre5)//4 -lre;
       -sre5*mre//2 (3+sre5)//4 -lre;
       -sre5*mre//2 -(3+sre5)//4 -lre;
       0 (1+sre5)//2 -lre;
       0 -(1+sre5)//2 -lre]
  res = convex_hull(ELF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J24: Gyroelongated pentagonal cupola")
  return res
end

function _johnson_solid(::Val{25})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-4 - 2*sr5 + (3*sr5+5)*mr)//8, "b")
  ELF, lre = Hecke.embedded_field(LF, real_embeddings(LF)[4])
  sre5 = ELF(sr5)
  mre = ELF(mr)
  V = [1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2+lre;
       -1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2+lre;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2+lre;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2+lre;
       0 mre mre*(sre5+1)//2+lre;
       (1+sre5)//4 mre*(sre5+3)//4 mre+lre;
       -(1+sre5)//4 mre*(sre5+3)//4 mre+lre;
       (3+sre5)//4 -mre//2 mre+lre;
       -(3+sre5)//4 -mre//2 mre+lre;
       0 -mre*(sre5+1)//2 mre+lre;
       1//2 mre*(sre5+5)//4 lre;
       1//2 -mre*(sre5+5)//4 lre;
       -1//2 mre*(sre5+5)//4 lre;
       -1//2 -mre*(sre5+5)//4 lre;
       (3+sre5)//4 sre5*mre//2 lre;
       (3+sre5)//4 -sre5*mre//2 lre;
       -(3+sre5)//4 sre5*mre//2 lre;
       -(3+sre5)//4 -sre5*mre//2 lre;
       (1+sre5)//2 0 lre;
       -(1+sre5)//2 0 lre;
       mre*(sre5+5)//4 1//2 -lre;
       mre*(sre5+5)//4 -1//2 -lre;
       -mre*(sre5+5)//4 1//2 -lre;
       -mre*(sre5+5)//4 -1//2 -lre;
       sre5*mre//2 (3+sre5)//4 -lre;
       sre5*mre//2 -(3+sre5)//4 -lre;
       -sre5*mre//2 (3+sre5)//4 -lre;
       -sre5*mre//2 -(3+sre5)//4 -lre;
       0 (1+sre5)//2 -lre;
       0 -(1+sre5)//2 -lre]
  res = convex_hull(ELF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J25: Gyroelongated pentagonal rotunda")
  return res
end

function _johnson_solid(::Val{30})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 -mre*(sre5+1)//4 mre*(sre5-1)//2;
       1//2 -mre*(sre5+1)//4 -mre*(sre5-1)//2;
       -1//2 -mre*(sre5+1)//4 mre*(sre5-1)//2;
       -1//2 -mre*(sre5+1)//4 -mre*(sre5-1)//2;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5-1)//2;
       (1+sre5)//4 mre*(sre5-1)//4 -mre*(sre5-1)//2;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5-1)//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -mre*(sre5-1)//2;
       0 mre mre*(sre5-1)//2;
       0 mre -mre*(sre5-1)//2;
       1//2 mre*(sre5+5)//4 0;
       1//2 -mre*(sre5+5)//4 0;
       -1//2 mre*(sre5+5)//4 0;
       -1//2 -mre*(sre5+5)//4 0;
       (3+sre5)//4 sre5*mre//2 0;
       (3+sre5)//4 -sre5*mre//2 0;
       -(3+sre5)//4 sre5*mre//2 0;
       -(3+sre5)//4 -sre5*mre//2 0;
       (1+sre5)//2 0 0;
       -(1+sre5)//2 0 0]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J30: Pentagonal orthobicupola")
  return res
end

function _johnson_solid(::Val{32})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [0 mre mre*(sre5+1)//2;
       1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2;
       -1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2;
       0 -mre*(sre5+1)//2 mre;
       (1+sre5)//4 mre*(sre5+3)//4 mre;
       -(1+sre5)//4 mre*(sre5+3)//4 mre;
       (3+sre5)//4 -mre//2 mre;
       -(3+sre5)//4 -mre//2 mre;
       1//2 mre*(sre5+5)//4 0;
       1//2 -mre*(sre5+5)//4 0;
       -1//2 mre*(sre5+5)//4 0;
       -1//2 -mre*(sre5+5)//4 0;
       (3+sre5)//4 sre5*mre//2 0;
       (3+sre5)//4 -sre5*mre//2 0;
       -(3+sre5)//4 sre5*mre//2 0;
       -(3+sre5)//4 -sre5*mre//2 0;
       (1+sre5)//2 0 0;
       -(1+sre5)//2 0 0;
       1//2 -mre*(sre5+1)//4 -mre*(sre5-1)//2;
       -1//2 -mre*(sre5+1)//4 -mre*(sre5-1)//2;
       (1+sre5)//4 mre*(sre5-1)//4 -mre*(sre5-1)//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -mre*(sre5-1)//2;
       0 mre -mre*(sre5-1)//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J32: Pentagonal orthocupolarotunda")
  return res
end

function _johnson_solid(::Val{33})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [0 mre mre*(sre5+1)//2;
       1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2;
       -1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2;
       0 -mre*(sre5+1)//2 mre;
       (1+sre5)//4 mre*(sre5+3)//4 mre;
       -(1+sre5)//4 mre*(sre5+3)//4 mre;
       (3+sre5)//4 -mre//2 mre;
       -(3+sre5)//4 -mre//2 mre;
       1//2 mre*(sre5+5)//4 0;
       1//2 -mre*(sre5+5)//4 0;
       -1//2 mre*(sre5+5)//4 0;
       -1//2 -mre*(sre5+5)//4 0;
       (3+sre5)//4 sre5*mre//2 0;
       (3+sre5)//4 -sre5*mre//2 0;
       -(3+sre5)//4 sre5*mre//2 0;
       -(3+sre5)//4 -sre5*mre//2 0;
       (1+sre5)//2 0 0;
       -(1+sre5)//2 0 0;
       1//2 mre*(sre5+1)//4 -mre*(sre5-1)//2;
       -1//2 mre*(sre5+1)//4 -mre*(sre5-1)//2;
       (1+sre5)//4 -mre*(sre5-1)//4 -mre*(sre5-1)//2;
       -(1+sre5)//4 -mre*(sre5-1)//4 -mre*(sre5-1)//2;
       0 -mre -mre*(sre5-1)//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J33: Pentagonal gyrocupolarotunda")
  return res
end

function _johnson_solid(::Val{34})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [0 mre mre*(sre5+1)//2;
       0 mre -mre*(sre5+1)//2;
       1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2;
       1//2 -mre*(sre5+1)//4 -mre*(sre5+1)//2;
       -1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2;
       -1//2 -mre*(sre5+1)//4 -mre*(sre5+1)//2;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2;
       (1+sre5)//4 mre*(sre5-1)//4 -mre*(sre5+1)//2;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -mre*(sre5+1)//2;
       0 -mre*(sre5+1)//2 mre;
       0 -mre*(sre5+1)//2 -mre;
       (1+sre5)//4 mre*(sre5+3)//4 mre;
       (1+sre5)//4 mre*(sre5+3)//4 -mre;
       -(1+sre5)//4 mre*(sre5+3)//4 mre;
       -(1+sre5)//4 mre*(sre5+3)//4 -mre;
       (3+sre5)//4 -mre//2 mre;
       (3+sre5)//4 -mre//2 -mre;
       -(3+sre5)//4 -mre//2 mre;
       -(3+sre5)//4 -mre//2 -mre;
       1//2 mre*(sre5+5)//4 0;
       1//2 -mre*(sre5+5)//4 0;
       -1//2 mre*(sre5+5)//4 0;
       -1//2 -mre*(sre5+5)//4 0;
       (3+sre5)//4 sre5*mre//2 0;
       (3+sre5)//4 -sre5*mre//2 0;
       -(3+sre5)//4 sre5*mre//2 0;
       -(3+sre5)//4 -sre5*mre//2 0;
       (1+sre5)//2 0 0;
       -(1+sre5)//2 0 0]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J34: Pentagonal orthobirotunda")
  return res
end

function _johnson_solid(::Val{35})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 sre3//2 1//2;
       1//2 sre3//2 -1//2;
       1//2 -sre3//2 1//2;
       1//2 -sre3//2 -1//2;
       -1//2 sre3//2 1//2;
       -1//2 sre3//2 -1//2;
       -1//2 -sre3//2 1//2;
       -1//2 -sre3//2 -1//2;
       1 0 1//2;
       1 0 -1//2;
       -1 0 1//2;
       -1 0 -1//2;
       1//2 -sre3//6 (3+2*sre2*sre3)//6;
       1//2 -sre3//6 -(3+2*sre2*sre3)//6;
       -1//2 -sre3//6 (3+2*sre2*sre3)//6;
       -1//2 -sre3//6 -(3+2*sre2*sre3)//6;
       0 sre3//3 (3+2*sre2*sre3)//6;
       0 sre3//3 -(3+2*sre2*sre3)//6]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J35: Elongated triangular orthobicupola")
  return res
end

function _johnson_solid(::Val{36})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 sre3//2 1//2;
       1//2 sre3//2 -1//2;
       1//2 -sre3//2 1//2;
       1//2 -sre3//2 -1//2;
       -1//2 sre3//2 1//2;
       -1//2 sre3//2 -1//2;
       -1//2 -sre3//2 1//2;
       -1//2 -sre3//2 -1//2;
       1 0 1//2;
       1 0 -1//2;
       -1 0 1//2;
       -1 0 -1//2;
       1//2 -sre3//6 (3+2*sre2*sre3)//6;
       1//2 sre3//6 -(3+2*sre2*sre3)//6;
       -1//2 -sre3//6 (3+2*sre2*sre3)//6;
       -1//2 sre3//6 -(3+2*sre2*sre3)//6;
       0 sre3//3 (3+2*sre2*sre3)//6;
       0 -sre3//3 -(3+2*sre2*sre3)//6]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J36: Elongated triangular gyrobicupola")
  return res
end

function _johnson_solid(::Val{38})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5-1))//2;
       1//2 -mre*(sre5+1)//4 -(1+mre*(sre5-1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5-1))//2;
       -1//2 -mre*(sre5+1)//4 -(1+mre*(sre5-1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5-1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 -(1+mre*(sre5-1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5-1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -(1+mre*(sre5-1))//2;
       0 mre (1+mre*(sre5-1))//2;
       0 mre -(1+mre*(sre5-1))//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J38: Elongated pentagonal orthobicupola")
  return res
end

function _johnson_solid(::Val{39})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5-1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5-1))//2;
       -1//2 mre*(sre5+1)//4 -(1+mre*(sre5-1))//2;
       1//2 mre*(sre5+1)//4 -(1+mre*(sre5-1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5-1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5-1))//2;
       -(1+sre5)//4 -mre*(sre5-1)//4 -(1+mre*(sre5-1))//2;
       (1+sre5)//4 -mre*(sre5-1)//4 -(1+mre*(sre5-1))//2;
       0 mre (1+mre*(sre5-1))//2;
       0 -mre -(1+mre*(sre5-1))//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J39: Elongated pentagonal gyrobicupola")
  return res
end

function _johnson_solid(::Val{40})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       0 mre (1+mre*(sre5+1))//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       0 -mre*(sre5+1)//2 mre+1//2;
       (1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       -(1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       (3+sre5)//4 -mre//2 mre+1//2;
       -(3+sre5)//4 -mre//2 mre+1//2;
       1//2 -mre*(sre5+1)//4 (-1-mre*(sre5-1))//2;
       -1//2 -mre*(sre5+1)//4 (-1-mre*(sre5-1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (-1-mre*(sre5-1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (-1-mre*(sre5-1))//2;
       0 mre (-1-mre*(sre5-1))//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J40: Elongated pentagonal orthocupolarotunda")
  return res
end

function _johnson_solid(::Val{41})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       0 mre (1+mre*(sre5+1))//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       0 -mre*(sre5+1)//2 mre+1//2;
       (1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       -(1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       (3+sre5)//4 -mre//2 mre+1//2;
       -(3+sre5)//4 -mre//2 mre+1//2;
       1//2 mre*(sre5+1)//4 (-1-mre*(sre5-1))//2;
       -1//2 mre*(sre5+1)//4 (-1-mre*(sre5-1))//2;
       (1+sre5)//4 -mre*(sre5-1)//4 (-1-mre*(sre5-1))//2;
       -(1+sre5)//4 -mre*(sre5-1)//4 (-1-mre*(sre5-1))//2;
       0 -mre (-1-mre*(sre5-1))//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J41: Elongated pentagonal gyrocupolarotunda")
  return res
end

function _johnson_solid(::Val{42})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       0 mre (1+mre*(sre5+1))//2;
       0 mre -(1+mre*(sre5+1))//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       1//2 -mre*(sre5+1)//4 -(1+mre*(sre5+1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       -1//2 -mre*(sre5+1)//4 -(1+mre*(sre5+1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 -(1+mre*(sre5+1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -(1+mre*(sre5+1))//2;
       0 -mre*(sre5+1)//2 mre+1//2;
       0 -mre*(sre5+1)//2 -mre-1//2;
       (1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       (1+sre5)//4 mre*(sre5+3)//4 -mre-1//2;
       -(1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       -(1+sre5)//4 mre*(sre5+3)//4 -mre-1//2;
       (3+sre5)//4 -mre//2 mre+1//2;
       (3+sre5)//4 -mre//2 -mre-1//2;
       -(3+sre5)//4 -mre//2 mre+1//2;
       -(3+sre5)//4 -mre//2 -mre-1//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J42: Elongated pentagonal orthobirotunda")
  return res
end

function _johnson_solid(::Val{43})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 mre*(sre5+5)//4 1//2;
       1//2 mre*(sre5+5)//4 -1//2;
       1//2 -mre*(sre5+5)//4 1//2;
       1//2 -mre*(sre5+5)//4 -1//2;
       -1//2 mre*(sre5+5)//4 1//2;
       -1//2 mre*(sre5+5)//4 -1//2;
       -1//2 -mre*(sre5+5)//4 1//2;
       -1//2 -mre*(sre5+5)//4 -1//2;
       (3+sre5)//4 sre5*mre//2 1//2;
       (3+sre5)//4 sre5*mre//2 -1//2;
       (3+sre5)//4 -sre5*mre//2 1//2;
       (3+sre5)//4 -sre5*mre//2 -1//2;
       -(3+sre5)//4 sre5*mre//2 1//2;
       -(3+sre5)//4 sre5*mre//2 -1//2;
       -(3+sre5)//4 -sre5*mre//2 1//2;
       -(3+sre5)//4 -sre5*mre//2 -1//2;
       (1+sre5)//2 0 1//2;
       (1+sre5)//2 0 -1//2;
       -(1+sre5)//2 0 1//2;
       -(1+sre5)//2 0 -1//2;
       0 mre (1+mre*(sre5+1))//2;
       0 -mre -(1+mre*(sre5+1))//2;
       1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       1//2 mre*(sre5+1)//4 -(1+mre*(sre5+1))//2;
       -1//2 -mre*(sre5+1)//4 (1+mre*(sre5+1))//2;
       -1//2 mre*(sre5+1)//4 -(1+mre*(sre5+1))//2;
       (1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       (1+sre5)//4 -mre*(sre5-1)//4 -(1+mre*(sre5+1))//2;
       -(1+sre5)//4 mre*(sre5-1)//4 (1+mre*(sre5+1))//2;
       -(1+sre5)//4 -mre*(sre5-1)//4 -(1+mre*(sre5+1))//2;
       0 -mre*(sre5+1)//2 mre+1//2;
       0 mre*(sre5+1)//2 -mre-1//2;
       (1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       (1+sre5)//4 -mre*(sre5+3)//4 -mre-1//2;
       -(1+sre5)//4 mre*(sre5+3)//4 mre+1//2;
       -(1+sre5)//4 -mre*(sre5+3)//4 -mre-1//2;
       (3+sre5)//4 -mre//2 mre+1//2;
       (3+sre5)//4 mre//2 -mre-1//2;
       -(3+sre5)//4 -mre//2 mre+1//2;
       -(3+sre5)//4 mre//2 -mre-1//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J43: Elongated pentagonal gyrobirotunda")
  return res
end

function _johnson_solid(::Val{44})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  sr2, sr3 = srv
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (sr3-1), "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre2 = EMF(sr2)
  sre3 = EMF(sr3)
  V = [1//2 -sre3//6 (2*sre2*sre3+3*mre)//6;
       -1//2 -sre3//6 (2*sre2*sre3+3*mre)//6;
       0 sre3//3 (2*sre2*sre3+3*mre)//6;
       1//2 sre3//2 mre//2;
       1//2 -sre3//2 mre//2;
       -1//2 sre3//2 mre//2;
       -1//2 -sre3//2 mre//2;
       1 0 mre//2;
       -1 0 mre//2;
       sre3//2 1//2 -mre//2;
       sre3//2 -1//2 -mre//2;
       -sre3//2 1//2 -mre//2;
       -sre3//2 -1//2 -mre//2;
       0 1 -mre//2;
       0 -1 -mre//2;
       -sre3//6 1//2 -(2*sre2*sre3+3*mre)//6;
       -sre3//6 -1//2 -(2*sre2*sre3+3*mre)//6;
       sre3//3 0 -(2*sre2*sre3+3*mre)//6]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J44: Gyroelongated triangular bicupola")
  return res
end

function _johnson_solid(::Val{45})
  NF, sr2 = quadratic_field(2)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (sr2+2), "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-2 - 2*sr2 + (sr2+2)*mr)//8, "b")
  ELF, lre = Hecke.embedded_field(LF, real_embeddings(LF)[6])
  sre2 = ELF(sr2)
  mre = ELF(mr)
  V = [1//2 1//2 (sre2+2*lre)//2;
       1//2 -1//2 (sre2+2*lre)//2;
       -1//2 1//2 (sre2+2*lre)//2;
       -1//2 -1//2 (sre2+2*lre)//2;
       1//2 (1+sre2)//2 lre;
       1//2 -(1+sre2)//2 lre;
       -1//2 (1+sre2)//2 lre;
       -1//2 -(1+sre2)//2 lre;
       (1+sre2)//2 1//2 lre;
       (1+sre2)//2 -1//2 lre;
       -(1+sre2)//2 1//2 lre;
       -(1+sre2)//2 -1//2 lre;
       0 mre//sre2 -lre;
       0 -mre//sre2 -lre;
       mre//sre2 0 -lre;
       -mre//sre2 0 -lre;
       mre//2 mre//2 -lre;
       mre//2 -mre//2 -lre;
       -mre//2 mre//2 -lre;
       -mre//2 -mre//2 -lre;
       mre//(2*sre2) mre*(1-sre2//2)//2 -(sre2+2*lre)//2;
       -mre//(2*sre2) -mre*(1-sre2//2)//2 -(sre2+2*lre)//2;
       -mre*(1-sre2//2)//2 mre//(2*sre2) -(sre2+2*lre)//2;
       mre*(1-sre2//2)//2 -mre//(2*sre2) -(sre2+2*lre)//2]
  res = convex_hull(ELF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J45: Gyroelongated square bicupola")
  return res
end

function _johnson_solid(::Val{46})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-4 - 2*sr5 + (3*sr5+5)*mr)//8, "b")
  ELF, lre = Hecke.embedded_field(LF, real_embeddings(LF)[4])
  sre5 = ELF(sr5)
  mre = ELF(mr)
  V = [1//2 -mre*(sre5+1)//4 mre*(sre5-1)//2+lre;
       -1//2 -mre*(sre5+1)//4 mre*(sre5-1)//2+lre;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5-1)//2+lre;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5-1)//2+lre;
       0 mre mre*(sre5-1)//2+lre;
       1//2 mre*(sre5+5)//4 lre;
       1//2 -mre*(sre5+5)//4 lre;
       -1//2 mre*(sre5+5)//4 lre;
       -1//2 -mre*(sre5+5)//4 lre;
       (3+sre5)//4 sre5*mre//2 lre;
       (3+sre5)//4 -sre5*mre//2 lre;
       -(3+sre5)//4 sre5*mre//2 lre;
       -(3+sre5)//4 -sre5*mre//2 lre;
       (1+sre5)//2 0 lre;
       -(1+sre5)//2 0 lre;
       mre*(sre5+5)//4 1//2 -lre;
       mre*(sre5+5)//4 -1//2 -lre;
       -mre*(sre5+5)//4 1//2 -lre;
       -mre*(sre5+5)//4 -1//2 -lre;
       sre5*mre//2 (3+sre5)//4 -lre;
       sre5*mre//2 -(3+sre5)//4 -lre;
       -sre5*mre//2 (3+sre5)//4 -lre;
       -sre5*mre//2 -(3+sre5)//4 -lre;
       0 (1+sre5)//2 -lre;
       0 -(1+sre5)//2 -lre;
       -mre*(sre5+1)//4 1//2 -mre*(sre5-1)//2-lre;
       -mre*(sre5+1)//4 -1//2 -mre*(sre5-1)//2-lre;
       mre*(sre5-1)//4 (1+sre5)//4 -mre*(sre5-1)//2-lre;
       mre*(sre5-1)//4 -(1+sre5)//4 -mre*(sre5-1)//2-lre;
       mre 0 -mre*(sre5-1)//2-lre]
  res = convex_hull(ELF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J46: Gyroelongated pentagonal bicupola")
  return res
end

function _johnson_solid(::Val{47})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-4 - 2*sr5 + (3*sr5+5)*mr)//8, "b")
  ELF, lre = Hecke.embedded_field(LF, real_embeddings(LF)[4])
  sre5 = ELF(sr5)
  mre = ELF(mr)
  V = [1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2+lre;
       -1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2+lre;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2+lre;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2+lre;
       0 mre mre*(sre5+1)//2+lre;
       (1+sre5)//4 mre*(sre5+3)//4 mre+lre;
       -(1+sre5)//4 mre*(sre5+3)//4 mre+lre;
       (3+sre5)//4 -mre//2 mre+lre;
       -(3+sre5)//4 -mre//2 mre+lre;
       0 -mre*(sre5+1)//2 mre+lre;
       1//2 mre*(sre5+5)//4 lre;
       1//2 -mre*(sre5+5)//4 lre;
       -1//2 mre*(sre5+5)//4 lre;
       -1//2 -mre*(sre5+5)//4 lre;
       (3+sre5)//4 sre5*mre//2 lre;
       (3+sre5)//4 -sre5*mre//2 lre;
       -(3+sre5)//4 sre5*mre//2 lre;
       -(3+sre5)//4 -sre5*mre//2 lre;
       (1+sre5)//2 0 lre;
       -(1+sre5)//2 0 lre;
       mre*(sre5+5)//4 1//2 -lre;
       mre*(sre5+5)//4 -1//2 -lre;
       -mre*(sre5+5)//4 1//2 -lre;
       -mre*(sre5+5)//4 -1//2 -lre;
       sre5*mre//2 (3+sre5)//4 -lre;
       sre5*mre//2 -(3+sre5)//4 -lre;
       -sre5*mre//2 (3+sre5)//4 -lre;
       -sre5*mre//2 -(3+sre5)//4 -lre;
       0 (1+sre5)//2 -lre;
       0 -(1+sre5)//2 -lre;
       -mre*(sre5+1)//4 1//2 -mre*(sre5-1)//2-lre;
       -mre*(sre5+1)//4 -1//2 -mre*(sre5-1)//2-lre;
       mre*(sre5-1)//4 (1+sre5)//4 -mre*(sre5-1)//2-lre;
       mre*(sre5-1)//4 -(1+sre5)//4 -mre*(sre5-1)//2-lre;
       mre 0 -mre*(sre5-1)//2-lre]
  res = convex_hull(ELF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J47: Gyroelongated pentagonal cupolarotunda")
  return res
end

function _johnson_solid(::Val{48})
  NF, sr5 = quadratic_field(5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-4 - 2*sr5 + (3*sr5+5)*mr)//8, "b")
  ELF, lre = Hecke.embedded_field(LF, real_embeddings(LF)[4])
  sre5 = ELF(sr5)
  mre = ELF(mr)
  V = [1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2+lre;
       -1//2 -mre*(sre5+1)//4 mre*(sre5+1)//2+lre;
       (1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2+lre;
       -(1+sre5)//4 mre*(sre5-1)//4 mre*(sre5+1)//2+lre;
       0 mre mre*(sre5+1)//2+lre;
       (1+sre5)//4 mre*(sre5+3)//4 mre+lre;
       -(1+sre5)//4 mre*(sre5+3)//4 mre+lre;
       (3+sre5)//4 -mre//2 mre+lre;
       -(3+sre5)//4 -mre//2 mre+lre;
       0 -mre*(sre5+1)//2 mre+lre;
       1//2 mre*(sre5+5)//4 lre;
       1//2 -mre*(sre5+5)//4 lre;
       -1//2 mre*(sre5+5)//4 lre;
       -1//2 -mre*(sre5+5)//4 lre;
       (3+sre5)//4 sre5*mre//2 lre;
       (3+sre5)//4 -sre5*mre//2 lre;
       -(3+sre5)//4 sre5*mre//2 lre;
       -(3+sre5)//4 -sre5*mre//2 lre;
       (1+sre5)//2 0 lre;
       -(1+sre5)//2 0 lre;
       mre*(sre5+5)//4 1//2 -lre;
       mre*(sre5+5)//4 -1//2 -lre;
       -mre*(sre5+5)//4 1//2 -lre;
       -mre*(sre5+5)//4 -1//2 -lre;
       sre5*mre//2 (3+sre5)//4 -lre;
       sre5*mre//2 -(3+sre5)//4 -lre;
       -sre5*mre//2 (3+sre5)//4 -lre;
       -sre5*mre//2 -(3+sre5)//4 -lre;
       0 (1+sre5)//2 -lre;
       0 -(1+sre5)//2 -lre;
       -mre*(sre5+1)//4 1//2 -mre*(sre5+1)//2-lre;
       -mre*(sre5+1)//4 -1//2 -mre*(sre5+1)//2-lre;
       mre*(sre5-1)//4 (1+sre5)//4 -mre*(sre5+1)//2-lre;
       mre*(sre5-1)//4 -(1+sre5)//4 -mre*(sre5+1)//2-lre;
       mre 0 -mre*(sre5+1)//2-lre;
       mre*(sre5+3)//4 (1+sre5)//4 -mre-lre;
       mre*(sre5+3)//4 -(1+sre5)//4 -mre-lre;
       -mre//2 (3+sre5)//4 -mre-lre;
       -mre//2 -(3+sre5)//4 -mre-lre;
       -mre*(sre5+1)//2 0 -mre-lre]
  res = convex_hull(ELF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J48: Gyroelongated pentagonal birotunda")
  return res
end

function _johnson_solid(::Val{49})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 -sre3//6 1//2;
       1//2 -sre3//6 -1//2;
       -1//2 -sre3//6 1//2;
       -1//2 -sre3//6 -1//2;
       0 sre3//3 1//2;
       0 sre3//3 -1//2;
       0 -(3*sre2+sre3)//6 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J49: Augmented triangular prism")
  return res
end

function _johnson_solid(::Val{50})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 -sre3//6 1//2;
       1//2 -sre3//6 -1//2;
       -1//2 -sre3//6 1//2;
       -1//2 -sre3//6 -1//2;
       0 sre3//3 1//2;
       0 sre3//3 -1//2;
       (1+sre2*sre3)//4 (3*sre2+sre3)//12 0;
       -(1+sre2*sre3)//4 (3*sre2+sre3)//12 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J50: Biaugmented triangular prism")
  return res
end

function _johnson_solid(::Val{51})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 -sre3//6 1//2;
       1//2 -sre3//6 -1//2;
       -1//2 -sre3//6 1//2;
       -1//2 -sre3//6 -1//2;
       0 sre3//3 1//2;
       0 sre3//3 -1//2;
       (1+sre2*sre3)//4 (3*sre2+sre3)//12 0;
       -(1+sre2*sre3)//4 (3*sre2+sre3)//12 0;
       0 -(3*sre2+sre3)//6 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J51: Triaugmented triangular prism")
  return res
end

function _johnson_solid(::Val{52})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 5, x^2 - 2], ["sqrt(5)", "sqrt(2)"])
  sr5, sr2 = srv
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[8])
  sre2 = EMF(sr2)
  sre5 = EMF(sr5)
  V = [1//2 -mre*(sre5+1)//4 1//2;
       1//2 -mre*(sre5+1)//4 -1//2;
       -1//2 -mre*(sre5+1)//4 1//2;
       -1//2 -mre*(sre5+1)//4 -1//2;
       (1+sre5)//4 mre*(sre5-1)//4 1//2;
       (1+sre5)//4 mre*(sre5-1)//4 -1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -1//2;
       0 mre 1//2;
       0 mre -1//2;
       0 -(sre2+mre*(sre5+1)//2)//2 0]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J52: Augmented pentagonal prism")
  return res
end

function _johnson_solid(::Val{53})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 5, x^2 - 2], ["sqrt(5)", "sqrt(2)"])
  sr5, sr2 = srv
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10, "a")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[8])
  sre2 = EMF(sr2)
  sre5 = EMF(sr5)
  V = [1//2 -mre*(sre5+1)//4 1//2;
       1//2 -mre*(sre5+1)//4 -1//2;
       -1//2 -mre*(sre5+1)//4 1//2;
       -1//2 -mre*(sre5+1)//4 -1//2;
       (1+sre5)//4 mre*(sre5-1)//4 1//2;
       (1+sre5)//4 mre*(sre5-1)//4 -1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 1//2;
       -(1+sre5)//4 mre*(sre5-1)//4 -1//2;
       0 mre 1//2;
       0 mre -1//2;
       (3+sre5+2*sre2*sre5*mre)//8 -(2*mre+sre2*sre5-sre2)//8 0;
       -(3+sre5+2*sre2*sre5*mre)//8 -(2*mre+sre2*sre5-sre2)//8 0]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J53: Biaugmented pentagonal prism")
  return res
end

function _johnson_solid(::Val{54})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 sre3//2 1//2;
       1//2 sre3//2 -1//2;
       1//2 -sre3//2 1//2;
       1//2 -sre3//2 -1//2;
       -1//2 sre3//2 1//2;
       -1//2 sre3//2 -1//2;
       -1//2 -sre3//2 1//2;
       -1//2 -sre3//2 -1//2;
       1 0 1//2;
       1 0 -1//2;
       -1 0 1//2;
       -1 0 -1//2;
       0 (sre2+sre3)//2 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J54: Augmented hexagonal prism")
  return res
end

function _johnson_solid(::Val{55})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 sre3//2 1//2;
       1//2 sre3//2 -1//2;
       1//2 -sre3//2 1//2;
       1//2 -sre3//2 -1//2;
       -1//2 sre3//2 1//2;
       -1//2 sre3//2 -1//2;
       -1//2 -sre3//2 1//2;
       -1//2 -sre3//2 -1//2;
       1 0 1//2;
       1 0 -1//2;
       -1 0 1//2;
       -1 0 -1//2;
       0 (sre2+sre3)//2 0;
       0 -(sre2+sre3)//2 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J55: Parabiaugmented hexagonal prism")
  return res
end

function _johnson_solid(::Val{56})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 sre3//2 1//2;
       1//2 sre3//2 -1//2;
       1//2 -sre3//2 1//2;
       1//2 -sre3//2 -1//2;
       -1//2 sre3//2 1//2;
       -1//2 sre3//2 -1//2;
       -1//2 -sre3//2 1//2;
       -1//2 -sre3//2 -1//2;
       1 0 1//2;
       1 0 -1//2;
       -1 0 1//2;
       -1 0 -1//2;
       (3+sre2*sre3)//4 (sre2+sre3)//4 0;
       -(3+sre2*sre3)//4 (sre2+sre3)//4 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J56: Metabiaugmented hexagonal prism")
  return res
end

function _johnson_solid(::Val{57})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3], ["sqrt(2)", "sqrt(3)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sre2, sre3 = srev
  V = [1//2 sre3//2 1//2;
       1//2 sre3//2 -1//2;
       1//2 -sre3//2 1//2;
       1//2 -sre3//2 -1//2;
       -1//2 sre3//2 1//2;
       -1//2 sre3//2 -1//2;
       -1//2 -sre3//2 1//2;
       -1//2 -sre3//2 -1//2;
       1 0 1//2;
       1 0 -1//2;
       -1 0 1//2;
       -1 0 -1//2;
       (3+sre2*sre3)//4 (sre2+sre3)//4 0;
       -(3+sre2*sre3)//4 (sre2+sre3)//4 0;
       0 -(sre2+sre3)//2 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J57: Triaugmented hexagonal prism")
  return res
end

function _johnson_solid(::Val{58})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(3+sre5)//4 1//2 0;
       (3+sre5)//4 -1//2 0;
       -(3+sre5)//4 1//2 0;
       -(3+sre5)//4 -1//2 0;
       1//2 0 (3+sre5)//4;
       1//2 0 -(3+sre5)//4;
       -1//2 0 (3+sre5)//4;
       -1//2 0 -(3+sre5)//4;
       0 (3+sre5)//4 1//2;
       0 (3+sre5)//4 -1//2;
       0 -(3+sre5)//4 1//2;
       0 -(3+sre5)//4 -1//2;
       (1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       0 (15+sre5)//20 (5+4*sre5)//10]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J58: Augmented dodecahedron")
  return res
end

function _johnson_solid(::Val{59})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(3+sre5)//4 1//2 0;
       (3+sre5)//4 -1//2 0;
       -(3+sre5)//4 1//2 0;
       -(3+sre5)//4 -1//2 0;
       1//2 0 (3+sre5)//4;
       1//2 0 -(3+sre5)//4;
       -1//2 0 (3+sre5)//4;
       -1//2 0 -(3+sre5)//4;
       0 (3+sre5)//4 1//2;
       0 (3+sre5)//4 -1//2;
       0 -(3+sre5)//4 1//2;
       0 -(3+sre5)//4 -1//2;
       (1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       0 (15+sre5)//20 (5+4*sre5)//10;
       0 -(15+sre5)//20 -(5+4*sre5)//10]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J59: Parabiaugmented dodecahedron")
  return res
end

function _johnson_solid(::Val{60})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(3+sre5)//4 1//2 0;
       (3+sre5)//4 -1//2 0;
       -(3+sre5)//4 1//2 0;
       -(3+sre5)//4 -1//2 0;
       1//2 0 (3+sre5)//4;
       1//2 0 -(3+sre5)//4;
       -1//2 0 (3+sre5)//4;
       -1//2 0 -(3+sre5)//4;
       0 (3+sre5)//4 1//2;
       0 (3+sre5)//4 -1//2;
       0 -(3+sre5)//4 1//2;
       0 -(3+sre5)//4 -1//2;
       (1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       0 (15+sre5)//20 (5+4*sre5)//10;
       0 (15+sre5)//20 -(5+4*sre5)//10]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J60: Metabiaugmented dodecahedron")
  return res
end

function _johnson_solid(::Val{61})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(3+sre5)//4 1//2 0;
       (3+sre5)//4 -1//2 0;
       -(3+sre5)//4 1//2 0;
       -(3+sre5)//4 -1//2 0;
       1//2 0 (3+sre5)//4;
       1//2 0 -(3+sre5)//4;
       -1//2 0 (3+sre5)//4;
       -1//2 0 -(3+sre5)//4;
       0 (3+sre5)//4 1//2;
       0 (3+sre5)//4 -1//2;
       0 -(3+sre5)//4 1//2;
       0 -(3+sre5)//4 -1//2;
       (1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 (1+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 (1+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//4 -(1+sre5)//4;
       0 (15+sre5)//20 (5+4*sre5)//10;
       0 (15+sre5)//20 -(5+4*sre5)//10;
       (15+sre5)//20 -(5+4*sre5)//10 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J61: Triaugmented dodecahedron")
  return res
end

function _johnson_solid(::Val{64})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3, x^2 - 5], ["sqrt(2)", "sqrt(3)", "sqrt(5)"])
  ENF, srev = Hecke.embedded_field(NF, real_embeddings(NF)[8])
  sre2, sre3, sre5 = srev
  V = [0 0 (sre3+2*sre2*sre3+sre3*sre5)//6;
       1//2 -sre3//6 (sre3+sre3*sre5)//6;
       -1//2 -sre3//6 (sre3+sre3*sre5)//6;
       0 sre3//3 (sre3+sre3*sre5)//6;
       (1+sre5)//4 -(sre3+sre3*sre5)//12 0;
       -(1+sre5)//4 -(sre3+sre3*sre5)//12 0;
       0 (sre3+sre3*sre5)//6 0;
       1//2 sre3//6 -sre3//3;
       -1//2 sre3//6 -sre3//3;
       0 -sre3//3 -sre3//3]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J64: Augmented tridiminished icosahedron")
  return res
end

function _johnson_solid(::Val{68})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [0 1//2 (5+3*sre5)//4;
       0 1//2 -(5+3*sre5)//4;
       0 -1//2 (5+3*sre5)//4;
       0 -1//2 -(5+3*sre5)//4;
       (5+3*sre5)//4 0 1//2;
       (5+3*sre5)//4 0 -1//2;
       -(5+3*sre5)//4 0 1//2;
       -(5+3*sre5)//4 0 -1//2;
       1//2 (5+3*sre5)//4 0;
       1//2 -(5+3*sre5)//4 0;
       -1//2 (5+3*sre5)//4 0;
       -1//2 -(5+3*sre5)//4 0;
       1//2 (3+sre5)//4 (3+sre5)//2;
       1//2 (3+sre5)//4 -(3+sre5)//2;
       1//2 -(3+sre5)//4 (3+sre5)//2;
       1//2 -(3+sre5)//4 -(3+sre5)//2;
       -1//2 (3+sre5)//4 (3+sre5)//2;
       -1//2 (3+sre5)//4 -(3+sre5)//2;
       -1//2 -(3+sre5)//4 (3+sre5)//2;
       -1//2 -(3+sre5)//4 -(3+sre5)//2;
       (3+sre5)//2 1//2 (3+sre5)//4;
       (3+sre5)//2 1//2 -(3+sre5)//4;
       (3+sre5)//2 -1//2 (3+sre5)//4;
       (3+sre5)//2 -1//2 -(3+sre5)//4;
       -(3+sre5)//2 1//2 (3+sre5)//4;
       -(3+sre5)//2 1//2 -(3+sre5)//4;
       -(3+sre5)//2 -1//2 (3+sre5)//4;
       -(3+sre5)//2 -1//2 -(3+sre5)//4;
       (3+sre5)//4 (3+sre5)//2 1//2;
       (3+sre5)//4 (3+sre5)//2 -1//2;
       (3+sre5)//4 -(3+sre5)//2 1//2;
       (3+sre5)//4 -(3+sre5)//2 -1//2;
       -(3+sre5)//4 (3+sre5)//2 1//2;
       -(3+sre5)//4 (3+sre5)//2 -1//2;
       -(3+sre5)//4 -(3+sre5)//2 1//2;
       -(3+sre5)//4 -(3+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       -1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       (1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       -(1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       0 (10+9*sre5)//10 (15+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J68: Augmented truncated dodecahedron")
  return res
end

function _johnson_solid(::Val{69})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [0 1//2 (5+3*sre5)//4;
       0 1//2 -(5+3*sre5)//4;
       0 -1//2 (5+3*sre5)//4;
       0 -1//2 -(5+3*sre5)//4;
       (5+3*sre5)//4 0 1//2;
       (5+3*sre5)//4 0 -1//2;
       -(5+3*sre5)//4 0 1//2;
       -(5+3*sre5)//4 0 -1//2;
       1//2 (5+3*sre5)//4 0;
       1//2 -(5+3*sre5)//4 0;
       -1//2 (5+3*sre5)//4 0;
       -1//2 -(5+3*sre5)//4 0;
       1//2 (3+sre5)//4 (3+sre5)//2;
       1//2 (3+sre5)//4 -(3+sre5)//2;
       1//2 -(3+sre5)//4 (3+sre5)//2;
       1//2 -(3+sre5)//4 -(3+sre5)//2;
       -1//2 (3+sre5)//4 (3+sre5)//2;
       -1//2 (3+sre5)//4 -(3+sre5)//2;
       -1//2 -(3+sre5)//4 (3+sre5)//2;
       -1//2 -(3+sre5)//4 -(3+sre5)//2;
       (3+sre5)//2 1//2 (3+sre5)//4;
       (3+sre5)//2 1//2 -(3+sre5)//4;
       (3+sre5)//2 -1//2 (3+sre5)//4;
       (3+sre5)//2 -1//2 -(3+sre5)//4;
       -(3+sre5)//2 1//2 (3+sre5)//4;
       -(3+sre5)//2 1//2 -(3+sre5)//4;
       -(3+sre5)//2 -1//2 (3+sre5)//4;
       -(3+sre5)//2 -1//2 -(3+sre5)//4;
       (3+sre5)//4 (3+sre5)//2 1//2;
       (3+sre5)//4 (3+sre5)//2 -1//2;
       (3+sre5)//4 -(3+sre5)//2 1//2;
       (3+sre5)//4 -(3+sre5)//2 -1//2;
       -(3+sre5)//4 (3+sre5)//2 1//2;
       -(3+sre5)//4 (3+sre5)//2 -1//2;
       -(3+sre5)//4 -(3+sre5)//2 1//2;
       -(3+sre5)//4 -(3+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       1//2 -(15+13*sre5)//20 -(15+3*sre5)//10;
       -1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       -1//2 -(15+13*sre5)//20 -(15+3*sre5)//10;
       (1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       (1+sre5)//4 -(25+13*sre5)//20 -(25+sre5)//20;
       -(1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       -(1+sre5)//4 -(25+13*sre5)//20 -(25+sre5)//20;
       0 (10+9*sre5)//10 (15+sre5)//20;
       0 -(10+9*sre5)//10 -(15+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J69: Parabiaugmented truncated dodecahedron")
  return res
end

function _johnson_solid(::Val{70})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [0 1//2 (5+3*sre5)//4;
       0 1//2 -(5+3*sre5)//4;
       0 -1//2 (5+3*sre5)//4;
       0 -1//2 -(5+3*sre5)//4;
       (5+3*sre5)//4 0 1//2;
       (5+3*sre5)//4 0 -1//2;
       -(5+3*sre5)//4 0 1//2;
       -(5+3*sre5)//4 0 -1//2;
       1//2 (5+3*sre5)//4 0;
       1//2 -(5+3*sre5)//4 0;
       -1//2 (5+3*sre5)//4 0;
       -1//2 -(5+3*sre5)//4 0;
       1//2 (3+sre5)//4 (3+sre5)//2;
       1//2 (3+sre5)//4 -(3+sre5)//2;
       1//2 -(3+sre5)//4 (3+sre5)//2;
       1//2 -(3+sre5)//4 -(3+sre5)//2;
       -1//2 (3+sre5)//4 (3+sre5)//2;
       -1//2 (3+sre5)//4 -(3+sre5)//2;
       -1//2 -(3+sre5)//4 (3+sre5)//2;
       -1//2 -(3+sre5)//4 -(3+sre5)//2;
       (3+sre5)//2 1//2 (3+sre5)//4;
       (3+sre5)//2 1//2 -(3+sre5)//4;
       (3+sre5)//2 -1//2 (3+sre5)//4;
       (3+sre5)//2 -1//2 -(3+sre5)//4;
       -(3+sre5)//2 1//2 (3+sre5)//4;
       -(3+sre5)//2 1//2 -(3+sre5)//4;
       -(3+sre5)//2 -1//2 (3+sre5)//4;
       -(3+sre5)//2 -1//2 -(3+sre5)//4;
       (3+sre5)//4 (3+sre5)//2 1//2;
       (3+sre5)//4 (3+sre5)//2 -1//2;
       (3+sre5)//4 -(3+sre5)//2 1//2;
       (3+sre5)//4 -(3+sre5)//2 -1//2;
       -(3+sre5)//4 (3+sre5)//2 1//2;
       -(3+sre5)//4 (3+sre5)//2 -1//2;
       -(3+sre5)//4 -(3+sre5)//2 1//2;
       -(3+sre5)//4 -(3+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       1//2 -(15+13*sre5)//20 (15+3*sre5)//10;
       -1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       -1//2 -(15+13*sre5)//20 (15+3*sre5)//10;
       (1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       (1+sre5)//4 -(25+13*sre5)//20 (25+sre5)//20;
       -(1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       -(1+sre5)//4 -(25+13*sre5)//20 (25+sre5)//20;
       0 (10+9*sre5)//10 (15+sre5)//20;
       0 -(10+9*sre5)//10 (15+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J70: Metabiaugmented truncated dodecahedron")
  return res
end

function _johnson_solid(::Val{71})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [0 1//2 (5+3*sre5)//4;
       0 1//2 -(5+3*sre5)//4;
       0 -1//2 (5+3*sre5)//4;
       0 -1//2 -(5+3*sre5)//4;
       (5+3*sre5)//4 0 1//2;
       (5+3*sre5)//4 0 -1//2;
       -(5+3*sre5)//4 0 1//2;
       -(5+3*sre5)//4 0 -1//2;
       1//2 (5+3*sre5)//4 0;
       1//2 -(5+3*sre5)//4 0;
       -1//2 (5+3*sre5)//4 0;
       -1//2 -(5+3*sre5)//4 0;
       1//2 (3+sre5)//4 (3+sre5)//2;
       1//2 (3+sre5)//4 -(3+sre5)//2;
       1//2 -(3+sre5)//4 (3+sre5)//2;
       1//2 -(3+sre5)//4 -(3+sre5)//2;
       -1//2 (3+sre5)//4 (3+sre5)//2;
       -1//2 (3+sre5)//4 -(3+sre5)//2;
       -1//2 -(3+sre5)//4 (3+sre5)//2;
       -1//2 -(3+sre5)//4 -(3+sre5)//2;
       (3+sre5)//2 1//2 (3+sre5)//4;
       (3+sre5)//2 1//2 -(3+sre5)//4;
       (3+sre5)//2 -1//2 (3+sre5)//4;
       (3+sre5)//2 -1//2 -(3+sre5)//4;
       -(3+sre5)//2 1//2 (3+sre5)//4;
       -(3+sre5)//2 1//2 -(3+sre5)//4;
       -(3+sre5)//2 -1//2 (3+sre5)//4;
       -(3+sre5)//2 -1//2 -(3+sre5)//4;
       (3+sre5)//4 (3+sre5)//2 1//2;
       (3+sre5)//4 (3+sre5)//2 -1//2;
       (3+sre5)//4 -(3+sre5)//2 1//2;
       (3+sre5)//4 -(3+sre5)//2 -1//2;
       -(3+sre5)//4 (3+sre5)//2 1//2;
       -(3+sre5)//4 (3+sre5)//2 -1//2;
       -(3+sre5)//4 -(3+sre5)//2 1//2;
       -(3+sre5)//4 -(3+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       (3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 (1+sre5)//2 -(2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 (2+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//2 -(2+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       (2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 (3+sre5)//4 -(1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 (1+sre5)//2;
       -(2+sre5)//2 -(3+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       (1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 (2+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 (3+sre5)//4;
       -(1+sre5)//2 -(2+sre5)//2 -(3+sre5)//4;
       1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       1//2 -(15+13*sre5)//20 (15+3*sre5)//10;
       -1//2 (15+13*sre5)//20 (15+3*sre5)//10;
       -1//2 -(15+13*sre5)//20 (15+3*sre5)//10;
       (1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       (1+sre5)//4 -(25+13*sre5)//20 (25+sre5)//20;
       -(1+sre5)//4 (25+13*sre5)//20 (25+sre5)//20;
       -(1+sre5)//4 -(25+13*sre5)//20 (25+sre5)//20;
       0 (10+9*sre5)//10 (15+sre5)//20;
       0 -(10+9*sre5)//10 (15+sre5)//20;
       -(15+3*sre5)//10 1//2 -(15+13*sre5)//20;
       -(15+3*sre5)//10 -1//2 -(15+13*sre5)//20;
       -(25+sre5)//20 (1+sre5)//4 -(25+13*sre5)//20;
       -(25+sre5)//20 -(1+sre5)//4 -(25+13*sre5)//20;
       -(15+sre5)//20 0 -(10+9*sre5)//10]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J71: Triaugmented truncated dodecahedron")
  return res
end

function _johnson_solid(::Val{72})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 -(3+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       0 -(3+sre5)//4 -(5+sre5)//4;
       0 -(3+sre5)//4 (5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 1//2 -(2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       -1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 (2+sre5)//2 -1//2;
       1//2 -(2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       -1//2 -(2+sre5)//2 -1//2;
       1//2 -(2+sre5)//2 1//2;
       -1//2 -(2+sre5)//2 1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 (3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 (3+sre5)//4;
       1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       (1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       0 (15+13*sre5)//20 (5+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J72: Gyrate rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{73})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 -(3+sre5)//4;
       0 -(3+sre5)//4 (5+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 1//2 -(2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       -1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 -(2+sre5)//2 1//2;
       -1//2 -(2+sre5)//2 1//2;
       1//2 (2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 (3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 (3+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       1//2 -(5+4*sre5)//10 -(10+3*sre5)//10;
       -1//2 -(5+4*sre5)//10 -(10+3*sre5)//10;
       (1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       (1+sre5)//4 -(5+2*sre5)//5 -(15+sre5)//20;
       -(1+sre5)//4 -(5+2*sre5)//5 -(15+sre5)//20;
       0 (15+13*sre5)//20 (5+sre5)//20;
       0 -(15+13*sre5)//20 -(5+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J73: Parabigyrate rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{74})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 -(3+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       0 -(3+sre5)//4 -(5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 1//2 -(2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       -1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 (2+sre5)//2 -1//2;
       1//2 -(2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       -1//2 -(2+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       (1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       (1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       0 (15+13*sre5)//20 (5+sre5)//20;
       0 -(15+13*sre5)//20 (5+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J74: Metabigyrate rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{75})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       0 -(3+sre5)//4 -(5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 (2+sre5)//2 -1//2;
       1//2 -(2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       -1//2 -(2+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       (1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       (1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       0 (15+13*sre5)//20 (5+sre5)//20;
       0 -(15+13*sre5)//20 (5+sre5)//20;
       -(10+3*sre5)//10 1//2 -(5+4*sre5)//10;
       -(10+3*sre5)//10 -1//2 -(5+4*sre5)//10;
       -(15+sre5)//20 (1+sre5)//4 -(5+2*sre5)//5;
       -(15+sre5)//20 -(1+sre5)//4 -(5+2*sre5)//5;
       -(5+sre5)//20 0 -(15+13*sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J75: Trigyrate rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{77})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 -(3+sre5)//4;
       0 -(3+sre5)//4 (5+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 1//2 -(2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       -1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 -(2+sre5)//2 1//2;
       -1//2 -(2+sre5)//2 1//2;
       1//2 (2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 (3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 (3+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       1//2 -(5+4*sre5)//10 -(10+3*sre5)//10;
       -1//2 -(5+4*sre5)//10 -(10+3*sre5)//10;
       (1+sre5)//4 -(5+2*sre5)//5 -(15+sre5)//20;
       -(1+sre5)//4 -(5+2*sre5)//5 -(15+sre5)//20;
       0 -(15+13*sre5)//20 -(5+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J77: Paragyrate diminished rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{78})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 -(3+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       0 -(3+sre5)//4 -(5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 1//2 -(2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       -1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 (2+sre5)//2 -1//2;
       1//2 -(2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       -1//2 -(2+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       (1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       0 -(15+13*sre5)//20 (5+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J78: Metagyrate diminished rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{79})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       0 -(3+sre5)//4 -(5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 (2+sre5)//2 -1//2;
       1//2 -(2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       -1//2 -(2+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 (5+4*sre5)//10 (10+3*sre5)//10;
       -1//2 -(5+4*sre5)//10 (10+3*sre5)//10;
       (1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       (1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 (5+2*sre5)//5 (15+sre5)//20;
       -(1+sre5)//4 -(5+2*sre5)//5 (15+sre5)//20;
       0 (15+13*sre5)//20 (5+sre5)//20;
       0 -(15+13*sre5)//20 (5+sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J79: Bigyrate diminished rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{82})
  NF, sr5 = quadratic_field(5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V = [(5+sre5)//4 0 (3+sre5)//4;
       -(5+sre5)//4 0 (3+sre5)//4;
       (5+sre5)//4 0 -(3+sre5)//4;
       0 (3+sre5)//4 -(5+sre5)//4;
       0 -(3+sre5)//4 -(5+sre5)//4;
       (3+sre5)//4 (5+sre5)//4 0;
       (3+sre5)//4 -(5+sre5)//4 0;
       -(3+sre5)//4 (5+sre5)//4 0;
       -(3+sre5)//4 -(5+sre5)//4 0;
       1//2 1//2 (2+sre5)//2;
       1//2 -1//2 (2+sre5)//2;
       -1//2 1//2 (2+sre5)//2;
       -1//2 -1//2 (2+sre5)//2;
       1//2 1//2 -(2+sre5)//2;
       1//2 -1//2 -(2+sre5)//2;
       (2+sre5)//2 1//2 1//2;
       (2+sre5)//2 1//2 -1//2;
       (2+sre5)//2 -1//2 1//2;
       (2+sre5)//2 -1//2 -1//2;
       -(2+sre5)//2 1//2 1//2;
       -(2+sre5)//2 1//2 -1//2;
       -(2+sre5)//2 -1//2 1//2;
       -(2+sre5)//2 -1//2 -1//2;
       1//2 (2+sre5)//2 -1//2;
       1//2 -(2+sre5)//2 -1//2;
       -1//2 (2+sre5)//2 -1//2;
       -1//2 -(2+sre5)//2 -1//2;
       (3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 (1+sre5)//4 (1+sre5)//2;
       -(3+sre5)//4 -(1+sre5)//4 (1+sre5)//2;
       (3+sre5)//4 (1+sre5)//4 -(1+sre5)//2;
       (3+sre5)//4 -(1+sre5)//4 -(1+sre5)//2;
       (1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       (1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 (3+sre5)//4 -(1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 (1+sre5)//4;
       -(1+sre5)//2 -(3+sre5)//4 -(1+sre5)//4;
       (1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       (1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 (1+sre5)//2 -(3+sre5)//4;
       -(1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4;
       -(10+3*sre5)//10 1//2 -(5+4*sre5)//10;
       -(10+3*sre5)//10 -1//2 -(5+4*sre5)//10;
       -(15+sre5)//20 (1+sre5)//4 -(5+2*sre5)//5;
       -(15+sre5)//20 -(1+sre5)//4 -(5+2*sre5)//5;
       -(5+sre5)//20 0 -(15+13*sre5)//20]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J82: Gyrate bidiminished rhombicosidodecahedron")
  return res
end

function _johnson_solid(::Val{84})
  Qx, x = QQ["x"]
  NF, q = number_field(2*x^3 + 11*x^2 + 4*x - 1, "a")
  NFy, y = NF["y"]
  MF, srq = number_field(y^2 - q, "sqrt(a)")
  EMF, sreq = Hecke.embedded_field(MF, real_embeddings(MF)[2])
  eq = EMF(q)
  r = sreq//2
  s = (eq^2//4 + 13*eq//8 + 13//8)*sreq
  t = eq^2//2 + 9*eq//4 + 1//4
  V = [t 0 -r;
       -t 0 -r;
       0 t r;
       0 -t r;
       1//2 0 s;
       -1//2 0 s;
       0 1//2 -s;
       0 -1//2 -s]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J84: Snub disphenoid")
  return res
end

function _johnson_solid(::Val{85})
  Qx, x = QQ["x"]
  NF, a = number_field(x^6 - 2*x^5 - 13*x^4 + 8*x^3 + 32*x^2 - 8*x - 4, "a")
  sr2 = -a^5//12 + 19*a^3//12 + a^2//2 - 25*a//6 + 1//3
  NFy, y = NF["y"]
  MF, b = number_field(y^2 - (1 - (1 - sr2//2)*a^2), "b")
  EMF, eb = Hecke.embedded_field(MF, real_embeddings(MF)[6])
  ea = EMF(a)
  ec = (7*ea^5//12 - 9*ea^4//4 - 43*ea^3//12 + 45*ea^2//4 + a//6 - 5//6)*eb
  sre2 = EMF(sr2)
  V = [1//2 1//2 ec//2;
       1//2 -1//2 ec//2;
       -1//2 1//2 ec//2;
       -1//2 -1//2 ec//2;
       sre2*ea//2 0 eb//2;
       -sre2*ea//2 0 eb//2;
       0 sre2*ea//2 eb//2;
       0 -sre2*ea//2 eb//2;
       ea//2 ea//2 -eb//2;
       ea//2 -ea//2 -eb//2;
       -ea//2 ea//2 -eb//2;
       -ea//2 -ea//2 -eb//2;
       0 sre2//2 -ec//2;
       0 -sre2//2 -ec//2;
       sre2//2 0 -ec//2;
       -sre2//2 0 -ec//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J85: Snub square antiprism")
  return res
end

function _johnson_solid(::Val{86})
  Qx, x = QQ["x"]
  NF, k = number_field(60*x^4 - 48*x^3 - 100*x^2 + 56*x + 23, "a")
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - 1 + k^2, "b")
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  ke = EMF(k)
  V = [0 1//2 mre;
       0 -1//2 mre;
       ke 1//2 0;
       ke -1//2 0;
       -ke 1//2 0;
       -ke -1//2 0;
       0 (2*ke^3//27+32*ke^2//135+334*ke//405-107//810) (1-2*ke^2)//(2*mre);
       0 -(2*ke^3//27+32*ke^2//135+334*ke//405-107//810) (1-2*ke^2)//(2*mre);
       1//2 0 (13*ke^3//27+208*ke^2//135-1733*ke//810-449//405)*mre;
       -1//2 0 (13*ke^3//27+208*ke^2//135-1733*ke//810-449//405)*mre]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J86: Sphenocorona")
  return res
end

function _johnson_solid(::Val{87})
  Qx, x = QQ["x"]
  NF, k = number_field(60*x^4 - 48*x^3 - 100*x^2 + 56*x + 23, "a")
  NFy, y = NF["y"]
  MF, mr = number_field([y^2 - 1 + k^2, y^2 - 2], [:b, Symbol("sqrt(2)")])
  EMF, mra = Hecke.embedded_field(MF, real_embeddings(MF)[8])
  ke = EMF(k)
  mre, sre2 = mra
  V = [0 1//2 mre;
       0 -1//2 mre;
       ke 1//2 0;
       ke -1//2 0;
       -ke 1//2 0;
       -ke -1//2 0;
       0 (2*ke^3//27+32*ke^2//135+334*ke//405-107//810) (1-2*ke^2)//(2*mre);
       0 -(2*ke^3//27+32*ke^2//135+334*ke//405-107//810) (1-2*ke^2)//(2*mre);
       1//2 0 (13*ke^3//27+208*ke^2//135-1733*ke//810-449//405)*mre;
       -1//2 0 (13*ke^3//27+208*ke^2//135-1733*ke//810-449//405)*mre;
       (ke+sre2*mre)//2 0 (ke*sre2+mre)//2]
  res = convex_hull(EMF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J87: Augmented sphenocorona")
  return res
end

function _johnson_solid(::Val{88})
  Qx, x = QQ["x"]
  NF, ks = number_field(1680*x^16- 4800*x^15 - 3712*x^14 + 17216*x^13+ 1568*x^12 - 24576*x^11 + 2464*x^10 + 17248*x^9- 3384*x^8 - 5584*x^7 + 2000*x^6+ 240*x^5- 776*x^4+ 304*x^3 + 200*x^2 - 56*x - 23, "a")
  NFy, y = NF["y"]
  UNF, us = number_field(y^2 - 1 + ks^2, "b")
  ENF, u = Hecke.embedded_field(UNF, real_embeddings(UNF)[4])
  k = ENF(ks) 
  v = (-779398396//3645*k^15 + 1527897016//3645*k^14 + 47481585556//54675*k^13 - 80358097444//54675*k^12 - 626829326//405*k^11 + 34715953556//18225*k^10 + 15225336298//10935*k^9 - 61415739374//54675*k^8 - 30047329289//54675*k^7 + 16348591376//54675*k^6 - 128888737//18225*k^5 - 282199511//6075*k^4 + 7240761619//109350*k^3 + 925235176//54675*k^2 - 1435429847//109350*k - 404623837//109350)*u
  w = (-4937536//3645*k^15 + 66072592//25515*k^14 + 1900033672//382725*k^13 - 429863104//54675*k^12 - 3078356//405*k^11 + 1002736772//127575*k^10 + 408983716//76545*k^9 - 1127158988//382725*k^8 - 481084868//382725*k^7 + 137102912//382725*k^6 - 37465594//127575*k^5 - 5910032//42525*k^4 + 59708339//382725*k^3 + 3939991//54675*k^2 - 4189357//382725*k - 2877197//382725)*u
  V = [0 1//2 u; 
      0 -1//2 u; 
      k 1//2 0;
      k -1//2 0;
      -k 1//2 0;
      -k -1//2 0;
      0 w//(2*u)+1//2 (1-2*k^2)//(2*u);
      0 -w//(2*u)-1//2 (1-2*k^2)//(2*u);
      1//2 0 -v//2;
      -1//2 0 -v//2;
      0 w*(2*k^2-1)//(2*(k^2-1)*u)+1//2 (2*k^4-1)//(2*u^3);
      0 -w*(2*k^2-1)//(2*(k^2-1)*u)-1//2 (2*k^4-1)//(2*u^3)]
  res = convex_hull(ENF, V; non_redundant = false)
  Polymake.setdescription!(pm_object(res), "Johnson solid J88: Sphenomegacorona")
  return res
end

function _johnson_solid(::Val{89})
  Qx, x = QQ["x"]
  NF, ks = number_field(26880*x^10 + 35328*x^9 - 25600*x^8 - 39680*x^7 + 6112x^6 + 13696*x^5 + 2128*x^4 - 1808*x^3 - 1119*x^2 + 494*x - 47, "a")
  NFy, y = NF["y"]
  UNF, us = number_field(y^2 - 1 + ks^2, "b")
  ENF, u = Hecke.embedded_field(UNF, real_embeddings(UNF)[6])
  k = ENF(ks) 
  v = (-174265//72*k^9 - 263299//72*k^8 + 354083//216*k^7 + 213419//54*k^6 + 333139//1728*k^5 - 2134337//1728*k^4 - 1505857//3456*k^3 + 147085//1728*k^2 + 6787789//55296*k - 1258781//55296)*u
  w = (-340207525//80208*k^9 - 524173255//80208*k^8 + 621900263//240624*k^7 + 411287453//60156*k^6 + 1019406967//1924992*k^5 - 3914171453//1924992*k^4 - 2949983917//3849984*k^3 + 217946413//1924992*k^2 + 12238010281//61599744*k - 2178441785//61599744)*u
  V = [1//2 1//2 u;
      1//2 -1//2 u;
      -1//2 1//2 u;
      -1//2 -1//2 u;
      1//2+k 1//2 0;
      1//2+k -1//2 0;
      -1//2-k 1//2 0;
      -1//2-k -1//2 0;
      1//2 0 -w//2;
      -1//2 0 -w//2;
      0 1//2+v//(2*u) v^2//(4*u);
      0 -1//2-v//(2*u) v^2//(4*u);
      0 (v*w+k+1)//(4*u^2) (2*k-1)*w//(4*(1-k))-v//(4*u^2);
      0 -(v*w+k+1)//(4*u^2) (2*k-1)*w//(4*(1-k))-v//(4*u^2)]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J89: Hebesphenomegacorona")
  return res
end

function _johnson_solid(::Val{90})
  Qx, x = QQ["x"]
  NF, ks = number_field(256*x^12 - 512*x^11 - 1664*x^10 + 3712*x^9 + 1552*x^8 - 6592*x^7 + 1248*x^6 + 4352*x^5 - 2024*x^4 - 944*x^3 + 672*x^2 - 24*x - 23, "a")
  NFy, y = NF["y"]
  UNF, us = number_field(y^2 - 1 + ks^2, "b")
  ENF, u = Hecke.embedded_field(UNF, real_embeddings(UNF)[6])
  k = ENF(ks) 
  v = (-12736//405*k^11 + 17728//405*k^10 + 32512//135*k^9 - 131488//405*k^8 - 184124//405*k^7 + 259436//405*k^6 + 125426//405*k^5 - 23542//45*k^4 - 13363//405*k^3 + 66277//405*k^2 - 17923//810*k - 5449//810)*u
  w = (-22592//405*k^11 + 28256//405*k^10 + 55664//135*k^9 - 200696//405*k^8 - 281488//405*k^7 + 360832//405*k^6 + 157612//405*k^5 - 28904//45*k^4 - 15476//405*k^3 + 72974//405*k^2 - 6898//405*k - 7073//810)*u
  V =  [1//2 0 u+v//4; 
        -1//2 0 u+v//4; 
        0 1//2 -u-v//4; 
        0 -1//2 -u-v//4;
        1//2+w//(2*u) 0 u-1//(2*u)+v//4;
        -1//2-w//(2*u) 0 u-1//(2*u)+v//4;
        0 1//2+w//(2*u) -u+1//(2*u)-v//4;
        0 -1//2-w//(2*u) -u+1//(2*u)-v//4;
        1//2 k v//4; 
        1//2 -k v//4;
        -1//2 k v//4; 
        -1//2 -k v//4;
        k 1//2 -v//4; 
        -k 1//2 -v//4;
        k -1//2 -v//4; 
        -k -1//2 -v//4]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J90: Disphenocingulum")
  return res
end

function _johnson_solid(::Val{92})
  Qx, x = QQ["x"]
  NF, srs = number_field([x^2 - 3, x^2 - 5], ["sqrt(3)", "sqrt(5)"])
  ENF, srse = Hecke.embedded_field(NF, real_embeddings(NF)[4])
  sr3, sr5 = srse
  V = [1//2 -sr3//6 sr3*(3+sr5)//6;
       -1//2 -sr3//6 sr3*(3+sr5)//6;
       0 sr3//3 sr3*(3+sr5)//6;
       1//2 sr3*(2+sr5)//6 sr3*(1+sr5)//6;
       -1//2 sr3*(2+sr5)//6 sr3*(1+sr5)//6;
       (3+sr5)//4 -sr3*(sr5-1)//12 sr3*(1+sr5)//6;
       -(3+sr5)//4 -sr3*(sr5-1)//12 sr3*(1+sr5)//6;
       (1+sr5)//4 -sr3*(5+sr5)//12 sr3*(1+sr5)//6;
       -(1+sr5)//4 -sr3*(5+sr5)//12 sr3*(1+sr5)//6;
       (3+sr5)//4 sr3*(3+sr5)//12 sr3//3;
       -(3+sr5)//4 sr3*(3+sr5)//12 sr3//3;
       0 -sr3*(3+sr5)//6 sr3//3; 
       1//2 sr3//2 0;
       1//2 -sr3//2 0;
       -1//2 sr3//2 0;
       -1//2 -sr3//2 0;
       1 0 0; 
       -1 0 0]
  res = convex_hull(ENF, V; non_redundant = true)
  Polymake.setdescription!(pm_object(res), "Johnson solid J92: Triangular hebesphenorotunda")
  return res
end

