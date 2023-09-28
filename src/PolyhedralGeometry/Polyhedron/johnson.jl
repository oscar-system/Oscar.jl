function _johnson_solid(::Val{9})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{10})
  Qx, x = QQ["x"]
  NF, qr8 = number_field(x^4 - 8)
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
  return convex_hull(ENF, V; non_redundant = true)
end

function _johnson_solid(::Val{13})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
  EMF, mre = Hecke.embedded_field(MF, real_embeddings(MF)[4])
  sre5 = EMF(sr5)
  V = [1//2 -mre*(sre5+1)//4 0;
       -1//2 -mre*(sre5+1)//4 0;
       (1+sre5)//4 mre*(sre5-1)//4 0;
       -(1+sre5)//4 mre*(sre5-1)//4 0;
       0 mre 0;
       0 0 mre*(sre5-1)//2;
       0 0 -mre*(sre5-1)//2]
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{16})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{17})
  Qx, x = QQ["x"]
  NF, qr8 = number_field(x^4 - 8)
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
  return convex_hull(ENF, V; non_redundant = true)
end

function _johnson_solid(::Val{18})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3])
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
  return convex_hull(ENF, V; non_redundant = true)
end

function _johnson_solid(::Val{20})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{21})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{22})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3])
  sr2, sr3 = srv
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (sr3-1))
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{23})
  Qx, x = QQ["x"]
  NF, sr2 = number_field(x^2 - 2)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (sr2+2))
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-2 - 2*sr2 + (sr2+2)*mr)//8)
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
  return convex_hull(ELF, V; non_redundant = true)
end

function _johnson_solid(::Val{24})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-4 - 2*sr5 + (3*sr5+5)*mr)//8)
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
  return convex_hull(ELF, V; non_redundant = true)
end

function _johnson_solid(::Val{25})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
  MFz, z = MF["z"]
  LF, lr = number_field(z^2 - (-4 - 2*sr5 + (3*sr5+5)*mr)//8)
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
  return convex_hull(ELF, V; non_redundant = true)
end

function _johnson_solid(::Val{30})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{32})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{33})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{34})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{35})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3])
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
  return convex_hull(ENF, V; non_redundant = true)
end

function _johnson_solid(::Val{36})
  Qx, x = QQ["x"]
  NF, srv = number_field([x^2 - 2, x^2 - 3])
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
  return convex_hull(ENF, V; non_redundant = true)
end

function _johnson_solid(::Val{38})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function _johnson_solid(::Val{39})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  NFy, y = NF["y"]
  MF, mr = number_field(y^2 - (5 + sr5)//10)
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
  return convex_hull(EMF, V; non_redundant = true)
end

function johnson_solid(::Val{83})
  Qx, x = QQ["x"]
  NF, sr5 = number_field(x^2 - 5)
  ENF, sre5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V =  [(5+sre5)//4 0 (3+sre5)//4; 
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
        -(1+sre5)//4 -(1+sre5)//2 -(3+sre5)//4]
end

function _johnson_solid(::Val{91})
  Qx, x = QQ["x"]
  NF, sr = number_field(x^2 - 5)
  ENF, sr5 = Hecke.embedded_field(NF, real_embeddings(NF)[2])
  V =  [1//2 0 (3+sr5)//4;
        1//2 0 -(3+sr5)//4;
        -1//2 0 (3+sr5)//4;
        -1//2 0 -(3+sr5)//4;
        (1+sr5)//4 1//2 1//2;
        (1+sr5)//4 1//2 -1//2;
        (1+sr5)//4 -1//2 1//2;
        (1+sr5)//4 -1//2 -1//2;
        -(1+sr5)//4 1//2 1//2;
        -(1+sr5)//4 1//2 -1//2;
        -(1+sr5)//4 -1//2 1//2;
        -(1+sr5)//4 -1//2 -1//2;
        0 (1+sr5)//4 0;
        0 -(1+sr5)//4 0]
  return convex_hull(ENF, V; non_redundant = true)
end

function _johnson_solid(::Val{92})
  Qx, x = QQ["x"]
  NF, srs = number_field([x^2 - 3, x^2 - 5])
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
  return convex_hull(ENF, V; non_redundant = true)
end

