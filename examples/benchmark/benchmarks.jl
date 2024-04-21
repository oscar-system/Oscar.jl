using Oscar
using GroebnerWalk
using BenchmarkTools

function benchmark(
  io,
  name::String,
  I::Ideal,
  target::MonomialOrdering,
  start::MonomialOrdering
)
  print(io, name, ",")
  t = @belapsed groebner_walk($I, $target, $start; algorithm=:standard)
  print(io, t, ",")
  t = @belapsed groebner_walk($I, $target, $start; algorithm=:generic)
  print(io, t, ",")
  t = @belapsed groebner_basis($I; ordering=$target)
  println(io, t)
end

p = 11863279
k = GF(p)
open("results.csv", "a") do io
  R, (x, y) = QQ[:x, :y]
  I = ideal([y^4 + x^3 - x^2 + x, x^4])
  benchmark(io, "simple", I, lex(R), default_ordering(R))

  R, (z0, z1, z2, z3, z4, z5, z6) = QQ[:z0, :z1, :z2, :z3, :z4, :z5, :z6]
  F = [
    z0 + z1 + z2 + z3 + z4 + z5 + z6,
    z0 * z1 + z1 * z2 + z2 * z3 + z3 * z4 + z4 * z5 + z5 * z6 + z6 * z0,
    z0 * z1 * z2 + z1 * z2 * z3 + z2 * z3 * z4 + z3 * z4 * z5 + z4 * z5 * z6 + z5 * z6 * z0 + z6 * z0 * z1,
    z0 * z1 * z2 * z3 + z1 * z2 * z3 * z4 + z2 * z3 * z4 * z5 + z3 * z4 * z5 * z6 + z4 * z5 * z6 * z0
    + z5 * z6 * z0 * z1 + z6 * z0 * z1 * z2,
    z0 * z1 * z2 * z3 * z4 + z1 * z2 * z3 * z4 * z5 + z2 * z3 * z4 * z5 * z6 + z3 * z4 * z5 * z6 * z0
    + z4 * z5 * z6 * z0 * z1 + z5 * z6 * z0 * z1 * z2 + z6 * z0 * z1 * z2 * z3,
    z0 * z1 * z2 * z3 * z4 * z5 + z1 * z2 * z3 * z4 * z5 * z6 + z2 * z3 * z4 * z5 * z6 * z0 + z3 * z4 * z5 * z6 * z0 * z1
    + z4 * z5 * z6 * z0 * z1 * z2 + z5 * z6 * z0 * z1 * z2 * z3 + z6 * z0 * z1 * z2 * z3 * z4,
    z0 * z1 * z2 * z3 * z4 * z5 * z6 - 1
  ]
  benchmark(io, "cyclic7-QQ", ideal(R, F), lex(R), default_ordering(R))

  R, (z0, z1, z2, z3, z4, z5, z6) = k[:z0, :z1, :z2, :z3, :z4, :z5, :z6]
  F = [
    z0 + z1 + z2 + z3 + z4 + z5 + z6,
    z0 * z1 + z1 * z2 + z2 * z3 + z3 * z4 + z4 * z5 + z5 * z6 + z6 * z0,
    z0 * z1 * z2 + z1 * z2 * z3 + z2 * z3 * z4 + z3 * z4 * z5 + z4 * z5 * z6 + z5 * z6 * z0 + z6 * z0 * z1,
    z0 * z1 * z2 * z3 + z1 * z2 * z3 * z4 + z2 * z3 * z4 * z5 + z3 * z4 * z5 * z6 + z4 * z5 * z6 * z0
    + z5 * z6 * z0 * z1 + z6 * z0 * z1 * z2,
    z0 * z1 * z2 * z3 * z4 + z1 * z2 * z3 * z4 * z5 + z2 * z3 * z4 * z5 * z6 + z3 * z4 * z5 * z6 * z0
    + z4 * z5 * z6 * z0 * z1 + z5 * z6 * z0 * z1 * z2 + z6 * z0 * z1 * z2 * z3,
    z0 * z1 * z2 * z3 * z4 * z5 + z1 * z2 * z3 * z4 * z5 * z6 + z2 * z3 * z4 * z5 * z6 * z0 + z3 * z4 * z5 * z6 * z0 * z1
    + z4 * z5 * z6 * z0 * z1 * z2 + z5 * z6 * z0 * z1 * z2 * z3 + z6 * z0 * z1 * z2 * z3 * z4,
    z0 * z1 * z2 * z3 * z4 * z5 * z6 - 1
  ]
  benchmark(io, "cyclic7-Fp", ideal(R, F), lex(R), default_ordering(R))

  # k = GF(11863279)
  # R, (x1, x2, x3, x4, x5, x6, x7) = k[:x1, :x2, :x3, :x4, :x5, :x6, :x7]
  # F = [
  #     1*x1+2*x2+2*x3+2*x4+2*x5+2*x6+2*x7-1,
  #     2*x4*x3+2*x5*x2+2*x6*x1+2*x7*x2-1*x6,
  #     1*x3^2+2*x4*x2+2*x5*x1+2*x6*x2+2*x7*x3-1*x5,
  #     2*x3*x2+2*x4*x1+2*x5*x2+2*x6*x3+2*x7*x4-1*x4,
  #     1*x2^2+2*x3*x1+2*x4*x2+2*x5*x3+2*x6*x4+2*x7*x5-1*x3,
  #     2*x2*x1+2*x3*x2+2*x4*x3+2*x5*x4+2*x6*x5+2*x7*x6-1*x2,
  #     1*x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2+2*x6^2+2*x7^2-1*x1
  # ];
  # benchmark(io, "katsura6-Fp", ideal(R, F), lex(R), default_ordering(R))

  # R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32) = QQ[
  #   :x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14, :x15, :x16, :x17, :x18, :x19, :x20, :x21, :x22, :x23, :x24, :x25, :x26, :x27, :x28, :x29, :x30, :x31, :x32
  # ]
  # F = [
  #   -x23*x32+x24*x31, -x22*x32+x24*x30, -x22*x31+x23*x30, -x21*x32+x24*x29, -x21*x31+x23*x29, -x21*x30+x22*x29, -x12*x32+x16*x28,  -x19*x28+x20*x27,
  #   -x11*x31+x15*x27, -x18*x28+x20*x26, -x18*x27+x19*x26, -x10*x30+x14*x26, -x17*x28+x20*x25, -x17*x27+x19*x25, -x17*x26+x18*x25, -x9*x29+x13*x25, x20*x8-x24*x4,
  #   -x17*x20-x17*x24-2*x17*x28-x17*x32+x18*x19+x18*x23+2*x18*x27+x18*x31+x19*x22+x19*x30-x20*x21-x20*x29-x21*x24-x21*x28-2*x21*x32+x22*x23+x22*x27+2*x22*x31+x23*x26-x24*x25-x25*x28-x25*x32+x26*x27+x26*x31+x27*x30-x28*x29-x29*x32+x30*x31,
  #   x19*x7-x23*x3, x18*x6-x2*x22, -x1*x21+x17*x5
  # ]
  # benchmark(io, "bayes148", ideal(R, F), lex(R), default_ordering(R))
end