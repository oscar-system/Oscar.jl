module benchmark_matrices

using BenchmarkTools
using Oscar

const reps = 10

# benchmark: matrix addition
function test_add(a)
   s = a
   for i = 1:reps
    s += a
   end
   return s
end


# benchmark: matrix multiplication
function test_mul(a)
   s = a
   for i = 1:reps
    s *= a
   end
   return s
end

# TODO: test_inv

# benchmark: min poly
function test_minpoly(a)
   s = minpoly(a)
   for i = 1:reps
     s += minpoly(a)
   end
   return s
end


# benchmark: char poly
function test_charpoly(a)
   s = charpoly(a)
   for i = 1:reps
     s += charpoly(a)
   end
   return s
end


# benchmark: order computation
# TODO: not even implemented for most of our matrix types
function test_order(a)
   s = order(a)
   for i = 1:reps
     s += order(a)
   end
   return s
end

# HACK HACK HACK
function Base.:+(x::Oscar.GAPGroupElem, y::Oscar.GAPGroupElem)
  return Oscar.group_element(parent(x), x.X+y.X)
end



# test data

# TODO: do something better than using identity matrix!!
# e.g. a "random" dense matrix; ideally invertible so that we
# can we can try inversion, too

n = 10  # 10 x 10 matrices: TODO: test with a range of different degrees

testdata = []

for (p,d) in [(2,1), (17,1), (17,2)]
  q = p^d
  F, _ = FiniteField(fmpz(p), d, "a")
  mats = []
  push!(mats, (one(GL(n,q)), "GAP GF("*string(p)*","*string(d)*")"))
  push!(mats, (one(GL(n,q)).X, "raw GAP GF("*string(p)*","*string(d)*")"))
  push!(mats, (identity_matrix(F, n), "FiniteField("*string(p)*","*string(d)*")"))
  if d == 1
    push!(mats, (identity_matrix(GF(p), n), "Oscar GF("*string(p)*","*string(d)*")"))
  end
  push!(testdata, (mats, "GF("*string(p)*","*string(d)*")"))
end

push!(testdata, ([
        (GAP.Globals.IdentityMat(n), "GAP Z"),
        (identity_matrix(ZZ, n), "Oscar ZZ"),
        (identity_matrix(QQ, n), "Oscar QQ"),
      ], "ZZ or QQ"))


suite = BenchmarkGroup()

for f in [test_add, test_mul]
#for f in [test_add, test_mul, test_minpoly, test_charpoly, test_order]
  println(nameof(f))
  suite[nameof(f)] = BenchmarkGroup()
  for (mats, ring) in testdata
    suite[nameof(f)][ring] = BenchmarkGroup()
    for (x, desc) in mats
      println("  input: ", desc)
      suite[nameof(f)][ring][desc] = @benchmarkable $(f)($x) # seconds=1
      f(x) # make sure it actually runs
    end
  end
end

end # module

# suite = benchmark_matrices.suite;
# tune!(suite)
# results = run(suite, verbose = true)
