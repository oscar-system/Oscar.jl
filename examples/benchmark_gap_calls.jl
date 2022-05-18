module benchmark_gap_calls

using Oscar, BenchmarkTools

const reps = 1000000

G = cyclic_group(10)

function test_isabelian_0()
  x = true
  for i = 1:reps
    x = x && is_abelian(G)
  end
end

function test_isabelian_1()
  x = true
  for i = 1:reps
    x = x && GAP.Globals.IsAbelian(G.X)
  end
end

function test_isabelian_1b()
  x = true
  g = G.X
  for i = 1:reps
    x = x && GAP.Globals.IsAbelian(g)
  end
end

function test_isabelian_2()
  f = GAP.Globals.IsAbelian
  x = true
  for i = 1:reps
    x = x && f(G.X)
  end
end

function test_isabelian_2b()
  f = GAP.Globals.IsAbelian
  x = true
  g = G.X
  for i = 1:reps
    x = x && f(g)
  end
end

function test_isabelian_2c_inner(f, g)
  x = true
  for i = 1:reps
    x = x && f(g)
  end
end
function test_isabelian_2c()
  f = GAP.Globals.IsAbelian
  g = G.X
  test_isabelian_2c_inner(f, g)
end
end # module

if false

@btime benchmark_gap_calls.test_isabelian_0()
@btime benchmark_gap_calls.test_isabelian_1()
@btime benchmark_gap_calls.test_isabelian_1b()
@btime benchmark_gap_calls.test_isabelian_2()
@btime benchmark_gap_calls.test_isabelian_2b()
@btime benchmark_gap_calls.test_isabelian_2c()

@btime test_isabelian_0()
@btime test_isabelian_1()
@btime test_isabelian_1b()
@btime test_isabelian_2()
@btime test_isabelian_2b()
@btime test_isabelian_2c()

end
