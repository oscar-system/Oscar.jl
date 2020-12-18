Read("bench.g");
BindConstant("reps", 10000000);
#reps := 1000000;

test_1 := function()
  local x, i;
  x := true;
  for i in [1..reps] do
    x := x and ReturnTrue();
  od;
end;

test_2 := function()
  local f, x, i;
  f := ReturnTrue;
  x := true;
  for i in [1..reps] do
    x := x and f();
  od;
end;

G := CyclicGroup(10);

test_isabelian_1 := function()
  local x, i;
  x := true;
  for i in [1..reps] do
    x := x and IsAbelian(G);
  od;
end;

test_isabelian_2 := function()
  local f, x, i;
  f := IsAbelian;
  x := true;
  for i in [1..reps] do
    x := x and f(G);
  od;
end;


if false then
Benchmark(test_1);
Benchmark(test_2);
Benchmark(test_isabelian_1);
Benchmark(test_isabelian_2);

Benchmark(test_2, rec(maxreps:=10));
fi;
