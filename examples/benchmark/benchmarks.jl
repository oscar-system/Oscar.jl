using Oscar
using GroebnerWalk
using BenchmarkTools

include("simple.jl")

open("results.csv", "a") do io
    simple(io)
end