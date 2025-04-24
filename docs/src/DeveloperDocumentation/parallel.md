```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Writing parallel methods

## Simple single-machine parallelization

In this file we explain and implement simple patterns for parallelization 
on a single machine with multiple cores. 

This can be used to deploy embarassingly parallel tasks on multiple 
cores and either wait for all of them to finish, or for the first 
successful computation on one of the cores. 

The pattern is the following.
   1. Choose a parallel task to be carried out.
   2. Wrap up the input data for the task in a concrete instance, say 
      `MyParallelTask`, of `ParallelTask` according to the rules below. 
   3. Implement `_compute(::MyParallelTask)` to carry out the task at hand.
   4. Use `parallel_all` and `parallel_any` on `Vector`s of `MyParallelTask` 
      to do things in parallel. 

## ParallelTask

In order for the generic code below to work, any concrete instance 
must be of the form 
```
struct MyParallelTask{T1, T2, ...} <: ParallelTask
  field1::T1
  field2::T2
  ...
end
```
that is with *concrete types* for the fields with *every one* of those 
types appearing in the same order as type parameters. 

It is also important to note that the serialization for each of the concrete
types in a `ParallelTask` must be implemented, see [serialization](@ref dev_serialization).

The following is a generic implementation which hopefully serves for 
most concrete tasks automatically. The methods might need to be overwritten, 
though!

## The `_compute` method

The method of `_compute` for the concrete task specifies what to do 
on the respective worker. The data will be extracted from the task's 
fields. The return value must be of the form `(success::Bool, result)` 
where `success` is to indicate whether we have an affirmative result 
in any reasonable sense. It is used to decide whether `all` or `any` 
of a given list of tasks to be done in parallel is achieved. 

The second value `result` can be any reasonable result of the computation 
which is to be returned to the main node. Note that you must not create 
new parents on the worker which are required for the contents of `return`, 
i.e. they need to use the parent-like objects sent from the main node. 

## Parallel methods
```@docs
parallel_all
parallel_any
```
