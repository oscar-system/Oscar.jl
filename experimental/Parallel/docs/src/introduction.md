```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Parallelization in Oscar

The functionality in the experimental package `Parallel` provides tools to do simple things
in parallel in Oscar. The main obstruction 
to using regular functions like `pmap` and friends is that most of the Oscar datatypes 
have `parent`s which need to be coherent throughout the worker pool and the main process. 
In order to achieve this, some more information on "type parameters" needs to be send 
up front when communicating with other processes; see the documentation of the 
[Oscar serialization](@ref dev_serialization) for details. 

Below we document specialized structs and methods which we overload so as to hopefully 
have a seamless integration of Oscar data types into `pmap` and other functions. 
Note that this is non-exhaustive work in progress and the user always remains 
responsible for serializing any type of object that they want to pass to workers! 

## `OscarWorkerPool`, `pmap`, and friends
```@docs 
OscarWorkerPool
oscar_worker_pool
```

## Slightly advanced distributed computing with centralized bookkeeping
```@docs
compute_distributed!
pop_task!
process_result!
```
