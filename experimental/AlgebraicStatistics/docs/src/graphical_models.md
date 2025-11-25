```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Abstract Graphical Model

Graphical models are a parametrized abstract type `GraphicalModel{T, L}` where `T` is the type of the graph and `L` is the type of labelings of the graph, this is either a `NamedTuple` type or `Nothing`. Being an abstract type there is no constructor for a graphical model. All graphical models are subtypes of `GraphicalModel`, and they inherit their type parameters from the graph that is passed to the constructor. See the section corresponding to the graphical model type you would like to work with. 

## Attributes

The following functions are attributes for all graphical models.
```@docs
model_ring(M::T) where T <: GraphicalModel
parameter_ring(M::T) where T <: GraphicalModel
parametrization(M::T) where T <: GraphicalModel
```

## Vanishing ideal

Ensuring that any concrete graphical model (one that is a subtype of `GraphicalModel`) have implementations for the above attributes ensures
that we can compute the vanishing ideal.
In can be the case that certain graphical models will benefit from certain algorithms, and in those cases the `vanishing_ideal` will be overloaded for that specific type.
In all other cases the vanishing ideal computation will fallback to this general function.
There is a keyword argument for setting algorithm which is set to `:eliminate` by default, other possibilities are `:f4` and `:markov`.

If you are only interested in generators up to a certain degree, considering using [`components_of_kernel`](@ref) with `parametrization`.
```@docs
vanishing_ideal(M::GraphicalModel; algorithm::Symbol = :eliminate)
```
# Extending functionality

One major design decision was to allow for the graphs and their labels to parametrize the graphical models. 
This allows users to create their own graphical model types and gives them the possibility to overload the appropriate attribute functions.
In this way, users can reuse the implementations for the functions defined for the abstract grpahical models, 
which avoids code duplication and also allows for the application of theorems for specific types of graphical models.

