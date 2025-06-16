```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Graphical Models

Graphical models are a parametrized abstract type `GraphicalModel{T, L}` where `T` is the type of graph (currently `Directed` or `Undirected`) and `L` is the type of labelings of the graph, this is either a `NamedTuple` type or `Nothing`. Being an abstract type there is no constructor for a graphical model. All graphical models are subtypes of `GraphicalModel`, and they inherit their type parameters from the graph that is passed to the constructor. See the section corresponding to the graphical model type you would like to work with. 



# Extending functionality

add in depth description of how to extend existing graphical models functionality.
