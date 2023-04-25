@doc raw"""
    literature_tate_model(arxiv_id::String, equ_nr::String)

Many Tate models have been created in the F-theory literature.
A significant number of them have even been given specific
names, for instance the "U(1)-restricted SU(5)-GUT model".
This method has access to a JSON-database, from which it can
look up such literature models. Provided that the model is in
our database, this method constructs the model including some
of the information available in the literature.

(Eventually, we hope to provide *all* information that is
available on such models in the literature.)

```jldoctest
julia> t = literature_tate_model("1109.3454", "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> v = toric_ambient_space(t)
Normal, simplicial toric variety

julia> a10,a21,a32,a43,a65,w,x,y,z = gens(cox_ring(v));

julia> I = ideal([x,y,w]);

julia> v2 = blow_up(v,I)
Normal toric variety

julia> cox_ring(v2)
Multivariate Polynomial Ring in 10 variables a10, a21, a32, a43, ..., e over Rational Field graded by
  a10 -> [0 0]
  a21 -> [0 0]
  a32 -> [0 0]
  a43 -> [0 0]
  a65 -> [0 0]
  w -> [1 0]
  x -> [1 2]
  y -> [1 3]
  z -> [0 1]
  e -> [-1 0]
```
"""
function literature_tate_model(arxiv_id::String, equ_nr::String)
  # Read all literature models
  values = JSON.parsefile(joinpath(@__DIR__, "database.json"))
  #return values
  
  # Find the entry that we are looking for
  fitting_entries = filter(v -> (v["arxiv_id"] == arxiv_id && v["equ_nr"] == equ_nr), values)
  
  # Check if multiple entries match the criteria
  @req length(fitting_entries) == 1 "Either no or more than one fitting entry found in database"
  
  # Construct the model in question
  model = fitting_entries[1]
  auxiliary_base_ring, vars = PolynomialRing(QQ, string.(model["variables"]))
  a1 = zero(auxiliary_base_ring)
  if length(model["a1"]) > 0
    a1 = prod([vars[p[1]]^p[2] for p in model["a1"]])
  end
  a2 = zero(auxiliary_base_ring)
  if length(model["a2"]) > 0
    a2 = prod([vars[p[1]]^p[2] for p in model["a2"]])
  end
  a3 = zero(auxiliary_base_ring)
  if length(model["a3"]) > 0
    a3 = prod([vars[p[1]]^p[2] for p in model["a3"]])
  end
  a4 = zero(auxiliary_base_ring)
  if length(model["a4"]) > 0
    a4 = prod([vars[p[1]]^p[2] for p in model["a4"]])
  end
  a6 = zero(auxiliary_base_ring)
  if length(model["a6"]) > 0
    a6 = prod([vars[p[1]]^p[2] for p in model["a6"]])
  end
  t = global_tate_model([a1, a2, a3, a4, a6], auxiliary_base_ring, Int.(model["dim"]))
  
  # Read in the resolutions and process them to the correct format
  resolutions = model["resolutions"];
  resolutions = [[string.(x) for x in r] for r in resolutions]
  set_attribute!(t, :resolutions => resolutions)
  
  # Remember saved attributes of this model, in particular the resolution sequences known
  set_attribute!(t, :arxiv_id => model["arxiv_id"])
  set_attribute!(t, :equ_nr => model["equ_nr"])
  set_attribute!(t, :description => model["description"])
  
  # Return the model
  return t
end
