###############################################################################
###############################################################################
### Definition and constructors
###############################################################################
###############################################################################

@doc Markdown.doc"""
    Cone(Rays)

     A polyhedral cone, not necessarily pointed, defined by the positive hull
     of the rays, Rays.
"""
struct Cone #a real polymake polyhedron
    pm_cone::Polymake.BigObject
end
function Cone(Rays::Union{Oscar.MatElem,AbstractMatrix})
    Cone(Polymake.polytope.Cone{Polymake.Rational}(
        INPUT_RAYS = matrix_for_polymake(Rays),
    ))
end
function Cone(Rays::Union{Oscar.MatElem,AbstractMatrix}, LS::Union{Oscar.MatElem,AbstractMatrix}; non_redundant::Bool=false)
   if non_redundant
       Cone(Polymake.polytope.Cone{Polymake.Rational}(
           RAYS = matrix_for_polymake(Rays),
           LINEALITY_SPACE = matrix_for_polymake(LS),
       ))
   else
       Cone(Polymake.polytope.Cone{Polymake.Rational}(
           INPUT_RAYS = matrix_for_polymake(Rays),
           INPUT_LINEALITY = matrix_for_polymake(LS),
       ))
   end
end

==(C0::Cone, C1::Cone) = Polymake.polytope.equal_polyhedra(pm_cone(C0), pm_cone(C1))


"""
    positive_hull(generators::Union{Oscar.MatElem,AbstractMatrix})

A polyhedral cone, not necessarily pointed, defined by the positive hull
of the `generators`. Redundant rays are allowed in the generators.
"""
function positive_hull(generators::Union{Oscar.MatElem,AbstractMatrix})
    # TODO: Filter out zero rows
    C=Polymake.polytope.Cone{Polymake.Rational}(INPUT_RAYS =
      matrix_for_polymake(remove_zero_rows(generators)))
    Cone(C)
end


"""
    pm_cone(C::Cone)

Get the underlying polymake `Cone`.
"""
pm_cone(C::Cone) = C.pm_cone


###############################################################################
###############################################################################
### Display
###############################################################################
###############################################################################

function Base.show(io::IO, C::Cone)
    print(io,"A polyhedral cone in ambient dimension $(ambient_dim(C))")
end

Polymake.visual(C::Cone; opts...) = Polymake.visual(pm_cone(C); opts...)

