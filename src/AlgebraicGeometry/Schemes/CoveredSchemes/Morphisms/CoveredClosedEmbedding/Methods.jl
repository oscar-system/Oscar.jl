### forwarding the essential getters
underlying_morphism(phi::CoveredClosedEmbedding) = phi.f

### additional functionality
@doc raw"""
    image_ideal(phi::CoveredClosedEmbedding)
    
For a closed embedding $\phi \colon X \to Y$
return the sheaf of ideals on $Y$ defining the image
of $\phi$.
"""
image_ideal(phi::CoveredClosedEmbedding) = phi.I
