export chow_ring, augmented_chow_ring, select

@doc Markdown.doc"""
    chow_ring(M::Matroid, ring::MPolyRing=nothing; extended=false)

Return the Chow ring of a matroid, optional also with the simplicial generators.

See [AHK18](@cite) and [BES21](@cite). 

# Examples
The following computes the chow ring of the Fano matroid.
```jldoctest
julia> M = fano_matroid();

julia> R = chow_ring(M);

julia> R[1]*R[8]
-x_{3,4,7}^2
```

# Examples
The following computes the Chow ring of the Fano matroid including variables for the simplicial generators.
```jldoctest
julia> M = fano_matroid();

julia> R = chow_ring(M, extended=true);

julia> f = R[22] + R[8] - R[29]
x_{1,2,3} + h_{1,2,3} - h_{1,2,3,4,5,6,7}

julia> f==0
true
```

The following computes the chow ring of the Fano matroid including variables for the simplicial generators.
```jldoctest
julia> M = uniform_matroid(3,3);

julia> GR, _ = GradedPolynomialRing(QQ,["a","b","c","d","e","f"]);

julia> R = chow_ring(M, GR);

julia> hilbert_series_reduced(R)
(t^2 + 4*t + 1, 1) 

```
"""
function chow_ring(M::Matroid, ring::Union{MPolyRing,Nothing}=nothing; extended::Bool=false)
        is_loopless(M) || error("Matroid has loops")
	Flats = flats(M)
	number_flats = length(Flats)
        number_flats >= 2 || error("matroid has to few flats")
        proper_flats = Flats[2:number_flats-1]

	if ring==nothing
		var_names = [replace(string("x_",S), "["=>"{", "]"=>"}", ", "=>",") for S in proper_flats] # create variable names, indexed by the proper flats of M
        	if extended
                	var_names = [var_names; [replace(string("h_",S), "["=>"{", "]"=>"}", ", "=>",") for S in [proper_flats;[Flats[number_flats]]]]] #add the variables for the simplicial generators
        	end
		ring, vars = PolynomialRing(QQ, var_names, cached=false)
	else
		if typeof(ring)<:MPolyRing_dec
			tR = ring.R
		else
			tR = ring
		end
		println(tR.nvars)
		println(length(proper_flats))
		tR.nvars == length(proper_flats) || error("the ring has the wrong number of variables")
		vars = [ring[i] for i in 1:tR.nvars]
	end

	I = linear_relations(ring, proper_flats, vars, M)
        J = quadratic_relations(ring, proper_flats, vars)
	Ex = elem_type(ring)[]
        if extended
                Ex = relations_extended_ring(ring, proper_flats, vars)
        end
        chow_modulus = ideal(ring, vcat(I, J, Ex))
        chow_ring, projection = quo(ring, chow_modulus)
        return chow_ring
end

function linear_relations(ring::MPolyRing, proper_flats::GroundsetType, vars::Vector, M::Matroid)
        alpha = zero(ring)
        relations = elem_type(ring)[]
        for i in M.groundset
                poly = ring(0)
                for index in findall(issubset([i],F) for F in proper_flats)
                        poly+= vars[index]
                end
                if i==M.groundset[1]
                        alpha = poly
                else
                        push!(relations, alpha-poly)
                end
        end
        return relations
end

function quadratic_relations(ring::MPolyRing, proper_flats::GroundsetType, vars::Vector)
	println(typeof(vars[1]))
        relations = elem_type(ring)[]
        for i in 1:length(proper_flats)
                F = proper_flats[i]
                for j in 1:i-1
                        G = proper_flats[j]
                        if !issubset(F,G) && !issubset(G,F)
                                push!(relations, vars[i]*vars[j])
                        end
                end
        end
        return relations
end

function relations_extended_ring(ring::MPolyRing, proper_flats::GroundsetType, vars::Vector)
        relations = elem_type(ring)[]
        s = length(proper_flats)
        # h_E = alpha = -x_E
        poly = ring(0)
        for index in findall(issubset([1],F) for F in proper_flats)
                poly+= vars[index]
        end
        push!(relations, poly-vars[2*s+1])

        # h_F = h_E - sum x_G where G is a proper flat containing F
        for i in 1:s
                F = proper_flats[i]
                poly = ring(0)
                for index in findall(issubset(F,G) for G in proper_flats)
                        poly+= vars[index]
                end
                push!(relations, poly+vars[s+i]-vars[2*s+1]) #add h_F and subtract h_E
        end
        return relations
end

@doc Markdown.doc"""
    augmented_chow_ring(M::Matroid)

Return an augmented Chow ring of a matroid. As described in BHMPW20(@cite). 

# Examples
```jldoctest
julia> M = fano_matroid();

julia> R = augmented_chow_ring(M);
```
"""
function augmented_chow_ring(M::Matroid)
        Flats = flats(M)
        sizeFlats = length(Flats)
        n = size_groundset(M)

        is_loopless(M) || error("Matroid has loops")
        sizeFlats>3 || error("Matroid has too few flats")

        proper_flats = Flats[1:sizeFlats-1]
        element_var_names = [string("y_", S) for S in M.groundset]
        flat_var_names = [replace(string("x_",S), "["=>"{", "]"=>"}", ", "=>",") for S in proper_flats]
        flat_var_names[1] = "x_{}" # Override "x_Any{}"
        var_names = vcat(element_var_names, flat_var_names)
        s = length(var_names)

        ring, vars = PolynomialRing(QQ, var_names)
        element_vars = vars[1:n]
        flat_vars = vars[n+1:s]

        I = augmented_linear_relations(ring, proper_flats, element_vars, flat_vars, M)
        J = augmented_quadratic_relations(ring, proper_flats, element_vars, flat_vars, M)

        chow_modulus = ideal(ring, vcat(I, J))
        chow_ring, projection = quo(ring, chow_modulus)

        return chow_ring
end

function augmented_linear_relations(ring::MPolyRing, proper_flats::GroundsetType, element_vars::Vector, flat_vars::Vector, M::Matroid)
        n = size_groundset(M)
        relations = elem_type(ring)[]
        i = 1
        for element in M.groundset
                relations[i] = element_vars[i]
                j = 1
                for proper_flat in proper_flats
                        if !(element in proper_flat)
                                relations[i] -= flat_vars[j]
                        end
                        j += 1
                end
                i += 1
        end
        return relations
end

function augmented_quadratic_relations(ring::MPolyRing, proper_flats::Vector, element_vars::Vector, flat_vars::Vector, M::Matroid)
        incomparable_polynomials = quadratic_relations(ring, proper_flats, flat_vars)
        xy_polynomials = elem_type(ring)[]

        i = 1
        for element in M.groundset
                j = 1
                for proper_flat in proper_flats
                        if !(element in proper_flat)
                                push!(xy_polynomials, element_vars[i] * flat_vars[j])
                        end
                        j += 1
                end
                i += 1
        end

        return vcat(incomparable_polynomials, xy_polynomials)
end

"""
A helper function to select indicies of a vector that do `include` elments of a given set and `exclude` anothers
"""
function select(include::Union{AbstractVector,Set},exclude::Union{AbstractVector,Set},set::Union{AbstractVector,Set})
        all = union(all...)
        return findall(s->issubset(include,s)&&issubset(s,setdiff(all,exclude)),set);
end
