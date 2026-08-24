# 5-sunlet + 5 sunlet glued along an edge
N = graph_from_edges(Directed, [[11,1],[11,12],[12,2],[12,8],[11,6],[6,7],[7,5],[8,7],[6,9],[8,10],[9,10],[9,3], [10,4]])
net = phylogenetic_network(N)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)

H = components_of_kernel(2, phi, show_progress = true)
#this ideal has dimension 13, and has degree 8
I = ideal(reduce(vcat, collect(values(H))))

# Another example, 3-sunlet + 6 sunlet glued along an edge
N2 = graph_from_edges(Directed, [[6,1],[12,6],[7,6],[12,7],[11,12],[11,5],[11,10],[10,4],[10,9],[9,3],[9,8],[8,2],[8,7]])
net2 = phylogenetic_network(N2)
PM2 = cavender_farris_neyman_model(net2)
R2, q2 = model_ring(PM2)
S2, a2 = parameter_ring(PM2)
phi2 = parametrization(PM2)

H2 = components_of_kernel(2, phi2, show_progress = true)
#this ideal has dimension 10
I2 = ideal(reduce(vcat, collect(values(H2))))

# One More example!
N = graph_from_edges(Directed, [[1,6],[6,7],[6,12],[7,12],[7,8],[8,2],[8,9],[9,3],[10,9],[10,4],[11,10],[11,5],[12,11]])
net = phylogenetic_network(N)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)

H = components_of_kernel(2, phi, show_progress = true)
I = ideal(reduce(vcat, collect(values(H))))

# Example: 6-sunlet + 4-sunlet
N = graph_from_edges(Directed, [[1,6],[6,7],[6,12],[7,12],[7,8],[8,2],[8,9],[9,3],[10,9],[10,4],[11,10],[11,5],[12,11]])
net = phylogenetic_network(N)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)

H = components_of_kernel(2, phi, show_progress = true)
I = ideal(reduce(vcat, collect(values(H))))

##########################################################
# Examples to see if choice of reticulation edges matter #
##########################################################

#The two networks below fill space (can verify by counting parameters)

# choice 1
N8 = graph_from_edges(Directed, [[1,4],[4,5],[4,8],[8,7],[5,6],[5,8],[7,6],[7,3],[6,2]])
net8 = phylogenetic_network(N8)
PM8 = cavender_farris_neyman_model(net8)
phi8 = parametrization(PM8)

#this is empty
H8 = components_of_kernel(2, phi8, show_progress = true)
I8 = ideal(reduce(vcat, collect(values(H8))))

#choice 2
N9 = graph_from_edges(Directed, [[1,4],[4,5],[4,8],[8,7],[5,6],[8,5],[7,6],[7,3],[6,2]])
net9 = phylogenetic_network(N9)
PM9 = cavender_farris_neyman_model(net9)
phi9 = parametrization(PM9)

#this is empty
H9 = components_of_kernel(2, phi9, show_progress = true)
I9 = ideal(reduce(vcat, collect(values(H9))))

examples = Dict()

function compute_stats_of_the_ideal(N)
    examples = load("examples")
    PM = cavender_farris_neyman_model(N)
    phi = parametrization(PM)
    H = components_of_kernel(2, phi, show_progress = true)
    I = ideal(reduce(vcat, collect(values(H))))
    examples[N] = (dim(I), degree(I), is_prime(I))
    save("examples", examples)
end

