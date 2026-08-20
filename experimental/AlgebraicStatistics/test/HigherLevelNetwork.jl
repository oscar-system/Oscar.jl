
# 5-sunlet + 4 sunlet glued along an edge
N1 = graph_from_edges(Directed, [[10, 9], [10, 11], [8, 9], [7, 8], [11, 7], [7, 6], [11, 12], [12, 6], [6, 1], [8, 2], [9, 3], [10, 4], [12, 5]])
net = phylogenetic_network(N1)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)

H = components_of_kernel(2, phi, show_progress = true)
I = ideal(reduce(vcat, collect(values(H))))


# two 6 sunlets glued along an edge
E = [[12,11], [12,13], [12,17],     
 [11,10], [11,3],
 [10,9],  [10,2],
 [13,14], [13,4],
 [14,15], [14,5],
 [15,16], [15,6],
 [16,17], [16,7],
 [17,18],                       
 [18,9],  [18,8],
 [9,1]]
N2 = graph_from_edges(Directed, E)
net = phylogenetic_network(N2)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM);
H = components_of_kernel(2, phi, show_progress = true)
I = ideal(reduce(vcat, collect(values(H))));
ngens(I)



N3 = graph_from_edges(Directed, [[3,1],[3,4],[4,2],[4,8],[3,6],[6,7],[7,5],[8,7],[6,9],[8,10],[9,10],[9,11], [10,12]])
net = phylogenetic_network(N3)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)
H = components_of_kernel(2, phi, show_progress = true)
I = ideal(reduce(vcat, collect(values(H))));
ngens(I)

N4 = graph_from_edges(Directed, [[6,1],[12,6],[7,6],[12,7],[11,12],[11,5],[11,10],[10,4],[10,9],[9,3],[9,8],[8,2],[8,7]])
net = phylogenetic_network(N4)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)
H = components_of_kernel(2, phi, show_progress = true)
I = ideal(reduce(vcat, collect(values(H))));
ngens(I)
