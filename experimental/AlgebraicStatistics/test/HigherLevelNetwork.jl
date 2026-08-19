
# 5-sunlet + 4 sunlet glued along an edge
N = graph_from_edges(Directed, [[10, 9], [10, 11], [8, 9], [7, 8], [11, 7], [7, 6], [11, 12], [12, 6], [6, 1], [8, 2], [9, 3], [10, 4], [12, 5]])
net = phylogenetic_network(N)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)

H = components_of_kernel(2, phi, show_progress = true)
<<<<<<< HEAD
I = ideal(reduce(vcat, collect(values(H))))
=======
I = ideal(reduce(vcat, collect(values(I))))
>>>>>>> origin


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

N = graph_from_edges(Directed, E)
net = phylogenetic_network(N)
PM = cavender_farris_neyman_model(net)
R, q = model_ring(PM)
S, a = parameter_ring(PM)
phi = parametrization(PM)
H = components_of_kernel(2, phi, show_progress = true)
I = ideal(reduce(vcat, collect(values(H))));
ngens(I)
