# A file with higher-level phylogenetic networks 
# leaves 1..n and interior nodes are ordered cyclically 

using Oscar
### Level 1 

# 4 sunlet

G_1_4_1 = graph_from_edges(Directed,[[5,1],[6,2], [7,3],[8,4],[5,6],[6,7],[7,8],[5,8]])


### Level 2 

# 2 leaves  

#diamond 

G_2_2 = graph_from_edges(Directed, [[3,1],[4,2],[5,4],[6,4],[5,3],[6,3],[7,5],[7,6]])

# 3 leaves 

#diamond with "interior leaf" 

#FIXME: doesn't have the right leaf set, apply a permutation to fix this! 
G_2_3_1 =  graph_from_edges(Directed, [[3,1],[4,2],[5,4],[6,4],[5,3],[6,3],[7,5],[7,6],[7,8]])

A = [[3,1],[4,2],[5,4],[6,4],[5,3],[6,3],[7,5],[7,6],[7,8]] 

perm = [1,2,3,7,5,6,4]
[permute!(a,perm) for a in A]
# pentagon 

G_2_3_2 = graph_from_edges(Directed, [[8,1],[4,8],[5,8],[5,4],[4,6],[5,7],[7,6],[6,2],[7,3]])

# 4 leaves 

#pentagon with "interior leaf" 

#FIXME: doesn't have the right leaf set, apply a permutation to fix this! 

G_2_4_1 = graph_from_edges(Directed, [[8,1],[4,8],[5,8],[5,9],[9,10],[9,4],[4,6],[5,7],[7,6],[6,2],[7,3]] ) 

# unbalanced hexagon 

G_2_4_2 = graph_from_edges(Directed,)
