# Ngraph

basically, input vertex and edge list and thats about it.

g=Ngraph.new

g.vertex=some_array

g.edge=array_of_vertex_pairs

Additionally, it has some typical graph generationg functions, for example

gn=Ngraph.cube(16) ## 16 dimensional hyper-cube graph. Obviously it has 65536 vertices, each has sixteen neighbours.

gn.bfs(0).map{|n|n.length} ## do a breadth-first search, and get each distanced vertices set size. the argument is breadth-first search start point. it will give you the same result, regardless of the starting point.

gn.vertex.length == gn.bfs(0).flatten.length ## check if the graph is connected, obviously true for this case.

gn.dbfs and gn.rdbfs does breadth first search with link direction.

gn.bfs takes a block argument which sets a terminal condition. It also takes direction as a second argument.

gn.scc(stt) returns a strong connected component of a directed graph

gn.update(vertice:vertices_list, edge:edge_list)
is added. it updates the existing network with new set of vertices and edges.
vertice given by name, while the edge must be given by the vertices index pairs.
this method returns a new network with node position inherited as much as possible from the original network.
must be useful for time-series network data.

gn.dpmds(plist)
is added. it is a new algorith developed to handle link direction within MDS framework.
the argment is a pivot node list given by vertices indices. 
