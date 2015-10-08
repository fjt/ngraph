# Ngraph

basically, input vertex and edge list and thats about it.

g=Ngraph.new

g.vertex=some_array

g.edge=array_of_vertex_pairs

Additionally, it has some typical graph generationg functions, for example

gn=Ngraph.cube(16) ## 16 dimensional hyper-cube graph. Obviously it has 65536 vertices, each has sixteen neighbours.

gn.bfs(0).map{|n|n.length} ## do a breadth-first search, and get each distanced vertices set size. the argument is breadth-first search start point. it will give you the same result, regardless of the starting point.

gn.vertex.length == gn.bfs(0).flatten.length ## check if the graph is connected, obviously true for this case.
