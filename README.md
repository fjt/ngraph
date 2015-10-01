# ngraph

basic usage

g=Ngraph.new
g.vertex=some_array
g.edge=array_of_vertex_pairs

gn=Ngraph.cube(16) ## 16 dimensional hyper-cube graph

gn.bfs(0).map{|n|n.length} ## do a breadth-first search, get each distanced vertices set size

gn.vertex.length == gn.bfs(0).flatten.length ## check if the graph is connected.

true ## true if connected
