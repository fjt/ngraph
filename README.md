# ngraph

basic usage

gn=Ngraph.cube(16) ## 16 dimensional hyper-cube graph
gn.bfs(0).map{|n|n.length} ## do a breadth-first search, get each distanced vertices set size
gn.vertex.length == gn.bfs(0).flatten.length 
true ## check if the graph is connected.
