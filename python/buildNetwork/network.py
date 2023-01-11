import osmnx as ox
import matplotlib.pyplot as plt

# Specify the name that is used to seach for the data
place_name = "Turin, Italy"

# Fetch OSM street network from the location
graph = ox.graph_from_place(place_name)
type(graph)
# fig, ax = ox.plot_graph(graph)
nodes, edges = ox.graph_to_gdfs(graph)

# Retrieve restaurants
#restaurants = ox.pois_from_place(place_name, amenities=['school'])

# Print some info about the road network that was retrieved.
print("Vertices: " + str(len(nodes)))
print("Edges: " + str(len(edges)))
# print("# POI: " + str(len(restaurants)))

# map vertices and edges
dictVertices = {}
vertices = []
id = 0
for el in zip(nodes["osmid"], nodes["x"], nodes["y"]):
	dictVertices[el[0]] = id
	newEl = (id, el[1], el[2])
	vertices.append(newEl)
	id = id+1

# Write out the set of vertices.
outfile = open("vertices.txt","w")
for v in vertices:
	outfile.write(str(v[0]) + " " + str(v[1]) + " " + str(v[2]) + "\n")
outfile.close()


# Create an undirected graph out of the graph we got.
dictEdges = {}
for el in zip(edges["u"], edges["v"]):
	origin = min(dictVertices[el[0]], dictVertices[el[1]])
	dest = max(dictVertices[el[0]], dictVertices[el[1]])
	if(dictEdges.get(origin) is None):
		dictEdges[origin] = {dest} # Create a singleton.
	else:
		dictEdges[origin].add(dest) # Expand the existing set.

# Write out the list of edges to a file.
outfile = open("edges.txt","w")
id = 0
for se in dictEdges:
	for e in dictEdges[se]:
		outfile.write(str(id) + " " + str(se) + " " + str(e) + "\n")
		id += 1
outfile.close()

