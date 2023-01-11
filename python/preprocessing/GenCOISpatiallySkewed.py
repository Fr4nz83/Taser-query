### IMPORTS ###
import numpy as np
import sys
import matplotlib.pyplot as plt
import argparse
from scipy import spatial


def associateEdgeToVertices(fileEdges):
	vertexDic = {}

	# Read the lines within the edge file.
	with open(fileEdges) as f:
    		lines = f.readlines()

	# Associate each vertex with the ID of the first edge they appear in.
	for l in lines:
		tokens = l.split()
		idEdge = int(tokens[0])
		idV1 = int(tokens[1])
		idV2 = int(tokens[2])
		if idV1 not in vertexDic:
			vertexDic[idV1] = (idEdge, 0.0)
		if idV2 not in vertexDic:
			vertexDic[idV2] = (idEdge, 1.0)

	return vertexDic


def getArrayVertices(fileVertices):

	with open(fileVertices) as f:
		lines = f.readlines()
	
	# Generate an array containing ONLY the geo-coordinates of the vertices. The position
	# of a vertex in the array corresponds to its ID.
	vertices = []
	for l in lines:
		tokens = l.split()
		vertices.append([float(tokens[1]), float(tokens[2])])

	return vertices, len(lines)



### MAIN ###
def main():
	
	# Setup args app.
	parser = argparse.ArgumentParser(description="This file generates a COI with POIs having skewed spatial distribution.")
	parser.add_argument("-v", "--vertices", required=True, type=str, help="File vertices")
	parser.add_argument("-e", "--edges", required=True, type=str, help="File edges")
	parser.add_argument("-o", "--output", required=True, type=str, help="COI output file")
	parser.add_argument("-np", "--numpoi", required=True, type=int, help="The number of POIs to insert in the COI")

	# Parse command line.
	args = parser.parse_args()

	# Create a dictionary where we take note of the first edge with which each vertex is associated.
	print("Creating dictionary \"Vertex => First Edge\"...")
	vertexDic = associateEdgeToVertices(args.edges)

	# DEBUG: for each vertex, print the associated edge.
	#for v in vertexDic:
	#	print("V: " + str(v) + " - E: " + str(vertexDic[v]))
		
	# Initialize an array containing the vertices' coordinates.
	print("Initializing vector geo-coordinates of vertices...")
	vertices, numVertices = getArrayVertices(args.vertices)
	
	# Print the pairs representing the vertices.
	#for v in vertices:
	#	print(v)

	# Generate the KD-Tree.
	print("Generating kd-tree...")
	tree = spatial.KDTree(vertices)

	# Generating the COI...	
	idVertexHotspot = np.random.randint(0, numVertices - 1)
	print("Generating COI centered on vertex with ID " + str(idVertexHotspot))
	distance, index = tree.query(vertices[idVertexHotspot], args.numpoi)


	# Print out the POIs of this COI in a file.
	print("Writing the COI in file " + args.output)
	outfile = open(args.output,"w")
	for i in range(len(index)):
		outfile.write(str(i) + " " + str(vertices[index[i]][0]) + " " + str(vertices[index[i]][1]) + " ")
		outfile.write(str(vertexDic[index[i]][0]) + " " + str(vertexDic[index[i]][1]) + "\r\n")
	outfile.close()		
		


if __name__ == "__main__":
	main()
