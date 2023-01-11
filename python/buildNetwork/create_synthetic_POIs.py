import sys
import numpy as np


### ENTRY POINT APP ###
def main():

	if(len(sys.argv) != 4):
		print("Error: the following arguments are expected: edges file, number of synthetic POIs to create, output file. Terminating...")
		return
	if (float(sys.argv[2]) <= 0.0):
		print("Error: the number of POIs to create must be > 0. Terminating...")
		return

	# Read the vertices from file
	with open(sys.argv[1]) as f:
    		edges = f.readlines()

	# Generate synthetic POIs from the set of edges.
	numPOIs = int(sys.argv[2])
	pois = []
	for id in range(numPOIs):
		idEdge = np.random.randint(0, len(edges))
		pois.append((idEdge, float(np.random.randint(0, 2))))
	
	# Write out the list of edges to a file.
	outfile = open(sys.argv[3],"w")
	id = 0
	for p in pois:
		outfile.write(str(id) + " " + "0.0" + " " + "0.0" + " " + str(p[0]) + " " + str(p[1]) + "\n")
		id += 1
	outfile.close()

if __name__ == "__main__":
	main()
