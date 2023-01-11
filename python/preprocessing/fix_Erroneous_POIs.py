### IMPORTS ###
import sys



### ENTRY POINT APP ###

def main():
	
	if(len(sys.argv) != 3): 
		print("Error. The following arguments are expected: edge file, POI file to correct. Terminating...")
		return


	# Read the edge file and find the largest ID.
	with open(sys.argv[1]) as f:
    		lines = f.readlines()

	maxIDEdge = 0;
	for l in lines:
		maxIDEdge = max(int(l.split(" ")[0]), maxIDEdge) 

	print("Reading Edge file: " + sys.argv[1])
	print("Lines edge file: " + str(len(lines)))
	print("Max ID edge: " + str(maxIDEdge))



	# Read the POI file and filter out the POIs associated with a wrong edge.
	with open(sys.argv[2]) as f:
    		lines = f.readlines()

	poi = []
	for l in lines:
		if(int(l.split(" ")[3]) <= maxIDEdge): poi.append(l)
		#else: print(l, end='')

	print("Reading POI file: " + sys.argv[2])
	print("Number of POIs in the file: " + str(len(lines)))
	print("Number of correct POIs: " + str(len(poi)))



	# Writing out the corrected file.
	print("Writing corrected POI file: " + sys.argv[2] + ".CORRECTED.txt")
	outfile = open(sys.argv[2] + ".CORRECTED.txt","w")
	for i in poi:
		outfile.write(i)
	outfile.close()
	

if __name__ == "__main__":
	main()
