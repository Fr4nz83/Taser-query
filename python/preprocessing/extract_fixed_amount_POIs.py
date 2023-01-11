### IMPORTS ###
import numpy as np
import random
import sys



### ENTRY POINT APP ###
def main():
	
	if(len(sys.argv) != 3):
		print("Error: the following arguments are expected: input file, number of POIs to extract. Terminating...")
		return


	# Getting the number of lines of the input file...
	with open(sys.argv[1]) as f:
    		lines = f.readlines()

	numExtract = min(int(sys.argv[2]), len(lines))
	print("Input file: " + sys.argv[1])
	print("POIs in input file: " + str(len(lines)))
	print("Number of POIs to extract: " + str(numExtract))


	# Extracting POIs...
	random.shuffle(lines)
	print("Shuffling and extracting " + str(numExtract) + " POIs...")
	outfile = open(sys.argv[1] + "." + str(numExtract) + ".txt","w")
	for i in range(0, numExtract):
		outfile.write(lines[i])
	outfile.close()

	return



if __name__ == "__main__":
	main()
