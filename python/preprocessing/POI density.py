### IMPORTS ###
import numpy as np
import random
import sys



### ENTRY POINT APP ###
def main():
	
	if(len(sys.argv) != 3):
		print("Error: the following arguments are expected: input file, fraction of POIs to select. Terminating...")
		return
	if (float(sys.argv[2]) <= 0.0) or (float(sys.argv[2]) > 1.0):
		print("Error: the fraction of selected POIs must fall in the (0,1] range. Terminating...")
		return

	# Getting the lines of the input file...
	with open(sys.argv[1]) as f:
    		lines = f.readlines()


	print("Input file: " + sys.argv[1])
	print("Fraction of POIs: " + sys.argv[2])
	print("Lines input file: " + str(len(lines)))


	# Extracting POIs...
	random.shuffle(lines)
	sizeFile = len(lines) * float(sys.argv[2])
	print("Shuffling and selecting " + str(int(sizeFile)) + " POIs...")
	outfile = open(sys.argv[1] + "." + sys.argv[2] + ".txt","w")
	for i in range(0, int(sizeFile)):
		outfile.write(lines[i])
	outfile.close()

	return



if __name__ == "__main__":
	main()
