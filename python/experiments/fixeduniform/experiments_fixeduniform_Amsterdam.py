from subprocess import call, Popen, PIPE
import sys
import glob

def main():

	# Some paths and names
	executable = "./bin/CASRQ"
	dataset = "Amsterdam"
	dir_roadnetwork = "./datasets/" + dataset + "/roadnetwork/"
	dir_poi = "./datasets/" + dataset + "/poi/fixed_uniform/"

	approach = int(sys.argv[1])
	num_rounds = int(sys.argv[2])
	fileVertices = glob.glob(dir_roadnetwork + "*Vertices*")[0]
	fileEdges = glob.glob(dir_roadnetwork + "*Edges*")[0]
	filesPOIs = sorted(glob.glob(dir_poi + "*.txt"), key=str.lower)


	# Build the line used to execute the application.
	invoke = []
	invoke.append(executable)
	invoke.append("-a" + str(approach))
	invoke.append("-r" + str(num_rounds))
	invoke.append("-v" + fileVertices)
	invoke.append("-e" + fileEdges)


	# Run the experiments.
	POIList = []
	outfile = open(dataset + "_fixed_uniform_" + str(approach) + "_results.txt","w")
	outfile.write("# col1: num. candidates - col2: OSR query time - col3: approach query time - col4: total query time\n\n")
	for i in range(0, 6):
		POIList.append("-p" + filesPOIs[i%len(filesPOIs)])
		p = Popen(invoke + POIList, stderr=PIPE)
		output, err = p.communicate()
		result = str(i+1) + " " + str(50**(i+1)) + " " + str(err).replace("b","").replace("'","") + "\n"
		outfile.write(str(result))
	outfile.close()

	return


if __name__ == "__main__":
	main()
