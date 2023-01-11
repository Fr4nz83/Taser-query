### IMPORTS ###
import numpy as np
import sys
import matplotlib.pyplot as plt



### FUNCTION DEFINITIONS ###

# Generates "size" integer values distributed according to a uniform distribution
# within the range [lb, ub].
def genUniformInt(lb, ub, size):
	return np.random.randint(lb,ub,size)

# Generates "size" float values distributed according to a uniform distribution
# within the range [lb,ub].
def genUniformFloat(lb, ub, size):
	return np.random.uniform(lb,ub,size)

# Generates "size" values distributed according to a Normal distribution
# having "mean = mu" and "std = sigma".
def genGaussian(mu, sigma, size):
	return np.random.normal(mu,sigma,size)

# Generates "size" values distributed according to a zipf distribution
# having "mean = mu" and "std = sigma".
def genZipf(a, size):
	return np.random.zipf(a, size)

# Returns the number of lines in "fname".
def getNumLinesFile(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1



### ENTRY POINT APP ###

def main():
	
	if(len(sys.argv) != 3): 
		print("Error. The following arguments are expected: input file, output file, distribution type. Terminating...")
		return

	with open(sys.argv[1]) as f:
    		lines = f.readlines()

	print("Input file: " + sys.argv[1])
	print("Lines input file: " + str(len(lines)))

	lb = 1
	ub = 100
	if(sys.argv[2] == "0"):
		outfile = open(sys.argv[1] + ".UNIFORM.txt","w")
		print("Generating uniform cost distribution.")
		for l in lines:
			outfile.write(l.replace("\r\n", "\n").replace("\n", " ") + str(np.random.randint(lb,ub)) + "\n")
		outfile.close()

	elif(sys.argv[2] == "1"):
		outfile = open(sys.argv[1] + ".NORMAL.txt","w")
		print("Generating Normal cost distribution.")
		mu = 0
		sigma = 1

		values = genGaussian(mu, sigma, len(lines))

		minval = min(values)
		maxval = max(values)
		for i in range(len(values)):
			values[i] = ((values[i] - minval) / (maxval - minval) * 99) + 1

		#plt.hist(values, normed=True, bins=20)
		#plt.show()

		for i in range(len(values)):
			outfile.write(lines[i].replace("\r\n", "\n").replace("\n", " ") + str(int(values[i])) + "\r\n")
		outfile.close()

	else:
		outfile = open(sys.argv[1] + ".ZIPF.txt","w")
		print("Generating zipf cost distribution.")
		a = 2.

		values = genZipf(a, len(lines))

		#plt.hist(values, normed=True, bins=100)
		#plt.show()

		for i in range(len(values)):
			outfile.write(lines[i].replace("\r\n", "\n").replace("\n", " ") + str(values[i]) + "\r\n")
		outfile.close()

if __name__ == "__main__":
	main()
