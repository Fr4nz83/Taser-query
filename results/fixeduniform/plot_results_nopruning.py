import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main(): 
	# Read data.
	dataNoPruning = pd.read_csv('nopruning.txt', sep=" ", skiprows=1, header=None, na_values='-')
	dataNoPruning.columns = ["length", "AMS CSE", "AMS Baseline", "OSLO CSE", "OSLO Baseline", "BER CSE", "BER Baseline"]
    
	# Configure the general properties of the plots.
	figNoPruningAMS = plt.figure(num=None, figsize=(8, 5))
	figNoPruningOSLO = plt.figure(num=None, figsize=(8, 5))
	figNoPruningBER = plt.figure(num=None, figsize=(8, 5))
	figs = [figNoPruningAMS.number, figNoPruningOSLO.number, figNoPruningBER.number]
    
    
	# Configure the general properties of the plots.
	SMALL_SIZE = 12
	MEDIUM_SIZE = 16
	BIGGER_SIZE = 18
	for i in figs:
		plt.figure(i)
        
		plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
		plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
		plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
		plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
		plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
		plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
		plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
	
		plt.autoscale(enable='True', axis='both')
		#plt.xscale("log")
		plt.yscale("log")
		plt.tick_params(top=False, right=True)
		plt.xticks(np.arange(min(dataNoPruning["length"]), max(dataNoPruning["length"]+1), 1.0)) 
		plt.xlabel('Length COI sequence')
		plt.ylabel('Execution time (sec.)')
    
    
	# Plot Amsterdam
	fig = plt.figure(figNoPruningAMS.number)
	ax = fig.add_subplot(111)
	figNoPruningAMS.suptitle('Amsterdam', y=0.96)
	ax.plot(dataNoPruning["length"], dataNoPruning["AMS CSE"], '-r^', label='LS-TASeR')
	ax.plot(dataNoPruning["length"], dataNoPruning["AMS Baseline"], '-go', label='Baseline')
	ax.legend()
	figNoPruningAMS.savefig("Amsterdam_pruning.pdf", bbox_inches = 'tight')


	# Plot Oslo
	fig = plt.figure(figNoPruningOSLO.number)
	ax = fig.add_subplot(111)
	figNoPruningOSLO.suptitle('Oslo', y=0.96)
	ax.plot(dataNoPruning["length"], dataNoPruning["OSLO CSE"], '-r^', label='LS-TASeR')
	ax.plot(dataNoPruning["length"], dataNoPruning["OSLO Baseline"], '-go', label='Baseline')
	ax.legend()
	figNoPruningOSLO.savefig("Oslo_pruning.pdf", bbox_inches = 'tight')


	# Plot Berlin
	fig = plt.figure(figNoPruningBER.number)
	ax = fig.add_subplot(111)
	figNoPruningBER.suptitle('Berlin', y=0.96)
	ax.plot(dataNoPruning["length"], dataNoPruning["BER CSE"], '-r^', label='LS-TASeR')
	ax.plot(dataNoPruning["length"], dataNoPruning["BER Baseline"], '-go', label='Baseline')
	ax.legend()
	figNoPruningBER.savefig("Berlin_pruning.pdf", bbox_inches = 'tight')

	return

if __name__ == "__main__":
    main()
