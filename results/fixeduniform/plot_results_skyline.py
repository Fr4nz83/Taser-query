import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main(): 
	# Read data.
	dataSkylines = pd.read_csv('skylines.txt', sep=" ", skiprows=1, header=None, na_values='-')
	dataSkylines.columns = ["length", "linear AMS", "conventional AMS", "linear OSLO", "conventional OSLO", "linear BER", "conventional BER"]
    
	# Configure the general properties of the plots.
	figSkyAMS = plt.figure(num=None, figsize=(8, 5))
	figSkyOSLO = plt.figure(num=None, figsize=(8, 5))
	figSkyBER = plt.figure(num=None, figsize=(8, 5))
	figs = [figSkyAMS.number, figSkyOSLO.number, figSkyBER.number]
    
    
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
		#plt.yscale("log")
		plt.tick_params(top=False, right=True)
		plt.xticks(np.arange(min(dataSkylines["length"]), max(dataSkylines["length"]+1), 1.0)) 
		plt.xlabel('Length COI sequence')
		plt.ylabel('Average number of results')
		plt.ylim(0, 90)
    
    
	# Plot Amsterdam
	fig = plt.figure(figSkyAMS.number)
	ax = fig.add_subplot(111)
	figSkyAMS.suptitle('Amsterdam', y=0.96)
	ax.plot(dataSkylines["length"], dataSkylines["linear AMS"], '-r^', label='Linear')
	ax.plot(dataSkylines["length"], dataSkylines["conventional AMS"], '-go', label='Conventional')
	ax.legend()
	figSkyAMS.savefig("Amsterdam_Skyline.pdf", bbox_inches = 'tight')


	# Plot Oslo
	fig = plt.figure(figSkyOSLO.number)
	ax = fig.add_subplot(111)
	figSkyOSLO.suptitle('Oslo', y=0.96)
	ax.plot(dataSkylines["length"], dataSkylines["linear OSLO"], '-r^', label='Linear')
	ax.plot(dataSkylines["length"], dataSkylines["conventional OSLO"], '-go', label='Conventional')
	ax.legend()
	figSkyOSLO.savefig("Oslo_Skyline.pdf", bbox_inches = 'tight')


	# Plot Berlin
	fig = plt.figure(figSkyBER.number)
	ax = fig.add_subplot(111)
	figSkyBER.suptitle('Berlin', y=0.96)
	ax.plot(dataSkylines["length"], dataSkylines["linear BER"], '-r^', label='Linear')
	ax.plot(dataSkylines["length"], dataSkylines["conventional BER"], '-go', label='Conventional')
	ax.legend()
	figSkyBER.savefig("Berlin_Skyline.pdf", bbox_inches = 'tight')

	return

if __name__ == "__main__":
    main()
