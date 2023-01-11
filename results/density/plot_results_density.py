import matplotlib.pyplot as plt
import pandas as pd

def main(): 
	# Read data.
	dataDensity = pd.read_csv('density_poi.txt', sep=" ", skiprows=1, header=None, na_values='-')
	dataDensity.columns = ["density", "amsterdam", "oslo", "berlin"]
    
	# Configure the general properties of the plots.
	figDensity = plt.figure(num=None, figsize=(8, 5))
	figs = [figDensity.number]
    
    
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
		plt.ylim(0,0.15)
		#plt.xscale("log")
		#plt.yscale("log")
		plt.xticks([5, 25, 50])
		plt.tick_params(top=False, right=True)    
		plt.xlabel('Number of POIs per COI')
		plt.ylabel('Execution time (sec.)')
    
    
	# Plot Density
	plt.figure(figDensity.number)
	plt.plot(dataDensity["density"], dataDensity["amsterdam"], '-ko', label='Amsterdam')
	plt.plot(dataDensity["density"], dataDensity["oslo"], '-rs', label='Oslo')
	plt.plot(dataDensity["density"], dataDensity["berlin"], '-bv', label='Berlin')
	plt.legend(loc='upper left')
	figDensity.savefig("POI_Density.pdf", bbox_inches = 'tight')

	return

if __name__ == "__main__":
    main()
