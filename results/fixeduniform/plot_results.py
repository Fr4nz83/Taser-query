import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def main(): 
    # First, read data of the baseline approach.
    columns = ["Length sequence", "CSR_PNE", "CSR_PNE_PRE"]

    dataModena = pd.read_csv('Modena2.txt', sep=" ", skiprows=1, header=None, na_values='-')
    dataModena.columns = columns

    dataAmsterdam = pd.read_csv('Amsterdam2.txt', sep=" ", skiprows=1, header=None, na_values='-')
    dataAmsterdam.columns = columns
    
    dataBerlin = pd.read_csv('Berlin2.txt', sep=" ", skiprows=1, header=None, na_values='-')
    dataBerlin.columns = columns

    dataOslo = pd.read_csv('Oslo2.txt', sep=" ", skiprows=1, header=None, na_values='-')
    dataOslo.columns = columns
    


    # Configure the general properties of the plots.
    figModena = plt.figure(num=None, figsize=(6, 5))
    figModenaSpeedup = plt.figure(num=None, figsize=(6, 5))
    figAmsterdam = plt.figure(num=None, figsize=(6, 5))
    figAmsterdamSpeedup = plt.figure(num=None, figsize=(6, 5))
    figBerlin = plt.figure(num=None, figsize=(6, 5))
    figBerlinSpeedup = plt.figure(num=None, figsize=(6, 5))
    figOslo = plt.figure(num=None, figsize=(6, 5))
    figOsloSpeedup = plt.figure(num=None, figsize=(6, 5))
    figs = [figModena.number, figModenaSpeedup.number, figAmsterdam.number, figAmsterdamSpeedup.number, figBerlin.number, figBerlinSpeedup.number, figOslo.number, figOsloSpeedup.number]
    
    
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
        plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
        plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
        
        plt.autoscale(enable='True', axis='both')
        #plt.xscale("log")
        #plt.yscale("log")
        plt.tick_params(top=False, right=True)    
        plt.xlabel('Length COI sequence')
        plt.ylabel('Execution time (sec.)')
        plt.ylim(0, 370)
    
    
    
    # Plot Modena
    fig = plt.figure(figModena.number)
    ax = fig.add_subplot(111)
    figModena.suptitle('Modena', y=0.96)
    ax.plot(dataModena["Length sequence"], dataModena["CSR_PNE"], '-r^', label='LS-TASeR')
    ax.plot(dataModena["Length sequence"], dataModena["CSR_PNE_PRE"], '-go', label='LS-TASeR (pre-computed)')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figModena.savefig("Modena_Fixed_Uniform.pdf", bbox_inches = 'tight')

    # Plot Modena (speedup)
    fig = plt.figure(figModenaSpeedup.number)
    ax = fig.add_subplot(111)
    figModenaSpeedup.suptitle('Modena', y=0.96)
    ax.plot(dataModena["Length sequence"], dataModena["CSR_PNE"]/dataModena["CSR_PNE_PRE"], '-bs', label='Speedup')
    ax.set_ylim(auto=True)
    ax.set_ylabel('Speedup')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figModenaSpeedup.savefig("Modena_Speedup_Fixed_Uniform.pdf", bbox_inches = 'tight')



    # Plot Amsterdam
    fig = plt.figure(figAmsterdam.number)
    ax = fig.add_subplot(111)
    figAmsterdam.suptitle('Amsterdam', y=0.96)
    ax.plot(dataAmsterdam["Length sequence"], dataAmsterdam["CSR_PNE"], '-r^', label='LS-TASeR')
    ax.plot(dataAmsterdam["Length sequence"], dataAmsterdam["CSR_PNE_PRE"], '-go', label='LS-TASeR (pre-computed)')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figAmsterdam.savefig("Amsterdam_Fixed_Uniform.pdf", bbox_inches = 'tight')

    # Plot Amsterdam (speedup)
    fig = plt.figure(figAmsterdamSpeedup.number)
    ax = fig.add_subplot(111)
    figAmsterdamSpeedup.suptitle('Amsterdam', y=0.96)
    ax.plot(dataAmsterdam["Length sequence"], dataAmsterdam["CSR_PNE"]/dataAmsterdam["CSR_PNE_PRE"], '-bs', label='Speedup')
    ax.set_ylim(auto=True)
    ax.set_ylabel('Speedup')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figAmsterdamSpeedup.savefig("Amsterdam_Speedup_Fixed_Uniform.pdf", bbox_inches = 'tight')
    
    
    # Plot Berlin.
    fig = plt.figure(figBerlin.number)
    ax = fig.add_subplot(111)
    figBerlin.suptitle('Berlin', y=0.96)
    ax.plot(dataBerlin["Length sequence"], dataBerlin["CSR_PNE"], '-r^', label='LS-TASeR')
    ax.plot(dataBerlin["Length sequence"], dataBerlin["CSR_PNE_PRE"], '-go', label='LS-TASeR (pre-computed)')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figBerlin.savefig("Berlin_Fixed_Uniform.pdf", bbox_inches = 'tight')

    # Plot Berlin (speedup).
    fig = plt.figure(figBerlinSpeedup.number)
    ax = fig.add_subplot(111)
    figBerlinSpeedup.suptitle('Berlin', y=0.96)
    ax.plot(dataBerlin["Length sequence"], dataBerlin["CSR_PNE"]/dataBerlin["CSR_PNE_PRE"], '-bs', label='Speedup')
    ax.set_ylim(auto=True)
    ax.set_ylabel('Speedup')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figBerlinSpeedup.savefig("Berlin_Speedup_Fixed_Uniform.pdf", bbox_inches = 'tight')


    # Plot Oslo.
    fig = plt.figure(figOslo.number)
    ax = fig.add_subplot(111)
    figOslo.suptitle('Oslo', y=0.96)
    ax.plot(dataOslo["Length sequence"], dataOslo["CSR_PNE"], '-r^', label='LS-TASeR')
    ax.plot(dataOslo["Length sequence"], dataOslo["CSR_PNE_PRE"], '-go', label='LS-TASeR (pre-computed)')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figOslo.savefig("Oslo_Fixed_Uniform.pdf", bbox_inches = 'tight')

    # Plot Oslo (speedup).
    fig = plt.figure(figOsloSpeedup.number)
    ax = fig.add_subplot(111)
    figOsloSpeedup.suptitle('Oslo', y=0.96)
    ax.plot(dataOslo["Length sequence"], dataOslo["CSR_PNE"]/dataOslo["CSR_PNE_PRE"], '-bs', label='Speedup')
    ax.set_ylim(auto=True)
    ax.set_ylabel('Speedup')
    ax.legend()
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1.0))
    figOsloSpeedup.savefig("Oslo_Speedup_Fixed_Uniform.pdf", bbox_inches = 'tight')




    return

if __name__ == "__main__":
    main()
