import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

import numpy as np

def block_limited(size):
    max_waves=32*104
    waves= (np.array(size)*2**10 /8 )/64
    return np.minimum(waves/max_waves,1) 

def plot_block_limited():
    x=2**np.arange(np.log2(min(data["size"])),np.log2(max(data["size"])), (np.log2(max(data["size"])) - np.log2(min(data["size"])))/n)
    
    plt.plot(x,block_limited(x)*100,"--",label="Theo. Block limited")


sns.set_style("whitegrid")
if __name__ == "__main__":

    theoretical_bandwidth= 32 * 32 * 1.6
    n=10000

    data= pd.read_csv("bandwidth_data.txt",sep='\s+')
    data= data[data["function"]=="Copy"]
    data["bw"]=data["bw"]/2**20 * 1e+3
    data["bw"]/= (theoretical_bandwidth)/100

    sns.scatterplot(data,x="size", y="bw",hue="function",s=100)
    plot_block_limited()

    plt.xscale("log")
    plt.ylim([0,100])
    plt.yscale("linear")
    plt.ylabel(r"Bandwidth [%]")
    plt.xlabel(r"Size [$KiB$]")
    plt.legend()
    plt.show()
