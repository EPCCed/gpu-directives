import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns

sns.set_style("whitegrid")
if __name__ == "__main__":

    data= pd.read_csv("bandwidth_data.txt",sep='\s+')
    data= data[data["function"]=="Copy"]
    data["bw"]/=2**10
    data["bw"]=data["bw"]/2**30 * 1e+9
    data["bw"]/= (32*32 * 1.6)
    sns.scatterplot(data,x="size", y="bw",hue="function",s=100)
    plt.xscale("log")
    plt.ylim([0,1])
    #plt.yscale("log")
    plt.ylabel(r"Bandwidth [%]")
    plt.xlabel(r"Size [$KiB$]")
    plt.legend()
    plt.show()
