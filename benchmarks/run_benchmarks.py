import numpy as np
import subprocess
import re
import pandas as pd



def get_performance_data(lines):

    regex_pattern=r"^([a-zA-Z]+)\W+([0-9]+(?:\.[0-9]+)?)\W+([0-9]+(?:\.[0-9]+)?)\W+([0-9]+(?:\.[0-9]+)?)\W+([0-9]+(?:\.[0-9]+)?)\W*$"

    data= {
        "function" : [] ,
        "bw" : [],
        "t_min" : [],
        "t_max" : [],
        "t_avg" : [],
    }

    for line in lines:
        m=re.match(regex_pattern,line)
        if m is not None:
            results=m.groups()
            data["function"].append(results[0])
            data["bw"].append(results[1])
            data["t_min"].append(results[2])
            data["t_max"].append(results[3])
            data["t_avg"].append(results[4])


    return pd.DataFrame(data)

if __name__ == "__main__":

    min_log= np.log2(8*2**10)
    max_log= np.log2( 16 * 2**30)
    n = 10
    sizes= 2**10 * np.int64( (2**(np.arange(min_log,max_log,(max_log-min_log)/n)) / 8  )/2**10)
    

    all_data=[]
    for size in sizes:
        try:
            output= subprocess.check_output(["./build/hip-stream" ,"-s",f"{size}"],stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as exc:
            print("Status : FAIL", exc.returncode, exc.output)
        else:        
            lines=output.decode('ascii').splitlines()
            data=get_performance_data(lines)
            data["size"]= size * 8 / 2**10

            all_data.append(data)

    print(pd.concat(all_data))