import numpy as np
import struct
import matplotlib
import matplotlib.pylab as plt
import argparse


import re




def sort_natural(l):
    """ 
    Sort the given list in the way that humans expect.
    """

    def tryint(s):
        try:
            return int(s)
        except:
            return s
        

    def alphanum_key(s):
        """ Turn a string into a list of string and number chunks.
                "z23a" -> ["z", 23, "a"]
        """
        
        return [ tryint(c) for c in re.split('([0-9]+)', s) ]

    l.sort(key=alphanum_key)


def get_size(f):
    return struct.unpack('1L', bytearray(f.read(8)))[0]

def get_array(f,n):
    return struct.unpack(f'{n}d', bytearray(f.read(8*n)))

def read_field(filename):
        with open(filename, mode='rb') as f: # b is important -> binary
            nx=get_size(f)
            ny=get_size(f)
            phi=np.array(get_array(f,(nx+2)*(ny+2)) ).reshape((nx+2,ny+2))[1:nx+1,1:ny+1]
            mesh=generate_mesh(nx,ny)
        return mesh,phi

def generate_mesh(nx,ny):
    x=-1 + np.arange(0,nx)*2/nx
    y=-1 + np.arange(0,ny)*2/ny
    X,Y=np.meshgrid(x,y)
    return X.transpose(),Y.transpose()


def make_plot(file,color):
    matplotlib.use('Qt5Agg')
    (X,Y),phi=read_field(file)
    r=np.sqrt(X**2 + Y**2)
    plt.plot( r.flatten(),phi.flatten(),"o",label="phi",color=color)

def get_diff(file1,file2):
    (X1,Y1),phi1=read_field(file1)
    (X2,Y2),phi2=read_field(file2)

    diff=np.sum(np.abs(phi1 - phi2))/len(phi1)

    return diff



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Simple helper program to check the output of the poisson solver.")

    parser.add_argument('files', type=str, nargs='+',help="Files containing the result of the simulations")
    
    parser.add_argument('--plot' ,help = "Whether to create a plot of the output data",action="store_true",default=False)
    parser.add_argument('--compare' ,help = "Compare the output between two different files",action="store_true",default=False)

    args = parser.parse_args()

    if args.compare:
        
        if len(args.files)==2:
            diff = get_diff(args.files[0],args.files[1])
            print(f"Distance between {args.files[0]} and {args.files[1]}:  { diff }")
        else:
             raise ValueError("You must provide exactly two files to compare.")
        
    if args.plot :
        

        files=sort_natural(args.files)
        colors = plt.cm.summer(np.linspace(0.8,0.1,len(args.files)))

        for i,file in enumerate(args.files):
                print(file)
                make_plot(file,color=colors[i])
        
        plt.xlabel("r")
        plt.ylabel(r"$\phi(r)$")       
        
        plt.show()
    
    