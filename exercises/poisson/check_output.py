import numpy as np
import struct
import matplotlib
import matplotlib.pylab as plt
import argparse

matplotlib.use('Qt5Agg')

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


def make_plot(file):
    (X,Y),phi=read_field(file)
    r=np.sqrt(X**2 + Y**2)
    plt.plot( r.flatten(),phi.flatten(),"o",label="phi")

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
        for file in args.files:
                make_plot(file)
        plt.show()
    
    