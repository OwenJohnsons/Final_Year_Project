import numpy as np
import glob 
from astropy.io import fits
from matplotlib import pyplot as plt 
from astropy.stats import *
from astropy.visualization import *



def SigmaCombine(filenames, nsig, iter, x1, x2, y1, y2): 
    '''
    Input 
    --- 
    filenames - list of path to fits file 
    nsig - standard deviation used in clipping 
    Iter - How many iterations of clipping that the function passes

    Output 
    --- 
    A single array equivalent to the inputted fits file. 
    '''

    datalist = []
    
    for file in filenames:
        with fits.open(file) as filename:
            datalist.append(filename[0].data[x1:x2, y1:y2])
        
    #converting this list to an array so that array operations can be performed
    loadeddata = np.asarray(datalist)
    #swapping around the dimensions of the array in order to obtain an array of every pixel value for each of the images. 
    Transposeddata=np.concatenate(np.swapaxes((loadeddata),0,2))

    i =0
    e = []
    while i < len(Transposeddata):
        e.append((np.ma.mean(sigma_clip(Transposeddata[i],nsig,iter))))
        i +=1
        #sigma clipping the values at each pixel
    
    sigmaclip_Transposed = np.asarray(e)

    #reversing our dimesion swaps to obtain an array of sigma clipped pixels
    SigmaClipped = np.asarray(np.array_split(sigmaclip_Transposed,len(loadeddata[0][2])))
    return SigmaClipped.T    

def quick_average(file_names, plot_cond): 
    master_data = [] 

    for i in range(0, len(file_names)): 
        data = fits.getdata(file_names[i])
        data = data[x1:x2, y1:y2]
        master_data.append(data)

    avg_data = np.average(master_data, axis=0)
    print(np.shape(avg_data))

    if plot_cond == True:
        ref = avg_data
        norm = ImageNormalize(ref, interval=ZScaleInterval(), stretch=SinhStretch())
        plt.imshow(ref, cmap=plt.cm.gray, norm=norm, interpolation='none')
        plt.tick_params(labelsize=16)
        # plt.gca().invert_yaxis()
        plt.show()

    return avg_data

def vmin_vmax(data):
    vmin = data.mean() - data.std()
    vmax = data.mean() + data.std()
    return vmin, vmax


def processed_fits(file_names, master_flat_norm, plot_cond): 
    master_data = [] 

    for i in range(0, len(file_names)): 
        data = fits.getdata(file_names[i])
        data = data[x1:x2, y1:y2]/master_flat_norm
        master_data.append(data)

    avg_data = (np.average(master_data, axis=0))
    print(np.shape(avg_data))

    if plot_cond == True:
        ref = avg_data
        norm = ImageNormalize(ref, interval=ZScaleInterval(), stretch=SinhStretch())
        plt.imshow(ref, cmap=plt.cm.gray, norm=norm, interpolation='none')
        plt.tick_params(labelsize=16)
        # plt.gca().invert_yaxis()
        plt.show()

    return avg_data

def find_nearest(array, value):
    array = np.asarray(array); idx = (np.abs(array - value)).argmin()
    return array[idx], idx