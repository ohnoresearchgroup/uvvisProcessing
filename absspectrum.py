import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class AbsSpectrum():
    def __init__(self,fullpath,horiba=False,bgrange=None,normrange=None):
        self.fullpath = fullpath

        if horiba is False:
            data = pd.read_csv(fullpath,skiprows=2, names =['wl','abs','extra'] )
            data = data.drop(columns='extra')
            data = data.dropna()
        else:
            data = pd.read_csv(fullpath,skiprows=12,sep='\t',names =['wl','abs','extra'])
            data = data.drop(columns='extra')
            data = data.dropna()


        self.data = data
        if bgrange is None:
            #set abscorr equal to abs
            self.data['abscorr'] = self.data['abs']
        else:
            #find indices of bgrange
            ridx = closest_index(self.data['wl'],bgrange[0])
            lidx = closest_index(self.data['wl'],bgrange[1])
            
            self.bg = np.mean(self.data['abs'][lidx:ridx])
            self.data['abscorr']=self.data['abs']-self.bg
        
        #normalized data
        if normrange is None:
            self.data['absnorm'] = self.data['abscorr']/np.max(self.data['abscorr'])
        else:
            #find indices of normrange
            ridx = closest_index(self.data['wl'],normrange[0])
            lidx = closest_index(self.data['wl'],normrange[1])
            self.data['absnorm'] = self.data['abscorr']/np.max(self.data['abscorr'][lidx:ridx])
        
        plt.figure()
        plt.plot(self.data['wl'],self.data['abscorr'])
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('Absorbance [a.u.]')
        
        print(data)
        
        
    def saveSpectrum(self,outputpath):
        self.data.to_csv(outputpath + "PROCESSED_" + self.name,index = False)

    def get_abs_at_point(self,wl):
        wavelengths = self.data['wl']
        absorbance = self.data['abscorr']
        idx = np.argmin(np.abs(wavelengths - wl))
        return absorbance[idx]
        
        
def closest_index(array, target):
    # Find the index of the element with the smallest difference to the target
    return min(range(len(array)), key=lambda i: abs(array[i] - target))